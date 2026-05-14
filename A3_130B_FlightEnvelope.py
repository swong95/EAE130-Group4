import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from ambiance import Atmosphere

# Aircraft Parameters
S = 600  # Wing area (ft^2)
b = 35.39616  # Wing span (ft)
#AR = b**2 / S  # Aspect ratio
AR = 2.08815  # Aspect ratio
MAC = 21.00147  # Mean Aerodynamic Chord (ft)

CL_max = 1.5  # Maximum Lift Coefficient (ASSUMED PLACEHOLDER VALUE)

CD_clean = 0.01070  # Clean configuration CD0

# Weight
Weight_Empty = 26615  # Empty weight (lb)
Weight_MTOW_Strike = 52028  # Maximum Takeoff Weight for Strike (lb)
Weight_MTOW_A2A = 48891  # Maximum Takeoff Weight for A2A (lb)

Weight_Strike = Weight_MTOW_Strike * 0.7   # Account for fuel burn for takeoff/climb
Weight_A2A = Weight_MTOW_A2A * 0.75        # Account for fuel burn for takeoff/climb

# Thrust (F135-PW-400 Engine)
Thrust_Dry = 28000  # Dry thrust (lbf)
Thrust_Afterburner = 43000  # Afterburner thrust (lbf)

# Structural / aerodynamic limits
q_max_psf = 2133  # Max dynamic pressure (psf) ~102 kPa, typical naval fighter structural limit
e_oswald = 0.8    # Oswald efficiency for swept wing

# Turn rate requirements (Raymer Ch. 17)
turn_alt_ft = 20000
turn_rate_req_deg = 8.0
turn_rate_des_deg = 10.0
Weight_mid_mission = (Weight_MTOW_Strike + Weight_Empty) / 2

X_MAX = 2500  # ft/s x-axis limit


# --- Atmosphere helpers ---

def density_at_altitude(h_ft):
    """Return ISA air density in slug/ft^3 at a given altitude in feet."""
    h_m = h_ft * 0.3048
    rho_si = Atmosphere(h_m).density[0]  # kg/m^3
    return rho_si * 0.00194032  # convert to slug/ft^3


rho_SL = density_at_altitude(0)


def speed_of_sound_at_altitude(h_ft):
    """Return speed of sound in ft/s at a given altitude in feet."""
    h_m = h_ft * 0.3048
    a_si = Atmosphere(h_m).speed_of_sound[0]  # m/s
    return a_si * 3.28084  # convert to ft/s


def thrust_at_altitude(T_SL, h_ft, n=1.0):
    """Return thrust scaled by density ratio with lapse exponent n."""
    return T_SL * (density_at_altitude(h_ft) / rho_SL) ** n


# --- Constraint functions ---

def stall_speed(Weight, rho):
    """V_stall = sqrt(2W / (rho * S * CL_max)) in ft/s."""
    return np.sqrt((2 * Weight) / (rho * S * CL_max))


def stall_speed_vs_altitude(Weight, h_min_ft=0, h_max_ft=58500, n=500):
    altitudes = np.linspace(h_min_ft, h_max_ft, n)
    v_stall = np.array([stall_speed(Weight, density_at_altitude(h)) for h in altitudes])
    return altitudes, v_stall


def corner_velocity_vs_altitude(Weight, g_load, h_min_ft=0, h_max_ft=58500, n=500):
    altitudes = np.linspace(h_min_ft, h_max_ft, n)
    v_star = np.array([np.sqrt((2 * g_load * Weight) / (density_at_altitude(h) * S * CL_max))
                       for h in altitudes])
    return altitudes, v_star


def max_speed(Thrust, rho, CD):
    """V_max = sqrt(2T / (rho * S * CD)) in ft/s (parasite drag limited)."""
    return np.sqrt((2 * Thrust) / (rho * S * CD))


def max_speed_vs_altitude(T_SL, CD, h_min_ft=0, h_max_ft=58500, n=500):
    altitudes = np.linspace(h_min_ft, h_max_ft, n)
    v_max = np.array([max_speed(thrust_at_altitude(T_SL, h), density_at_altitude(h), CD)
                      for h in altitudes])
    return altitudes, v_max


def mach_limit_vs_altitude(mach, h_min_ft=0, h_max_ft=58500, n=500):
    altitudes = np.linspace(h_min_ft, h_max_ft, n)
    v_mach = np.array([mach * speed_of_sound_at_altitude(h) for h in altitudes])
    return altitudes, v_mach


def q_limit_vs_altitude(q_max, h_min_ft=0, h_max_ft=58500, n=500):
    altitudes = np.linspace(h_min_ft, h_max_ft, n)
    v_q = np.array([np.sqrt(2 * q_max / density_at_altitude(h)) for h in altitudes])
    return altitudes, v_q


def aero_ceiling(T_SL, Weight, h_min_ft=0, h_max_ft=80000, n=4000):
    """
    Find aerodynamic ceiling: lowest altitude where thrust < D_min.
    D_min = 2 * W * sqrt(CD0 * k), where k = 1 / (pi * AR * e).
    """
    k = 1.0 / (np.pi * AR * e_oswald)
    D_min = 2 * Weight * np.sqrt(CD_clean * k)
    altitudes = np.linspace(h_min_ft, h_max_ft, n)
    for h in altitudes:
        if thrust_at_altitude(T_SL, h) < D_min:
            return h
    return h_max_ft


def sustained_turn_speed(turn_rate_deg_s, T_SL, Weight, h_ft, n_v=2000):
    """
    Find speeds sustaining the required turn rate at h_ft.
    Ref: Raymer, "Aircraft Design: A Conceptual Approach", Ch. 17.
    """
    g = 32.174  # ft/s^2
    omega = np.radians(turn_rate_deg_s)
    rho = density_at_altitude(h_ft)
    T_h = thrust_at_altitude(T_SL, h_ft)
    k = 1.0 / (np.pi * AR * e_oswald)

    v_stall_1g = stall_speed(Weight, rho)
    speeds = np.linspace(v_stall_1g, max_speed(T_h, rho, CD_clean), n_v)

    results = []
    for V in speeds:
        q = 0.5 * rho * V**2
        n_required = np.sqrt(1 + (omega * V / g)**2)
        n_clmax = (q * S * CL_max) / Weight
        drag_term = (T_h / (q * S)) - CD_clean
        if drag_term <= 0:
            continue
        n_thrust = np.sqrt(drag_term / k) * (q * S) / Weight
        if min(n_clmax, n_thrust) >= n_required:
            results.append((V, min(n_clmax, n_thrust)))
    return results


# --- Diagnostics ---
print("=== Max Speed Diagnostic ===")
for label, T in [("Dry", Thrust_Dry), ("Afterburner", Thrust_Afterburner)]:
    for h_ft in [0, 15000, 30000, 40000, 50000]:
        rho = density_at_altitude(h_ft)
        T_h = thrust_at_altitude(T, h_ft)
        v = max_speed(T_h, rho, CD_clean)
        mach = v / speed_of_sound_at_altitude(h_ft)
        h_label = "SL" if h_ft == 0 else f"{h_ft // 1000}k ft"
        print(f"  {label:12s} @ {h_label:6s}: {v:7.1f} ft/s  (M{mach:.2f})")

print()
print("=== Aerodynamic Ceilings ===")
for label, T in [("Dry", Thrust_Dry), ("Afterburner", Thrust_Afterburner)]:
    for wlabel, W in [("MTOW Strike", Weight_MTOW_Strike), ("MTOW A2A", Weight_MTOW_A2A)]:
        h_ceil = aero_ceiling(T, W)
        print(f"  {label:12s} / {wlabel}: {h_ceil / 1000:.1f}k ft")

print()
print("=== Sustained Turn Rate @ 20,000 ft (mid-mission weight) ===")
for label, T in [("Dry", Thrust_Dry), ("Afterburner", Thrust_Afterburner)]:
    for rate_label, rate in [("8 deg/s (req)", turn_rate_req_deg), ("10 deg/s (des)", turn_rate_des_deg)]:
        pts = sustained_turn_speed(rate, T, Weight_mid_mission, turn_alt_ft)
        if pts:
            v_min, _ = pts[0]
            v_max_turn, _ = pts[-1]
            print(f"  {label:12s} {rate_label}: achievable {v_min:.0f} - {v_max_turn:.0f} ft/s")
        else:
            print(f"  {label:12s} {rate_label}: NOT achievable")
print()


# --- Pre-compute all boundary arrays on shared altitude grid ---
N = 500
h_ceil_AB_A2A = aero_ceiling(Thrust_Afterburner, Weight_MTOW_A2A)
alts = np.linspace(0, h_ceil_AB_A2A, N)
alts_k = alts / 1000  # altitude in 1000 ft for plotting

_, vs_stall_strike = stall_speed_vs_altitude(Weight_MTOW_Strike, h_max_ft=h_ceil_AB_A2A, n=N)
_, vs_stall_a2a    = stall_speed_vs_altitude(Weight_MTOW_A2A,    h_max_ft=h_ceil_AB_A2A, n=N)
_, vm_dry          = max_speed_vs_altitude(Thrust_Dry,         CD_clean, h_max_ft=h_ceil_AB_A2A, n=N)
_, vm_ab           = max_speed_vs_altitude(Thrust_Afterburner, CD_clean, h_max_ft=h_ceil_AB_A2A, n=N)
_, v_mach23        = mach_limit_vs_altitude(2.3, h_max_ft=h_ceil_AB_A2A, n=N)
_, v_q             = q_limit_vs_altitude(q_max_psf, h_max_ft=h_ceil_AB_A2A, n=N)

_, vs_1g_strike = stall_speed_vs_altitude(Weight_MTOW_Strike, h_max_ft=h_ceil_AB_A2A, n=N)
_, vs_1g_a2a    = stall_speed_vs_altitude(Weight_MTOW_A2A,    h_max_ft=h_ceil_AB_A2A, n=N)
_, v_star_7g    = corner_velocity_vs_altitude(Weight_MTOW_Strike, g_load=7, h_max_ft=h_ceil_AB_A2A, n=N)
_, v_star_8g    = corner_velocity_vs_altitude(Weight_MTOW_Strike, g_load=8, h_max_ft=h_ceil_AB_A2A, n=N)

# Ceiling altitudes for internal horizontal lines
h_ceil_dry_strike = aero_ceiling(Thrust_Dry,         Weight_MTOW_Strike)
h_ceil_ab_strike  = aero_ceiling(Thrust_Afterburner, Weight_MTOW_Strike)
h_ceil_dry_a2a    = aero_ceiling(Thrust_Dry,         Weight_MTOW_A2A)

# Envelope boundaries:
#   Left  = A2A stall (lift limit)
#   Right = min(AB max, Mach 2.3, q_max) — whichever is most restrictive at each altitude
# This produces the classic "thumb" shape where q_max cuts off the bottom-right
# (high-speed at low altitude is structurally limited).
v_right = np.minimum(np.minimum(vm_ab, v_mach23), v_q)

# Walk the boundary CCW: up the left edge, then back down the right edge.
env_v = np.concatenate([vs_stall_a2a, v_right[::-1]])
env_h = np.concatenate([alts_k,       alts_k[::-1]])


# --- Plot ---
fig, ax = plt.subplots(figsize=(8, 8), dpi=150)

# Build envelope clip path (in data coordinates)
clip_verts = list(zip(env_v, env_h))
clip_path = Path(clip_verts)
clip_patch = PathPatch(clip_path, transform=ax.transData, facecolor="none", edgecolor="none")
ax.add_patch(clip_patch)

def clipped(artist):
    """Clip an artist (or list of artists) to the envelope boundary."""
    if isinstance(artist, (list, tuple)):
        for a in artist:
            a.set_clip_path(clip_patch)
    else:
        artist.set_clip_path(clip_patch)
    return artist

# Fill main envelope interior (light blue)
ax.fill(env_v, env_h, color="lightsteelblue", alpha=0.25, zorder=0)

# Shade region between Mach 2.3 and afterburner max speed (where AB exceeds Mach limit)
mask = vm_ab > v_mach23
if mask.any():
    ax.fill_betweenx(alts_k, v_mach23, np.minimum(vm_ab, v_q),
                     where=mask, color="orange", alpha=0.20, zorder=1,
                     label="AB > M2.3")

# Outer boundary lines (these define the envelope, no clipping needed)
ax.plot(vs_stall_a2a, alts_k, color="C1", linestyle="-",  lw=1.5, label="Stall (A2A)")
ax.plot(vm_ab,    alts_k, color="C4", linestyle="-",  lw=1.5, label="V_max AB")
ax.plot(v_mach23, alts_k, color="C7", linestyle="-.", lw=1.5, label="M2.3")
ax.plot(v_q,      alts_k, color="C8", linestyle="-",  lw=1.5, label=f"q_max ({q_max_psf} psf)")

# Top ceiling line (forms top of envelope, no clipping)
ax.axhline(h_ceil_AB_A2A / 1000, color="C4", linestyle="-",  lw=1.5, alpha=0.8,
           label=f"Ceiling AB/A2A ({h_ceil_AB_A2A/1000:.0f}k)")

# Internal lines — all clipped to envelope
clipped(ax.plot(vs_stall_strike, alts_k, color="C0", linestyle="--", lw=1.5, label="Stall (Strike)"))
clipped(ax.plot(vm_dry, alts_k, color="C3", linestyle="-", lw=1.5, label="V_max Dry"))

# Internal ceiling lines (clipped to envelope)
for label, h_c, color, ls in [
    ("Ceil Dry/Strike", h_ceil_dry_strike, "C3", ":"),
    ("Ceil AB/Strike",  h_ceil_ab_strike,  "C4", ":"),
    ("Ceil Dry/A2A",    h_ceil_dry_a2a,    "C3", ":"),
]:
    line = ax.axhline(h_c / 1000, color=color, linestyle=ls, lw=1.0, alpha=0.5,
                      label=f"{label} ({h_c/1000:.0f}k)")
    clipped(line)

# Maneuvering stall curves (clipped)
for W, color, label in [(Weight_MTOW_Strike, "C0", "Strike"),
                         (Weight_MTOW_A2A,    "C1", "A2A")]:
    _, vs_1g = stall_speed_vs_altitude(W, h_max_ft=h_ceil_AB_A2A, n=N)
    for g_load, ls in [(7, "--"), (8, ":")]:
        clipped(ax.plot(vs_1g * np.sqrt(g_load), alts_k,
                        color=color, linestyle=ls, lw=1.0, alpha=0.6,
                        label=f"{g_load}g Stall ({label})"))

# Corner velocity curves and shading (clipped)
clipped(ax.plot(v_star_7g, alts_k, color="C5", linestyle="-", lw=1.2, label="V* 7g"))
clipped(ax.plot(v_star_8g, alts_k, color="C6", linestyle="-", lw=1.2, label="V* 8g"))
margin_fill = ax.fill_betweenx(alts_k, v_star_7g, v_star_8g, alpha=0.15, color="C5", label="7-8g band")
clipped(margin_fill)

# Sustained turn rate markers at 20,000 ft (clipped)
turn_colors = {"Dry": "C3", "Afterburner": "C4"}
for t_label, T in [("Dry", Thrust_Dry), ("AB", Thrust_Afterburner)]:
    thrust_key = "Afterburner" if t_label == "AB" else t_label
    for rate_label, rate, marker in [("8°/s", turn_rate_req_deg, "^"),
                                     ("10°/s", turn_rate_des_deg, "D")]:
        pts = sustained_turn_speed(rate, T, Weight_mid_mission, turn_alt_ft)
        if pts:
            v_lo, _ = pts[0]
            v_hi, _ = pts[-1]
            color = turn_colors[thrust_key]
            clipped(ax.plot([v_lo, v_hi], [turn_alt_ft / 1000, turn_alt_ft / 1000],
                            color=color, linewidth=4, alpha=0.5, solid_capstyle="butt"))
            clipped(ax.scatter(v_lo, turn_alt_ft / 1000, marker=marker, s=60, color=color, zorder=6,
                               label=f"{rate_label} {t_label}"))

# Design flight conditions (clipped)
design_points = [
    ("Cruise",      0.85, 35000),
    ("A2A Dash",    1.6,  30000),
    ("Strike Dash", 0.85, 0),
]
for dp_label, mach, h_ft in design_points:
    v = mach * speed_of_sound_at_altitude(h_ft)
    h_plot = max(h_ft, 1000)  # nudge sea-level points up so the star isn't clipped
    clipped(ax.scatter(v, h_plot / 1000, marker="*", s=120, zorder=6))
    h_label = "SL" if h_ft == 0 else f"{h_ft // 1000}k ft"
    ax.annotate(f"{dp_label}\nM{mach} / {h_label}",
                xy=(v, h_plot / 1000), xytext=(8, 4),
                textcoords="offset points", fontsize=7)

ax.set_xlabel("Airspeed (ft/s)")
ax.set_ylabel("Altitude (1000 ft)")
ax.set_xlim(0, X_MAX)
ax.set_ylim(0, h_ceil_AB_A2A / 1000 * 1.05)
ax.set_title("F/A-XX Flight Envelope")
ax.legend(fontsize=6, loc="upper left")

# Secondary x-axis on top showing knots (1 ft/s = 0.592484 kts)
ax_kts = ax.twiny()
ax_kts.set_xlim(0, X_MAX * 0.592484)
ax_kts.set_xlabel("Airspeed (kts)")

plt.tight_layout()
plt.show()