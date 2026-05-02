import numpy as np
import matplotlib.pyplot as plt
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
turn_alt_ft = 20000          # Altitude for turn rate constraint (ft)
turn_rate_req_deg  = 8.0     # Minimum required sustained turn rate (deg/s)
turn_rate_des_deg  = 10.0    # Desired sustained turn rate (deg/s)
# Mid-mission weight: average of strike takeoff and empty
Weight_mid_mission = (Weight_MTOW_Strike + Weight_Empty) / 2


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
    """Return arrays of altitude (ft) and 1g stall speed (ft/s) across altitude range."""
    altitudes = np.linspace(h_min_ft, h_max_ft, n)
    v_stall = np.array([stall_speed(Weight, density_at_altitude(h)) for h in altitudes])
    return altitudes, v_stall


def corner_velocity_vs_altitude(Weight, g_load, h_min_ft=0, h_max_ft=58500, n=500):
    """Return arrays of altitude (ft) and corner velocity (ft/s) for a given g-load."""
    altitudes = np.linspace(h_min_ft, h_max_ft, n)
    v_star = np.array([np.sqrt((2 * g_load * Weight) / (density_at_altitude(h) * S * CL_max))
                       for h in altitudes])
    return altitudes, v_star


def max_speed(Thrust, rho, CD):
    """V_max = sqrt(2T / (rho * S * CD)) in ft/s (parasite drag limited)."""
    return np.sqrt((2 * Thrust) / (rho * S * CD))


def max_speed_vs_altitude(T_SL, CD, h_min_ft=0, h_max_ft=58500, n=500):
    """Return arrays of altitude (ft) and max speed (ft/s) across altitude range."""
    altitudes = np.linspace(h_min_ft, h_max_ft, n)
    v_max = np.array([max_speed(thrust_at_altitude(T_SL, h), density_at_altitude(h), CD)
                      for h in altitudes])
    return altitudes, v_max


def mach_limit_vs_altitude(mach, h_min_ft=0, h_max_ft=58500, n=500):
    """Return arrays of altitude (ft) and speed (ft/s) for a constant Mach limit."""
    altitudes = np.linspace(h_min_ft, h_max_ft, n)
    v_mach = np.array([mach * speed_of_sound_at_altitude(h) for h in altitudes])
    return altitudes, v_mach


def q_limit_vs_altitude(q_max, h_min_ft=0, h_max_ft=58500, n=500):
    """Return arrays of altitude (ft) and speed (ft/s) at the dynamic pressure limit."""
    altitudes = np.linspace(h_min_ft, h_max_ft, n)
    v_q = np.array([np.sqrt(2 * q_max / density_at_altitude(h)) for h in altitudes])
    return altitudes, v_q


def aero_ceiling(T_SL, Weight, h_min_ft=0, h_max_ft=80000, n=4000):
    """
    Find aerodynamic ceiling: lowest altitude where thrust < D_min.
    D_min = 2 * W * sqrt(CD0 * k), where k = 1 / (pi * AR * e).
    Returns ceiling altitude in feet.
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
    Find the speed(s) at which the aircraft sustains the required turn rate at h_ft.
    Uses Raymer Ch. 17: omega = g * sqrt(n_z^2 - 1) / V, solved for n_z, then
    checks whether thrust can sustain that load factor at each candidate speed.

    Load factor from turn rate: n_z = sqrt(1 + (omega * V / g)^2)
    Thrust-limited n_z: n_z_thrust = T(h) / W * (V / V_stall_1g) ... via
        T = D = q*S*(CD0 + k*CL^2), CL = n_z*W / (q*S)
        Solving: n_z_max_thrust = sqrt((T/W * q*S - (q*S)^2 * CD0 / W) / (k * W))
    CLmax-limited n_z: n_z_clmax = q * S * CL_max / W

    Returns list of (V ft/s, n_z) pairs where sustained turn rate >= target.
    Ref: Raymer, "Aircraft Design: A Conceptual Approach", Ch. 17.
    """
    g = 32.174  # ft/s^2
    omega = np.radians(turn_rate_deg_s)  # rad/s
    rho = density_at_altitude(h_ft)
    T_h = thrust_at_altitude(T_SL, h_ft)
    k = 1.0 / (np.pi * AR * e_oswald)

    v_stall_1g = stall_speed(Weight, rho)
    speeds = np.linspace(v_stall_1g, max_speed(T_h, rho, CD_clean), n_v)

    results = []
    for V in speeds:
        q = 0.5 * rho * V**2

        # Load factor required to achieve the target turn rate at this speed
        n_required = np.sqrt(1 + (omega * V / g)**2)

        # CLmax limit on load factor
        n_clmax = (q * S * CL_max) / Weight

        # Thrust limit on load factor: T = D, D = q*S*(CD0 + k*(n*W/(q*S))^2)
        # Rearranging: k*(n*W)^2/(q*S) = T/S*q - CD0*q  => solve for n
        drag_term = (T_h / (q * S)) - CD_clean
        if drag_term <= 0:
            continue  # thrust insufficient to overcome parasite drag at this speed
        n_thrust = np.sqrt(drag_term / k) * (q * S) / Weight

        n_sustained = min(n_clmax, n_thrust)

        if n_sustained >= n_required:
            results.append((V, n_sustained))

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


# --- Combined flight envelope ---
fig, ax = plt.subplots()

# 1g stall boundaries (left edge)
stall_configs = [
    ("Stall - MTOW Strike", Weight_MTOW_Strike, "C0"),
    ("Stall - MTOW A2A",    Weight_MTOW_A2A,    "C1"),
]
for label, W, color in stall_configs:
    alts, vs = stall_speed_vs_altitude(W)
    ax.plot(vs, alts / 1000, label=label, color=color, linestyle="--")

# Max speed boundaries (right edge)
max_configs = [
    ("Max Speed - Dry",         Thrust_Dry,         "C3"),
    ("Max Speed - Afterburner", Thrust_Afterburner, "C4"),
]
for label, T, color in max_configs:
    alts, vm = max_speed_vs_altitude(T, CD_clean)
    ax.plot(vm, alts / 1000, label=label, color=color, linestyle="-")

# Mach 2.3 limit
alts_mach, v_mach = mach_limit_vs_altitude(2.3)
ax.plot(v_mach, alts_mach / 1000, color="C7", linestyle="-.", label="Mach 2.3 Limit")

# Dynamic pressure limit
alts_q, v_q = q_limit_vs_altitude(q_max_psf)
ax.plot(v_q, alts_q / 1000, color="C8", linestyle="-.", label=f"q_max = {q_max_psf} psf")

# Aerodynamic ceilings as horizontal lines
ceil_configs = [
    ("Ceil - Dry/Strike",  Thrust_Dry,         Weight_MTOW_Strike, "C3"),
    ("Ceil - AB/Strike",   Thrust_Afterburner, Weight_MTOW_Strike, "C4"),
    ("Ceil - Dry/A2A",     Thrust_Dry,         Weight_MTOW_A2A,    "C3"),
    ("Ceil - AB/A2A",      Thrust_Afterburner, Weight_MTOW_A2A,    "C4"),
]
for label, T, W, color in ceil_configs:
    h_ceil = aero_ceiling(T, W)
    ax.axhline(h_ceil / 1000, color=color, linestyle=":", alpha=0.6,
               label=f"{label}: {h_ceil / 1000:.1f}k ft")

# Maneuvering stall curves: V_stall_n = V_stall_1g * sqrt(n)
for label, W, color in stall_configs:
    alts, vs_1g = stall_speed_vs_altitude(W)
    for g_load, ls in [(7, "--"), (8, ":")]:
        ax.plot(vs_1g * np.sqrt(g_load), alts / 1000,
                label=f"Stall {g_load}g - {label.split('- ')[1]}",
                color=color, linestyle=ls, alpha=0.6)

# Corner velocity curves and shading (MTOW Strike as design weight)
W_design = Weight_MTOW_Strike
alts_7g, v_star_7g = corner_velocity_vs_altitude(W_design, g_load=7)
alts_8g, v_star_8g = corner_velocity_vs_altitude(W_design, g_load=8)

ax.plot(v_star_7g, alts_7g / 1000, color="C5", linestyle="-", label="Corner V* 7g (MTOW Strike)")
ax.plot(v_star_8g, alts_8g / 1000, color="C6", linestyle="-", label="Corner V* 8g (MTOW Strike)")
ax.fill_betweenx(alts_7g / 1000, v_star_7g, v_star_8g, alpha=0.15, color="C5", label="7g-8g margin")


# Sustained turn rate markers at 20,000 ft
turn_colors = {"Dry": "C3", "Afterburner": "C4"}
for t_label, T in [("Dry", Thrust_Dry), ("Afterburner", Thrust_Afterburner)]:
    for rate_label, rate, marker in [("8 deg/s", turn_rate_req_deg, "^"),
                                     ("10 deg/s", turn_rate_des_deg, "D")]:
        pts = sustained_turn_speed(rate, T, Weight_mid_mission, turn_alt_ft)
        if pts:
            # Plot the speed range as a horizontal bar at 20k ft
            v_lo, _ = pts[0]
            v_hi, _ = pts[-1]
            color = turn_colors[t_label]
            ax.plot([v_lo, v_hi], [turn_alt_ft / 1000, turn_alt_ft / 1000],
                    color=color, linewidth=4, alpha=0.5, solid_capstyle="butt")
            ax.scatter(v_lo, turn_alt_ft / 1000, marker=marker, s=60,
                       color=color, zorder=6,
                       label=f"Turn {rate_label} {t_label} ({v_lo:.0f}-{v_hi:.0f} ft/s)")

# Design flight conditions
design_points = [
    ("Cruise",      0.85, 35000),
    ("A2A Dash",    1.6,  30000),
    ("Strike Dash", 0.85, 0),
]
for dp_label, mach, h_ft in design_points:
    v = mach * speed_of_sound_at_altitude(h_ft)
    ax.scatter(v, h_ft / 1000, marker="*", s=120, zorder=6)
    h_label = "SL" if h_ft == 0 else f"{h_ft // 1000}k ft"
    ax.annotate(f"{dp_label}\nM{mach} / {h_label}",
                xy=(v, h_ft / 1000), xytext=(8, 4),
                textcoords="offset points", fontsize=7)

ax.set_xlabel("Speed (ft/s)")
ax.set_ylabel("Altitude (1000 ft)")
ax.set_title("Flight Envelope")
ax.legend(fontsize=7)
plt.tight_layout()
plt.show()