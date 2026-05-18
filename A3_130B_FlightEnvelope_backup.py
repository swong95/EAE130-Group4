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

CD_clean = 0.01070  # Clean configuration

# Weight
Weight_Empty = 26615  # Empty weight (lb)
Weight_MTOW_Strike = 52028  # Maximum Takeoff Weight for Strike (lb)
Weight_MTOW_A2A = 48891  # Maximum Takeoff Weight for A2A (lb)

Weight_Strike = Weight_MTOW_Strike * 0.7 # Account for fuel burn for takeoff/climb
Weight_A2A = Weight_MTOW_A2A * 0.75 # Account for fuel burn for takeoff/climb

# Thrust (F135-PW-400 Engine)
Thrust_Dry = 28000  # Dry thrust (lbf)
Thrust_Afterburner = 43000  # Afterburner thrust (lbf)


# Stall Condition Constraint (For Stall Speed, vary rho with altitude)
def stall_speed(Weight, rho):
    """
    Calculate stall speed (V_stall) in ft/s given weight (lb) and air density (slug/ft^3).
    V_stall = sqrt((2 * Weight) / (rho * S * CL_max))
    """
    V_stall = np.sqrt((2 * Weight) / (rho * S * CL_max))
    return V_stall


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


def mach_limit_vs_altitude(mach, h_min_ft=0, h_max_ft=50000, n=500):
    """Return arrays of altitude (ft) and speed (ft/s) for a constant Mach limit."""
    altitudes = np.linspace(h_min_ft, h_max_ft, n)
    v_mach = np.array([mach * speed_of_sound_at_altitude(h) for h in altitudes])
    return altitudes, v_mach


def thrust_at_altitude(T_SL, h_ft, n=1.0):
    """Return thrust scaled by density ratio with lapse exponent n."""
    return T_SL * (density_at_altitude(h_ft) / rho_SL) ** n


def stall_speed_vs_altitude(Weight, h_min_ft=0, h_max_ft=50000, n=500):
    """Return arrays of altitude (ft) and stall speed (ft/s) across altitude range."""
    altitudes = np.linspace(h_min_ft, h_max_ft, n)
    v_stall = np.array([stall_speed(Weight, density_at_altitude(h)) for h in altitudes])
    return altitudes, v_stall


def corner_velocity_vs_altitude(Weight, g_load, h_min_ft=0, h_max_ft=50000, n=500):
    """Return arrays of altitude (ft) and corner velocity (ft/s) for a given g-load."""
    altitudes = np.linspace(h_min_ft, h_max_ft, n)
    v_star = np.array([np.sqrt((2 * g_load * Weight) / (density_at_altitude(h) * S * CL_max)) for h in altitudes])
    return altitudes, v_star


# Max Speed Constraint
def max_speed(Thrust, rho, CD):
    """
    Calculate max speed (V_max) in ft/s given thrust (lbf) and air density (slug/ft^3).
    V_max = sqrt((2 * Thrust) / (rho * S * CD))
    """
    V_max = np.sqrt((2 * Thrust) / (rho * S * CD))
    return V_max


def max_speed_vs_altitude(T_SL, CD, h_min_ft=0, h_max_ft=50000, n=500):
    """Return arrays of altitude (ft) and max speed (ft/s) across altitude range."""
    altitudes = np.linspace(h_min_ft, h_max_ft, n)
    v_max = np.array([max_speed(thrust_at_altitude(T_SL, h), density_at_altitude(h), CD) for h in altitudes])
    return altitudes, v_max


# Stall speed vs altitude for each mission weight (for testing purposes)
#fig, ax = plt.subplots()
#for label, W in [("MTOW Strike", Weight_MTOW_Strike), ("MTOW A2A", Weight_MTOW_A2A)]:
#    alts, vs = stall_speed_vs_altitude(W)
#    ax.plot(vs, alts / 1000, label=label)
#ax.set_xlabel("Stall Speed (ft/s)")
#ax.set_ylabel("Altitude (1000 ft)")
#ax.set_title("Stall Speed vs Altitude")
#ax.legend()
#plt.tight_layout()
#plt.show()


# Max speed vs altitude (for testing purposes)
#fig, ax = plt.subplots()
#for label, T in [("Dry", Thrust_Dry), ("Afterburner", Thrust_Afterburner)]:
#    alts, vm = max_speed_vs_altitude(T, CD_clean)
#    ax.plot(vm, alts / 1000, label=label)
#ax.set_xlabel("Max Speed (ft/s)")
#ax.set_ylabel("Altitude (1000 ft)")
#ax.set_title("Max Speed vs Altitude")
#ax.legend()
#plt.tight_layout()
#plt.show()


# Max speed diagnostic output
print("=== Max Speed Diagnostic ===")
for label, T in [("Dry", Thrust_Dry), ("Afterburner", Thrust_Afterburner)]:
    for h_ft in [0, 15000, 30000, 40000, 50000]:
        rho = density_at_altitude(h_ft)
        T_h = thrust_at_altitude(T, h_ft)
        v = max_speed(T_h, rho, CD_clean)
        mach = v / speed_of_sound_at_altitude(h_ft)
        h_label = "SL" if h_ft == 0 else f"{h_ft//1000}k ft"
        print(f"  {label:12s} @ {h_label:6s}: {v:7.1f} ft/s  (M{mach:.2f})")
print()


# Combined flight envelope (stall boundary + max speed boundary)
fig, ax = plt.subplots()

# Stall boundaries (left edge)
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

# Maneuvering stall curves: V_stall_n = V_stall_1g * sqrt(n)
for label, W, color in stall_configs:
    alts, vs_1g = stall_speed_vs_altitude(W)
    for g_load, ls in [(7, "--"), (8, ":")]:
        ax.plot(vs_1g * np.sqrt(g_load), alts / 1000,
                label=f"Stall {g_load}g - {label.split('- ')[1]}",
                color=color, linestyle=ls, alpha=0.6)

# Corner velocity curves and shading (use MTOW Strike as design weight)
W_design = Weight_MTOW_Strike
alts_7g, v_star_7g = corner_velocity_vs_altitude(W_design, g_load=7)
alts_8g, v_star_8g = corner_velocity_vs_altitude(W_design, g_load=8)

ax.plot(v_star_7g, alts_7g / 1000, color="C5", linestyle="-",  label="Corner V* 7g (MTOW Strike)")
ax.plot(v_star_8g, alts_8g / 1000, color="C6", linestyle="-",  label="Corner V* 8g (MTOW Strike)")

# Shade margin between 7g and 8g corner velocity
ax.fill_betweenx(alts_7g / 1000, v_star_7g, v_star_8g, alpha=0.15, color="C5", label="7gâ€“8g margin")

# Mark sea level and 15,000 ft corner velocity points
for g_load, v_star, color in [(7, v_star_7g, "C5"), (8, v_star_8g, "C6")]:
    for h_target, h_label in [(0, "SL"), (15000, "15k ft")]:
        idx = np.argmin(np.abs(alts_7g - h_target))
        ax.scatter(v_star[idx], alts_7g[idx] / 1000, color=color, zorder=5)
        ax.annotate(f"V* {g_load}g @ {h_label}\n{v_star[idx]:.0f} ft/s",
                    xy=(v_star[idx], alts_7g[idx] / 1000),
                    xytext=(8, 6), textcoords="offset points", fontsize=7)

# Design flight conditions
design_points = [
    ("Cruise",       0.85, 35000),
    ("A2A Dash",     1.6,  30000),
    ("Strike Dash",  0.85, 0),
]
for dp_label, mach, h_ft in design_points:
    v = mach * speed_of_sound_at_altitude(h_ft)
    ax.scatter(v, h_ft / 1000, marker="*", s=120, zorder=6)
    h_label = "SL" if h_ft == 0 else f"{h_ft//1000}k ft"
    ax.annotate(f"{dp_label}\nM{mach} / {h_label}",
                xy=(v, h_ft / 1000), xytext=(8, 4),
                textcoords="offset points", fontsize=7)

ax.set_xlabel("Speed (ft/s)")
ax.set_ylabel("Altitude (1000 ft)")
ax.set_title("Flight Envelope (Stall & Max Speed Boundaries)")
ax.legend(fontsize=7)
plt.tight_layout()
plt.show()

