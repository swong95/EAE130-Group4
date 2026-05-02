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

CD_clean = 0.04378  # Clean configuration

# Weight
Weight_Empty = 26615  # Empty weight (lb)
Weight_MTOW_Strike = 52028  # Maximum Takeoff Weight for Strike (lb)
Weight_MTOW_A2A = 48891  # Maximum Takeoff Weight for A2A (lb)

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


def thrust_at_altitude(T_SL, h_ft, n=1.0):
    """Return thrust scaled by density ratio with lapse exponent n."""
    return T_SL * (density_at_altitude(h_ft) / rho_SL) ** n


def stall_speed_vs_altitude(Weight, h_min_ft=0, h_max_ft=50000, n=500):
    """Return arrays of altitude (ft) and stall speed (ft/s) across altitude range."""
    altitudes = np.linspace(h_min_ft, h_max_ft, n)
    v_stall = np.array([stall_speed(Weight, density_at_altitude(h)) for h in altitudes])
    return altitudes, v_stall


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
fig, ax = plt.subplots()
for label, W in [("MTOW Strike", Weight_MTOW_Strike), ("MTOW A2A", Weight_MTOW_A2A)]:
    alts, vs = stall_speed_vs_altitude(W)
    ax.plot(vs, alts / 1000, label=label)
ax.set_xlabel("Stall Speed (ft/s)")
ax.set_ylabel("Altitude (1000 ft)")
ax.set_title("Stall Speed vs Altitude")
ax.legend()
plt.tight_layout()
plt.show()


# Max speed vs altitude (for testing purposes)
fig, ax = plt.subplots()
for label, T in [("Dry", Thrust_Dry), ("Afterburner", Thrust_Afterburner)]:
    alts, vm = max_speed_vs_altitude(T, CD_clean)
    ax.plot(vm, alts / 1000, label=label)
ax.set_xlabel("Max Speed (ft/s)")
ax.set_ylabel("Altitude (1000 ft)")
ax.set_title("Max Speed vs Altitude")
ax.legend()
plt.tight_layout()
plt.show()


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

ax.set_xlabel("Speed (ft/s)")
ax.set_ylabel("Altitude (1000 ft)")
ax.set_title("Flight Envelope (Stall & Max Speed Boundaries)")
ax.legend()
plt.tight_layout()
plt.show()
