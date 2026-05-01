import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# FIXED MISSION PARAMETERS
# =============================================================================
S_ref = 600

CL = 0.6
CD0 = 0.02
e = 0.85

payload = 15500
fuel = 18000
OEW_base = 26000  # baseline empty weight

# structural sensitivity coefficients (conceptual scaling)
k_struct_AR = 0.12
k_struct_tc = 0.25


# =============================================================================
# MTOW MODEL (conceptual sizing closure)
# =============================================================================
def mtow_model(AR, tc):

    # aerodynamic effects
    induced_drag = CL**2 / (np.pi * AR * e)
    CD = CD0 + induced_drag * (1 + 0.2/tc)

    LD = CL / CD

    # structural scaling (proxy for OEW growth)
    OEW = OEW_base * (1 + k_struct_AR*(AR - 4)) * (1 - k_struct_tc*(tc - 0.12))

    MTOW = OEW + fuel + payload

    return MTOW, LD


# =============================================================================
# DESIGN SPACE
# =============================================================================
AR_vals = np.linspace(2.0, 8.0, 40)
tc_vals = np.linspace(0.08, 0.18, 40)

MTOW_grid = np.zeros((len(tc_vals), len(AR_vals)))
LD_grid   = np.zeros((len(tc_vals), len(AR_vals)))


# =============================================================================
# EVALUATION
# =============================================================================
for i, tc in enumerate(tc_vals):
    for j, AR in enumerate(AR_vals):

        MTOW, LD = mtow_model(AR, tc)

        MTOW_grid[i, j] = MTOW
        LD_grid[i, j] = LD


# =============================================================================
# CARPET PLOT
# =============================================================================
fig, ax = plt.subplots(figsize=(11, 8))

# constant AR lines
for j in range(len(AR_vals)):
    ax.plot(MTOW_grid[:, j], LD_grid[:, j],
            color='steelblue', alpha=0.7, linewidth=0.8)

# constant t/c lines
for i in range(len(tc_vals)):
    ax.plot(MTOW_grid[i, :], LD_grid[i, :],
            color='darkorange', alpha=0.7, linewidth=0.8)

# labels (optional)
for j in range(0, len(AR_vals), 8):
    ax.text(MTOW_grid[-1, j], LD_grid[-1, j],
            f'AR={AR_vals[j]:.1f}',
            fontsize=8, color='steelblue')

for i in range(0, len(tc_vals), 8):
    ax.text(MTOW_grid[i, 0], LD_grid[i, 0],
            f't/c={tc_vals[i]:.2f}',
            fontsize=8, color='darkorange')

# design point
AR_d = 5.0
tc_d = 0.12
MTOW_d, LD_d = mtow_model(AR_d, tc_d)

ax.scatter(MTOW_d, LD_d, color='red', s=90)

ax.set_xlabel("MTOW (lb)")
ax.set_ylabel("Cruise L/D")
ax.set_title("Trade Study 1: Aero-Structural Efficiency\nAR vs t/c → System-Level MTOW")

ax.grid(True, linestyle='--', alpha=0.4)

plt.tight_layout()
plt.show()