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

    # induced drag
    induced_drag = CL**2 / (np.pi * AR * e)

    # === FORM FACTOR MODEL ===
    Cf = 0.004 # from previous cad 
    Swet_Sref = 1440.5/660
    xc = 0.3 # max thickness location (~30% chord)

    # Using form factor model to account for drag increase due to thickness 
    FF = 1 + 0.6*xc*tc + 100*(tc**4)

    CD0_eff = Cf * FF * Swet_Sref

    CD = CD0_eff + induced_drag

    LD = CL / CD

    # structural scaling
    OEW = OEW_base * (1 + k_struct_AR*(AR - 4)) * (1 - k_struct_tc*(tc - 0.12))

    MTOW = OEW + fuel + payload

    return MTOW, LD


# =============================================================================
# DESIGN SPACE
# =============================================================================
AR_vals = np.linspace(1.5, 4, 8)
tc_vals = np.linspace(0.04, 0.15, 5)

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
    ax.plot(LD_grid[:, j], MTOW_grid[:, j],   # <-- swapped
            color='steelblue', alpha=0.7, linewidth=0.8)

# constant t/c lines
for i in range(len(tc_vals)):
    ax.plot(LD_grid[i, :], MTOW_grid[i, :],   # <-- swapped
            color='darkorange', alpha=0.7, linewidth=0.8)

# labels
for j in range(0, len(AR_vals)):
    ax.text(LD_grid[-1, j], MTOW_grid[-1, j]-100,   # <-- swapped
            f'AR={AR_vals[j]:.1f}',
            fontsize=8, color='steelblue')

for i in range(0, len(tc_vals), 2):
    ax.text(LD_grid[i, -1]+0.1, MTOW_grid[i, -1],   # <-- moved to end
            f't/c={tc_vals[i]:.2f}',
            fontsize=8, color='darkorange')

# design point
AR_d = 2.02
tc_d = 0.09
MTOW_d, LD_d = mtow_model(AR_d, tc_d)

ax.scatter(LD_d, MTOW_d, color='red', s=90)   # <-- swapped

# design point label
ax.text(LD_d-.75, MTOW_d+100,
        f'  Design Point\nAR={AR_d}, t/c={tc_d}',
        fontsize=10, color='red')

# axis labels (swapped)
ax.set_xlabel("Cruise L/D")
ax.set_ylabel("MTOW (lb)")

ax.set_title("Trade Study 1: \n Cruise L/D, AR, t/c → MTOW")

ax.grid(True, linestyle='--', alpha=0.4)

plt.tight_layout()
plt.show()