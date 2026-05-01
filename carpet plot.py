import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# CARPET PLOT: Flyaway Cost vs Wing Loading
# Sweeps Aspect Ratio (AR) and Wing Area (S_w)
# Weight model based on Raymer Chapter 15 equations
# =============================================================================

# -----------------------------------------------------------------------------
# GEOMETRY & STRUCTURAL PARAMETERS (from OpenVSP model / Raymer)
# -----------------------------------------------------------------------------
K_dw   = 0.768          # delta wing correction factor (Raymer)
K_vs   = 1.0            # non-variable-sweep wing
tc_root = 0.125         # root thickness-to-chord ratio
biglamb = np.deg2rad(40)  # wing sweep at 25% MAC [rad]
lamb    = 3.28876 / 32.88764  # wing taper ratio (c_tip / c_root)
S_csw   = 200           # wing control surface area [ft^2]
wing_fudge = 0.85       # composite material weight reduction factor

# Vertical tail geometry
K_rht    = 1.047        # rolling horizontal tail factor (2 surfaces move independently)
H_tH_v   = 0            # no horizontal tail above fuselage
S_vt     = 100          # vertical tail area [ft^2]
M_val    = 2.0          # Mach number
L_t      = 8.72490      # tail moment arm [ft]
S_r      = 100          # rudder area [ft^2]
AR_vt    = 1.00176      # vertical tail aspect ratio
lamb_vt  = 2.29229 / 9.16916   # vertical tail taper ratio
biglamb_vt = np.deg2rad(45)    # vertical tail sweep at 25% MAC [rad]
tail_fudge = 0.83       # composite material weight reduction factor

# Fuselage geometry
K_dwf     = 0.774       # delta wing + fuselage correction (Raymer)
L_fus     = 40.28139    # fuselage structural length [ft]
D_fus     = 5.4         # fuselage depth/height at max point [ft]
W_fus     = 13.232      # fuselage structural width [ft]
fuselage_fudge = 0.90   # composite material weight reduction factor

# -----------------------------------------------------------------------------
# LANDING GEAR PARAMETERS
# -----------------------------------------------------------------------------
K_cb  = 1.0             # crossbeam factor (conventional)
K_tpg = 1.0             # tripod gear factor
W_l   = 51908 * 0.78    # design landing gross weight [lb]
N_l   = 3 * 1.5         # ultimate landing load factor
L_m   = 60              # main gear length [in] (placeholder)
L_n   = 60 + 2.704*12   # nose gear length [in] (placeholder)
N_nw  = 2               # number of nose wheels

# -----------------------------------------------------------------------------
# PROPULSION PARAMETERS
# -----------------------------------------------------------------------------
N_en   = 1              # number of engines
T      = 43000          # total installed thrust [lbf]
T_e    = 43000          # thrust per engine [lbf]
W_en   = 6422           # engine weight [lb] (F135 estimate)
K_vg   = 1.62           # variable geometry air induction factor
K_d    = 2.75           # duct shape constant (more elliptical)
L_d    = 5.809          # duct length [ft]
L_s    = 5.809          # single duct length [ft]
D_e    = 43 / 12        # engine inlet diameter [ft] (43 in -> ft)
L_sh   = 220 / 12       # engine shroud length [ft] (guess)
L_ec   = 10.503         # length from engine front to cockpit [ft]
AIS_fudge = 0.85        # composite weight reduction for air induction system

# Fuel system
V_i = 120 * 7.48052     # integral tank volume [gal]  (120 ft^3 -> gal)
V_t = V_i               # total fuel volume [gal]
V_p = V_i               # self-sealing "protected" tank volume [gal]
N_t = 1                 # number of fuel tanks
SFC = 0.886             # specific fuel consumption [lb/lbf/hr]

# -----------------------------------------------------------------------------
# SYSTEMS & EQUIPMENT PARAMETERS
# -----------------------------------------------------------------------------
M       = 2             # Mach number (used in flight controls formula)
S_cs    = S_r + S_csw   # total control surface area [ft^2] (rudder + wing CS)
N_s     = 1             # number of flight control systems
N_c     = 1             # number of crew
N_ci    = 1.0           # pilot config factor (1.0 = single pilot)
K_vsh   = 1.0           # variable sweep hydraulics factor
N_u     = 10            # number of hydraulic utility functions (5-15 typical)
K_mc    = 1.45          # electronics mission completion factor
R_kva   = 160           # electrical system rating [kVA] (110-160 for fighters)
L_a     = L_ec          # electrical routing distance, generators to avionics [ft]
N_gen   = N_en          # number of generators
W_uav   = 2500          # uninstalled avionics weight [lb] (from RFP)
S_fw    = 400           # firewall surface area [ft^2] (placeholder)

# Fixed non-structural weights
W_canopy    = 200       # canopy weight [lb] (F-22 rough estimate)
W_pilot     = 200       # pilot + flight gear [lb]

# Landing/catapult gear fractions of W_dg
W_arrgear_frac = 0.008
W_catgear_frac = 0.003

# -----------------------------------------------------------------------------
# PAYLOAD (A2A mission — primary cost driver)
# -----------------------------------------------------------------------------
N_AIM120C = 6;  W_AIM120C = 358    # AIM-120C AMRAAM [lb]
N_AIM9X   = 2;  W_AIM9X   = 186    # AIM-9X Sidewinder [lb]
N_JDAM    = 4;  W_JDAM    = 1015   # JDAM [lb]
W_py_AIM120C = 42.96               # pylon weight per AIM-120C [lb]
W_py_AIM9X   = 22.32               # pylon weight per AIM-9X [lb]
W_py_JDAM    = 121.8               # pylon weight per JDAM [lb]

# -----------------------------------------------------------------------------
# COST MODEL PARAMETERS (Raymer Chapter 18 DAPCA IV)
# -----------------------------------------------------------------------------
Q     = 500             # production quantity
CPI   = 1.43            # cost index (inflation adjustment to base year)
V_vel = 1040            # max velocity [KTAS] for cost equations
CPI_av = 2.96           # avionics cost index

# DAPCA IV coefficients, exponents, and labour rates
# Indices: 0=engineering, 1=tooling, 2=manufacturing, 3=quality control,
#          4=development support, 5=flight test, 6=manufacturing materials
R      = [115, 118, 98, 108, 1, 1, 1]      # labour rates [$/hr] or multipliers
coeff  = [4.86, 5.99, 7.37, 0.076, 91.3, 2498, 22.1]
W_exp  = [0.777, 0.777, 0.82, 0, 0.630, 0.325, 0.921]
V_exp  = [0.894, 0.696, 0.484, 0, 1.3, 0.822, 0.621]
Q_exp  = [0.163, 0.263, 0.641, 0, 0, 0, 0.799]

rdte = [0, 4, 5]        # indices for R&D/test/evaluation cost terms
prod = [1, 2, 3, 6]     # indices for production cost terms

C_engine = 20_400_000   # flyaway engine cost per aircraft [$] (P&W F135 estimate)

# Fixed design gross weight (W_dg).
# NOTE: For a fully converged carpet plot, W_dg should be iterated to
# convergence with W_total inside compute_outputs(). Here it is held fixed
# at the converged design-point value, giving correct trends but not
# fully self-consistent off-design points.
W_dg = 51908.0          # design gross weight [lb]
N_z  = 1.5 * 7          # ultimate load factor (1.5 * limit load factor)

# Fixed fuel and payload weights used in W_total
W_fuel_a2a   = 18967.3  # A2A mission fuel [lb]
W_payload_a2a = N_AIM120C*W_AIM120C + N_AIM9X*W_AIM9X  # A2A payload [lb]
W_useful      = W_pilot + N_AIM120C*W_py_AIM120C + N_AIM9X*W_py_AIM9X + N_JDAM*W_py_JDAM


# =============================================================================
# WEIGHT + COST FUNCTION
# Inputs : AR  — wing aspect ratio
#          S_w — trapezoidal wing area [ft^2]
# Outputs: C_flyaway — flyaway unit cost [$]
#          WS        — wing loading [lb/ft^2]
# =============================================================================
def compute_outputs(AR, S_w):

    # ---- Wing weight (Raymer 15.25, delta-wing form) ----
    W_wing = (0.0103 * K_dw * K_vs
              * (W_dg * N_z)**0.5
              * S_w**0.622
              * AR**0.785
              * tc_root**-0.4
              * (1 + lamb)**0.05
              * np.cos(biglamb)**-1.0
              * S_csw**0.04
              * wing_fudge)

    # ---- Vertical tail weight (Raymer 15.27) ----
    W_vt = (0.452 * K_rht
            * (1 + H_tH_v)**0.5
            * (W_dg * N_z)**0.488
            * S_vt**0.718
            * M_val**0.341
            * L_t**-1.0
            * (1 + S_r/S_vt)**0.348
            * AR_vt**0.223
            * (1 + lamb_vt)**0.25
            * np.cos(biglamb_vt)**-0.323
            * tail_fudge)

    # ---- Fuselage weight (Raymer 15.28, delta-wing form) ----
    W_fuselage = (0.499 * K_dwf
                  * W_dg**0.35
                  * N_z**0.25
                  * L_fus**0.5
                  * D_fus**0.849
                  * (W_fus * 0.685)   # note: Raymer has W^0.685; coded as W*0.685 per original
                  * fuselage_fudge)

    # ---- Landing gear weights (Raymer 15.29–30) ----
    W_mlg = K_cb * K_tpg * (W_l * N_l)**0.25 * L_m**0.973
    W_nlg = (W_l * N_l)**0.290 * L_n**0.5 * N_nw**0.525

    # ---- Propulsion group (Raymer 15.31–39) ----
    W_emounts    = 0.013  * N_en**0.795 * T**0.579  * N_z
    W_firewall   = 1.13   * S_fw
    W_ensect     = 0.01   * W_en**0.717 * N_en      * N_z
    W_ais        = (13.29 * K_vg
                    * L_d**0.643 * K_d**0.182
                    * N_en**1.498
                    * (L_s/L_d)**-0.373
                    * D_e * AIS_fudge)
    W_tailpipe   = 3.5  * D_e * L_sh * N_en
    W_encool     = 4.55 * D_e * L_sh * N_en
    W_oilcool    = 37.82 * N_en**1.023
    W_encontrols = 10.5  * N_en**1.008 * L_ec**0.222
    W_starter    = 0.025 * T_e**0.760  * N_en**0.72
    W_fuelsystanks = (7.45
                      * V_t**0.47
                      * (1 + V_i/V_t)**-0.095
                      * (1 + V_p/V_t)
                      * N_t**0.066
                      * N_en**0.052
                      * (T * SFC / 1000)**0.249)

    # ---- Systems / equipment group (Raymer 15.40–49) ----
    W_flightcont = 36.28  * M**0.003   * S_cs**0.489 * N_s**0.484 * N_c**0.127
    W_instr      = 8.0 + 36.37 * N_en**0.676 * N_t**0.237 + 26.4 * (1 + N_ci)**1.356
    W_hydraulics = 37.23  * K_vsh      * N_u**0.664
    W_elec       = (172.2 * K_mc
                    * R_kva**0.152
                    * N_c**0.10
                    * L_a**0.10
                    * N_gen**0.091)
    W_avionics   = 2.117  * W_uav**0.933
    W_furn       = 217.6  * N_c
    W_acai       = 201.6  * ((W_uav + 200*N_c) / 1000)**0.735
    W_handgear   = 3.2e-4 * W_dg
    W_arrgear    = W_arrgear_frac * W_dg
    W_catgear    = W_catgear_frac * W_dg

    # ---- Group weight totals ----
    W_structures = (W_wing + W_vt + W_fuselage + W_ais + W_canopy
                    + W_mlg + W_nlg + W_ensect + W_handgear + W_arrgear + W_catgear)

    W_propulsion = (W_en + W_fuelsystanks + W_emounts + W_encontrols
                    + W_encool + W_oilcool + W_tailpipe + W_starter)

    W_equip_static = (W_avionics + W_firewall + W_flightcont + W_instr
                      + W_hydraulics + W_elec + W_acai + W_furn)

    # Total weight (A2A mission, full fuel + payload)
    W_total = (W_structures + W_propulsion + W_equip_static
               + W_payload_a2a + W_useful + W_fuel_a2a)

    # ---- Wing loading ----
    WS = W_total / S_w

    # ---- DAPCA IV cost model ----
    avionics_cost = 2000 * W_avionics * CPI_av   # installed avionics cost [$]

    # Production cost terms
    C_prod = []
    for i in prod:
        if i == 3:
            # Quality control hours = 13% of manufacturing hours (index 2)
            mfg_hours = coeff[2] * W_total**W_exp[2] * V_vel**V_exp[2] * Q**Q_exp[2]
            C_prod.append(coeff[i] * R[i] * mfg_hours * CPI)
        else:
            C_prod.append(coeff[i] * W_total**W_exp[i] * V_vel**V_exp[i] * Q**Q_exp[i] * R[i] * CPI)

    # Flyaway cost per aircraft = (total production cost / quantity) + engine + avionics
    C_flyaway = (sum(C_prod) / Q) + C_engine + avionics_cost

    return C_flyaway, WS


# =============================================================================
# GRID SWEEP
# =============================================================================
n_AR = 25
n_S  = 25
AR_vals = np.linspace(1.5, 4.0, n_AR)
S_vals  = np.linspace(400, 900, n_S)

C_grid  = np.zeros((n_AR, n_S))
WS_grid = np.zeros((n_AR, n_S))

for i, A in enumerate(AR_vals):
    for j, S in enumerate(S_vals):
        C_grid[i, j], WS_grid[i, j] = compute_outputs(A, S)


# =============================================================================
# DESIGN POINT
# =============================================================================
AR_design = 2.08815
S_design  = 600.0
C_design, WS_design = compute_outputs(AR_design, S_design)


# =============================================================================
# CARPET PLOT
# =============================================================================
fig, ax = plt.subplots(figsize=(11, 8))

# --- Constant AR contour lines ---
for i in range(n_AR):
    ax.plot(C_grid[i, :], WS_grid[i, :], color='steelblue', linewidth=0.9, alpha=0.7)

# --- Constant S_w contour lines ---
for j in range(n_S):
    ax.plot(C_grid[:, j], WS_grid[:, j], color='darkorange', linewidth=0.9, alpha=0.7)

# --- Labels: AR lines (annotate at the high-S end) ---
x_offset = 0.008 * (C_grid.max() - C_grid.min())
y_offset = 0.008 * (WS_grid.max() - WS_grid.min())
for i in range(0, n_AR, 5):
    ax.text(C_grid[i, -1] + x_offset, WS_grid[i, -1],
            f'AR={AR_vals[i]:.1f}', fontsize=8.5,
            color='steelblue', ha='left', va='center')

# --- Labels: S_w lines (annotate at the low-AR end) ---
for j in range(0, n_S, 5):
    ax.text(C_grid[0, j], WS_grid[0, j] - y_offset,
            f'S={S_vals[j]:.0f}', fontsize=8.5,
            color='darkorange', ha='center', va='top')

# --- Design point marker ---
ax.scatter(C_design, WS_design, color='red', s=90, zorder=5)
ax.annotate(
    f'Design point\nAR={AR_design:.2f}, S={S_design:.0f} ft²\n'
    f'Cost=${C_design/1e6:.2f}M  W/S={WS_design:.1f} lb/ft²',
    xy=(C_design, WS_design),
    xytext=(20, 10), textcoords='offset points',
    fontsize=9, color='red', ha='left', va='bottom',
    arrowprops=dict(arrowstyle='->', color='red', lw=1.2),
    bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='none', alpha=0.85)
)

# --- Legend proxies for carpet line colours ---
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color='steelblue',  lw=1.5, label='Constant AR'),
    Line2D([0], [0], color='darkorange', lw=1.5, label='Constant S_w'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='red',
           markersize=8, label='Design point'),
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=9)

# --- Formatting ---
ax.xaxis.set_major_formatter(
    plt.FuncFormatter(lambda x, _: f'${x/1e6:.0f}M'))
ax.set_xlabel('Flyaway Cost ($ per aircraft)', fontsize=13)
ax.set_ylabel('Wing Loading  W/S  [lb/ft²]', fontsize=13)
ax.set_title('Trade Study: Flyaway Cost vs Wing Loading\n'
             '(Raymer Ch.15 weight model — A2A mission)', fontsize=14)
ax.grid(True, linestyle='--', alpha=0.4)
plt.tight_layout()
plt.savefig('carpet_plot.png', dpi=150)
plt.show()