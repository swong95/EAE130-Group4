import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# FIXED GLOBAL SYSTEM CONSTANTS (no name collisions)
# =============================================================================
M_sys = 2.0  # Mach number used ONLY in systems equations

# =============================================================================
# CORE PARAMETERS (simplified + unchanged structure)
# =============================================================================
K_dw   = 0.768
K_vs   = 1.0
tc_root = 0.125
biglamb = np.deg2rad(40)
lamb    = 3.28876 / 32.88764
S_csw   = 200
wing_fudge = 0.85

K_rht    = 1.047
H_tH_v   = 0
S_vt     = 100
M_val    = 2.0
L_t      = 8.72490
S_r      = 100
AR_vt    = 1.00176
lamb_vt  = 2.29229 / 9.16916
biglamb_vt = np.deg2rad(45)
tail_fudge = 0.83

K_dwf     = 0.774
L_fus     = 40.28139
D_fus     = 5.4
W_fus     = 13.232
fuselage_fudge = 0.90

# Landing gear
K_cb  = 1.0
K_tpg = 1.0
W_l   = 51908 * 0.78
N_l   = 3 * 1.5
L_m   = 60
L_n   = 60 + 2.704*12
N_nw  = 2

# Propulsion
N_en   = 1
T      = 43000
T_e    = 43000
W_en   = 6422
K_vg   = 1.62
K_d    = 2.75
L_d    = 5.809
L_s    = 5.809
D_e    = 43 / 12
L_sh   = 220 / 12
L_ec   = 10.503
AIS_fudge = 0.85

# Fuel system
V_i = 120 * 7.48052
V_t = V_i
V_p = V_i
N_t = 1
SFC = 0.886

# Systems
S_csw = 200
S_cs  = S_csw
N_s   = 1
N_c   = 1
N_ci  = 1.0
K_vsh = 1.0
N_u   = 10
K_mc  = 1.45
R_kva = 160
L_a   = L_ec
N_gen = N_en
W_uav = 2500

# Payload
N_AIM120C = 6; W_AIM120C = 358
N_AIM9X   = 2; W_AIM9X   = 186

W_fuel_a2a = 18967.3

W_dg = 51908.0
N_z  = 1.5 * 7

# =============================================================================
# MAIN MODEL
# =============================================================================
def compute_outputs(AR, S_w):

    # wing span
    b = np.sqrt(AR * S_w)

    # wing weight
    W_wing = (0.0103 * K_dw * K_vs
              * (W_dg * N_z)**0.5
              * S_w**0.622
              * AR**0.785
              * tc_root**-0.4
              * (1 + lamb)**0.05
              * np.cos(biglamb)**-1.0
              * S_csw**0.04
              * wing_fudge)

    # vertical tail
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

    # fuselage
    W_fuselage = (0.499 * K_dwf
                  * W_dg**0.35
                  * N_z**0.25
                  * L_fus**0.5
                  * D_fus**0.849
                  * (W_fus * 0.685)
                  * fuselage_fudge)

    # propulsion group
    W_emounts  = 0.013 * N_en**0.795 * T**0.579 * N_z
    W_firewall = 1.13 * 400
    W_ensect   = 0.01 * W_en**0.717 * N_en * N_z

    W_ais      = (13.29 * K_vg * L_d**0.643 * K_d**0.182 *
                 N_en**1.498 * (L_s/L_d)**-0.373 * D_e * AIS_fudge)

    W_tailpipe = 3.5 * D_e * L_sh * N_en
    W_encool   = 4.55 * D_e * L_sh * N_en
    W_oilcool  = 37.82 * N_en**1.023
    W_starter  = 0.025 * T_e**0.760 * N_en**0.72

    W_fuelsystanks = (7.45 * V_t**0.47 *
                      (1 + V_i/V_t)**-0.095 *
                      (1 + V_p/V_t) *
                      N_t**0.066 *
                      N_en**0.052 *
                      (T * SFC / 1000)**0.249)

    # SYSTEMS (FIXED M BUG HERE)
    W_flightcont = 36.28 * M_sys**0.003 * S_cs**0.489 * N_s**0.484 * N_c**0.127
    W_instr      = 8.0 + 36.37 * N_en**0.676 * N_t**0.237 + 26.4 * (1 + N_ci)**1.356
    W_hydraulics = 37.23 * K_vsh * N_u**0.664
    W_elec       = 172.2 * K_mc * R_kva**0.152 * N_c**0.1 * L_a**0.1 * N_gen**0.091
    W_avionics   = 2.117 * W_uav**0.933

    # totals
    W_struct = W_wing + W_vt + W_fuselage
    W_prop   = (W_emounts + W_firewall + W_ensect + W_tailpipe +
                W_encool + W_oilcool + W_starter + W_fuelsystanks)
    W_sys    = W_flightcont + W_instr + W_hydraulics + W_elec + W_avionics

    W_payload = N_AIM120C*W_AIM120C + N_AIM9X*W_AIM9X

    W_total = W_struct + W_prop + W_sys + W_payload + W_fuel_a2a

    WS = W_total / S_w

    # -------------------------
    # FUEL METRIC (CARPET Z)
    # -------------------------
    fuel_burn = W_fuel_a2a * (W_total / W_dg)

    return fuel_burn, WS


# =============================================================================
# GRID
# =============================================================================
n_AR = 10
n_S  = 10

AR_vals = np.linspace(1.5, 4.0, n_AR)
S_vals  = np.linspace(400, 900, n_S)

F_grid  = np.zeros((n_AR, n_S))
WS_grid = np.zeros((n_AR, n_S))

for i, A in enumerate(AR_vals):
    for j, S in enumerate(S_vals):
        F_grid[i, j], WS_grid[i, j] = compute_outputs(A, S)


# =============================================================================
# DESIGN POINT
# =============================================================================
AR_design = 2.09
S_design  = 600

F_design, WS_design = compute_outputs(AR_design, S_design)


# =============================================================================
# PLOT
# =============================================================================
fig, ax = plt.subplots(figsize=(11, 8))

# --- Constant AR lines ---
for i in range(n_AR):
    ax.plot(F_grid[i, :], WS_grid[i, :],
            color='steelblue', lw=0.9, alpha=0.7)

    # label AR lines (place near right edge)
    mid_idx = -1
    ax.text(F_grid[i, mid_idx],
            WS_grid[i, mid_idx],
            f"AR={AR_vals[i]:.1f}",
            fontsize=8,
            color='steelblue',
            ha='left',
            va='center')

# --- Constant S_w lines ---
for j in range(n_S):
    ax.plot(F_grid[:, j], WS_grid[:, j],
            color='darkorange', lw=0.9, alpha=0.7)

    # label S lines (place near bottom edge)
    ax.text(F_grid[0, j],
            WS_grid[0, j],
            f"S={S_vals[j]:.0f}",
            fontsize=8,
            color='darkorange',
            ha='center',
            va='top')

# --- Design point ---
ax.scatter(F_design, WS_design, color='red', s=90)

ax.annotate(
    f'AR={AR_design:.2f}, S={S_design:.0f}',
    (F_design, WS_design),
    xytext=(15, 10),
    textcoords='offset points',
    color='red',
    bbox=dict(fc='white', alpha=0.85)
)

# --- Labels ---
ax.set_xlabel('Fuel Burn [lb]')
ax.set_ylabel('Wing Loading W/S [lb/ft²]')
ax.set_title('Carpet Plot: Fuel Burn vs Wing Loading (AR × S_w)')
ax.grid(True, linestyle='--', alpha=0.4)

plt.tight_layout()
plt.show()