import numpy as np

# ----------------------
# Aircraft & Mission Inputs
# ----------------------
S_wing = 600     # ft^2
S_ht = 0        # ft^2
S_vt = 112.693       # ft^2
S_wet_fuselage = 339.275  # ft^2
TOGW_guess = 70000       # lbs
T_0 = 25000              # lbf per engine
num_engines = 1
num_pilot = 1
avg_wt_person = 200
W_crew = num_pilot * avg_wt_person

# Payload (example)
W_AIM120C = 358
W_AIM9X = 186
W_JDAM = 1015
W_avionics = 2500
num_AIM120C = 6
num_AIM9X = 2
num_JDAM = 4
W_payload = num_AIM120C*W_AIM120C + num_AIM9X*W_AIM9X + W_avionics

# Mission parameters
L_D_max = 10
R = 700*2          # nm
E = 20/60          # hr
c = 0.8            # lb/(lbf*hr)
V = 589*0.85       # knots

# ----------------------
# Functions for Weight/Fuel
# ----------------------
def calculate_engine_weight(T_0):
    W_eng_dry = 0.521 * T_0**0.9
    W_eng_oil = 0.082 * T_0**0.65
    W_eng_rev = 0.034 * T_0
    W_eng_control = 0.26 * T_0**0.5
    W_eng_start = 9.33 * (W_eng_dry/1000)**1.078
    return W_eng_dry + W_eng_oil + W_eng_rev + W_eng_control + W_eng_start

def calculate_empty_weight(S_wing, S_ht, S_vt, S_wet_fuselage, TOGW, T_0, num_engines):
    W_wing = S_wing * 9
    W_ht = S_ht * 4
    W_vt = S_vt * 5.3
    W_fuselage = S_wet_fuselage * 4.8
    W_landing_gear = 0.045 * TOGW
    W_engines = calculate_engine_weight(T_0) * num_engines * 1.3
    W_all_else = 0.17 * TOGW
    return W_wing + W_ht + W_vt + W_fuselage + W_landing_gear + W_engines + W_all_else

def calculate_weight_fraction(L_D_max, R, E, c, V):
    L_D = 0.866 * L_D_max
    W3_W2 = np.exp((-R*c) / (V*L_D))
    W4_W3 = np.exp((-E*c) / L_D)
    W1_W0 = 0.970
    W2_W1 = 0.985
    W5_W4 = 0.995
    W5_W0 = W5_W4 * W4_W3 * W3_W2 * W2_W1 * W1_W0
    Wf_W0 = (1 - W5_W0) * 1.06
    return Wf_W0

# ----------------------
# Weight/Fuel Calculations
# ----------------------
W_empty = calculate_empty_weight(S_wing, S_ht, S_vt, S_wet_fuselage, TOGW_guess, T_0, num_engines)
Wf_W0 = calculate_weight_fraction(L_D_max, R, E, c, V)
W_fuel = Wf_W0 * TOGW_guess
W0 = W_empty + W_crew + W_payload + W_fuel

# ----------------------
# Cost Inputs
# ----------------------      # lb, for cost equations
V_max = 936.8              # knots
Quant_p = 500
T_max = 43000              # lbs thrust (afterburner)
Mach_max = 1.6
Turb_inlet = 4059.67       # deg R
N_eng_per_ac = 1
N_eng = N_eng_per_ac * Quant_p
FTA = 2
C_avionics = 30e6

# Cost Escalation Factor (CEF)
base_year = 1986 
then_year = 2026 
b_CEF = 5.17053 + 0.104981*(base_year - 2006) 
t_CEF = 5.17053 + 0.104981*(then_year - 2006) 
CEF = t_CEF / b_CEF

# Program Hours
H_E = 4.86 * W_empty**0.777 * V_max**0.894 * Quant_p**0.163 
H_T = 5.99 * W_empty**0.777 * V_max**0.696 * Quant_p**0.263 
H_M = 7.37 * W_empty**0.820 * V_max**0.484 * Quant_p**0.641 
H_Q = 0.133 * H_M 

# Costs
C_delv = (45.42 * W_empty**0.630 * V_max**1.3) * CEF 
C_flt = (1243.03 * W_empty**0.325 * V_max**0.822 * FTA**1.21) * CEF 
C_mfg = (11.0 * W_empty**0.921 * V_max**0.621 * Quant_p**0.799) * CEF 
C_eng = (1548 * (0.043*T_max + 243.25 * Mach_max + 0.969*Turb_inlet-2228)) * CEF

# Rates
R_E = 2.576 * then_year - 5058 
R_T = 2.833 * then_year - 5666 
R_M = 2.316 * then_year - 4552 
R_Q = 2.600 * then_year - 5112 

# RDT&E & Flyaway Costs
C_RDTE = H_E * R_E + H_T * R_T + H_M * R_M + H_Q * R_Q + C_delv + C_flt
C_flyaway = ((H_M * R_M + H_Q * R_Q + C_mfg) / Quant_p) + (C_eng * N_eng_per_ac) + C_avionics
Total_cost = C_RDTE + Quant_p * C_flyaway

# ----------------------
# Print Results
# ----------------------
print(f"Estimated Empty Weight: {W_empty:,.2f} lb")
print(f"Fuel Fraction Wf/W0: {Wf_W0:.3f}")
print(f"Fuel Weight: {W_fuel:,.2f} lb")
print(f"Takeoff Gross Weight (W0): {W0:,.2f} lb\n")

print(f"Total RDT&E Costs: ${C_RDTE:,.2f}")
print(f"Unit Flyaway Cost: ${C_flyaway:,.2f}")
print(f"Total Costs ({Quant_p} units): ${Total_cost:,.2f}")