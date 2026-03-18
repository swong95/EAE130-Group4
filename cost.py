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
L_D_max = 12
L_D_norm = 0.866*L_D_max
R = 700*2          # nm
E = 20/60          # hr
c = 0.8            # lb/(lbf*hr)
V = 589*0.85       # knots
FH = (R / V)+ 0.25 # flight hours = range / cruise velocity + standard tolerance allowance time 
BH = FH + 0.3 # block hours. = flight hours + ground ops 

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

# fuel costs
SFC = 0.8 # from A4 code, specific fuel consumption lb/lbfhr 
W_midflight = W0 - 0.5*W_fuel
T_cruise = W_midflight / L_D_norm
rho_fuel = 6.7 # lb/gal gives range of  6.47-7.01 lb/U.S. gallon for JP-* fuel  
fuel_burn = SFC * T_cruise / rho_fuel # average fuel burn rate
FH_YR_AC = 500 # upper range for fighters, flight hrs/ year/ annual costs
fuel_price = 2.14 # dollars/gallon as of march 17.... to be changed due to irl circumstances.
C_annual_fuel = fuel_burn * FH_YR_AC * fuel_price

# maintainence costs
MMH_FH = 20 # upper range for fighter 
Maint_labor_costs = 150 # 150 per hor
C_annual_maint_labor = MMH_FH * FH_YR_AC * Maint_labor_costs
Crew_ratio = 1.1 # from raymer table 18.1
C_less_engine = C_flyaway - C_eng 
C_mat = 3.3*CEF*Crew_ratio*(C_less_engine/10**6)+7.04+(58*((C_eng/10**6)-13)) * FH_YR_AC  # annual material costs
C_annual_maint = C_annual_maint_labor + C_mat

# crew costs 
C_crew_BH = ((35*(V_max*W0/10**5)**0.3)+84) * CEF # dollar per block hour pg 512 raymer eqn 18.10 for 2 man crew
C_annual_crew = C_crew_BH * (BH/FH) * FH_YR_AC

DOC = C_annual_crew + C_annual_fuel + C_annual_maint

Total_DOC = Total_cost + DOC
# Print Results
# ----------------------
print(f"Estimated Empty Weight: {W_empty:,.2f} lb")
print(f"Fuel Fraction Wf/W0: {Wf_W0:.3f}")
print(f"Fuel Weight: {W_fuel:,.2f} lb")
print(f"Takeoff Gross Weight (W0): {W0:,.2f} lb\n")

print(f"Total RDT&E Costs: ${C_RDTE:,.2f}")
print(f"Unit Flyaway Cost: ${C_flyaway:,.2f}")
print(f"Total Costs, RDT&E and Flyaway for ({Quant_p} units): ${Total_cost:,.2f}")
print(f"Annual Fuel Costs per unit: ${C_annual_fuel:,.2f}")
print(f"Annual Crew Costs per unit: ${C_annual_crew:,.2f}")
print(f"Annual Maintainence Labor Costs per unit: ${C_annual_maint:,.2f}")
print(f"Annual Direct Operating Cost per unit: ${DOC:,.2f}")
print(f"Total Costs, RDT&E, Flyaway and DOC for ({Quant_p} units): ${Total_DOC:,.2f}")
