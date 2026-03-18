import math
import numpy as np

# inputs 
W_empty = 25000 # lb
V_max = 936.8 # knots, use that "A Mach 1.6 dash speed at 30,000 ft. is required."
Quant_p = 500 # production quantity
T_max = 43000 # lbs max thurst including afterburner of F-135 
Mach_max = 1.6 # max engine mach
Turb_inlet = 4059.67 # deg R given that it is 3600 F in wikipedia
N_eng_per_ac = 1 # Number of engines per aircraft
N_eng = N_eng_per_ac * Quant_p # total production quantity of engines
FTA = 2 # number of flight test aircraft
C_avionics = 30e6 # cost per unit

# CEF (cost escalation factor) values
base_year = 1986 
then_year = 2026 
b_CEF = 5.17053 + 0.104981*(base_year - 2006) 
t_CEF = 5.17053 + 0.104981*(then_year - 2006) 
CEF = t_CEF / b_CEF

# hours (Total for program)
H_E = 4.86 * W_empty**0.777 * V_max**0.894 * Quant_p**0.163 
H_T = 5.99 * W_empty**0.777 * V_max**0.696 * Quant_p**0.263 
H_M = 7.37 * W_empty**0.820 * V_max**0.484 * Quant_p**0.641 
H_Q = 0.133 * H_M 

# costs (Total for program)
C_delv = (45.42 * W_empty**0.630 * V_max**1.3) * CEF 
C_flt = (1243.03 * W_empty**0.325 * V_max**0.822 * FTA**1.21) * CEF 
C_mfg = (11.0 * W_empty**0.921 * V_max**0.621 * Quant_p**0.799) * CEF 
C_eng = (1548 * (0.043*T_max + 243.25 * Mach_max + 0.969*Turb_inlet-2228)) * CEF # cost PER ENGINE

# rates, idk where i got these in 
R_E = 2.576 * then_year - 5058 
R_T = 2.833 * then_year - 5666 
R_M = 2.316 * then_year - 4552 
R_Q = 2.600 * then_year - 5112 

# RDT&E for total program
C_RDTE = H_E * R_E + H_T * R_T + H_M * R_M + H_Q * R_Q + C_delv + C_flt
print(f"Total RDT&E Costs: {C_RDTE:,.2f}")

# Flyaway is unit cost. 
# Divide total manufacturing labor, quality labor, and materials by Quant_p.
# Add engine cost (C_eng) and avionics.
C_flyaway = ((H_M * R_M + H_Q * R_Q + C_mfg) / Quant_p) + (C_eng * N_eng_per_ac) + C_avionics
print(f"Unit Flyaway Cost: {C_flyaway:,.2f}")

Total_cost = C_RDTE + 500 * C_flyaway
print(f"Total Costs: {Total_cost:,.2f}")