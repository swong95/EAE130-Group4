import numpy as np
import matplotlib.pyplot as plt

# --- Given Constants ---
AR = 3.36462 # Aspect Ratio
s = 59.88512 # Span [ft]
s_ref = 1053.94 # Reference Area [ft^2]
S_wet = 2255.877 # Wetted Area [ft^2]
c_f = 0.0040 # Skin Friction Coefficient 

# --- Other Constants ---
rho_SL = 0.002377  # slugs/ft^3 (air density at sea level)
rho_landing_rhoSL = 1 # density ratio (rho/rho_SL) for landing constraint calculation
rho_20000 = 0.001267 # slugs/ft^3 (air density at 20,000 ft)
g = 32.174 # ft/s^2 (gravity)

# Calculate base C_D_0 from wetted area
C_D_0_calc = c_f * (S_wet / s_ref)


# See drag_polars_fin.py for drag polars and drag polar calculations

# Requirements and Regulations 

# Takeoff Constraint 
# Using catapult constraints for takeoff to determine takeoff velocity -> determine upper bound of wing loading for takeoff constraint

C_L_max_takeoff = 2  # CL_max for takeoff, estimate from Roskam (Table from Slide 11 in 06-PreliminarySizing_Part2.pdf)
V_takeoff = 130 # knots, estimated from catapult diagram assuming MTOW = 70,000 lb and catapult energy level of 200
V_takeoff = V_takeoff * 1.6878 # convert to ft/s
# Using stall condition to determine upper bound of wingloading for takeoff constraint
takeoff_W_S = 0.5 * rho_SL * V_takeoff**2 * C_L_max_takeoff/1.21
print("Takeoff Wing Loading Constraint (W/S): " + str(round(takeoff_W_S, 2)) + " lb/ft^2")

# Landing Field Length Constraint
V_approach = 244.732/1.2  # ft/s (145 knots)
V_arrest = 1.05 * V_approach # ft/s, estimated from RFP as 5% higher than approach speed
C_L_max_landing = 2.1  # CL_max for landing, estimate from Roskam (Table from Slide 11 in 06-PreliminarySizing_Part2.pdf)
landing_W_S = 0.5 * rho_SL * V_arrest**2 * C_L_max_landing

# old assumptions
# s_land = (V_approach**2) / (2 * 32.2 * (0.3))  # landing distance requirement from FAR (ft)
# # Using .3 as deceleration factor estimation
# s_a = 1000  # ft, allowance distance 
# Wl_Wto = 0.65  # Max landing to take off weight fraction

# Climb 
# Took values from example code, might want to double check these later

k_s_to = 1.2
k_s_approach = 1.1 
C_D_0 = 0.01288
W_to = 70000.0 # lbs
SEROC_to = 200 # ft/min
SEROC_approach = 500 # ft/min
Vs_to = np.sqrt(2*W_to / (rho_SL * s_ref * C_L_max_takeoff))
Vs_approach = np.sqrt(2*W_to / (rho_SL * s_ref * C_L_max_landing))
G_to = (SEROC_to / 60) / (k_s_to * Vs_to)
G_approach = (SEROC_approach / 60) / (k_s_approach * Vs_approach)
e = 0.8 # (taking high end of oswald eff factor for clean config as estimate) 
k = 1/(np.pi * e * AR)
coef_1_climb = ((k_s_to**2 * C_D_0/C_L_max_takeoff) + k*(C_L_max_takeoff/k_s_to**2) + G_to) # Climb constraint for takeoff
coef_2_climb = ((k_s_approach**2 * C_D_0/C_L_max_landing) + k*(C_L_max_landing/k_s_approach**2) + G_approach) # Climb constraint for approach
print("Climb gradient coefficient:", coef_1_climb)

# Cruise / Dash
# Using dash speed for air to air mission since assignment asks for dash speed
V_dash_a2a = 589* 1.6878*1.6 ## ft/s (Speed) [Using Ma = 1.6 at 30,000 ft dash speed for Air to Air Combat]
                     ## Speed of soud pulled from Engineers Edge Table for 30,000 
V_dash_strike = 589* 1.6878*.85 ## ft/s (Speed) [Using Ma = .85 dash speed for Strike]
rho_a2a = 0.000891  # slugs/ft^3 (air density at 30,000 ft)
q = 1/2 * rho_a2a * V_dash_a2a**2
e_dash = 0.8  # assuming clean config for dash

# Sustained Turn Constraint
V_turn = 548.538 # ft/s (325 knots) From Raymer of estimated sustained turn speed for fighter aircraft
phi = 8 # deg/s turn rate, RFP
n = np.sqrt(((phi*(np.pi/180))*V_turn/g)**2 + 1) # load factor for sustained turn 
q_turn = 0.5 * rho_20000 * V_turn**2 # turning dynamic pressure at 20,000 ft
turn_coeff1 = q_turn * C_D_0_calc 
turn_coeff2 = n**2 / (q_turn * np.pi * AR * e_dash)

# WS Plot
WS = np.linspace(1,150,300)

TW_takeoff = takeoff_W_S
# TW_landing = (rho_landing_rhoSL * C_L_max_takeoff) / (80 * Wl_Wto) * (s_land + s_a) * np.ones(30)
TW_landing = landing_W_S
print("TW_Landing:", TW_landing)
TW_takeoff_climb = coef_1_climb * np.ones(len(WS))
TW_approach_climb = coef_2_climb * np.ones(len(WS))
TW_turn = turn_coeff1*(1/WS) + turn_coeff2 * WS
TW_cruise_a2a = (q * C_D_0_calc) / WS + (WS) / (q * np.pi * AR * e_dash)
TW_cruise_strike = (V_dash_strike * C_D_0_calc) / WS + (WS) / (V_dash_strike * np.pi * AR * e_dash)

plt.figure(figsize=(8,4))
plt.title('T/W - W/S')
plt.xlabel("W/S $(lb/ft^2)$")
plt.ylabel("T/W")
plt.plot(WS, TW_takeoff_climb, label='Takeoff Climb Constraint', linestyle='-', linewidth=2)
plt.plot(WS, TW_approach_climb, label='Approach Climb Constraint', linestyle='-', linewidth=2)
plt.plot(WS, TW_cruise_a2a, label='Cruise Constraint (Air to Air)', linestyle='-', linewidth=2)
plt.plot(WS, TW_cruise_strike, label='Cruise Constraint (Strike)', linestyle='-', linewidth=2)
plt.plot(WS, TW_turn, label='Sustained Turn Constraint', linestyle='-', linewidth=2)
ymin, ymax = plt.ylim()
plt.vlines(TW_takeoff, ymin, ymax, label='Takeoff Constraint', colors='blue', linewidth=2)
plt.vlines(TW_landing, ymin, ymax, label='Landing Constraint', colors='orange', linewidth=2)
plt.ylim(0, 1.5)
plt.legend(loc='best')

# ----- Shading feasible region -----
WS_max = min(takeoff_W_S, landing_W_S)     # most restrictive vertical constraint
TW_req = np.maximum.reduce([
    TW_takeoff_climb,
    TW_approach_climb,
    TW_cruise_a2a,
    TW_cruise_strike,
    TW_turn
])

# Make sure ylim is set before shading so we know the top fill boundary
plt.ylim(0, 1.5)
y_top = plt.ylim()[1]

feasible_mask = (WS <= WS_max) & (TW_req <= y_top)

plt.fill_between(
    WS, TW_req, y_top,
    where=feasible_mask,
    alpha=0.18,
    label='Feasible Design Space'
)

# ----- F/A-18 Super Hornet design point -----
WS_f18 = 127.0   # lb/ft^2
TW_f18 = 0.93

plt.scatter(
    WS_f18,
    TW_f18,
    s=120,
    marker='*',
    color='black',
    zorder=5,
    label='F/A-18E/F Super Hornet'
)

plt.annotate(
    'F/A-18E/F',
    (WS_f18, TW_f18),
    textcoords="offset points",
    xytext=(8,8),
    fontsize=9
)
# ----- Our chosen design point -----
WS_design = 72.3
TW_design = 0.33

plt.scatter(
    WS_design,
    TW_design,
    s=110,
    marker='o',
    color='red',
    zorder=6,
    label='Chosen Design'
)

plt.annotate(
    'Chosen Design Point',
    (WS_design, TW_design),
    textcoords="offset points",
    xytext=(8,-12),
    fontsize=9
)

plt.legend(loc='best')

plt.legend(loc='best')

plt.show()
plt.close()

# Design point
TW_design = 0.328
WS_design = 72.30 # lbf/ft^2

# Takeoff Distance Check



# Landing Distance Check
# Arresting Gear
C_L_max = 1.8
V_stall= np.sqrt(2*W_to / (rho_SL * s_ref * C_L_max))
design_point_takeoff_weight = WS_design * s_ref
landing_weight = 0.78 * design_point_takeoff_weight # Landing to Takeoff weight ratio from lecture 6 slide 39
# landing weight roughly ~= 51500 lbs
V_arrest_engage = V_stall*1.15
print("Stall Speed (ft/s): ", V_stall)
print("Landing Engagement Speed (ft/s): ", V_arrest_engage)
print("Landing Engagement Speed (kts): ", V_arrest_engage * 0.592484)

Force_hook = 150000
Force_hook_avg = 0.8 * Force_hook # Constant Force Assumption

# Distance Calculation
landing_distance = ((landing_weight/g)*V_arrest_engage**2)/(2*Force_hook_avg)
print("Landing Distance (ft): ", landing_distance)
nimitz_landing_field_length = 700 #ft
if landing_distance < nimitz_landing_field_length:
    print("Arresting gear succesfully stops aircraft")
else:
    print("Aircraft falls off carrier during landing")