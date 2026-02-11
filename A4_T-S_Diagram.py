import numpy as np
import matplotlib.pyplot as plt

# constants and cd0 from openvsp
AR = 3.36462
CD0 = 0.01359

# efficiency factors for each configuration
e_clean = 0.80
e_to    = 0.75
e_lf    = 0.70
e_lg    = e_clean 

# delta CD0 for each configuration
dCD0_to   = 0.020
dCD0_lf   = 0.075
dCD0_gear = 0.025

# lift coeff ranges
cl_clean   = np.linspace(-0.9, 0.9, 100)
cl_takeoff_no_gear = np.linspace(-2.0, 2.0, 100)
cl_takeoff_gear = np.linspace(-2.0, 2.0, 100)
cl_landing_no_gear = np.linspace(-2.6, 2.6, 100)
cl_landing_gear = np.linspace(-2.6, 2.6, 100)

# using drag polar eqn CD = (CD0 + dCD0) + (CL^2 / (pi * AR * e))

cd_clean = CD0 + (cl_clean**2 / (np.pi * AR * e_clean))

cd_takeoff_no_gear = (CD0 + dCD0_to) + (cl_takeoff_no_gear**2 / (np.pi * AR * e_to))

cd_takeoff_gear = (CD0 + dCD0_to+dCD0_gear) + (cl_takeoff_gear**2 / (np.pi * AR * e_to))

cd_landing_no_gear = (CD0 + dCD0_lf) + (cl_landing_no_gear**2 / (np.pi * AR * e_lf))

cd_landing_gear = (CD0 + dCD0_lf+dCD0_gear) + (cl_landing_gear**2 / (np.pi * AR * e_lf))

# plots
plt.figure(figsize=(10, 6))
plt.title('Fighter Aircraft Drag Polars', fontsize=14)
plt.xlabel("$C_D$", fontsize=12)
plt.ylabel("$C_L$", fontsize=12)

plt.plot(cd_clean,   cl_clean,   label=f'Clean', linewidth=2)
plt.plot(cd_takeoff_no_gear, cl_takeoff_no_gear, label=f'Takeoff Flaps With Gears Up', linewidth=2)
plt.plot(cd_takeoff_gear, cl_takeoff_gear, label=f'Takeoff Flaps With Gears Down', linewidth=2)
plt.plot(cd_landing_no_gear, cl_landing_no_gear, label=f'Landing Flaps With Gears Up', linewidth=2)
plt.plot(cd_landing_gear, cl_landing_gear, label=f'Landing Flaps With Gears Down', linewidth=2)

plt.grid(True, linestyle=':', alpha=0.6)
plt.legend(loc='best')
plt.show()

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

plt.figure(figsize=(14,8))

plt.rcParams.update({
    "font.size": 14,
    "axes.labelsize": 16,
    "axes.titlesize": 18,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 12
})

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

# ***************************A4 CODE START *********************************
# Based on A3_final_submission.py (All code above)


# Inner Loop: For fixed S and T, find W_0
"""
Inputs:
S: Wing Area [ft^2]
T: Thrust [lbf] or per Engine Thrust: T_0 = T/n_eng
W_payload: payload weight [lbm or slug]
W_crew: crew weight [lbm or slug]

Outputs:
W_0_converged: Converged gross takeoff weight
"""
def calculate_engine_weight(T_0):
    """
    Given: single-engine thrust (T_0) [lbf]
    Find: single engine weight estimate (W_eng) [lbm]
    """
    W_eng_dry = = 0.521 * T_0**0.9
    W_eng_oil = 0.082 * T_0**0.65
    W_eng_rev = 0.034 * T_0
    W_eng_control = 0.26 * T_0**0.5
    W_eng_start = 9.33 * (W_eng_dry/1000) ** 1.078
    W_eng = W_eng_dry + W_eng_oil + W_eng_rev + W_eng_control + W_eng_start
    return W_eng

def calculate_empty_weight(S_wing, S_ht, S_vt, S_wet_fuselage, TOGW, T_0, num_eng):
    W_wing = S_wing * 9 # From Roskam Weight Model Table - 9x multiplier
    W_ht = S_ht * 4
    W_vt = S_vt * 5.3
    W_fuselage = S_wet_fuselage * 4.8
    W_landing_gear = TOGW *0.045 # Landing Gear - Navy
    single_engine_weight = calculate_engine_weight(T_0)
    W_engines = single_engine_weight * num_eng * 1.3
    W_all_else = TOGW * 0.17
    W_empty = W_wing + W_ht + W_vt + W_fuselage + W_landing_gear + W_engines + W_all_else
    return W_empty

def calculate_weight_fraction(L_D_max, R, E, c, V):
    """This function calculates the weight fractions for cruise and loiter/descent phases based on the Breguet range and endurance equations, and also other terms.
    Args:
        L_D_max (float): Maximum lift-to-drag ratio of the aircraft.
        R (float): Range in nautical miles.
        E (float): Endurance in hours.
        c (float): Specific fuel consumption in lb/(lbf hr).
        V (float): Velocity in knots."""
    L_D = 0.94 * L_D_max

    W3_W2 = np.exp((-R*c) / (V*L_D))  # cruise
    # print("Cruise Fuel Fraction (W3/W2): " + str(round(W3_W2, 3)))

    W4_W3 = np.exp((-E*c) / (L_D))    # loiter/descent
    # print("Loiter Fuel Fraction (W4/W3): " + str(round(W4_W3, 3)))

    W1_W0 = 0.970   # engine start & takeoff
    W2_W1 = 0.985   # climb
    W5_W4 = 0.995   # landing

    W5_W0 = W5_W4 * W4_W3 * W3_W2 * W2_W1 * W1_W0
    # print("Final Fuel Fraction (W5/W0): " + str(round(W5_W0, 3)))

    Wf_W0 = (1 - W5_W0) * 1.06    # compute fuel fraction
    # print("Total Fuel Fraction Wf/W0: {:.3f}".format(Wf_W0))

    return Wf_W0

def inner_loop_weight(TOGW_guess, S_wing, S_ht, S_vt, 
                      S_wet_fuselage, num_eng, W_crew, 
                      W_payload, T_0, err=1e-6, max_iter=200):
    
    W0_history = []
    delta = np.inf
    it = 0
    while delta > err and it < max_iter:
        # 1) fuel fraction (could be constant or updated)
        Wf_W0 = calculate_weight_fraction(L_D_max, R, E, c, V)
        W_fuel = Wf_W0 * TOGW_guess

        # 2) empty weight based on current TOGW guess + geometry + thrust
        W_empty = calculate_empty_weight(
            S_wing, S_ht, S_vt, S_wet_fuselage,
            TOGW_guess, T_0, num_engines
        )

        # 3) new gross weight
        W0_new = W_empty + w_crew + w_payload + W_fuel
        W0_history.append(W0_new)

        # 4) convergence check
        delta = abs(W0_new - TOGW_guess) / max(abs(W0_new), 1e-9)

        # 5) update
        TOGW_guess = W0_new
        it += 1

    converged = (delta <= err)
    return TOGW_guess, converged, it, np.array(W0_history)
# End Inner Loop
