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

# Drag Polar
# plt.figure(figsize=(10, 6))
# plt.title('Fighter Aircraft Drag Polars', fontsize=14)
# plt.xlabel("$C_D$", fontsize=12)
# plt.ylabel("$C_L$", fontsize=12)

# plt.plot(cd_clean,   cl_clean,   label=f'Clean', linewidth=2)
# plt.plot(cd_takeoff_no_gear, cl_takeoff_no_gear, label=f'Takeoff Flaps With Gears Up', linewidth=2)
# plt.plot(cd_takeoff_gear, cl_takeoff_gear, label=f'Takeoff Flaps With Gears Down', linewidth=2)
# plt.plot(cd_landing_no_gear, cl_landing_no_gear, label=f'Landing Flaps With Gears Up', linewidth=2)
# plt.plot(cd_landing_gear, cl_landing_gear, label=f'Landing Flaps With Gears Down', linewidth=2)

# plt.grid(True, linestyle=':', alpha=0.6)
# plt.legend(loc='best')
# plt.show()

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


# Cruise / Dash
# Using dash speed for air to air mission since assignment asks for dash speed
V_dash_a2a = 589* 1.6878*1.6 ## ft/s (Speed) [Using Ma = 1.6 at 30,000 ft dash speed for Air to Air Combat]
                     ## Speed of soud pulled from Engineers Edge Table for 30,000 ft
V_dash_strike = 661 * 1.6878 * .85 ## ft/s (Speed) [Using Ma = .85 dash speed for Strike]
                     ## Speed of sound pulled from Engineers Edge Table for sea level
rho_a2a = 0.000891  # slugs/ft^3 (air density at 30,000 ft)
q = 1/2 * rho_a2a * V_dash_a2a**2
q_strike = 1/2 * rho_SL * V_dash_strike**2
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
TW_cruise_strike = (q_strike * C_D_0_calc) / WS + (WS) / (q_strike * np.pi * AR * e_dash)
a2a_cruise_coef_1 = q * C_D_0_calc
a2a_cruise_coef_2 = 1 / (q * np.pi * AR * e_dash)
strike_cruise_coef_1 = q_strike * C_D_0_calc
strike_cruise_coef_2 = 1 / (q_strike * np.pi * AR * e_dash)
print("a2a coefs 1 & 2: ", a2a_cruise_coef_1, " & ", a2a_cruise_coef_2)
print("strike coefs 1 & 2: ", strike_cruise_coef_1, " & ", strike_cruise_coef_2)

plt.figure(figsize=(14,8))

plt.rcParams.update({
    "font.size": 14,
    "axes.labelsize": 16,
    "axes.titlesize": 18,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 12
})

# plt.title('T/W - W/S')
# plt.xlabel("W/S $(lb/ft^2)$")
# plt.ylabel("T/W")
# plt.plot(WS, TW_takeoff_climb, label='Takeoff Climb Constraint', linestyle='-', linewidth=2)
# plt.plot(WS, TW_approach_climb, label='Approach Climb Constraint', linestyle='-', linewidth=2)
# plt.plot(WS, TW_cruise_a2a, label='Cruise Constraint (Air to Air)', linestyle='-', linewidth=2)
# plt.plot(WS, TW_cruise_strike, label='Cruise Constraint (Strike)', linestyle='-', linewidth=2)
# plt.plot(WS, TW_turn, label='Sustained Turn Constraint', linestyle='-', linewidth=2)
# ymin, ymax = plt.ylim()
# plt.vlines(TW_takeoff, ymin, ymax, label='Takeoff Constraint', colors='blue', linewidth=2)
# plt.vlines(TW_landing, ymin, ymax, label='Landing Constraint', colors='orange', linewidth=2)
# plt.ylim(0, 1.5)
# plt.legend(loc='best')

# # ----- Shading feasible region -----
# WS_max = min(takeoff_W_S, landing_W_S)     # most restrictive vertical constraint
# TW_req = np.maximum.reduce([
#     TW_takeoff_climb,
#     TW_approach_climb,
#     TW_cruise_a2a,
#     TW_cruise_strike,
#     TW_turn
# ])

# # Make sure ylim is set before shading so we know the top fill boundary
# plt.ylim(0, 1.5)
# y_top = plt.ylim()[1]

# feasible_mask = (WS <= WS_max) & (TW_req <= y_top)

# plt.fill_between(
#     WS, TW_req, y_top,
#     where=feasible_mask,
#     alpha=0.18,
#     label='Feasible Design Space'
# )

# # ----- F/A-18 Super Hornet design point -----
# WS_f18 = 127.0   # lb/ft^2
# TW_f18 = 0.93

# plt.scatter(
#     WS_f18,
#     TW_f18,
#     s=120,
#     marker='*',
#     color='black',
#     zorder=5,
#     label='F/A-18E/F Super Hornet'
# )

# plt.annotate(
#     'F/A-18E/F',
#     (WS_f18, TW_f18),
#     textcoords="offset points",
#     xytext=(8,8),
#     fontsize=9
# )
# # ----- Our chosen design point -----
# WS_design = 72.3
# TW_design = 0.33

# plt.scatter(
#     WS_design,
#     TW_design,
#     s=110,
#     marker='o',
#     color='red',
#     zorder=6,
#     label='Chosen Design'
# )

# plt.annotate(
#     'Chosen Design Point',
#     (WS_design, TW_design),
#     textcoords="offset points",
#     xytext=(8,-12),
#     fontsize=9
# )

# plt.legend(loc='best')

# plt.legend(loc='best')

# plt.show()
# plt.close()

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
# print("Stall Speed (ft/s): ", V_stall)
# print("Landing Engagement Speed (ft/s): ", V_arrest_engage)
# print("Landing Engagement Speed (kts): ", V_arrest_engage * 0.592484)

Force_hook = 150000
Force_hook_avg = 0.8 * Force_hook # Constant Force Assumption

# Distance Calculation
landing_distance = ((landing_weight/g)*V_arrest_engage**2)/(2*Force_hook_avg)
# print("Landing Distance (ft): ", landing_distance)
# nimitz_landing_field_length = 700 #ft
# if landing_distance < nimitz_landing_field_length:
#     print("Arresting gear succesfully stops aircraft")
# else:
#     print("Aircraft falls off carrier during landing")

# ***************************A4 CODE START *********************************
# Based on A3_final_submission.py (All code above)


# Inner Loop: For fixed S and T, find W_0
def calculate_engine_weight(T_0):
    """Calculate the single engine weight based on the given thrust using empirical relationships.
    Args:
        T_0 (float): Thrust in pounds-force (lbf).
    Returns:
        float: Estimated engine weight in pounds (lb).
    """
    W_eng_dry = 0.521 * T_0**0.9
    W_eng_oil = 0.082 * T_0**0.65
    W_eng_rev = 0.034 * T_0
    W_eng_control = 0.26 * T_0**0.5
    W_eng_start = 9.33 * (W_eng_dry/1000) ** 1.078
    W_eng = W_eng_dry + W_eng_oil + W_eng_rev + W_eng_control + W_eng_start
    return W_eng

# # For Test the function with a sample thrust value: GE90 has max thrust of around 85000 lbf and the weight is around 17400 lb.
# T_0 = 85000
# Engine_weight = calculate_engine_weight(T_0)
# print("Engine weight:", Engine_weight, "lb")

    
def calculate_empty_weight(S_wing, S_ht, S_vt, S_wet_fuselage, TOGW, T_0 , num_engines):
    W_wing = S_wing * 9
    W_ht = S_ht * 4
    W_vt = S_vt * 5.3
    W_fuselage = S_wet_fuselage * 4.8
    W_landing_gear = 0.045 * TOGW
    Engine_weight = calculate_engine_weight(T_0)
    W_engines = Engine_weight * num_engines * 1.3
    W_all_else = 0.17 * TOGW
    W_empty = W_wing + W_ht + W_vt + W_fuselage + W_landing_gear + W_engines + W_all_else
    return W_empty

# # For Test the function with sample values: for boeing 777-200 er, it is around 297300 lb
# S_wing = 4605
# S_ht = 1097
# S_vt = 709.9
# S_wet_fuselage = 11354
# TOGW = 580000  # Example value for Takeoff Gross Weight
# T_0 = 95000  # Example value for thrust per engine
# num_engines = 2  # Example number of engines
# estimated_empty_weight = calculate_empty_weight(S_wing, S_ht, S_vt, S_wet_fuselage, TOGW, T_0 , num_engines)
# print("Estimated empty weight:", estimated_empty_weight, "lb")



def calculate_weight_fraction(L_D_max, R, E, c, V):
    """This function calculates the weight fractions for cruise and loiter/descent phases based on the Breguet range and endurance equations, and also other terms.
    Args:
        L_D_max (float): Maximum lift-to-drag ratio of the aircraft.
        R (float): Range in nautical miles.
        E (float): Endurance in hours.
        c (float): Specific fuel consumption in lb/(lbf hr).
        V (float): Velocity in knots."""
    
    L_D = 0.866 * L_D_max

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


def inner_loop_weight(
    TOGW_guess,
    S_wing, S_ht, S_vt, S_wet_fuselage,
    num_engines, w_crew, w_payload, T_0,
    err=1e-6,
    max_iter=200
):
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

def outer_loop_thrust_for_one_constraint(
    S_wing_grid,
    TOGW_guess_init,
    T_total_guess_init,      # total thrust guess (all engines), lbf
    num_engines,
    S_ht, S_vt, S_wet_fuselage,
    W_crew, W_payload,
    type, T_W_coef_1, T_W_coef_2,
    tol_T_rel=1e-3,          
    max_iter_T=100,
    relax=1.0                # optional damping: 0.3~1.0 (use <1 if oscillation)
):
    
    T_total_converged = []
    W0_converged = []
    iter_counts = []
    iter_w_counts = []
    T_total_history_allS = []  # list of arrays (one per S)

    for S_wing in S_wing_grid:

        # Initialize outer loop for this S
        T_total = T_total_guess_init
        T_hist = []

        for k in range(max_iter_T):
            # Convert total thrust to per-engine thrust for the weight model
            T_0 = T_total / num_engines

            # Inner loop: converge weight for (S, T_0)
            W0, wconv, it_w, W0_hist = inner_loop_weight(
                TOGW_guess_init,
                S_wing, S_ht, S_vt, S_wet_fuselage,
                num_engines, W_crew, W_payload, T_0
            )

            # Wing loading from converged weight
            WS = W0 / S_wing

            # Constraint equations: compute required T/W from W/S
            # Type 1 refers to T/W driven constraints like the climb constraints. These are the horizontal 
            # lines on the T/W - W/S constraint diagram and are simply given as constants in terms of T/W
            # Example for takeoff climb: TW_req = T_W_takeoff_climb

            # Type 2 refers to W/S driven constraints like takeoff and landing. These are the vertical lines
            # on the T/W - W/S constraint diagram and are given as constants in terms of W/S. This likely requires
            # a rewritten outer loop to find wing surface area, S instead of T like we already have.
            
            # Type 3 refers to T/W driven constraints like cruise/dash and turn. These are in the form of
            # TW_req = coef_1 / WS + coef_2 * WS.

            if type == 1:
                TW_req = T_W_coef_1
            elif type == 2:
                TW_req = T_W_coef_1 / WS
            elif type == 3:
                TW_req = T_W_coef_1 / WS + T_W_coef_2 * WS
            
            # Required total thrust
            T_req = TW_req * W0

            # Store history
            T_hist.append(T_total)

            # Check outer convergence
            if abs(T_req - T_total) / max(abs(T_total), 1e-9) < tol_T_rel:
                T_total = T_req
                break

            # Update thrust (optionally relaxed damping)
            T_total = (1 - relax) * T_total + relax * T_req

        # Save results for this S
        T_total_converged.append(T_total)
        W0_converged.append(W0)
        iter_counts.append(k+1)
        T_total_history_allS.append(np.array(T_hist))
        iter_w_counts.append(it_w)
        

    return (np.array(T_total_converged),
            np.array(W0_converged),
            np.array(iter_counts),
            T_total_history_allS,
            W0, wconv, iter_w_counts, W0_hist)



# Pilot Weight
num_pilot = 1
avg_wt_person = 200  # lb

W_crew = num_pilot * avg_wt_person  # lb

# Payload Weight
W_AIM120C = 358         # lb (from: https://www.navair.navy.mil/product/AMRAAM)
W_AIM9X  = 186          # lb (from: https://www.navair.navy.mil/product/AIM-9X-Sidewinder)
W_JDAM   = 1015         # lb (from: https://en.wikipedia.org/wiki/Mark_83_bomb)
W_avionics = 2500       # lb (2500 lb avionics/sensors converted to kg)

num_AIM120C = 6
num_AIM9X  = 2
num_JDAM   = 4

# Calculate payload weights for both missions
W_payload_a2a = (num_AIM120C*W_AIM120C + num_AIM9X*W_AIM9X + W_avionics)  # lb
print("W_payload_air2air: " + str(W_payload_a2a) + " lb")

W_payload_strike = (num_JDAM*W_JDAM + num_AIM9X*W_AIM9X + W_avionics)  # lb
print("W_payload_strike: " + str(W_payload_strike) + " lb")

L_D_max = 9
R = 700 * 2
E = 20 / 60
c = 0.8
V = 589 * 0.85

S_ht = 0
S_vt = 136
S_wet_fuselage = 687
num_engines = 2

# Set grid of wing areas to analyze
S_wing_grid = list(range(0, 6000, 1))  # Example range of wing areas to analyze

TOGW_guess_init = 70000  # Initial guess for Takeoff Gross Weight in pounds
T_total_guess_init = 15000 * num_engines  # Initial guess for total thrust in pounds-force

T_W_approach_climb = coef_2_climb

T_approach_climb_curve, W0_curve, n_iter_T, T_hist_allS, W0_final, wconv_final, it_w_final, W0_hist_final = outer_loop_thrust_for_one_constraint(
    S_wing_grid=S_wing_grid,
    TOGW_guess_init=TOGW_guess_init,
    T_total_guess_init=T_total_guess_init,
    num_engines=num_engines,
    S_ht=S_ht, S_vt=S_vt, S_wet_fuselage=S_wet_fuselage,
    W_crew=W_crew, W_payload=W_payload_a2a,
    type = 1,
    T_W_coef_1 = T_W_approach_climb,
    T_W_coef_2 = 0,
    tol_T_rel=1e-6,
    max_iter_T=200,
    relax=1
)

T_W_takeoff_climb = coef_1_climb

T_takeoff_climb_curve, W0_curve, n_iter_T, T_hist_allS, W0_final, wconv_final, it_w_final, W0_hist_final = outer_loop_thrust_for_one_constraint(
    S_wing_grid=S_wing_grid,
    TOGW_guess_init=TOGW_guess_init,
    T_total_guess_init=T_total_guess_init,
    num_engines=num_engines,
    S_ht=S_ht, S_vt=S_vt, S_wet_fuselage=S_wet_fuselage,
    W_crew=W_crew, W_payload=W_payload_a2a,
    type = 1,
    T_W_coef_1 = T_W_takeoff_climb,
    T_W_coef_2 = 0,
    tol_T_rel=1e-6,
    max_iter_T=200,
    relax=1
)

coeff_1_turn = turn_coeff1
coeff_2_turn = turn_coeff2

T_turn_curve, W0_curve, n_iter_T, T_hist_allS, W0_final, wconv_final, it_w_final, W0_hist_final = outer_loop_thrust_for_one_constraint(
    S_wing_grid=S_wing_grid,
    TOGW_guess_init=TOGW_guess_init,
    T_total_guess_init=T_total_guess_init,
    num_engines=num_engines,
    S_ht=S_ht, S_vt=S_vt, S_wet_fuselage=S_wet_fuselage,
    W_crew=W_crew, W_payload=W_payload_a2a,
    type = 3,
    T_W_coef_1 = coeff_1_turn,
    T_W_coef_2 = coeff_2_turn,
    tol_T_rel=1e-6,
    max_iter_T=200,
    relax=1
)

a2a_coeff_1_cruise = a2a_cruise_coef_1
a2a_coeff_2_cruise = a2a_cruise_coef_2
T_a2a_cruise_curve, W0_curve, n_iter_T, T_hist_allS, W0_final, wconv_final, it_w_final, W0_hist_final = outer_loop_thrust_for_one_constraint(
    S_wing_grid=S_wing_grid,
    TOGW_guess_init=TOGW_guess_init,
    T_total_guess_init=T_total_guess_init,
    num_engines=num_engines,
    S_ht=S_ht, S_vt=S_vt, S_wet_fuselage=S_wet_fuselage,
    W_crew=W_crew, W_payload=W_payload_a2a,
    type = 3,
    T_W_coef_1 = a2a_coeff_1_cruise,
    T_W_coef_2 = a2a_coeff_2_cruise,
    tol_T_rel=1e-6,
    max_iter_T=200,
    relax=1
)

strike_coeff_1_cruise = strike_cruise_coef_1
strike_coeff_2_cruise = strike_cruise_coef_2  
T_strike_cruise_curve, W0_curve, n_iter_T, T_hist_allS, W0_final, wconv_final, it_w_final, W0_hist_final = outer_loop_thrust_for_one_constraint(
    S_wing_grid=S_wing_grid,
    TOGW_guess_init=TOGW_guess_init,
    T_total_guess_init=T_total_guess_init,
    num_engines=num_engines,
    S_ht=S_ht, S_vt=S_vt, S_wet_fuselage=S_wet_fuselage,
    W_crew=W_crew, W_payload=W_payload_a2a,
    type = 3,
    T_W_coef_1 = strike_coeff_1_cruise,
    T_W_coef_2 = strike_coeff_2_cruise,
    tol_T_rel=1e-6,
    max_iter_T=200,
    relax=1
)

T_takeoff_curve, W0_curve, n_iter_T, T_hist_allS, W0_final, wconv_final, it_w_final, W0_hist_final = outer_loop_thrust_for_one_constraint(
    S_wing_grid=S_wing_grid,
    TOGW_guess_init=TOGW_guess_init,
    T_total_guess_init=T_total_guess_init,
    num_engines=num_engines,
    S_ht=S_ht, S_vt=S_vt, S_wet_fuselage=S_wet_fuselage,
    W_crew=W_crew, W_payload=W_payload_a2a,
    type = 2,
    T_W_coef_1 = T_W_takeoff_climb,
    T_W_coef_2 = 0,
    tol_T_rel=1e-6,
    max_iter_T=200,
    relax=1
)


# ============================================================
# LANDING LENGTH CONSTRAINT  (Kevin — fix this block)
# ============================================================
# The landing constraint limits wing loading at touchdown:
#   W_land / S  <=  landing_W_S
# where W_land = Wl_W0 * W0 (landing weight as a fraction of TOGW).
#
# Unlike the other constraints, landing does NOT constrain T/W —
# it constrains S for a given W0.  So we sweep over thrust values
# (T_landing_grid) and, for each T, find the minimum S that keeps
# the landing wing loading within limits.
#
# Because W0 depends on S (through wing weight in the empty-weight
# model), we iterate a small inner loop until S and W0 are mutually
# consistent on the landing boundary: S = W0 * Wl_W0 / landing_W_S

Wl_W0 = 0.78        # landing-to-TOGW weight ratio (matches line 265)
max_land_iter = 100  # maximum iterations for the S–W0 consistency loop
land_s_tol = 1e-5    # convergence tolerance on S (relative)

# Range of thrust values to sweep — adjust endpoints if the curve
# doesn't cover the region of interest on the plot
T_landing_grid = np.linspace(5000, 120000, 500)

S_landing_curve = []  # minimum S at each thrust level
T_landing_valid  = []  # thrust values for which the loop converged

for T_total_land in T_landing_grid:
    T_0_land = T_total_land / num_engines

    # Start the S guess at the reference wing area
    S_guess = s_ref

    wconv_land = False
    for _ in range(max_land_iter):
        # Evaluate W0 for this (S_guess, T) using the existing inner loop
        W0_land, wconv_land, _, _ = inner_loop_weight(
            TOGW_guess_init,
            S_guess, S_ht, S_vt, S_wet_fuselage,
            num_engines, W_crew, W_payload_a2a, T_0_land
        )

        # Landing boundary condition: S such that W_land / S = landing_W_S
        S_new = W0_land * Wl_W0 / landing_W_S

        # Check convergence on S
        if abs(S_new - S_guess) / max(S_guess, 1e-9) < land_s_tol:
            S_guess = S_new
            break

        S_guess = S_new  # update and repeat

    # Only keep points where the weight inner loop also converged
    if wconv_land:
        S_landing_curve.append(S_guess)
        T_landing_valid.append(T_total_land)

S_landing_curve = np.array(S_landing_curve)
T_landing_valid  = np.array(T_landing_valid)
# ============================================================
# END LANDING LENGTH CONSTRAINT
# ============================================================

# A3 Design Point
TW_design = 0.328
WS_design = 72.30 # lbf/ft^2

T_design = TW_design * 66000
S_design = 66000 / WS_design


# Deciding if loops converged or not
print('Inner loop never iterated more than ',max(it_w_final), ' times, which is less than the chosen max of 200 meaning the loop converged.')
print('Outer loop never iterated more than ',max(n_iter_T), ' times, which is less than the chosen max of 200 meaning the loop converged.')

# Plot the resulting T vs S curve from the outer loop convergence
T_actual_F18 = .93*127*500
S_actual_F18 = 500
print(f'Actual T for F-18: {T_actual_F18} lbf, Actual S for F-18: {S_actual_F18} ft^2')

T_design_twinF414 = T_design # 22,000 lbf per GE F-414-400
S_design_twinF414 = S_design

plt.figure(figsize=(16,9))
plt.title('Converged T vs S for Approach Climbing Constraint')
plt.xlabel("Wing Area S (ft^2)")
plt.ylabel("Total Thrust T (lbf)")
plt.plot(S_actual_F18, T_actual_F18, label='Actual F/A-18 E/F Super Hornet', marker='x', markersize=10, color='red')
plt.plot(S_design_twinF414, T_design_twinF414, label='F/A-XX Design Point', marker='x', markersize=10, color='red')
plt.plot(S_wing_grid, T_approach_climb_curve, label='Approach Climb Constraint')
plt.plot(S_wing_grid, T_takeoff_climb_curve, label = 'Takeoff Climb Constraint')
plt.plot(S_wing_grid, T_turn_curve, label = 'Turn Constraint')
plt.plot(S_wing_grid, T_a2a_cruise_curve, label = 'Air-to-Air Dash Constraint')
plt.plot(S_wing_grid, T_strike_cruise_curve, label = 'Strike Dash Constraint')
plt.plot(S_landing_curve, T_landing_valid, label='Landing Constraint', linewidth=2)  # landing line
plt.xlim(0, 2000)   # realistic wing area range for a fighter (ft^2); F/A-18 is 500 ft^2
plt.ylim(0, 80000) # realistic total thrust range (lbf); adjust if curves are cut off
plt.legend(loc='best')
plt.grid()
plt.show()
