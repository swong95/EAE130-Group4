import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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
            # -----------------------------------------
            # For cruise as example:
            if type == 1:
                TW_req = T_W_coef_1
            elif type == 2:
                TW_req = T_W_coef_1 / WS
            elif type == 3:
                TW_req = T_W_coef_1 / WS + T_W_coef_2 * WS
            # For takeoff as example:
            # TW_req = coef_takeoff_constraint*WS

            # # For landing
            # TW_req = coef_landing_constraint

            # # For climb
            # TW_req = coef_1_climb_constraint
            # -----------------------------------------
            
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

T_W_approach_climb = 0.25909968

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

T_W_takeoff_climb = 0.19013610119245802

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

coeff_1_turn = 1.6319978296085682
coeff_2_turn = 0.004136010310430455

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

a2a_coeff_1_cruise = 9.649831911798403
a2a_coeff_2_cruise = 0.0001049213548068663
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

strike_coeff_1_cruise = 9.150428708110605
strike_coeff_2_cruise = 0.00011064765052450515
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





# Deciding if loops converged or not
print('Inner loop never iterated more than ',max(it_w_final), ' times, which is less than the chosen max of 200 meaning the loop converged.')
print('Outer loop never iterated more than ',max(n_iter_T), ' times, which is less than the chosen max of 200 meaning the loop converged.')

# Plot the resulting T vs S curve from the outer loop convergence
T_actual_F18 = .93*127*500
S_actual_F18 = 500
print(f'Actual T for F-18: {T_actual_F18} lbf, Actual S for F-18: {S_actual_F18} ft^2')

plt.figure(figsize=(16,9))
plt.title('Converged T vs S for Approach Climbing Constraint')
plt.xlabel("Wing Area S (ft^2)")
plt.ylabel("Total Thrust T (lbf)")
plt.plot(S_actual_F18, T_actual_F18, label='Actual F/A-18 E/F Super Hornet', marker='x', markersize=10, color='red')
plt.plot(S_wing_grid, T_approach_climb_curve, label='Converged T for Approach Climb Constraint')
plt.plot(S_wing_grid, T_takeoff_climb_curve, label = 'Converged T for Takeoff Climb Constraint')
plt.plot(S_wing_grid, T_turn_curve, label = 'Converged T for Turn Constraint')
plt.plot(S_wing_grid, T_a2a_cruise_curve, label = 'Converged T for Air-to-Air Dash Constraint')
plt.plot(S_wing_grid, T_strike_cruise_curve, label = 'Converged T for Strike Dash Constraint')
plt.legend(loc='best')
plt.grid()
plt.show()

