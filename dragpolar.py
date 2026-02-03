import numpy as np
import matplotlib.pyplot as plt

# --- Given Constants ---
AR = 2.80672
s = 59.88512 
s_ref = 1053.94
S_wet = 2255.877
c_f = 0.0040

# Calculate base C_D_0 from wetted area
C_D_0_calc = c_f * (S_wet / s_ref)
print(f"Calculated base Zero-lift drag coefficient (C_D_0): {C_D_0_calc:.5f}")

# --- Drag Polar Equations Implementation ---
# Define Lift Coefficient ranges for each configuration
cL_clean = np.linspace(-0.9, 0.9, 100)
cL_takeoff = np.linspace(-2.0, 2.0, 100)
cL_landing = np.linspace(-2.6, 2.6, 100)

# Sean - modified drag polors to follow eq C_D = (C_D_0 + Delta_C_D_0) + (1/pi e AR) * C_L^2

# 1. Clean Configuration: CD = 0.01597 + 0.0406 * CL^2 (oringinal eq before modifications below)
clean_drag = C_D_0_calc + 1/(np.pi*0.8*AR) * cL_clean**2

# 2. Takeoff Flaps: CD = 0.02597 + 0.0433 * CL^2 (oringinal eq before modifications below)
takeoff_drag = (C_D_0_calc + 0.01) + 1/(np.pi*0.75*AR) * cL_takeoff**2

# 3. Landing Flaps: CD = 0.07097 + 0.0464 * CL^2 (oringinal eq before modifications below)
landing_flaps_drag = (C_D_0_calc + 0.055) + 1/(np.pi*0.75*AR) * cL_landing**2

# 4. Landing Gear: CD = 0.03097 + 0.0406 * CL^2 (oringinal eq before modifications below)
# (Calculated using the landing CL range)
landing_gear_drag = (C_D_0_calc + 0.015) + 1/(np.pi*AR) * cL_landing**2

# --- Plotting ---
plt.figure(figsize=(10, 6))
plt.title('Aircraft Drag Polars - Multi-Configuration')
plt.xlabel("Drag Coefficient ($C_D$)")
plt.ylabel("Lift Coefficient ($C_L$)")

plt.plot(clean_drag, cL_clean, label='Clean', linestyle='-', linewidth=2)
plt.plot(takeoff_drag, cL_takeoff, label='Takeoff Flaps', linestyle='-', linewidth=2)
plt.plot(landing_flaps_drag, cL_landing, label='Landing Flaps', linestyle='-', linewidth=2)
plt.plot(landing_gear_drag, cL_landing, label='Landing Gear', linestyle='-', linewidth=2)

# Grid and Legend adjustments
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.legend(loc='lower right')

# Set x-limits to focus on the flight envelope
plt.xlim(0, 0.45) 
plt.tight_layout()

plt.show()
plt.close(fig1)

# Requirements and Regulations 
# Takeoff Field Length
BFL_eq = 32 # ft (approximation from catapult stroke length from RFP, not sure if it makes sense but going to run with this for now)
TOP = BFL_eq/37.5  # Takeoff Parameter (TOP) in ft/knot^2
rho_rhoSL = 1.0  # density ratio at sea level
C_L_max_takeoff = 2.2  # max lift coefficient with takeoff flaps

# Landing Field Length
s_land = BFL_eq * .6  # FAR Requirements
s_a = 0 # Allowance distance (set to zero since not sure how this applies to catapult)
V_approach = 244.732  # ft/s (145 knots)
Wl_Wto = 0.65  # Max landing to take off weight fraction

# Climb 
# Took values from example code, might want to double check these later
N_eng = 2
k_s = 1.2 
C_L_max_climb = 2.2
G = 0.012
e = 0.8 # (taking high end of oswald eff factor for clean config as estimate) 
coef_1_climb = (1/0.8) * (N_eng / (N_eng - 1)) * ((k_s**2) / C_L_max_climb * C_D_0_calc + C_L_max_climb / (np.pi * AR * e * k_s**2) + G)
print("Climb gradient coefficient:", coef_1_climb)

# Cruise / Dash
# Using dash speed for air to air mission since assignment asks for dash speed
V_dash_a2a = 589*1.6 ## knots (Speed) [Using Ma = 1.6 at 30,000 ft dash speed for Air to Air Combat]
                     ## Speed of soud pulled from Engineers Edge Table for 30,000 ft
rho_a2a = 0.000891  # slugs/ft^3 (air density at 30,000 ft)
q = 1/2 * rho_a2a * V_dash_a2a**2
e_dash = 0.8  # assuming clean config for dash

# WS Plot
WS = np.linspace(1,300,30)

TW_takeoff = (WS) / (TOP * rho_rhoSL * C_L_max_takeoff)
TW_landing = (rho_rhoSL * C_L_max_takeoff) / (80 * Wl_Wto) * (s_land + s_a) * np.ones(len(WS))
TW_climb = coef_1_climb * np.ones(len(WS))
TW_cruise = (q * C_D_0_calc) / WS + (WS) / (q * np.pi * AR * e_dash)

plt.figure(figsize=(8,4))
plt.title('T/W - W/S')
plt.xlabel("W/S $(lb/ft^2)$")
plt.ylabel("T/W")
plt.plot(WS, TW_takeoff, label='Takeoff field length', linestyle='-', linewidth=2)
plt.plot(TW_landing, np.linspace(0,1,30), label='Landing field length', linestyle='-', linewidth=2)
plt.plot(WS, TW_climb, label='Takeoff climb', linestyle='-', linewidth=2)
plt.plot(WS, TW_cruise, label='Cruise', linestyle='-', linewidth=2)
plt.ylim(0, 0.5)
plt.legend(loc='best')
plt.show()
plt.close(fig2)

print("TOP:", TOP)