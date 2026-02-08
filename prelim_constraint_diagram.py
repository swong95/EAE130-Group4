import numpy as np
import matplotlib.pyplot as plt

# --- 1. Constants and Atmospheric Data ---
AR = 2.80672
s_ref = 1053.94
S_wet = 2255.877
c_f = 0.0040
C_D_0 = c_f * (S_wet / s_ref)  
e = 0.8                        
rho_SL = 0.002377              # slugs/ft^3
rho_30k = 0.000891             
sigma_30k = rho_30k / rho_SL    
g = 32.174
k = 1/(np.pi * e * AR)                     

# --- 2. Mission Requirements ---
W_to = 31255.0                 
V_dash_knots = 589 * 1.6       
V_dash = V_dash_knots * 1.6878 
q_dash = 0.5 * rho_30k * V_dash**2

# Takeoff/Landing Parameters
s_to_req = 2500                # ft
s_land_req = 3500              # ft
CL_max_to = 2.0                
CL_max_land = 2.4              
Wl_Wto = 0.85

#Climb Requirements
k_s = 1.2
Vs_to = np.sqrt(2*W_to / (rho_SL * s_ref * CL_max_to))

SEROC_to = 200                 # ft/min
SEROC_approach = 500           # ft/min
G_to = (SEROC_to / 60) / (k_s * Vs_to)                     
# Come back to this one
# G_approach = (SEROC_approach / 60) / (k_s * Vs_to)                 

# --- 3. Baseline Design Point (Dry Thrust) ---
# Total dry thrust for 2x F404-GE-402 engines (~11,000 lbf each)
thrust_max_dry = 11000.0 * 2   

ws_baseline = W_to / s_ref
tw_baseline_dry = thrust_max_dry / W_to

# --- 4. Constraint Calculations ---
WS = np.linspace(20, 250, 100)

# A. Takeoff Field Length (TOFL)
K_to = 37.5
TW_tofl = (K_to * WS) / (1.0 * CL_max_to * s_to_req)

# B. Landing Field Length (LFL)
WS_max_land = (s_land_req * 1.0 * CL_max_land) / (80.0 * Wl_Wto) 

# C. Cruise / Dash (M1.6 at 30,000 ft)
thrust_lapse = sigma_30k**0.7 
TW_dash_alt = (q_dash * C_D_0 / WS) + (WS / (q_dash * np.pi * AR * e))
TW_dash_sl = TW_dash_alt / thrust_lapse

# D. Takeoff Climb
L_D_climb = 0.5 * np.sqrt(np.pi * AR * e / C_D_0) 
TW_climb = ((k_s**2 * C_D_0/CL_max_to) + k*(CL_max_to/k_s**2) + G_to) * np.ones_like(WS)

# --- 5. Plotting ---
plt.figure(figsize=(10, 7))

# Explicitly requested title and labels
plt.title('T/W - W/S', fontsize=14)
plt.xlabel('W/S ($lb/ft^2$)', fontsize=12)
plt.ylabel('T/W', fontsize=12)

# Plot Constraints
plt.plot(WS, TW_tofl, label=f'Takeoff Field Length ({s_to_req} ft)', color='blue', lw=2)
plt.plot(WS, TW_dash_sl, label='Cruise', color='purple', lw=2)
plt.plot(WS, TW_climb, label='Takeoff Climb', color='orange', ls='-', lw=2)
plt.axvline(x=WS_max_land, color='red', label=f'Landing Field Length ({s_land_req} ft)', lw=2)

# Plot Baseline Design Point (Dry Power Only)
plt.scatter(ws_baseline, tw_baseline_dry, color='red', edgecolor='black', s=120, zorder=5, 
            label=f'Baseline: T/W={tw_baseline_dry:.2f}, W/S={ws_baseline:.1f}')

# Shading Valid Region
y_min_valid = np.maximum(TW_tofl, TW_dash_sl)
y_min_valid = np.maximum(y_min_valid, TW_climb)
plt.fill_between(WS, y_min_valid, 1.0, where=(WS <= WS_max_land) & (y_min_valid <= 1.0), 
                 color='green', alpha=0.15, label='Design Window')

plt.grid(True, alpha=0.3)

# --- AXIS LIMITS ---
plt.ylim(0, 1.0) 
plt.xlim(0, 250)

# --- LEGEND ---
plt.legend(loc='upper right', frameon=True, shadow=True, fontsize='small')

plt.tight_layout()
plt.show()
