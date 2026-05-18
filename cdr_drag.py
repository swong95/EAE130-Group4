import csv
import math
import os

import numpy as np
import matplotlib
if os.environ.get("MPL_AGG"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Geometry from OpenVSP
cf_c = 0.25             # Flap chord : Wing chord -- both flaps share this ratio
Sf_begin1 = 0.602941    # Flap 1 local span start w.r.t. total span
Sf_end1 = 0.774510      # Flap 1 local span end w.r.t. total span
c1_root = 10.33612      # Xsec with flap 1 root chord (ft)
c1_tip = 3.28876        # Xsec with flap 1 tip chord (ft)
Sf_begin2 = 0.450980    # Flap 2 local span start w.r.t total span 
Sf_end2 = 0.573529      # Flap 2 local span end w.r.t. total span
c2_root = 25.55608      # Xsec with flap 2 root chord (ft)
c2_tip = 10.33612       # Xsec with flap 2 tip chord (ft)
Sw = 600                # Wing area (ft^2)
b = 35.39616            # Wing total span (ft)
c = 17.99347            # Wing total chord (ft)

# Reynolds Number -- Conditions for each config
V_takeoff = 130 * 1.6878        # knots --> ft/s, estimated from catapult diagram assuming MTOW = 70,000 lb and catapult energy level of 200
V_approach = 244.732/1.2        # ft/s (145 knots)
M = 0.85                        # Cruise mach (average of Superhornet's)
V_cruise = M * 994.7            # V = M*a -- a at 30k ft (ft/s)
rho_SL = 0.002377               # slug/ft^3 (air density at sea level)
rho_30k = 0.0008893             # slug/ft^3 (air density at 30,000 ft.)
mu_SL = 0.0000003737            # lbf*s/ft^2 (viscosity at sea level)
mu_30k = 0.0000003106           # lbf*s/ft^2 (viscosity at 30,000 ft.)
Re_SL_takeoff = (rho_SL * V_takeoff * c) / mu_SL
Re_SL_landing = (rho_SL * V_approach * c) / mu_SL
Re_30k = (rho_30k * V_cruise * c) / mu_30k

# Mach Numbers & CL at each config 
M_takeoff = V_takeoff/1116.4    # M = V/a -- a at SL
M_landing = V_approach/1116.4   # M = V/a -- a at SL
CL_cruise = 70000 / (0.5 * rho_30k * V_cruise**2 * Sw)      # CL at cruise to find wave drag
CL_takeoff = 70000 / (0.5 * rho_SL * V_takeoff**2 * Sw)     # CL at takeoff to find wave drag
CL_landing = 70000 / (0.5 * rho_SL * V_approach**2 * Sw)    # CL at landing to find wave drag

# Flapped Wing Area (ft^2)
delta_Sf1 = (Sf_end1 - Sf_begin1) * b/2     # Flap 1 span (ft)
delta_Sf2 = (Sf_end2 - Sf_begin2) * b/2     # Flap 2 span (ft)
cf1 = cf_c * (c1_root - Sf_end1 * (c1_root - c1_tip))       # Flapped wing root chord at its spanwise location and tapered -- Flap 1
cf2 = cf_c * (c1_root - Sf_begin1 * (c1_root - c1_tip))     # Flapped wing tip chord at its spanwise location and tapered -- Flap 1
cf3 = cf_c * (c2_root - Sf_end2 * (c2_root - c2_tip))       # Flapped wing root chord at its spanwise location and tapered -- Flap 2
cf4 = cf_c * (c2_root - Sf_begin2 * (c2_root - c2_tip))     # Flapped wing tip chord at its spanwise location and tapered -- Flap 2
Sf1 = (cf1 + cf2)/2 * delta_Sf1     # Trapezoid area = (a+b)/2 * h -- Flap 1 area
Sf2 = (cf3 + cf4)/2 * delta_Sf2     # Flap 2 area (ft^2)
Sf = 2* (Sf1 + Sf2)                 # Total flap area on whole wing (ft^2)

# Zero-lift Drag from OpenVSP
CD0_clean = 0.00336
CD0_takeoff_no_gears = 0.005
CD0_takeoff_gears = 0.007
CD0_landing_no_gears = 0.01103
CD0_landing_gears = 0.0122

# Wave Drag from OpenVSP -- CDwave = CDtot - CD0 - CDi from OpenVSP at CL config
CDwave_clean = -0.16-0.003-(-0.17)
CDwave_takeoff_no_gears = 0     # M = 0.2 < 0.8 so no shock --> no wave drag
CDwave_takeoff_gears = 0        # M = 0.2 < 0.8 so no shock --> no wave drag
CDwave_landing_no_gears = 0     # M = 0.18 < 0.8 so no shock --> no wave drag
CDwave_landing_gears = 0        # M = 0.18 < 0.8 so no shock --> no wave drag

# Lift-induced Drag -- tabulated CDi vs CL from AVL alpha sweep (alpha_sweep.py).
def _load_cdi_table(csv_path):
    cls, cdis = [], []
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            cls.append(float(row["CL"]))
            cdis.append(float(row["CDi"]))
    order = np.argsort(cls)
    return np.asarray(cls)[order], np.asarray(cdis)[order]

_HERE = os.path.dirname(os.path.abspath(__file__))
CL_cruise_tab,  CDi_cruise_tab  = _load_cdi_table(os.path.join(_HERE, "alpha_sweep_cruise.csv"))
CL_takeoff_tab, CDi_takeoff_tab = _load_cdi_table(os.path.join(_HERE, "alpha_sweep_takeoff.csv"))
CL_landing_tab, CDi_landing_tab = _load_cdi_table(os.path.join(_HERE, "alpha_sweep_landing.csv"))

def CDi_of(CL, CL_tab, CDi_tab):
    return np.interp(CL, CL_tab, CDi_tab)

# Extrapolates the CL and CD that are outside of VSPAero sweeps
#def fit_cdi(CL_tab, CDi_tab):
#    return np.polyfit(CL_tab, CDi_tab, 2)

#def CDi_of(CL, coeffs):
#    return np.polyval(coeffs, CL)

#coeffs_clean   = fit_cdi(CL_cruise_tab,  CDi_cruise_tab)
#coeffs_takeoff = fit_cdi(CL_takeoff_tab, CDi_takeoff_tab)
#coeffs_landing = fit_cdi(CL_landing_tab, CDi_landing_tab)

# CDi at the design flight condition
CDi_clean   = float(CDi_of(CL_cruise,   CL_cruise_tab,  CDi_cruise_tab))
CDi_takeoff = float(CDi_of(CL_takeoff,  CL_takeoff_tab, CDi_takeoff_tab))
CDi_landing = float(CDi_of(CL_landing,  CL_landing_tab, CDi_landing_tab))

#CDi_clean   = float(CDi_of(CL_cruise, coeffs_clean))
#CDi_takeoff = float(CDi_of(CL_takeoff, coeffs_takeoff))
#CDi_landing = float(CDi_of(CL_landing, coeffs_landing))

# Trimmed Drag from AVL -- CDtrim = CDtot,trimmed - CDtot,untrimmed
CD_trim_clean = 0.03342 - 0.00594
CD_trim_takeoff = 0.16855 - 0.09049
CD_trim_landing = 0.06248 - 0.00798         # Defl 32.5
# CD_trim_landing = 0.01318 - 0.00164       # Defl 15

# Flap Deflection Angle for Takeoff -- from Raymar
deltaf_takeoff = 30
deltaf_takeoff_rad = math.radians(deltaf_takeoff)

# Flap Deflection Angle for Landing -- from Raymar 
deltaf_landing = 65
deltaf_landing_rad = math.radians(deltaf_landing)

# Flap Drag
CD_flaps_takeoff = 1.7 * (cf_c)**1.38 * (Sf/Sw) * (math.sin(deltaf_takeoff_rad))**2 # Flap drag takeoff -- plain flaps
CD_flaps_landing = 1.7 * (cf_c)**1.38 * (Sf/Sw) * (math.sin(deltaf_landing_rad))**2 # Flap drag landing -- plain flaps
CD_flaps_clean = 0  # No flap deflection in clean config so no flap drag
print(CD_flaps_takeoff)
print(CD_flaps_landing)

# Drag Polar
CD_clean = (CD0_clean + CD_flaps_clean + CD_trim_clean + CDwave_clean) + CDi_clean
CD_takeoff_no_gears = (CD0_takeoff_no_gears + CD_flaps_takeoff + CD_trim_takeoff + CDwave_takeoff_no_gears) + CDi_takeoff
CD_takeoff_gears = (CD0_takeoff_gears + CD_flaps_takeoff + CD_trim_takeoff + CDwave_takeoff_gears) + CDi_takeoff
CD_landing_no_gears = (CD0_landing_no_gears + CD_flaps_landing + CD_trim_landing + CDwave_landing_no_gears) + CDi_landing
CD_landing_gears = (CD0_landing_gears + CD_flaps_landing + CD_trim_landing + CDwave_landing_gears) + CDi_landing

# Print Results
print("Clean CD: ", CD_clean)
print("Takeoff, gears up CD: ", CD_takeoff_no_gears)
print("Takeoff, gears down CD: ", CD_takeoff_gears)
print("Landing, gears up CD: ", CD_landing_no_gears)
print("Landing, gears down CD: ", CD_landing_gears)

# CL range
CL_plot = np.linspace(-0.5, 1.5, 350)

def _polar(CL_grid, CD0, CD_flaps, CD_trim, CDwave, CL_tab, CDi_tab):
    in_range = (CL_grid >= CL_tab[0]) & (CL_grid <= CL_tab[-1])
    CDi = CDi_of(CL_grid, CL_tab, CDi_tab)
    CD = (CD0 + CD_flaps + CD_trim + CDwave) + CDi
    return CL_grid[in_range], CD[in_range]

#camber_clean = -0.005
#camber_takeoff = -0.008
#camber_landing = -0.01

#def _polar(CL_grid, CD0, camber_slope, CD_flaps, CD_trim, CDwave, coeffs):

    # induced drag from AVL fit
#    CDi = CDi_of(CL_grid, coeffs)

    # camber-dependent profile drag shift
#    CD0_camber = CD0 + camber_slope * CL_grid

    # total drag
#    CD = CD0_camber + CD_flaps + CD_trim + CDwave + CDi

#    return CL_grid, CD

# Drag polars -- inviscid term comes from the AVL alpha sweep
CL_clean_p,            CD_clean_polar           = _polar(CL_plot, CD0_clean,            CD_flaps_clean,   CD_trim_clean,   CDwave_clean,            CL_cruise_tab,  CDi_cruise_tab)
CL_takeoff_no_gears_p, CD_takeoff_no_gears_polar = _polar(CL_plot, CD0_takeoff_no_gears, CD_flaps_takeoff, CD_trim_takeoff, CDwave_takeoff_no_gears, CL_takeoff_tab, CDi_takeoff_tab)
CL_takeoff_gears_p,    CD_takeoff_gears_polar    = _polar(CL_plot, CD0_takeoff_gears,    CD_flaps_takeoff, CD_trim_takeoff, CDwave_takeoff_gears,    CL_takeoff_tab, CDi_takeoff_tab)
CL_landing_no_gears_p, CD_landing_no_gears_polar = _polar(CL_plot, CD0_landing_no_gears, CD_flaps_landing, CD_trim_landing, CDwave_landing_no_gears, CL_landing_tab, CDi_landing_tab)
CL_landing_gears_p,    CD_landing_gears_polar    = _polar(CL_plot, CD0_landing_gears,    CD_flaps_landing, CD_trim_landing, CDwave_landing_gears,    CL_landing_tab, CDi_landing_tab)

#CL_clean_p, CD_clean_polar = _polar(CL_plot, CD0_clean, camber_clean, CD_flaps_clean, CD_trim_clean, CDwave_clean, coeffs_clean)
#CL_takeoff_no_gears_p, CD_takeoff_no_gears_polar = _polar(CL_plot, CD0_takeoff_no_gears, camber_takeoff, CD_flaps_takeoff, CD_trim_takeoff, CDwave_takeoff_no_gears, coeffs_takeoff)
#CL_takeoff_gears_p, CD_takeoff_gears_polar = _polar(CL_plot, CD0_takeoff_gears, camber_takeoff, CD_flaps_takeoff, CD_trim_takeoff, CDwave_takeoff_gears, coeffs_takeoff)
#CL_landing_no_gears_p, CD_landing_no_gears_polar = _polar(CL_plot, CD0_landing_no_gears, camber_landing, CD_flaps_landing, CD_trim_landing, CDwave_landing_no_gears, coeffs_landing)
#CL_landing_gears_p, CD_landing_gears_polar = _polar( CL_plot, CD0_landing_gears, camber_landing, CD_flaps_landing, CD_trim_landing, CDwave_landing_gears, coeffs_landing)

# Plot
plt.figure(figsize=(8,6))
plt.plot(CD_clean_polar,           CL_clean_p,            label="Clean")
plt.plot(CD_takeoff_no_gears_polar, CL_takeoff_no_gears_p, label="Takeoff (Gear Up)")
plt.plot(CD_takeoff_gears_polar,    CL_takeoff_gears_p,    label="Takeoff (Gear Down)")
plt.plot(CD_landing_no_gears_polar, CL_landing_no_gears_p, label="Landing (Gear Up)")
plt.plot(CD_landing_gears_polar,    CL_landing_gears_p,    label="Landing (Gear Down)")

plt.xlabel("$C_D$")
plt.ylabel("$C_L$")
plt.title("Fighter Aircraft Concept Drag Polars")
plt.legend()
plt.grid(True)
plt.xlim(0,1.2)
plt.ylim(-1,2)

plt.savefig(os.path.join(_HERE, "drag_polar_big_text.png"), dpi=140, bbox_inches="tight")
plt.show()
