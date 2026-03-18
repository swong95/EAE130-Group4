import math
import numpy as np


def m_to_ft(m):
    return m * 3.28084
def sqm_to_sqft(sq_m):
    return sq_m * 10.7639

# Empennage sizing #
# Given values related to Vertical Tail #
c_vt = 0.09               # vertical tail volume coefficient  
L_vt = 26.4431975      # vertical tail moment arm (m) --> (ft)

# Given values related to Wing #
b_w  = 59.88512      # wing span tip-to-tip (m) --> (ft)
c_w  = 5.91142     # wing chord (ft^2)
S_w  = 1059.85985 # wing area (m^2) --> (ft^2)



# Given values related to Horizontal Tail #
c_ht = 0             # horizontal tail volume coefficient
L_ht = 0     # horizontal tail moment arm (m) --> (ft)


S_vt = (c_vt*b_w*S_w)/L_vt # Estimated Vertical Tail Area (ft^2)
S_ht = 0 # Estimated Horizontal Tail Area (ft^2)


S_vt_actual = 145.24352   # Actual vertical tail area (m^2) --> (ft^2)
S_ht_actual = 0  # Actual horizontal tail area (m^2) --> (ft^2)

print("Estimated Vertical Tail Area = {} ft^2".format(S_vt))
print("Actual Vertical Tail Area = {} ft^2".format(S_vt_actual))
print("Estimated Horizontal Tail Area = {} ft^2".format(S_ht))
print("Actual Horizontal Tail Area = {} ft^2".format(S_ht_actual))

AR_w  = 3.36462# Aspect ratio of wing
lambda_w = math.radians(31.6) # Sweep angle of wing (radians)

AR_h  = 0  # Aspect ratio of horizontal stabilizer
lambda_h = math.radians(35)   # Sweep angle of horizontal stabilizer (radians)

eta_w = 0.97   # Difference factor between the theoretical section lift curve slope for the wing
# eta_h =    # Difference factor between the theoretical section lift curve slope for the horizontal tail

M     = 0.85 # Mach number

CL_a_w  = (2*np.pi*AR_w)/(((2)+(np.sqrt((((AR_w/eta_w)**2)*(1+(np.tan(lambda_w))**2-M**2))+(4))))) # Lift curve slope of wing, / radian
# CL_a_h0 = (2*np.pi*AR_h)/(((2)+(np.sqrt((((AR_h/eta_h)**2)*(1+(np.tan(lambda_h))**2-M**2))+(4))))) # Lift curve slope of horizantle tail, / radian

print("Lift curve slope of wing = {} / radian".format(CL_a_w))
# print("Lift curve slope of horizontal tail = {} / radian".format(CL_a_h0))

de_dalpha = 2*CL_a_w / (np.pi * AR_w)
print("Downwash: %.3f / radian" %de_dalpha)
# CL_a_h = CL_a_h0 / (1 - de_dalpha)
# print("Lift curve slope of horizontal tail corrected for downwash = {} / radian".format(CL_a_h))
Kf  = 0.344          # Empirical factor (Assumed)
Lf  = 40.095 # Fuselage length (m) --> (ft)
Wf  = 8.19   # Maximum width of fuselage (m) --> (ft)

dCmf_dCL = (Kf * (Wf ** 2) * Lf) / (S_w * c_w * CL_a_w)
print("dCmf_dCL = {}".format(dCmf_dCL))
x_cg = 24.78           # Aircraft center of gravity (ft) assumed
x_25MAC = 7.057         # Distance from nose to 25% MAC (ft) assumed

SM = (x_cg-x_25MAC) / (c_w) - 0 / (CL_a_w * S_w * c_w) + dCmf_dCL
print("SM = {}".format(-SM))