import numpy as np
import matplotlib.pyplot as plt

# Aircraft Parameters
S = 600  # Wing area (ft^2)
b = 35.39616  # Wing span (ft)
#AR = b**2 / S  # Aspect ratio
AR = 2.08815  # Aspect ratio
MAC = 21.00147  # Mean Aerodynamic Chord (ft)

# Weight
Weight_Empty = 26615  # Empty weight (lb)
Weight_MTOW_Strike = 52028  # Maximum Takeoff Weight for Strike (lb)
Weight_MTOW_A2A = 48891  # Maximum Takeoff Weight for A2A (lb)

# Thrust (F135-PW-400 Engine)
Thrust_Dry = 28000  # Dry thrust (lbf)
Thrust_Afterburner = 43000  # Afterburner thrust (lbf)


# Stall Condition (For Stall Speed, vary rho with altitude)
