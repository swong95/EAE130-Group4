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
cl_takeoff = np.linspace(-2.0, 2.0, 100)
cl_landing = np.linspace(-2.6, 2.6, 100)

# using drag polar eqn CD = (CD0 + dCD0) + (CL^2 / (pi * AR * e))

cd_clean = CD0 + (cl_clean**2 / (np.pi * AR * e_clean))

cd_takeoff = (CD0 + dCD0_to) + (cl_takeoff**2 / (np.pi * AR * e_to))

cd_landing = (CD0 + dCD0_lf) + (cl_landing**2 / (np.pi * AR * e_lf))

cd_gear = (CD0 + dCD0_gear) + (cl_landing**2 / (np.pi * AR * e_lg))

# plots
plt.figure(figsize=(10, 6))
plt.title('Fighter Aircraft Drag Polars', fontsize=14)
plt.xlabel("$C_D$", fontsize=12)
plt.ylabel("$C_L$", fontsize=12)

plt.plot(cd_clean,   cl_clean,   label=f'Clean', linewidth=2)
plt.plot(cd_takeoff, cl_takeoff, label=f'Takeoff Flaps', linewidth=2, linestyle='dotted')
plt.plot(cd_landing, cl_landing, label=f'Landing Flaps', linewidth=2)
plt.plot(cd_gear,    cl_landing, label=f'Landing Gear', linewidth=2)

plt.grid(True, linestyle=':', alpha=0.6)
plt.legend(loc='best')
plt.show()