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