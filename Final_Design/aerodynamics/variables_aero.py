import numpy as np

A = 10 # [-] Aspect ratio
S = 15 # [m2] Wing surface
V = 50 # [m/s]

taper = 0.9 # [-] Taper ratio
twist = np.radians(2)
gamma = 0
# Clmax_r = 1.6 # Clmax root airfoil
# Clmax_t = 1.5 # Clmax tip airfoil

# Cla_r = 0.7/np.radians(7)
# Cla_t = 1/np.radians(10)
# Cd0_r = 0.006
# Cd0_t = 0.004
# a0_r = np.radians(-2)
# a0_t = np.radians(-4)


Clmax_t = 1.6 # Clmax root airfoil
Clmax_r = 1.5 # Clmax tip airfoil

Cla_t = 6.016
Cla_r = 1/np.radians(10)
Cd0_t = 0.004
Cd0_r = 0.004
a0_t = np.radians(-3)
a0_r = np.radians(-4)

"""
VALIDATION DATA
Clmax_t = 1.6 # Clmax root airfoil
Clmax_r = 1.5 # Clmax tip airfoil

Cla_t = 0.7/np.radians(7)
Cla_r = 1/np.radians(10)
Cd0_t = 0.006
Cd0_r = 0.004
a0_t = np.radians(-2)
a0_r = np.radians(-4)
"""