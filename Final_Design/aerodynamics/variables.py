import numpy as np

A = 10.1 # [-] Aspect ratio
S = 15.6 # [m2] Wing surface

taper = 0.9 # [-] Taper ratio
twist = np.radians(2)
gamma = 0

w_fuselage = 1.2
CD0 = 0.0055*3.7

hwl = 0.4
kwl = 2.1

# NACA 4415
Clmax_r = 1.6
Cla_r = 1/np.radians(10)
Cd0_r = 0.007
a0_r = np.radians(-4)
deltaAlphaStall_r = np.radians(8)

Clmax_t = Clmax_r
Cla_t = Cla_r
Cd0_t = Cd0_r
a0_t = a0_r
deltaAlphaStall_t = deltaAlphaStall_r


# NACA 651412
# Clmax_t = 1.5
# Cla_t = 6.016
# Cd0_t = 0.0055
# a0_t = np.radians(-3)
# deltaAlphaStall_t = np.radians(2.5)
