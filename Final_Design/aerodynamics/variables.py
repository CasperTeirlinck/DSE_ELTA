import numpy as np

A = 10.1 # [-] Aspect ratio
# S = 15.6 # [m2] Wing surface
S = 14.738


taper = 0.9 # [-] Taper ratio
twist = np.radians(2)
gamma = 0

hwl = 0.3 # [m] winglet height
kwl = 2.1 # [-] winglet factor

S_wet_fus = 17.507
l_fus = 9.420
fus_A_max = 1.218
w_fuselage = 1.05

BLturbratio_fus = 1
BLturbratio_wing = 0.65
l_gear = 0.8
w_gear = 0.175
dCD_gear = 0.15
flap_area_ratio = 0.3

# Empennage
S_h = 3.83
S_v = 1.3
MAC_emp = 1.5
BLturbratio_emp = 0.65
tc_emp=0.12
xc_emp=0.3
V_stall=23.15
rho_cruise=1.04
visc=1.8e-5

# NACA 4415
Clmax_r = 1.6
Cla_r = 1/np.radians(10)
Cd0_r = 0.007
a0_r = np.radians(-4)
tc_airfoil = 0.15
xc_airfoil = 0.3

Clmax_t = Clmax_r
Cla_t = Cla_r
Cd0_t = Cd0_r
a0_t = a0_r








# NACA 651412
# Clmax_t = 1.5
# Cla_t = 6.016
# Cd0_t = 0.0055
# a0_t = np.radians(-3)
# deltaAlphaStall_t = np.radians(2.5)

# w_fuselage = 1.06
# # CD0 = 0.0055*3.7
# CD0 = 0.0050*3.7

# w_fuselage = 1.2
# S_wetted_fus = 19.182
# BLturb_ratio_fus = 1
# r_tail = 0.05

# CD0 = 0.

# MAC = 1
# CD_misc = 0.
