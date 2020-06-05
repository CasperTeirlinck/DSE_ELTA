'''
This script performs the sizing of the control surfaces.
Author: Bob
'''

import numpy as np
from math import pi

# Iteration parameters
sizing = True
i = 0                   # [-]       Number of iterations
step = 0.01             # [m]       Increase/Decrease of aileron length at every iteration

# Inputs
S = 23.0                # [m2]      Wing surface area
b = 16.6                # [m]       Wing span
taper = 0.4             # [-]       Wing taper ratio
cr = 1.98               # [m]       Wing root chord
cla = 2*pi              # [/rad]    Lift curve slope
cd0 = 0.03              # [-]       Zero lift drag coefficient

dphi = 60*pi/180        # [rad]     Bank angle
dt = 5                  # [s]       Bank time
VTO = 30                # [m/s]     Take-off speed
VL = 25.2*1.1           # [m/s]     Landing speed

b1 = 6.36               # [m]       Aileron start
b2 = 7.8                # [m]       Aileron end
clear_tip = 0.5         # [m]       Distance from the tip that should be clear of control surfaces
da = 20*pi/180          # [rad]     Aileron deflection angle
clda = 0.046825*180/pi  # [/rad]    Change in the airfoil’s lift coefficient with aileron deflection

# Parameter calculations
# Required roll rate
p_req = dphi / dt

# Check initial b2
if b2>(b/2-clear_tip):
    print('Initial value for b2 is too large! Change either b2 or the clearance from the wing tip.')
else:
    pass

# Create history lists
b1lst = [b1]
b2lst = [b2]

# Perform iterations
while sizing and i<100:
    # Calculate roll damping
    Clp = -(cla + cd0)*cr*b/(24*S) * (1 + 3*taper)

    # Calculate roll authority
    Clda = clda*cr/(S*b) * ((b2**2-b1**2)+4*(taper-1)/(3*b)*(b2**3-b1**3))

    # Calculate roll rate for take-off and landing
    p_TO = -Clda/Clp*da * 2*VTO/b
    p_L = -Clda/Clp*da * 2*VL/b

    # Check whether both p are larger than required
    # If one of the p is smaller than required
    if p_TO<p_req or p_L<p_req:
        b2 += step/2
        # Check whether b2 can become larger
        if b2>(b/2-clear_tip):
            b2 -=step/2
            b1 -= step
        else:
            b1 -= step/2

    # If both ps are larger than required
    else:
        if abs(b1-b1lst[-1])<(step/2) or abs(b2-b2lst[-1])<(step/2):
            sizing = False
        else:
            b2 -= step/2
            b1 += step/2

    # Add values to history lists
    b1lst.append(b1)
    b2lst.append(b2)

    # Add iteration
    i += 1

print('Done in',i,'iteration(s)')
