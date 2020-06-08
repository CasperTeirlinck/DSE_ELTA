'''
This script performs the sizing of the control surfaces.
Author: Bob
'''

import numpy as np
from math import pi

'''
aileron_sizing()        Sizing of the ailerons

Inputs:
    variables[class]:   The class should contain:
                         - b1 [float]:          Aileron start
                         - b2 [float]:          Aileron end

Outputs:
    variables [class]:  The class contains updated values of:
                         - b1 [float]:          Aileron start
                         - b2 [float]:          Aileron end
'''

def aileron_sizing(variables):
    # Iteration parameters
    sizing = True
    i = 0                       # [-]       Number of iterations
    max_i = 1000                # [-]       Maximum number of iterations
    step = 0.01                 # [m]       Increase/Decrease of aileron length at every iteration

    # Inputs
    dphi = 60*pi/180            # [rad]     Bank angle
    dt = 5                      # [s]       Bank time
    VTO = 1.05*25.2             # [m/s]     Take-off speed
    VL = 1.1*25.2               # [m/s]     Landing speed

    bf = 1.5                    # [m]       Fuselage width

    S = 23.0                    # [m2]      Wing surface area
    b = 16.6                    # [m]       Wing span
    taper = 0.4                 # [-]       Wing taper ratio
    cr = 1.98                   # [m]       Wing root chord
    cla_TO = 2*pi               # [/rad]    Take-off configuration lift curve slope
    cla_L = 2*pi                # [/rad]    Landing configuration lift curve slope
    cd0_TO = 0.03               # [-]       Take-off configuration zero lift drag coefficient
    cd0_L = 0.03                # [-]       Landing configuration zero lift drag coefficient

    b1 = variables.b1           # [m]       Aileron start
    b2 = variables.b2           # [m]       Aileron end
    clear_tip = 0.5             # [m]       Distance from the tip that should be clear of control surfaces
    da = 20*pi/180              # [rad]     Aileron deflection angle
    clda_TO = 0.046825*180/pi   # [/rad]    Take-off configuration change in the airfoil’s lift coefficient with aileron deflection
    clda_L = 0.046825*180/pi    # [/rad]    Landing configuration Change in the airfoil’s lift coefficient with aileron deflection

    sm = 0.1                    # [-]       Safety margin

    # Parameter calculations
    # Required roll rate
    p_req = dphi/dt * (1+sm)

    # Check initial b1
    if b1<bf/2:
        print('Initial value for b1 is too small! Change b1.')
    # Check initial b2
    elif b2>(b/2-clear_tip):
        print('Initial value for b2 is too large! Change either b2 or the clearance from the wing tip.')
    else:
        pass

    # Create history lists
    b1lst = [b1]
    b2lst = [b2]

    # Roll rate calculation
    def roll_rate(V,cla,cd0,clda):
        # Calculate roll damping
        Clp = -(cla + cd0)*cr*b/(24*S) * (1 + 3*taper)

        # Calculate roll authority
        Clda = clda*cr/(S*b) * ((b2**2-b1**2)+4*(taper-1)/(3*b)*(b2**3-b1**3))

        # Calculate roll rate for take-off and landing
        p = -Clda/Clp*da * 2*V/b

        return p

    # Perform iterations
    while sizing and i<100:
        # Calculate roll rate
        p_TO = roll_rate(VTO,cla_TO,cd0_TO,clda_TO)
        p_L = roll_rate(VL,cla_L,cd0_L,clda_L)

        # Check whether both p are larger than required
        # If one of the p is smaller than required
        if p_TO<p_req or p_L<p_req:
            b2 += step/2
            # Check whether b2 can be shifted further
            if b2>(b/2-clear_tip):
                b2 -=step/2
                b1 -= step
            else:
                b1 -= step/2

        # If both ps are larger than required
        else:
            # Check for convergence
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

        # Check aileron size
        if b1<bf/2:
            print('Aileron became too large.')
            sizing = False
        # Check for maximum number of iterations
        elif i>=max_i:
            sizing = False
            print('Aileron sizing did not converge within',i,'iterations.')
        else:
            pass

    # Transform to numpy arrays
    #b1lst = np.array(b1lst)
    #b2lst = np.array(b2lst)

    # Update values in variables class
    variables.b1 = b1
    variables.b2 = b2

    return variables

# Test
class Test_variables_sc:
    def __init__(self):
        self.b1 = 6.36          # [m]       Aileron start
        self.b2 = 7.8           # [m]       Aileron end

if __name__ ==  "__main__":
    test_v = Test_variables_sc()

    test_v = aileron_sizing(test_v)