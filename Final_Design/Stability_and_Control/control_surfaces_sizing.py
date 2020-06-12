'''
This script performs the sizing of the control surfaces.
Author: Bob
'''

import numpy as np
from math import pi,cos,tan,atan,sqrt
from wing_properties import chord,sweep
import matplotlib.pyplot as plt

'''
aileron_sizing()        Sizing of the ailerons

Inputs:
    variables[class]:   The class should contain:
                         - b1 [float]:          Aileron start

Outputs:
    variables [class]:  The class contains updated values of:
                         - b1 [float]:          Aileron start
                         - b2 [float]:          Aileron end

V&V: Verified
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

    S = 21.09                   # [m2]      Wing surface area
    b = sqrt(S*10.1)            # [m]       Wing span
    taper = 0.4                 # [-]       Wing taper ratio
    cr = 2*S/((1+taper)*b)      # [m]       Wing root chord
    cla = 4.87                  # [/rad]    Take-off configuration lift curve slope
    cd0_TO = 0.02               # [-]       Take-off configuration zero lift drag coefficient
    cd0_L = 0.02                # [-]       Landing configuration zero lift drag coefficient

    b1 = variables.b1           # [m]       Aileron start
    b2 = variables.b2           # [m]       Aileron end
    clear_tip = 0.05*(b/2)      # [m]       Distance from the tip that should be clear of control surfaces
    da = 20*pi/180              # [rad]     Aileron deflection angle
    clda_TO = 0.046825*180/pi   # [/rad]    Take-off configuration change in the airfoil’s lift coefficient with aileron deflection
    clda_L = 0.046825*180/pi    # [/rad]    Landing configuration Change in the airfoil’s lift coefficient with aileron deflection

    sm = 0.1                    # [-]       Safety margin

    # Parameter calculations
    # Required roll rate
    p_req = dphi/dt * (1+sm)

    # Aileron end
    b2 = b/2 - clear_tip

    # Check initial b1
    if b1<bf/2:
        print('Initial value for b1 is too small! Change b1.')
    else:
        pass

    # Create history list
    b1lst = [b1]

    # Roll rate calculation
    def roll_rate(V,cd0,clda):
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
        p_TO = roll_rate(VTO,cd0_TO,clda_TO)
        p_L = roll_rate(VL,cd0_L,clda_L)

        # Most critical roll rate
        p_crit = max(p_TO,p_L)

        # Check whether p is larger than required
        # If p is smaller than required
        if p_crit<p_req:
            b1 -= step

        # If p is larger than required
        else:
            # Check for convergence
            if abs(b1-b1lst[-1])<(step/2):
                sizing = False
            else:
                b1 += step

        # Add values to history lists
        b1lst.append(b1)

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


'''
flap_sizing()               Sizing of the ailerons

Inputs:
    variables [class]:      The class should contain:
                             (- b1 [float]:          Aileron start)
                             - f1 [float]:          Flap start
                             - f2 [float]:          Flap end
    fix_position [string]:  Choose 'fuselage end' or 'aileron start'

Outputs:
    variables [class]:  The class contains updated values of:
                         - f1 [float]:          Flap start
                         - f2 [float]:          Flap end

V&V:    Verified
'''

def flap_sizing(variables,fix_position='fuselage end'):
    # Check input
    fix_positionlst = ['fuselage end','aileron start']
    if not fix_position in fix_positionlst:
        print("Wrong fix_position input ("+fix_position+"). Choose 'fuselage end' or 'aileron start'")
    else: pass

    # Inputs
    S = 21.09                   # [m2]      Wing surface area
    b = sqrt(S*10.1)            # [m]       Wing span
    sweepc4 = 0                 # [rad]     Wing quarter chord sweep angle
    taper = 0.4                 # [rad]     Wing taper ratio
    cr = 2*S/((1+taper)*b)      # [m]       Wing root chord

    CLmax_req = 1.8             # [-]       Required maximum lift coefficient
    CLmax_wing = 1.4            # [-]       Wing maximum lift coefficient
    CLa = 2*pi                  # [/rad]    Wing lift curve slope

    dClmax = 1.13               # [-]
    da0l_airfoil = -15*pi/180   # [rad]

    cfc = 0.8                   # [-]       Start of the flap as percentage of the chord

    sm = 0.1                    # [-]       Safety margin

    # Parameter calculations
    # Flap star/end location
    # Chord at flap start/end location
    if fix_position == 'fuselage end':
        bf = 1.5                # [m]       Fuselage width
        d_ff = 0.05             # [m]       Spacing between fuselage and flap
        f1 = bf/2 + d_ff
        cf1 = chord(f1, b, sweepc4, taper, cr)
    else:
        b1 = variables.b1       # [m]       Aileron start
        d_af = 0.05             # [m]       Spacing between flap and aileron
        f2 = b1 - d_af
        cf2 = chord(f2, b, sweepc4, taper, cr)
    # Leading edge sweep angle
    sweepLE = sweep(0,b,sweepc4,taper,cr)

    # Hinge line sweep angle
    sweep_hinge = sweep(cfc,b,sweepc4,taper,cr)

    # Trailing edge sweep angles
    sweepTE = sweep(1,b,sweepc4,taper,cr)

    # Increase in lift coefficient
    dCLmax = CLmax_req - CLmax_wing

    # Required flapped surface
    SwfS = dCLmax/(0.9*dClmax*cos(sweep_hinge)) * (1+sm)
    Swf = SwfS*S

    # Shift in zero lift angle of attack
    da0L = da0l_airfoil*SwfS*cos(sweep_hinge)

    # Change in lift curve slope
    CLa_flapped = CLa

    # Flap span calculation
    # Solving the equation:
    # 'fuselage end': cf1*bfl - 0.5*bf^2*tan(sweepLE) + 0.5*bf^2*tan(sweepTE) = Swf/2
    # 'aileron start': cf2*bfl + 0.5*bf^2*tan(sweepLE) - 0.5*bf^2*tan(sweepTE) = Swf/2
    # a*bfl^2 + b*bfl + c = 0

    # Fuselage end
    if fix_position == 'fuselage end':
        A = 0.5*(-tan(sweepLE)+tan(sweepTE))
        B = cf1
        C = -Swf/2

        D = B**2 - 4*A*C

        if not D >= 0:
            print('There is a problem with the flap sizing! (control_surfaces_sizing.py)')
            return
        else:
            bfllst = [0,0]
            bfllst[0] = (-B+sqrt(D))/(2*A)
            bfllst[1] = (-B-sqrt(D))/(2*A)
            if bfllst[0]>0 and (f1+bfllst[0])<(b/2):
                bfl = bfllst[0]
            elif bfllst[1]>0 and (f1+bfllst[1])<(b/2):
                bfl = bfllst[1]
            else:
                print("Flap is too large, it doesn't fit on the wing!")

    # Aileron start
    else:
        A = 0.5*(tan(sweepLE)-tan(sweepTE))
        B = cf2
        C = -Swf / 2

        D = B**2 - 4*A*C

        if not D >= 0:
            print('There is a problem with the flap sizing! (control_surfaces_sizing.py)')
            return
        else:
            bfllst = [0,0]
            bfllst[0] = (-B+sqrt(D))/(2*A)
            bfllst[1] = (-B-sqrt(D))/(2*A)
            if bfllst[0]>0 and (f2-bfllst[0])>0:
                bfl = bfllst[0]
            elif bfllst[1]>0 and (f2-bfllst[1])>0:
                bfl = bfllst[1]
            else:
                print("Flap is too large, it doesn't fit on the wing!")

    # Flap start/end location
    if fix_position == 'fuselage end':
        f2 = f1 + bfl
    else:
        f1 = f2 - bfl

    # Update variables class
    variables.Swf = Swf
    variables.f1 = f1
    variables.f2 = f2

    return variables


# Test
class Test_variables_sc:
    def __init__(self):
        self.b1 = 6.36          # [m]       Aileron start
        self.b2 = None          # [m]       Aileron end
        self.Swf = None         # [m2]      Flapped area
        self.f1 = None          # [m]       Flap start
        self.f2 = None          # [m]       Flap end

if __name__ ==  "__main__":
    test_v = Test_variables_sc()

    test_v = aileron_sizing(test_v)
    test_v = flap_sizing(test_v,fix_position='fuselage end')
    print(test_v.b1)
    print(test_v.b2)
    print(test_v.f1)
    print(test_v.f2)
