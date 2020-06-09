'''
This script contains functions to calculate wing properties
Author: Bob
'''

from math import tan,atan

'''
sweep() :               Calculates the sweep angle at c'/c

Inputs:
    pc [float]:         Fraction of the chord (between 0 and 1) [-]
    b [float]:          Wing span [m]
    sweepc4 [float]:    Quarter chord sweep angle [rad]
    taper [float]:      Taper ratio [-]
    cr [float]:         Root chord [m]

Outputs:
    Sweep angle at c'/c [rad]

V&V:    Verified
'''

def sweep(pc,b,sweepc4,taper,cr):
    sweepLE = atan(tan(sweepc4) - cr/(2*b)*(taper-1))
    return atan(tan(sweepLE) + 2*pc*cr*(taper-1)/b)


'''
chord() :               Calculates the chord at location y

Inputs:
    y [float]:          y location [m]
    b [float]:          Wing span [m]
    sweepc4 [float]:    Quarter chord sweep angle [rad]
    taper [float]:      Taper ratio [-]
    cr [float]:         Root chord [m]

Outputs:
    Chord at location y [rad]

V&V:    Verified
'''

def chord(y,b,sweepc4,taper,cr):
    sweepLE = sweep(0,b,sweepc4,taper,cr)
    sweepTE = sweep(1,b,sweepc4,taper,cr)
    return cr-y*tan(sweepLE) + y*tan(sweepTE)


'''
XLEMAC() :              Calculates the distance nose - leading edge of the MAC

Inputs:
    lfn [float]:        Distance nose - wing leading edge root chord [m]
    sweepc4 [float]:    Quarter chord sweep angle [rad]
    cr [float]:         Root chord [m]
    b [float]:          Span [m]
    taper [float]:      Taper ratio [-]

Outputs:
    xmac [float]:       Distance LE root chord - LE of the MAC [m]
    xlemac [float]:     Distance nose - leading edge of the MAC [m]

V&V:    Verified
'''

def XLEMAC(lfn,b,sweepc4,taper,cr):
    sweepLE = sweep(0,b,sweepc4,taper,cr)
    ymac = (b/6)*((1+2*taper)/(1+taper))
    xmac = ymac*tan(sweepLE)
    xlemac = lfn + xmac

    return xlemac,xmac