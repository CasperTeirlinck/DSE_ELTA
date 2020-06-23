'''
This script contains functions to calculate wing properties
Author: Bob
'''

from math import tan,atan,cos,pi,sqrt

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
XMAC() :                Calculates the distance wing leading edge root chord - leading edge of the MAC

Inputs:
    b [float]:          Span [m]
    sweepc4 [float]:    Quarter chord sweep angle [rad]
    taper [float]:      Taper ratio [-]
    cr [float]:         Root chord [m]

Outputs:
    xmac [float]:       Distance LE root chord - LE of the MAC [m]

V&V:    Verified
'''

def XMAC(b,sweepc4,taper,cr):
    sweepLE = sweep(0,b,sweepc4,taper,cr)
    ymac = (b/6)*((1+2*taper)/(1+taper))
    return ymac*tan(sweepLE)


'''
XLEMAC() :              Calculates the distance nose - leading edge of the MAC

Inputs:
    lfn [float]:        Distance nose - wing leading edge root chord [m]
    sweepc4 [float]:    Quarter chord sweep angle [rad]
    cr [float]:         Root chord [m]
    b [float]:          Span [m]
    taper [float]:      Taper ratio [-]

Outputs:
    xlemac [float]:     Distance nose - leading edge of the MAC [m]

V&V:    Verified
'''

def XLEMAC(lfn,b,sweepc4,taper,cr):
    xmac = XMAC(b,sweepc4,taper,cr)
    return lfn + xmac


'''
downwash() :            Calculates the downwash gradient

Inputs:
    lh [float]:         Tail arm [m]
    b [float]:          Wing span [m]
    zh [float]:         Horizontal tail height [m]
    zw [float]:         Wing height [m]
    twist [float]:      Wing twist [rad]
    sweep [float]:      Wing quarter chord sweep angel [rad]
    CLa [float]:        Wing lift rate coefficient [-]
    A [float]:          Wing aspect ratio [-]

Outputs:
    deda [float]:       Downwash gradient [-]
'''

def downwash(lh,b,zh,zw,twist,sweep,CLa,A):
    r = lh*2/b
    mtv = 2/b * ((zh-zw)+lh*tan(twist)) * cos(twist)
    KeLambda = (0.1124 + 0.1265*sweep + sweep**2)/(r**2) + 0.1025/r + 2
    KeLambda0 = 0.1124/(r**2) + 0.1024/r + 2
    deda = KeLambda/KeLambda0 * (r/(r**2 + mtv**2)*0.4876/sqrt(r**2+0.6319+mtv**2)+
                                 (1+(r**2/(r**2+0.7915+5.0734*mtv**2))**0.3113)
                                 *(1-sqrt(mtv**2/(1+mtv**2)))) * CLa/(pi*A)
    return deda


'''
percMAC() :             Calculates the postion w.r.t. the MAC

Inputs:
    xcg [float]:        Centre of gravity [m]
    xlemac [float]:     Distanc nose - LE MAC [m]
    MAC [float]:        Mean Aerodynamic Chord [m]

Outputs:
    xcg [float]:        Centre of gravity [%MAC]
'''

def percMAC(xcg,xlemac,MAC):
    return (xcg-xlemac)/MAC