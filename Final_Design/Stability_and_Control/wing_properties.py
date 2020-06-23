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
    

'''

def downwash(lh,bw,zh,zw,twistwr,sweepw,CLaw,Aw):
    r = lh*2/bw
    mtv = 2/bw * ((zh-zw)+lh*tan(twistwr)) * cos(twistwr)
    KeLambda = (0.1124 + 0.1265*sweepw + sweepw**2)/(r**2) + 0.1025/r + 2
    KeLambda0 = 0.1124/(r**2) + 0.1024/r + 2
    deda = KeLambda/KeLambda0 * (r/(r**2 + mtv**2)*0.4876/sqrt(r**2+0.6319+mtv**2)+
                                 (1+(r**2/(r**2+0.7915+5.0734*mtv**2))**0.3113)
                                 *(1-sqrt(mtv**2/(1+mtv**2)))) * CLaw/(pi*Aw)
    return deda