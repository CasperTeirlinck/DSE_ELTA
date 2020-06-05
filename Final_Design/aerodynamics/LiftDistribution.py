
import numpy as np

import variables_aero as v


def TransformSpan(y,wingspan):
    return -np.arccos(2*y/wingspan)


def TransformTheta(theta,wingspan):
    return -.5*wingspan*np.cos(theta)


def CalculateChord(theta,taper,S,span):
    y = TransformTheta(theta,wingspan)
    return 2*S/(wingspan+taper*wingspan)*(1-(1/wingspan-taper/wingspan)*abs(2*y))


def CalculateTwistAngle(theta,wingspan,twist,gamma):
    y = TransformTheta(theta,wingspan)
    if gamma != 0:
        C1 = twist/(1-exp(gamma*wingspan))
        C2 = C1*exp(awingspan)
        return C1 - C2*exp(-gamma*y)
    else:
        print("NOT IMPLEMENTED YET")
        return None


def CalculateLiftSlope(theta,wingspan,a_airfoil1,a_airfoil2):
    return 2*np.pi


def CalculateDistributionCoefficients(wingspan,sweep,liftslope_1,liftslope2,zeroliftangle_1,zeroliftangle_2,N,Fuselage=False):
    b = wingspan







alpha = np.radians(5)
V = 50 # [m/s]

Acoeff = np.array([
    [0.2316, 0], # A_n, aL=0 _n
    [0.0277, 0], 
    [0.0040, 0]
])
A = [ A[0]*alpha + A[1] for A in Acoeff ]

CL = np.pi * v.A * A[0]

def circ(theha):
    for n in range(1, )
    return 2 * v.b * V * np.sum(A * np.sin(n*theta))

