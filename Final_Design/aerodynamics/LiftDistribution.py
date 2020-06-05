
import numpy as np
from numpy import linalg as la
import variables_aero as v

# Verified
def TransformSpan(y,wingspan):
    return -np.arccos(2*y/wingspan)

# Verified
def TransformTheta(theta,wingspan):
    return -.5*wingspan*np.cos(theta)

# NOT verified
def CalculateChord(theta,taper,wingsurface,wingspan):
    y = TransformTheta(theta,wingspan)
    return 2*wingsurface/(wingspan+taper*wingspan)*(1-(1/wingspan-taper/wingspan)*abs(2*y))

# NOT verified
def CalculateTwistAngle(theta,wingspan,twist,gamma):
    y = TransformTheta(theta,wingspan)
    if gamma != 0:
        C1 = twist/(1-exp(gamma*wingspan))
        C2 = C1*exp(gamma*wingspan)
        return C1 - C2*exp(-gamma*y)
    else:
        print("NOT IMPLEMENTED YET")
        return None

# NOT Verified
def CalculateLiftSlope(theta,wingspan,a_airfoil1,a_airfoil2):
    return 2*np.pi

# Verified without lift slope implementation
def CalculateDistributionCoefficients(wingspan,taper,wingsurface,liftslope_1,liftslope_2,zeroliftangle_1,zeroliftangle_2,N,Fuselage=False):
    matrix = np.ndarray((N,N))
    samplepoints = np.arange(0,wingspan/2,wingspan/(2*N))
    samplepoints = TransformSpan(samplepoints,wingspan)
    
    for i in range(N):
        for j in range(N):
            theta_sample = samplepoints[i]
            
            a0 = CalculateLiftSlope(theta_sample,wingspan,liftslope_1,liftslope_2)
            c = CalculateChord(theta_sample,taper,wingsurface,wingspan)

            element = np.sin(theta_sample*(2*j+1))*(4*wingspan/(a0*c) + (2*j+1)/np.sin(theta_sample))
            if abs(element) > 1e-5:
                matrix[i,j] = element

            else:
                matrix[i,j] = 0.




def calculateLiftDistribution(alpha, Acoeff):
    A = [ A[0]*alpha + A[1] for A in Acoeff ]

    CL = np.pi * v.A * A[0]

    # def circ(theta):
    #     for n in range(1, )
    #     return 2 * v.b * V * np.sum(A * np.sin(n*theta))


if __name__ == "__main__":
    M = CalculateDistributionCoefficients(16,0.4,24,6,6,0,0,30)
    M_inverse = la.inv(M)
    
    Acoeff = np.array([
        [0.2316, 0], # A_n, aL=0 _n
        [0.0277, 0], 
        [0.0040, 0]
    ])
    calculateLiftDistribution(np.radians(5), Acoeff)