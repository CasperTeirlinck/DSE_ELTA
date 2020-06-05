
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

def CalculateZeroLiftAngle(theta,wingspan,twist,zeroliftangle_1,zeroliftangle_2,gamma):
    return 0

# NOT Verified
def CalculateLiftSlope(theta,wingspan,a_airfoil1,a_airfoil2):
    return 2*np.pi

def CalculateFuselageContribution():
    return 0


# Verified without lift slope & twist implementation
def CalculateDistributionCoefficients(wingsurface,wingspan,taper,twist,liftslope_1,liftslope_2,zeroliftangle_1,zeroliftangle_2,gamma,N,FuselageIncluded=False):
    matrix = np.ndarray((N,N)) # Create sample matrix
    column2 = np.zeros((N,1)) # Create column for twist and fuselage contributions

    samplepoints = np.arange(0,wingspan/2,wingspan/(2*N)) # Create uniform sampling distribution
    samplepoints = TransformSpan(samplepoints,wingspan)   # Convert sample points to theta-coordinates
    
    for i in range(N):
        theta_sample = samplepoints[i] # Use sample point i
        a0 = CalculateLiftSlope(theta_sample,wingspan,liftslope_1,liftslope_2) # Calculate the lift slope of sample point i
        c = CalculateChord(theta_sample,taper,wingsurface,wingspan)            # Calculate the chord of sample point i
        zeroliftangle = CalculateZeroLiftAngle(theta_sample,wingspan,twist,zeroliftangle_1,zeroliftangle_2,gamma) # Calculate the zero lift angle
        fuselageangle = FuselageIncluded*CalculateFuselageContribution() # Calculate the fuselage contribution to the angle of attack.
        
        column2[i] = - zeroliftangle - fuselageangle # Calculate element (i,1) of the coefficient matrix

        for j in range(N): 
            element = np.sin(theta_sample*(2*j+1))*(4*wingspan/(a0*c) + (2*j+1)/np.sin(theta_sample)) # Calculate element (i,j) of the matrix
            
            # Add element to matrix; if element is close to zero, add zero.
            if abs(element) > 1e-6:
                matrix[i,j] = element

            else:
                matrix[i,j] = 0.
    
    # Calculate inverse of the matrix
    matrix_inverse = la.inv(matrix)


    # Calculate the columns of the coefficient matrix
    column1 = np.matmul(matrix_inverse, np.ones((N,1)))
    column2 = np.matmul(matrix_inverse, column2)

    print("COLUMN 1: \n",column1)
    print("COLUMN 2: \n",column2)

    # Merge columns into one matrix
    coefficientmatrix = np.concatenate((column1,column2),axis=1)

    return coefficientmatrix

def calculateLiftDistribution(alpha, Acoeff):
    A = [ A[0]*alpha + A[1] for A in Acoeff ]

    CL = np.pi * v.A * A[0]

    # def circ(theta):
    #     for n in range(1, )
    #     return 2 * v.b * V * np.sum(A * np.sin(n*theta))


if __name__ == "__main__":
    

    Acoeff = np.array([
        [0.2316, 0], # A_n, aL=0 _n
        [0.0277, 0], 
        [0.0040, 0]
    ])
    calculateLiftDistribution(np.radians(5), Acoeff)
