
import numpy as np
from numpy import linalg as la
import variables_aero as v
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

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
    
    matrix_inverse = la.inv(matrix) # Calculate inverse of the matrix

    # Calculate the columns of the coefficient matrix
    column1 = np.matmul(matrix_inverse, np.ones((N,1)))
    column2 = np.matmul(matrix_inverse, column2)

    coefficientmatrix = np.concatenate((column1,column2),axis=1) # Merge columns into one matrix

    return coefficientmatrix

def calcLiftDistribution(alphaRange, Acoeff, b, taper, S, V, rho, AR):
    A = np.array([ A[0]*alpha + A[1] for A in Acoeff ])
    
    CL = np.pi * AR * A[0]

    def circ(theta):
        nsum = 0
        for n in range(A.size):
            i = n
            n = 2*n+1
            nsum += A[i]*np.sin(n*theta)
        return 2 * b * V * nsum

    Cl_distr = []
    yPnts = np.linspace(-b/2, b/2, 100)
    for y in yPnts:
        theta = -TransformSpan(y, b)
        Cl = (2 * circ(theta))/(V * CalculateChord(theta, taper, S, b))
        Cl_distr.append(Cl)

    return Cl_distr, yPnts, CL

def calcCLmax(alphaRange, yPnts, Cl_distr, ClmaxDistr):
    ClmaxDistr = ClmaxDistr[yPnts]

    for alpha in range(-15, 15):
        alpha = np.radians(alpha)
        Cl_distr, y_plot, CL = calcLiftDistribution(alpha, Acoeff, b=v.b, taper=1, S=v.S, V=v.V, rho=v.rho, AR=v.A)

def plotLiftDistribution(y, Cl_range, ClmaxDistr):
    fig = plt.figure(figsize=(10, 4.5))
    ax1 = fig.add_subplot(111)

    for distr in Cl_range:
        ax1.plot(y, distr, linewidth=2, color='green', marker='o', fillstyle='none', markevery=5)
        ax1.plot(y, ClmaxDistr(y), linewidth=1, color='red')

    ax1.axvline(x=0, linewidth=2, color='black')
    ax1.axhline(y=0, linewidth=2, color='black')
    ax1.set_xlabel('Wingspan [m]')
    ax1.set_ylabel('')
    ax1.xaxis.grid(color='black', linestyle='--')
    ax1.yaxis.grid(color='black', linestyle='--')
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.show()

if __name__ == "__main__":
    # M = CalculateDistributionCoefficients(16,0.4,24,6,6,0,0,30)
    # M_inverse = la.inv(M)
    
    Acoeff = np.array([
        [0.2316, 0], # A_n, aL=0 _n
        [0.0277, 0], 
        [0.0040, 0]
    ])

    Cl_distr_range = []
    for alpha in range(5, 6):
        alpha = np.radians(alpha)
        Cl_distr, y_plot, CL = calcLiftDistribution(alpha, Acoeff, b=v.b, taper=1, S=v.S, V=v.V, rho=v.rho, AR=v.A)
        Cl_distr_range.append(Cl_distr)
    
    ClmaxDistr = lambda y: (v.Clmax_t - v.Clmax_r)/(v.b/2) * abs(y) + v.Clmax_r

    print(f'CL = {CL}')
    plotLiftDistribution(y_plot, Cl_distr_range, ClmaxDistr)
