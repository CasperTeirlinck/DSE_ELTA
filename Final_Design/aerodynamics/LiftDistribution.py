
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
        C1 = twist/(1-np.exp(gamma*wingspan))
        C2 = C1*np.exp(gamma*wingspan)
        return C1 - C2*np.exp(-gamma*y)
    else:
        
        return twist*(1-abs(y)/wingspan)

# NOT verified
def CalculateZeroLiftAngle(theta,wingspan,twist,zeroliftangle_root,zeroliftangle_tip,gamma):
    alpha_geometric = CalculateTwistAngle(theta,wingspan,twist,gamma)
    y = TransformTheta(theta,wingspan)
    #return 2./wingspan(zeroliftangle_tip - zeroliftangle_root)*abs(y) - zeroliftangle_root - alpha_geometric
    return 0

# NOT Verified
def CalculateLiftSlope(theta,wingspan,liftslope_root,liftslope_tip):
    y = TransformTheta(theta,wingspan)
    return 2./wingspan*(liftslope_tip - liftslope_root)*abs(y) + liftslope_root
    

def CalculateFuselageContribution():
    return 0


# Verified without lift slope & twist implementation
def CalculateDistributionCoefficients(wingsurface,wingspan,taper,twist,liftslope_root,liftslope_tip,zeroliftangle_root,zeroliftangle_tip,gamma,N,FuselageIncluded=False):
    matrix = np.ndarray((N,N)) # Create sample matrix
    column2 = np.zeros((N,1)) # Create column for twist and fuselage contributions

    #samplepoints = np.arange(0,wingspan/2,wingspan/(2*N)) # Create uniform sampling distribution
    #samplepoints = TransformSpan(samplepoints,wingspan)   # Convert sample points to theta-coordinates
    samplepoints = np.linspace((wingspan/2-np.pi/2)/N,np.pi/2,N)
    # print(len(samplepoints)-N)
    for i in range(N):
        theta_sample = samplepoints[i] # Use sample point i
        a0 = CalculateLiftSlope(theta_sample,wingspan,liftslope_root,liftslope_tip) # Calculate the lift slope of sample point i
        c = CalculateChord(theta_sample,taper,wingsurface,wingspan)            # Calculate the chord of sample point i
        zeroliftangle = CalculateZeroLiftAngle(theta_sample,wingspan,twist,zeroliftangle_root,zeroliftangle_tip,gamma) # Calculate the zero lift angle
        fuselageangle = FuselageIncluded*CalculateFuselageContribution() # Calculate the fuselage contribution to the angle of attack.
        
        column2[i] = - zeroliftangle - fuselageangle # Calculate element (i,1) of the coefficient matrix

        for j in range(N): 
            element = np.sin(theta_sample*(2*j+1))*(4*wingspan/(a0*c) + (2*j+1)/np.sin(theta_sample)) # Calculate element (i,j) of the matrix
            
            # Add element to matrix; if element is close to zero, add zero.
            if abs(element) > 1e-4:
                matrix[i,j] = element

            else:
                matrix[i,j] = 0.
    
    matrix_inverse = la.inv(matrix) # Calculate inverse of the matrix
    # print(la.det(matrix))
    # Calculate the columns of the coefficient matrix
    column1 = np.matmul(matrix_inverse, np.ones((N,1)))
    column2 = np.matmul(matrix_inverse, column2)

    coefficientmatrix = np.concatenate((column1,column2),axis=1) # Merge columns into one matrix
    
    return coefficientmatrix

def calcLiftDistribution(alpha, Acoeff, b, taper, S, V, rho, AR):
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

def calcCLmax(alphaRange, ClmaxDistr):

    diff = []
    alphaRange = np.arange(*alphaRange, 1)

    for alpha in alphaRange:
        alpha = np.radians(alpha)
        Cl_distr, yPnts, CL = calcLiftDistribution(alpha, Acoeff, b=v.b, taper=1, S=v.S, V=v.V, rho=v.rho, AR=v.A)

        Clmax_idx = np.argmax(Cl_distr)
        diff.append( np.abs( Cl_distr[Clmax_idx] - ClmaxDistr((Clmax_idx/yPnts.size - 0.5)*np.max(yPnts)*2) ) )

    alphaMax = alphaRange[np.argmin(diff)]
    Cl_distr, yPnts, CL = calcLiftDistribution(np.radians(alphaMax), Acoeff, b=v.b, taper=1, S=v.S, V=v.V, rho=v.rho, AR=v.A)
    
    return CL, alphaMax

def plotLiftDistribution(y, Cl_range, ClmaxDistr=None):
    fig = plt.figure(figsize=(10, 4.5))
    ax1 = fig.add_subplot(111)

    for distr in Cl_range:
        ax1.plot(y, distr, linewidth=2, color='green', marker='o', fillstyle='none', markevery=5)
        if ClmaxDistr: ax1.plot(y, ClmaxDistr(y), linewidth=1, color='red')

    ax1.axvline(x=0, linewidth=2, color='black')
    ax1.axhline(y=0, linewidth=2, color='black')
    ax1.set_xlabel('Wingspan [m]')
    ax1.set_ylabel('')
    ax1.xaxis.grid(color='black', linestyle='--')
    ax1.yaxis.grid(color='black', linestyle='--')
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.show()

if __name__ == "__main__":
    M = CalculateDistributionCoefficients(v.S, v.b, 0.4, 0, 1.85*np.pi, 2.1*np.pi, 0, 0, 0, 20)
    
    # print(M)
    
    Acoeff = np.array([
        [0.2316, 0], # A_n, aL=0 _n
        [0.0277, 0], 
        [0.0040, 0]
    ])

    Acoeff = M

    # CLrange:
    Cl_distr_range = []
    for alpha in range(13,14):
        alpha = np.radians(alpha)
        Cl_distr, yPnts, CL = calcLiftDistribution(alpha, Acoeff, b=v.b, taper=1, S=v.S, V=v.V, rho=v.rho, AR=v.A)
        Cl_distr_range.append(Cl_distr)
    
    print(f'CL = {CL}')
    plotLiftDistribution(yPnts, Cl_distr_range)

    # CLmax:
    ClmaxDistr = lambda y: (v.Clmax_t - v.Clmax_r)/(v.b/2) * abs(y) + v.Clmax_r
    CLmax, alphaMax = calcCLmax([0, 20], ClmaxDistr)
    print(f'CLmax = {round(CLmax, 2)} @ a = {alphaMax}')
