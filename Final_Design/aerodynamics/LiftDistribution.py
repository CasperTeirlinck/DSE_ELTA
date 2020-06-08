
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

# Verified
def CalculateChord(theta,taper,wingsurface,wingspan):
    y = TransformTheta(theta,wingspan)
    return 2*wingsurface/(wingspan+taper*wingspan)*(1-(1/wingspan-taper/wingspan)*abs(2*y))

# Verified
def CalculateTwistAngle(theta,wingspan,twist,gamma):
    y = TransformTheta(theta,wingspan)
    b_half = .5*wingspan
    if gamma != 0:
        C1 = twist/(1-np.exp(gamma*b_half))
        C2 = C1*np.exp(gamma*b_half)
        return C1 - C2*np.exp(-gamma*abs(y))
    else:
        
        return twist*(1-abs(y)/wingspan)

# Verified
def CalculateZeroLiftAngle(theta,wingspan,twist,zeroliftangle_root,zeroliftangle_tip,gamma):
    alpha_geometric = CalculateTwistAngle(theta,wingspan,twist,gamma)
    y = TransformTheta(theta,wingspan)
    return 2./wingspan*(zeroliftangle_tip - zeroliftangle_root)*abs(y) + zeroliftangle_root - alpha_geometric

# Verified
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


def calcCLmax(alphaRange, ClmaxDistr, Acoeff, b=v.b, taper=v.taper, S=v.S, V=v.V, rho=v.rho, AR=v.A):

    alphaRange = np.arange(*alphaRange, 0.1)

    alphaMax = None
    Cl_distrMax = None
    yPntsMax = None
    CLmax = None

    for alpha in alphaRange:
        if alphaMax: break

        alpha = np.radians(alpha)
        Cl_distr, yPnts, CL = calcLiftDistribution(alpha, Acoeff, b, taper, S, V, rho, AR)

        for Cl, y in zip(Cl_distr, yPnts):
            Clmax = ClmaxDistr((y/yPnts.size - 0.5)*np.max(yPnts)*2)

            if np.abs(Cl - ClmaxDistr(y)) <= 0.01:
                alphaMax = alpha
                Cl_distrMax = Cl_distr
                yPntsMax = yPnts
                CLmax = CL
                break

    if not alphaMax: return None
    return CLmax, alphaMax, yPntsMax, Cl_distrMax

def calcCD0wing(thicknesstochord,frictioncoefficient=0.0055):
    formfactor = 1 + 2.7*thicknesstochord + 100*thicknesstochord**4
    return 2*frictioncoefficient*formfactor

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

    # Test convergence
    if False:
        CLlist = []
        alpha = np.radians(13)
        Nmax = 100
        for N in range(1, Nmax):
            Acoeff = CalculateDistributionCoefficients(v.S, v.b, v.taper, v.twist, v.Cla_r, v.Cla_t, v.a0_r, v.a0_t, 0, N)
            A1 = Acoeff[0][0]*alpha + Acoeff[0][1]
            CL = np.pi * v.A * A1
            print(f'N: {N}')
            CLlist.append(CL)
        plt.plot(range(1, Nmax), CLlist)
        plt.show()
    
    # wingsurface,wingspan,taper,twist,liftslope_root,liftslope_tip,zeroliftangle_root,zeroliftangle_tip,gamma,N,FuselageIncluded=False
    M = CalculateDistributionCoefficients(v.S, v.b, v.taper, v.twist, v.Cla_r, v.Cla_t, v.a0_r, v.a0_t, 0, 1000)
    
    Acoeff = M

    # CLrange:
    if True:
        Cl_distr_range = []
        alphaRange = [5]
        # alphaRange = np.arange(0, 15, 1)
        for alpha in alphaRange:
            alpha = np.radians(alpha)
            Cl_distr, yPnts, CL = calcLiftDistribution(alpha, Acoeff, b=v.b, taper=v.taper, S=v.S, V=v.V, rho=v.rho, AR=v.A)
            Cl_distr_range.append(Cl_distr)
        
        c_tip = CalculateChord(0, v.taper, v.S, v.b)
        c_root = CalculateChord(np.pi/2, v.taper, v.S, v.b)
        print(f'CL = {CL} @ a = {alphaRange[-1]}    c_root={c_root}, c_tip={c_tip}, b={v.b}')
        plotLiftDistribution(yPnts, Cl_distr_range)

    # CLa:
    if True:
        CLa = np.pi*v.A*Acoeff[0][0]
        print(f'CLa = {CLa}')

    # CLmax:
    if False:
        ClmaxDistr = lambda y: (v.Clmax_t - v.Clmax_r)/(v.b/2) * abs(y) + v.Clmax_r
        CLmax, alphaMax, yPnts, Cl_distr = calcCLmax([0, 20], ClmaxDistr, Acoeff)
        print(f'CLmax = {round(CLmax, 2)} @ a = {round(np.degrees(alphaMax), 2)}')
        plotLiftDistribution(yPnts, [Cl_distr], ClmaxDistr)
