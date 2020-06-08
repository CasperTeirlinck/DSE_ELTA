
import numpy as np
from numpy import linalg as la
import variables_aero as v
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

class WingPlanform:
    def __init__(self, S, A, taper, twist, gamma):
        self.S = S
        self.A = A
        self.taper = taper
        self.twist = twist
        self.gamma = gamma
        self.b = np.sqrt(A*S)

        self.Clmax_r = None
        self.Clmax_t = None
        self.Cla_r = None
        self.Cla_t = None
        self.a0_r = None
        self.a0_t = None

        self.coeff = None

    def setAirfoils(self, Clmax_r, Clmax_t, Cla_r, Cla_t, a0_r, a0_t, Cd0_r, Cd0_t):
        self.Clmax_r = Clmax_r
        self.Clmax_t = Clmax_t
        self.Cd0_r = Cd0_r
        self.Cd0_t = Cd0_t
        self.Cla_r = Cla_r
        self.Cla_t = Cla_t
        self.a0_r = a0_r
        self.a0_t = a0_t

    def transformTheta(self, theta, b): # Verified
        return -.5*b*np.cos(theta)

    def transformSpan(self, y, b): # Verified
        return -np.arccos(2*y/b)

    def calculateChord(self, theta, taper, S, b): # Verified
        y = self.transformTheta(theta, b)
        return 2*S/(b + taper*b)*(1-(1/b - taper/b)*abs(2*y))

    def calcCoefficients(self, N, FuselageIncluded=False): # Verified without lift slope & twist implementation
        
        def _calcLiftSlope(theta, b, Cla_r, Cla_t): # Verified
            y = self.transformTheta(theta, b)
            return 2./self.b*(Cla_t - Cla_r)*abs(y) + Cla_r
        
        def _calculateTwistAngle(theta, b, twist, gamma): # Verified
            y = self.transformTheta(theta, b)
            b_half = .5*b
            if gamma != 0:
                C1 = twist/(1-np.exp(gamma*b_half))
                C2 = C1*np.exp(gamma*b_half)
                return C1 - C2*np.exp(-gamma*abs(y))
            else:
                return twist*(1-abs(y)/b)

        def _calculateZeroLiftAngle(theta, b, twist, a0_r, a0_t, gamma):
            alpha_geometric = _calculateTwistAngle(theta, b, twist, gamma)
            y = self.transformTheta(theta, b)
            return 2./b*(a0_t - a0_r)*abs(y) + a0_r - alpha_geometric

        def _calculateFuselageContribution():
            return 0

        matrix = np.ndarray((N,N)) # Create sample matrix
        column2 = np.zeros((N,1)) # Create column for twist and fuselage contributions

        samplepoints = np.linspace((self.b/2-np.pi/2)/N, np.pi/2, N)

        for i in range(N):
            theta_sample = samplepoints[i] # Use sample point i
            
            a0 = _calcLiftSlope(theta_sample, self.b, self.Cla_r, self.Cla_t) # Calculate the lift slope of sample point i
            c = self.calculateChord(theta_sample, self.taper, self.S, self.b) # Calculate the chord of sample point i
        
            zeroliftangle = _calculateZeroLiftAngle(theta_sample, self.b, self.twist, self.a0_r , self.a0_t, self.gamma) # Calculate the zero lift angle
            fuselageangle = FuselageIncluded*_calculateFuselageContribution() # Calculate the fuselage contribution to the angle of attack.
            column2[i] = - zeroliftangle - fuselageangle # Calculate element (i,1) of the coefficient matrix

            for j in range(N): 
                element = np.sin(theta_sample*(2*j+1))*(4*self.b/(a0*c) + (2*j+1)/np.sin(theta_sample)) # Calculate element (i,j) of the matrix
                
                # Add element to matrix; if element is close to zero, add zero.
                if abs(element) > 1e-4:
                    matrix[i,j] = element
                else:
                    matrix[i,j] = 0.
        
        matrix_inverse = la.inv(matrix) # Calculate inverse of the matrix
        column1 = np.matmul(matrix_inverse, np.ones((N,1)))
        column2 = np.matmul(matrix_inverse, column2)

        coefficientmatrix = np.concatenate((column1,column2),axis=1) # Merge columns into one matrix
        
        self.coeff = coefficientmatrix

        return coefficientmatrix

    def calcCL(self, alpha):

        A1 = self.coeff[0][0]*alpha + self.coeff[0][1]
        return np.pi * self.A * A1
        
    def calcLiftDistribution(self, alpha, N, showCircComponents=False):

        A = np.array([ A[0]*alpha + A[1] for A in self.coeff ])

        def _showCircComponents():
            testfig = plt.figure(figsize=(10, 4.5))
            testax = testfig.add_subplot(111)
            testy = np.linspace(-self.b/2, self.b/2, 1000)
            print(testy)
            for n in range(A.size):
                i = n
                n = 2*n+1

                term = lambda theta: 2 * self.b * A[i]*np.sin(n*theta)
                testax.plot(testy, term(self.transformSpan(testy, self.b)))

            plt.show()    
        if showCircComponents: _showCircComponents()
        
        def _circ(theta):
            nsum = 0
            for n in range(A.size):
                i = n
                n = 2*n+1
                nsum += A[i]*np.sin(n*theta)
            return 2 * self.b * nsum

        Cl_distr = []
        yPnts = np.linspace(-self.b/2, self.b/2, N)
        for y in yPnts:
            theta = -self.transformSpan(y, self.b)
            Cl = (2 * _circ(theta))/(self.calculateChord(theta, self.taper, self.S, self.b))
            Cl_distr.append(Cl)

        return Cl_distr, yPnts

    def calcCLa(self):
        CLa = np.pi*self.A*self.coeff[0][0]
        return CLa
    
    def calcCLmax(self):

        alphaRange = np.radians(np.arange(0, 20, 0.1))
        ClmaxDistr = lambda y: (self.Clmax_t - self.Clmax_r)/(self.b/2) * abs(y) + self.Clmax_r

        alphaMax = None

        for alpha in alphaRange:
            if alphaMax: break

            Cl_distr, yPnts = self.calcLiftDistribution(alpha, 100)

            for Cl, y in zip(Cl_distr, yPnts):
                if np.abs(Cl - ClmaxDistr(y)) <= 0.01:
                    alphaMax = alpha
                    Cl_distrMax = Cl_distr
                    yPntsMax = yPnts
                    CLmax = self.calcCL(alphaMax)
                    break

        if not alphaMax: return None
        return CLmax, alphaMax, Cl_distrMax, yPntsMax, ClmaxDistr

    def calcCD0wing(self, thicknesstochord, frictioncoefficient=0.0055):
        formfactor = 1 + 2.7*thicknesstochord + 100*thicknesstochord**4
        return 2*frictioncoefficient*formfactor

    def calcCD0wing(S, b, taper):
        # DO NOT USE WITH LOW TAPER RATIOS
        croot = self.calculateChord(np.pi/2, taper, S, b)
        ctip = self.calculateChord(0, taper, S, b)
        cmean = np.mean([croot, ctip])
        return croot*self.Cd0_r/(2*cmean) + ctip*self.Cd0_t/(2*cmean)

    def calcCDi(self, alpha):
        
        def _delta(A):
            nsum = 0
            for n in range(A.size):
                i = n
                n = 2*n+1
                nsum += n*(A[i]/A[0])**2
            return nsum
        
        A = np.array([ A[0]*alpha + A[1] for A in self.coeff ])
        CL = self.calcCL(alpha)
        return CL**2/(np.pi*self.A) * (1 + _delta(A))

def plotLiftDistribution(y, Cl_range, ClmaxDistr=None):
    fig = plt.figure(figsize=(10, 4.5))
    ax1 = fig.add_subplot(111)

    for distr in Cl_range:
        ax1.plot(y, distr, linewidth=2, color='green', marker='o', fillstyle='none', markevery=5)
        if ClmaxDistr: ax1.plot(y, ClmaxDistr(y), linewidth=1.5, linestyle='-.', color='black')

    ax1.axvline(x=0, linewidth=2, color='black')
    ax1.axhline(y=0, linewidth=2, color='black')
    ax1.set_xlabel('Wingspan [m]')
    ax1.set_ylabel('')
    ax1.xaxis.grid(color='black', linestyle='--')
    ax1.yaxis.grid(color='black', linestyle='--')
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.show()


if __name__ == "__main__":

    wing = WingPlanform(v.S, v.A, v.taper, v.twist, v.gamma)
    wing.setAirfoils(v.Clmax_r, v.Clmax_t, v.Cla_r, v.Cla_t, v.a0_r, v.a0_t, v.Cd0_r, v.Cd0_t)
    wing.calcCoefficients(100)

    """ Test convergence """
    if False:
        CLlist = []
        alpha = np.radians(13)
        Nmax = 100
        for N in range(1, Nmax):
            coeff = wing.calcCoefficients(N)
            A1 = coeff[0][0]*alpha + coeff[0][1]
            CL = np.pi * v.A * A1
            print(f'N: {N}')
            CLlist.append(CL)
        plt.plot(range(1, Nmax), CLlist)
        plt.show()

    """ Lift Distribution """
    if True:
        alpha = np.radians(5)
        CL = wing.calcCL(alpha)
        print(f'CL = {round(CL, 2)} @ a = {round(np.degrees(alpha), 2)}')

        CDi = wing.calcCDi(alpha)
        print(f'CDi = {CDi} @ a = {round(np.degrees(alpha), 2)}')

        Cl_distr, yPnts = wing.calcLiftDistribution(alpha, 100)
        plotLiftDistribution(yPnts, [Cl_distr])

    # Cl_distr_range = []
    # # alphaRange = [5]
    # alphaRange = np.arange(-15, 15, 1)
    # for alpha in alphaRange:
    #     alpha = np.radians(alpha)
    #     Cl_distr, yPnts = wing.calcLiftDistribution(alpha)
    #     Cl_distr_range.append(Cl_distr)
    # plotLiftDistribution(yPnts, Cl_distr_range)

    """ CLa """
    CLa = wing.calcCLa()
    print(f'CLa = {CLa}')

    """ CLmax """
    if False:
        CLmax, alphaMax, Cl_distrMax, yPntsMax, ClmaxDistr = wing.calcCLmax()
        print(f'CLmax = {round(CLmax, 2)} @ a = {round(np.degrees(alphaMax), 2)}')
        plotLiftDistribution(yPntsMax, [Cl_distrMax], ClmaxDistr)