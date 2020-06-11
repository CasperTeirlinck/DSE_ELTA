import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os

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

    def setAirfoils(self, Clmax_r, Clmax_t, Cla_r, Cla_t, a0_r, a0_t, Cd0_r, Cd0_t, deltaAlphaStall_r=0, deltaAlphaStall_t=0):
        self.Clmax_r = Clmax_r
        self.Clmax_t = Clmax_t
        self.Cd0_r = Cd0_r
        self.Cd0_t = Cd0_t
        self.Cla_r = Cla_r
        self.Cla_t = Cla_t
        self.a0_r = a0_r
        self.a0_t = a0_t
        self.deltaAlphaStall_r = deltaAlphaStall_r
        self.deltaAlphaStall_t = deltaAlphaStall_t

    def transformTheta(self, theta, b): # Verified
        return -.5*b*np.cos(theta)

    def transformSpan(self, y, b): # Verified
        return -np.arccos(2*y/b)
    
    def calculateChord(self, theta, taper, S, b): # Verified
        y = self.transformTheta(theta, b)
        cr = (2*S)/(b*(1+taper))
        ct = taper*cr
        return 2*(ct-cr)/b * abs(y) + cr

    def calcCoefficients(self, N, tipCutoff=0.9, FuselageIncluded=False): # Verified without lift slope & twist implementation
        
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
            
            return 0 #Rdu - 1

        matrix = np.ndarray((N,N)) # Create sample matrix
        column2 = np.zeros((N,1)) # Create column for twist and fuselage contributions

        samplepoints = np.linspace(( self.b/2*tipCutoff - np.pi/2)/N, np.pi/2, N)

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

        def _alphai(theta):
            nsum = 0
            for n in range(A.size):
                i = n
                n = 2*n+1
                nsum += n*A[i]* np.sin(n*theta)/np.sin(theta)
            return nsum

        Cl_distr = []
        CDi_distr = []
        yPnts = np.linspace(-self.b/2, self.b/2, N)
        # yPntsCDi = []
        for y in yPnts:
            theta = -self.transformSpan(y, self.b)
            Cl = (2 * _circ(theta))/(self.calculateChord(theta, self.taper, self.S, self.b))
            Cl_distr.append(Cl)
            # if abs(y) < self.b/2*0.9:
                # CDi_distr.append(Cl*_alphai(theta))
                # yPntsCDi.append(y)

        return Cl_distr, yPnts#, CDi_distr, yPntsCDi

    def calcCLa(self):
        CLa = np.pi*self.A*self.coeff[0][0]
        return CLa
    
    def calcCLmax(self, plotProgression=False, printMaxLoc=False):

        alphaStep = 0.2
        alphaRange = np.radians(np.arange(8, 20, alphaStep))
        ClmaxDistr = lambda y: (self.Clmax_t - self.Clmax_r)/(self.b/2) * abs(y) + self.Clmax_r

        alphaMax = None
        alphaMaxLoc = None
        for alpha in alphaRange:
            if alphaMax: break

            Cl_distr, yPnts = self.calcLiftDistribution(alpha, 100)
            Cl_distr = np.array_split(Cl_distr, 2)[1]
            yPnts = np.array_split(yPnts, 2)[1]

            for Cl, y in zip(Cl_distr, yPnts):
                if np.abs(Cl - ClmaxDistr(y)) <= 0.01:
                    alphaMax = alpha
                    Cl_distrMax = Cl_distr
                    yPntsMax = yPnts
                    CLmax = self.calcCL(alphaMax)
                    alphaMaxLoc = y
                    if printMaxLoc: print(f'CLmax location @ y = {round(alphaMaxLoc, 2)}')
                    break

        if plotProgression:
            stallAlphaDistrSmooth = lambda y: (self.deltaAlphaStall_t - self.deltaAlphaStall_r)/(self.b/2) * abs(y) + self.deltaAlphaStall_r
            stallProgression = []
            stallProgressionSmooth = []

            for alpha in np.arange(alphaMax, alphaRange[-1], np.radians(alphaStep)):
                Cl_distr, yPnts = self.calcLiftDistribution(alpha, 100)
                Cl_distr = np.array_split(Cl_distr, 2)[1]
                yPnts = np.array_split(yPnts, 2)[1]

                idxs = np.argwhere(np.diff(np.sign(Cl_distr - ClmaxDistr(yPnts)))).flatten()

                for idx in idxs:
                    stallProgression.append([yPnts[idx], alpha])
                    stallProgressionSmooth.append([yPnts[idx], alpha + stallAlphaDistrSmooth(yPnts[idx])])

            stallProgression = sorted(stallProgression, key=lambda dataPnt: dataPnt[0])
            stallProgression = np.array(stallProgression)
            stallProgressionSmooth = sorted(stallProgressionSmooth, key=lambda dataPnt: dataPnt[0])
            stallProgressionSmooth = np.array(stallProgressionSmooth)

            fig = plt.figure(figsize=(10, 4.5))
            ax1 = fig.add_subplot(111)

            ax1.plot(stallProgression[:,0], np.degrees(stallProgression[:,1]), linewidth=2, color='red', marker='', fillstyle='none', markevery=4, label='stall onset')
            ax1.plot(-1*stallProgression[:,0], np.degrees(stallProgression[:,1]), linewidth=2, color='red', marker='', fillstyle='none', markevery=4)

            ax1.plot(stallProgressionSmooth[:,0], np.degrees(stallProgressionSmooth[:,1]), linewidth=2, color='blue', linestyle='-.', label='full stall')
            ax1.plot(-1*stallProgressionSmooth[:,0], np.degrees(stallProgressionSmooth[:,1]), linewidth=2, color='blue', linestyle='-.')

            ax1.axvline(x=0, linewidth=2, color='black')
            ax1.axhline(y=0, linewidth=2, color='black')
            ax1.set_xlabel('Wingspan [m]')
            ax1.set_ylabel('alpha [deg]')
            ax1.xaxis.grid(color='black', linestyle='--')
            ax1.yaxis.grid(color='black', linestyle='--')
            plt.legend(loc='lower right')
            fig.suptitle('Spanwise Stall Progression', fontsize=16, y=0.97)
            plt.tight_layout(rect=[0, 0, 1, 0.93])
            plt.show()

        if not alphaMax:
            print('alphaMax not found')
            return None
        return CLmax, alphaMax, Cl_distrMax, yPntsMax, ClmaxDistr, alphaMaxLoc

    def calcCD0wing(S, b, taper):
        # DO NOT USE WITH LOW TAPER RATIOS
        croot = self.calculateChord(np.pi/2, taper, S, b)
        ctip = self.calculateChord(0, taper, S, b)
        cmean = np.mean([croot, ctip])
        return croot*self.Cd0_r/(2*cmean) + ctip*self.Cd0_t/(2*cmean)

    def calcCDi(self, alpha):
        ###!!! ONLY RUN FOR TIPCUTOFF < 0.8 !!!###

        def _delta(A):
            deltasum = 0
            for n in range(1, A.size):
                i = n
                n = 2*n+1
                deltasum += n*(A[i]/A[0])**2
            return deltasum

        def _nsum(A):
            nsum = 0
            for n in range(A.size):
                i = n
                n = 2*n+1
                nsum += n*A[i]**2
            return nsum
        
        A = np.array([ A[0]*alpha + A[1] for A in self.coeff ])

        CL = self.calcCL(alpha)
        CDi1 = (CL**2)/(np.pi*self.A) * (1 + _delta(A))
        CDi2 = np.pi*self.A * _nsum(A)

        e = 1/(1+_delta(A))
        return CDi1, e   

    def calcAlphai(self, alpha, N):

        def _alphai(A, theta):
            nsum = 0
            for n in range(A.size):
                i = n
                n = 2*n+1
                nsum += n*A[i]* np.sin(n*theta)/np.sin(theta)
            return nsum

        A = np.array([ A[0]*alpha + A[1] for A in self.coeff ])

        alphai_distr = []
        tipCutoff = 0.9
        yPnts = np.linspace(-self.b/2 * tipCutoff, self.b/2 * tipCutoff, N)
        for y in yPnts:
            theta = self.transformSpan(y, self.b)
            alphai_distr.append(-_alphai(A, theta))

        return alphai_distr, yPnts

