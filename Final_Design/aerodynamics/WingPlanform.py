import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os
from math import pi,tan,atan,cos,sqrt

class WingPlanform:
    def __init__(self, S, A, taper, twist, gamma, CD0):
        self.S = S
        self.A = A
        self.taper = taper
        self.twist = twist
        self.gamma = gamma

        self.b = np.sqrt(A*S)
        self.c_r = (2*S)/(self.b*(1 + taper))
        self.c_t = taper*self.c_r
        self.MAC = (2/3)*self.c_r*((1 + taper + taper**2)/(1 + taper))
        self.sweepLE = np.arctan(-self.c_r/(2*self.b)*(taper - 1))
        self.YMAC = self.b/6 * (1 + 2*taper)/(1 + taper)
        self.XMAC = self.YMAC*np.tan(self.sweepLE)

        self.Clmax_r = None
        self.Clmax_t = None
        self.Cla_r = None
        self.Cla_t = None
        self.a0_r = None
        self.a0_t = None

        self.CD0 = CD0
        self.hwl = None
        self.kwl = None

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

    def setWinglets(self, hwl, kwl):
        self.hwl = hwl
        self.kwl = kwl

    def transformTheta(self, theta, b): # Verified
        return -.5*b*np.cos(theta)

    def transformSpan(self, y, b): # Verified
        return -np.arccos(2*y/b)
    
    def calculateChord(self, theta, taper, S, b): # Verified
        y = self.transformTheta(theta, b)
        return 2*(self.c_t - self.c_r)/b * abs(y) + self.c_r

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
        yPnts = np.linspace(-self.b/2, self.b/2, N)
        for y in yPnts:
            theta = -self.transformSpan(y, self.b)
            Cl = (2 * _circ(theta))/(self.calculateChord(theta, self.taper, self.S, self.b))
            Cl_distr.append(Cl)

        return Cl_distr, yPnts

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
            stallProgression = []

            for alpha in np.arange(alphaMax, alphaRange[-1], np.radians(alphaStep)):
                Cl_distr, yPnts = self.calcLiftDistribution(alpha, 100)
                Cl_distr = np.array_split(Cl_distr, 2)[1]
                yPnts = np.array_split(yPnts, 2)[1]

                idxs = np.argwhere(np.diff(np.sign(Cl_distr - ClmaxDistr(yPnts)))).flatten()

                for idx in idxs:
                    stallProgression.append([yPnts[idx], alpha])

            stallProgression = sorted(stallProgression, key=lambda dataPnt: dataPnt[0])
            stallProgression = np.array(stallProgression)

            fig = plt.figure(figsize=(10, 4.5))
            ax1 = fig.add_subplot(111)

            ax1.plot(stallProgression[:,0], np.degrees(stallProgression[:,1]), linewidth=2, color='red', marker='', fillstyle='none', markevery=4, label='stall onset')
            ax1.plot(-1*stallProgression[:,0], np.degrees(stallProgression[:,1]), linewidth=2, color='red', marker='', fillstyle='none', markevery=4)

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

    def calcCD0wing(self, S, b, taper):
        # DO NOT USE WITH LOW TAPER RATIOS
        croot = self.calculateChord(np.pi/2, taper, S, b)
        ctip = self.calculateChord(0, taper, S, b)
        cmean = np.mean([croot, ctip])
        return croot*self.Cd0_r/(2*cmean) + ctip*self.Cd0_t/(2*cmean)

    def calcDelta(self, alpha):
        A = np.array([ A[0]*alpha + A[1] for A in self.coeff ])
        deltasum = 0
        for n in range(1, A.size):
            i = n
            n = 2*n+1
            deltasum += n*(A[i]/A[0])**2
        return deltasum

    def calcCDi(self, alpha):
        ###!!! ONLY RUN FOR TIPCUTOFF < 0.8 !!!###

        A = np.array([ A[0]*alpha + A[1] for A in self.coeff ])

        def _nsum():
            nsum = 0
            for n in range(A.size):
                i = n
                n = 2*n+1
                nsum += n*A[i]**2
            return nsum

        CDi = np.pi*self.A * _nsum()

        return CDi

    def calcespan(self):
        elist = []
        for alpha in np.radians(np.linspace(2, 10, 3)):
            elist.append( 1/(1 + self.calcDelta(alpha)) )

        return sum(elist)/len(elist)

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

    def calcOswald(self, fuselagewidth, hasWinglets=False):
        k_fuselage = 1-2*(fuselagewidth/self.b)**2
        Q = 1/(self.calcespan()*k_fuselage)
        P = 0.38*self.CD0
        k_winglet = (1+2*self.hwl/(self.kwl*self.b))**2
        
        if not hasWinglets:
            return 1/(Q+P*np.pi*self.A)

        else:
            return k_winglet/(Q+P*np.pi*self.A)

    def calcsweep(self,pc, b, sweepc4, taper, cr):
        sweepLE = atan(tan(sweepc4) - cr / (2 * b) * (taper - 1))
        return atan(tan(sweepLE) + 2 * pc * cr * (taper - 1) / b)

    def calcchord(self,y, b, sweepc4, taper, cr):
        sweepLE = self.calcsweep(0, b, sweepc4, taper, cr)
        sweepTE = self.calcsweep(1, b, sweepc4, taper, cr)
        return cr - y * tan(sweepLE) + y * tan(sweepTE)

    def aileron_sizing(self,VTO,VL,fuselagewidth,clear_tip):
        # Iteration parameters
        sizing = True
        i = 0                       # [-]       Number of iterations
        max_i = 1000                # [-]       Maximum number of iterations
        step = 0.01                 # [m]       Increase/Decrease of aileron length at every iteration

        # Inputs
        dphi = 60*pi/180            # [rad]     Bank angle
        dtTO = 5                    # [s]       Bank time take-off
        dtL = 4                     # [s]       Bank time landing
        #VTO = 1.05*25.2             # [m/s]     Take-off speed
        #VL = 1.1*25.2               # [m/s]     Landing speed

        bf = fuselagewidth          # [m]       Fuselage width

        S = self.S                  # [m2]      Wing surface area
        b = self.b                  # [m]       Wing span
        taper = self.taper          # [-]       Wing taper ratio
        cr = self.c_r               # [m]       Wing root chord
        cla = self.calcCLa()        # [/rad]    Lift curve slope
        cd0_TO = 0.02               # [-]       Take-off configuration zero lift drag coefficient
        cd0_L = 0.02                # [-]       Landing configuration zero lift drag coefficient

        b1 = 6.36                   # [m]       Aileron start
        #clear_tip = 0.05*(b/2)      # [m]       Distance from the tip that should be clear of control surfaces
        da = 20*pi/180              # [rad]     Aileron deflection angle
        clda = 0.046825*180/pi      # [/rad]    Change in the airfoil’s lift coefficient with aileron deflection

        sm = 0.1  # [-]       Safety margin

        # Parameter calculations
        # Required roll rate
        p_reqTO = dphi/dtTO * (1+sm)
        p_reqL = dphi/dtL * (1+sm)

        # Aileron end
        b2 = b/2 - clear_tip

        # Check initial b1
        if b1 < bf / 2:
            print('Initial value for b1 is too small! Change b1.')
        else:
            pass

        # Create history list
        b1lst = [b1]

        # Roll rate calculation
        def roll_rate(V, cd0):
            # Calculate roll damping
            Clp = -(cla + cd0) * cr * b / (24 * S) * (1 + 3 * taper)

            # Calculate roll authority
            Clda = clda * cr / (S * b) * ((b2 ** 2 - b1 ** 2) + 4 * (taper - 1) / (3 * b) * (b2 ** 3 - b1 ** 3))

            # Calculate roll rate for take-off and landing
            p = -Clda / Clp * da * 2 * V / b

            return p

        # Perform iterations
        while sizing and i < 100:
            # Calculate roll rate
            p_TO = roll_rate(VTO, cd0_TO)
            p_L = roll_rate(VL, cd0_L)

            # Check whether p is larger than required
            # If p is smaller than required
            if p_TO < p_reqTO or p_L < p_reqL:
                b1 -= step

            # If p is larger than required
            else:
                # Check for convergence
                if abs(b1 - b1lst[-1]) < (step / 2):
                    sizing = False
                else:
                    b1 += step

            # Add values to history lists
            b1lst.append(b1)

            # Add iteration
            i += 1

            # Check aileron size
            if b1 < bf / 2:
                print('Aileron became too large.')
                sizing = False
            # Check for maximum number of iterations
            elif i >= max_i:
                sizing = False
                print('Aileron sizing did not converge within', i, 'iterations.')
            else:
                pass

        return b1,b2

    def flap_sizing(self, fix_position='fuselage end'):
        # Check input
        fix_positionlst = ['fuselage end', 'aileron start']
        if not fix_position in fix_positionlst:
            print("Wrong fix_position input (" + fix_position + "). Choose 'fuselage end' or 'aileron start'")
        else:
            pass

        # Inputs
        S = self.S                  # [m2]      Wing surface area
        b = self.b                  # [m]       Wing span
        sweepc4 = 0                 # [rad]     Wing quarter chord sweep angle
        taper = self.taper          # [rad]     Wing taper ratio
        cr = self.c_r               # [m]       Wing root chord

        CLmax_req = 2               # [-]       Required maximum lift coefficient
        CLmax_wing = 1.51           # [-]       Wing maximum lift coefficient
        CLa = 2 * pi                # [/rad]    Wing lift curve slope

        dClmax = 1.25 * 0.96        # [-]
        da0l_airfoil = -15*pi/180   # [rad]

        cfc = 0.8                   # [-]       Start of the flap as percentage of the chord

        sm = 0.1                    # [-]       Safety margin

        # Parameter calculations
        # Flap star/end location
        # Chord at flap start/end location
        if fix_position == 'fuselage end':
            bf = 1.5                # [m]       Fuselage width
            d_ff = 0.05             # [m]       Spacing between fuselage and flap
            f1 = bf / 2 + d_ff
            cf1 = self.calcchord(f1, b, sweepc4, taper, cr)
        #else:
        #    b1 = variables.b1       # [m]       Aileron start
        #    d_af = 0.05             # [m]       Spacing between flap and aileron
        #    f2 = b1 - d_af
        #    cf2 = chord(f2, b, sweepc4, taper, cr)
        # Leading edge sweep angle
        sweepLE = self.calcsweep(0, b, sweepc4, taper, cr)

        # Hinge line sweep angle
        sweep_hinge = self.calcsweep(cfc, b, sweepc4, taper, cr)

        # Trailing edge sweep angles
        sweepTE = self.calcsweep(1, b, sweepc4, taper, cr)

        # Increase in lift coefficient
        dCLmax = (CLmax_req - CLmax_wing) * (1 + sm)

        # Required flapped surface
        SwfS = dCLmax / (0.9 * dClmax * cos(sweep_hinge))
        Swf = SwfS * S

        # Shift in zero lift angle of attack
        da0L = da0l_airfoil * SwfS * cos(sweep_hinge)

        # Change in lift curve slope
        CLa_flapped = CLa

        # Flap span calculation
        # Solving the equation:
        # 'fuselage end': cf1*bfl - 0.5*bf^2*tan(sweepLE) + 0.5*bf^2*tan(sweepTE) = Swf/2
        # 'aileron start': cf2*bfl + 0.5*bf^2*tan(sweepLE) - 0.5*bf^2*tan(sweepTE) = Swf/2
        # a*bfl^2 + b*bfl + c = 0

        # Fuselage end
        if fix_position == 'fuselage end':
            A = 0.5 * (-tan(sweepLE) + tan(sweepTE))
            B = cf1
            C = -Swf / 2

            D = B ** 2 - 4 * A * C

            if not D >= 0:
                print('There is a problem with the flap sizing!')
                return
            else:
                bfllst = [0, 0]
                bfllst[0] = (-B + sqrt(D)) / (2 * A)
                bfllst[1] = (-B - sqrt(D)) / (2 * A)
                if bfllst[0] > 0 and (f1 + bfllst[0]) < (b / 2):
                    bfl = bfllst[0]
                elif bfllst[1] > 0 and (f1 + bfllst[1]) < (b / 2):
                    bfl = bfllst[1]
                else:
                    print("Flap is too large, it doesn't fit on the wing!")

        # Aileron start
        #else:
        #    A = 0.5 * (tan(sweepLE) - tan(sweepTE))
        #    B = cf2
        #    C = -Swf / 2

        #    D = B ** 2 - 4 * A * C

        #    if not D >= 0:
        #        print('There is a problem with the flap sizing!')
        #        return
        #    else:
        #        bfllst = [0, 0]
        #        bfllst[0] = (-B + sqrt(D)) / (2 * A)
        #        bfllst[1] = (-B - sqrt(D)) / (2 * A)
        #        if bfllst[0] > 0 and (f2 - bfllst[0]) > 0:
        #            bfl = bfllst[0]
        #        elif bfllst[1] > 0 and (f2 - bfllst[1]) > 0:
        #            bfl = bfllst[1]
        #        else:
        #            print("Flap is too large, it doesn't fit on the wing!")
        else:
            pass

        # Flap start/end location
        if fix_position == 'fuselage end':
            f2 = f1 + bfl
        #else:
        #    f1 = f2 - bfl
        else:
            pass

        return f1,f2,Swf