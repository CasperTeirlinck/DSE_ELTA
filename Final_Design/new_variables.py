import numpy as np
from math import sqrt,cos,pi
import matplotlib.pyplot as plt


class NewVariables:
    def __init__(self,haswinglets,wingletheight):
        self.init_general()
        self.init_aerodynamics(haswinglets,wingletheight)
        self.init_weight()
        self.init_propulsion()
        self.init_sc()

    def init_general(self):
        self.WTO = None
        self.Woew_classII = None
        self.WPL = 1961.33
        self.rho0 = 1.225
        self.g0 = 9.80665
        self.R = 287.05
        self.gamma = 1.4 
        self.lmbda = -0.0065

        self.Vcruise = 50           # [m/s]         Cruise speed
        self.hcruise = 914.4        # [m]           Cruise altitude

        self.fuselagelength = 9.420
        self.fuselagewidth = 1.05
        self.fuselageheight = 1.211
        self.fuselagefrontalarea = 1.218
        self.fuselagewettedarea = 17.507

        self.h_landinggear = None
        self.w_landinggear = None

        self.h_htail = None

        self.W_wing = None
        self.W_fgroup = None

        self.cg_wing = None         # [m]           Distance LE root chord - wing cg

        self.xcg_fgroup = None
        self.xcg_fuselage = None
        self.xcgPL = None
        self.xcgbat = None
        self.xcg_min = None
        self.xcg_max = None

    def init_aerodynamics(self,haswinglets,wingletheight):
        # Conditions
        self.clean_config = None
        self.hasWinglets = haswinglets
        
        # General wing geometry
        self._S = 300
        self._A = 2
        self._taper = 0.467
        self._twist = np.radians(6)
        self._gammaL = 0

        # Complementary wing geometry
        self._b = np.sqrt(self.A * self.S)
        self._c_r = (2 * self.S) / (self.b * (1 + self.taper))
        self.c_t = self.taper * self.c_r
        self.MAC = (2 / 3) * self.c_r * ((1 + self.taper + self.taper ** 2) / (1 + self.taper))
        self._sweepLE = np.arctan(-self.c_r / (2 * self.b) * (self.taper - 1))
        self._YMAC = self.b / 6 * (1 + 2 * self.taper) / (1 + self.taper)
        self.XMAC = self.YMAC * np.tan(self.sweepLE)
        self.Snet = self.S - self.calculateChord(self.transformSpan(.25*self.fuselagewidth,self.b),self.taper,self.S,self.b)*self.fuselagewidth

        # Flap geometry
        self.flapstart = None
        self.flapend = None
        self.flapspan = None
        self.flapaffectedarea = None

        # Exterior design inputs
        self.BLturbratio_wing = 0.65
        self.BLturbratio_fus = 1
        self.BLturbratio_emp = 0.65


        # Constants & statistical variables
        self.kwl = 2.1
        self.dCD_landinggear = 35

        # Airfoil properties
        self.Clmax_r = 1.6
        self.Clmax_t = self.Clmax_r
        self.Cla_r = 1/np.radians(10)
        self.Cla_t = self.Cla_r
        self.a0_r = np.radians(-4)
        self.a0_t = self.a0_r

        # Freely adjustable variables
        self.hwl = wingletheight

        # Flap variables
        self.CLmax_req = 2                  # [-]       Required maximum lift coefficient
        self.dClmax = 1.25 * 0.96           # [-]
        self.da0l_airfoil = -15*pi/180      # [rad]
        self.cfc = 0.8                      # [-]       Start of the flap as percentage of the chord
        self.d_ff = 0.05                    # [m]       Spacing between fuselage and flap

        # Aileron variables
        self.dphi = 60*pi/180               # [rad]     Bank angle
        self.dtTO = 5                       # [s]       Bank time take-off
        self.dtL = 4                        # [s]       Bank time landing
        self.clear_tip = 0.05*(b/2)         # [m]       Distance from the tip that should be clear of control surfaces
        self.da = 20*pi/180                 # [rad]     Aileron deflection angle
        self.clda = 0.046825*180/pi         # [/rad]    Take-off configuration change in the airfoil’s lift coefficient with aileron deflection


        # Output variables
        self.wing_CL_alpha = None
        self.wing_CL_max = None
        self.wing_alpha_max = None
        self.CD0clean = 0.01912
        self.CD0flap = None
        if haswinglets:
            self.eclean = 0.84
        if not haswinglets:
            self.eclean = 0.79
        self.eflaps = None

        # Working variables
        self.coeff = None

    def init_weight(self):
        pass

    def init_propulsion(self):
        # Sizing
        self.W_batt       = None           # Battery weight in Newtons
        self.v_batt       = None           # Battery volume in liters
        self.W_motor      = 30 * 9.80665   # Motor weight in Newtons
        self.W_shaft      = 4.48 * 9.80665 # Engine shaft weight in Newtons
        self.W_prop       = 12 * 9.80665   # Propeller weight in Newtons
        self.P_max        = 65 * 1000      # Maximum power produced by the engine in W

        # Battery characteristics
        self.batt_E_prod  = None           # Energy produced by the battery in Wh
        self.batt_Ns      = None           # Number of battery cells in series
        self.batt_Np      = None           # Number of battery cells in parallel
        self.batt_N       = None           # Number of battery cells
        self.eff_batt     = 0.95           # Battery efficiency
        self.eff_prop     = 0.85           # Propeller efficiency
        self.eff_tot_prop = 0.95 * 0.88    # Total propulsive efficiency
        self.V_req_batt   = 400            # Required voltage in Volts
        self.I_req_batt   = 189.75         # Required current in Amps
        self.DoD          = 90             # Depth of discharge of the battery

        # Requirements
        self.range_m      = 250*1000       # Range requirement in meters
        self.endurance_s  = 2.5*3600       # Endurance requirement in seconds

        # Battery cell data of the Samsung 21700 50E
        self.batt_cell_diameter   = 0.0211                                                          # Cell diameter [m]
        self.batt_cell_length     = 0.0707                                                          # Cell length [m]
        self.batt_cell_volume     = np.pi * (self.batt_cell_diameter / 2) ** 2 * self.batt_cell_length # [m^3]
        self.batt_cell_mass       = 0.0687                                                          # Cell mass [kg]
        self.batt_cell_V_nom      = 3.6                                                             # Nominal voltage [V]
        self.batt_cell_V_max      = 4.2                                                             # Maximum voltage [V]
        self.batt_cell_V_cut      = 2.5                                                             # Cut-off voltage [V]
        self.batt_cell_I_max      = 9.8                                                             # Maximum discharge current [A]
        self.batt_cell_C_Ah       = 5.0                                                             # Capacity [Ah]
        self.batt_cell_C_Wh       = self.batt_cell_C_Ah * self.batt_cell_V_nom                      # Capacity [Wh]
        self.batt_cell_E_spec     = self.batt_cell_C_Wh / self.batt_cell_mass                       # Capacity per kg of weight [Wh/kg]
        self.batt_cell_E_vol_spec = self.batt_cell_C_Wh / self.batt_cell_volume                     # Capacity per unit volume [Wh/m^3]
        self.batt_cell_P          = self.batt_cell_I_max * self.batt_cell_V_nom                     # Maximum power [W]

    def init_sc(self):
        # Fuselage variables
        self.Sfs = 5.258            # [m2]          Fuselage lateral surface area
        self.hf1 = 1.146            # [m]           Fuselage nose height
        self.hf2 = 0.306            # [m]           Fuselage tail height
        self.bf1 = 0.960            # [m]           Fuselage nose width
        self.bf2 = 0.243            # [m]           Fuselage tail width


        # Propeller geometry parameters
        self.Bp = 2.5               # [-]           Number of blades per propeller
        self.xp1 = 2                # [m]           1st propeller location
        self.xp2 = 2                # [m]           2nd propeller location
        self.Dp1 = 2                # [m2]          1st propeller disk diameter
        self.Dp2 = 2                # [m2]          2nd propeller disk diameter


        # Wing variables
        self.xlemac = None          # [m]           Distance nose - leading edge mean aerodynamic chord
        self.Snet = None            # [m2]          Net wing surface area

        self.VhV = sqrt(0.85)       # [-]           Tail/wing speed ratio
        self.eta = 0.95             # [-]           Airfoil efficiency coefficient
        self.Cm0af = -0.1           # [-]           Airfoil zero lift pitching moment coefficient
        self.mu1 = 0.3              # [-]           Flap coefficient 1
        self.deda = None            # [-]           Downwash gradient
        self.cc = 1                 # [-]           Chord ratio (extended flap/clean)
        self.CnBi = 0.024           # [-]           Wing configuration lateral stability component


        # Horizontal tail variables
        self.lh = 6.9               # [m]           Horizontal tail arm
        self.Sh = None              # [m2]          Minimum required horizontal tail surface
        self.bh = None              # [m]           Horizontal tail span
        self.Ah = None              # [-]           Horizontal tail aspect ratio
        self.sweeph = 0             # [rad]         Horizontal tail half chord sweep
        self.ch_r = None            # [m]           Horizontal tail root chord
        self.ih = 0                 # [rad]         Horizontal tail incidence angle

        self.CLh_L = -0.8           # [-]           Horizontal tail landing configuration lift coefficient


        # Vertical tail variables
        self.lv = 6.9               # [m]           Vertical tail arm
        self.Sv = None              # [m2]          Vertical tail surface

        self.CnB = None             # [-]           Directional stability coefficient

        '''
        # Flight performance parameters
        self.a_pitch = 12*pi/180    # [rad/s]       Take-off pitch angular velocity
        self.VTO = 1.05 * 25.2      # [m/s]         Take-off velocity
        self.rhoTO = 1.225          # [kg/m3]       Take-off density
        self.VL = 1.1 * 25.2
        self.mu = 0.05              # [-]           Take-off friction factor

        # Centre of gravity
        self.zcg = 1                # [m]           Centre of gravity height
        self.xcg_f = 2.3            # [m]           Fuselage centre of gravity location

        # Component locations
        self.xmg = 2                # [m]           Main gear location
        self.xacw = 1               # [m]           Wing/fuselage aerodynamic centre location
        self.xach = 6               # [m]           Horizontal tail aerodynamic centre location

        self.zmg = 0                # [m]           Main gear height
        self.zT = 1                 # [m]           Thrust vector height
        self.zD = 0.5               # [m]           Drag vector heigh

        # Elevator geometry parameters
        self.bebh = 1               # [-]           Elevator span
        self.de_max = 20*pi/180     # [rad]         Maximum elevator deflection
        self.de_min = -25*pi/180    # [rad]         Minimum elevator deflection (maximum negative)
        
        # Horizontal tail aerodynamic parameters
        self.CLh_TO = None          # [-]           Horizontal tail take-off lift coefficient
        self.CLah = 4               # [/rad]        Horizontal lift curve slope
        '''

    @property
    def S(self):
        return self._S

    @S.setter
    def S(self, val):
        self._S = val
        self.b = np.sqrt(self.A * self.S)
        self.Snet = self.S - self.calculateChord(self.transformSpan(.25*w_fuselage,self.b),self.taper,self.S,self.b)*w_fuselage
        # self.calcCoefficients()       # Maybe enable? Disabled for performance reasons.

    @property
    def A(self):
        return self._A

    @A.setter
    def A(self, val):
        self._A = val
        self.b = np.sqrt(self.A * self.S)

    @property
    def b(self):
        return self._b

    @b.setter
    def b(self, val):
        self._b = val
        self.c_r = (2 * self.S) / (self.b * (1 + self.taper))
        self.YMAC = self.b / 6 * (1 + 2 * self.taper) / (1 + self.taper)
        self.Snet = self.S - self.calculateChord(self.transformSpan(.25*w_fuselage,self.b),self.taper,self.S,self.b)*w_fuselage
        # self.calcCoefficients()       # Maybe enable? Disabled for performance reasons.

    @property
    def taper(self):
        return self._taper

    @taper.setter
    def taper(self, val):
        self._taper = val
        self.c_r = (2 * self.S) / (self.b * (1 + self.taper))
        self.MAC = (2 / 3) * self.c_r * ((1 + self.taper + self.taper ** 2) / (1 + self.taper))
        self.sweepLE = np.arctan(-self.c_r / (2 * self.b) * (self.taper - 1))
        self.YMAC = self.b / 6 * (1 + 2 * self.taper) / (1 + self.taper)
        self.Snet = self.S - self.calculateChord(self.transformSpan(.25*w_fuselage,self.b),self.taper,self.S,self.b)*w_fuselage
        # self.calcCoefficients()       # Maybe enable? Disabled for performance reasons.

    @property
    def twist(self):
        return self._twist

    @twist.setter
    def twist(self, val):
        self._twist = val
        # self.calcCoefficients()       # Maybe enable? Disabled for performance reasons.

    @property
    def gammaL(self):
        return self._gammaL

    @gammaL.setter
    def gammaL(self, val):
        self._gammaL = val
        # self.calcCoefficients()       # Maybe enable? Disabled for performance reasons.

    @property
    def c_r(self):
        return self._c_r

    @c_r.setter
    def c_r(self, val):
        self._c_r = val
        self.c_t = self.taper * self.c_r
        self.MAC = (2 / 3) * self.c_r * ((1 + self.taper + self.taper ** 2) / (1 + self.taper))
        self.sweepLE = np.arctan(-self.c_r / (2 * self.b) * (self.taper - 1))

    @property
    def sweepLE(self):
        return self._sweepLE

    @sweepLE.setter
    def sweepLE(self, val):
        self._sweepLE = val
        self.XMAC = self.YMAC * np.tan(self.sweepLE)

    @property
    def YMAC(self):
        return self._YMAC

    @YMAC.setter
    def YMAC(self, val):
        self._YMAC = val
        self.XMAC = self.YMAC * np.tan(self.sweepLE)


    ###################################

    def update_WTO(self):
        self.WTO = self.Woew_classII + self.Wbat + self.WPL

    def setAirfoils(self, Clmax_r, Clmax_t, Cla_r, Cla_t, a0_r, a0_t, deltaAlphaStall_r=0,
                    deltaAlphaStall_t=0):
        self.Clmax_r = Clmax_r
        self.Clmax_t = Clmax_t
        self.Cla_r = Cla_r
        self.Cla_t = Cla_t
        self.a0_r = a0_r
        self.a0_t = a0_t
        self.deltaAlphaStall_r = deltaAlphaStall_r
        self.deltaAlphaStall_t = deltaAlphaStall_t

    def setWinglets(self, hwl, kwl):
        self.hwl = hwl
        self.kwl = kwl

    def transformTheta(self, theta, b):  # Verified
        return -.5 * b * np.cos(theta)

    def transformSpan(self, y, b):  # Verified
        return -np.arccos(2 * y / b)

    def calculateChord(self, theta, taper, S, b):  # Verified
        y = self.transformTheta(theta, b)
        return 2 * (self.c_t - self.c_r) / b * abs(y) + self.c_r

    def calcCoefficients(self, N=150, tipCutoff=0.9,
                         FuselageIncluded=False):  # Verified without lift slope & twist implementation

        def _calcLiftSlope(theta, b, Cla_r, Cla_t):  # Verified
            y = self.transformTheta(theta, b)
            return 2. / self.b * (Cla_t - Cla_r) * abs(y) + Cla_r

        def _calculateTwistAngle(theta, b, twist, gammaL):  # Verified
            y = self.transformTheta(theta, b)
            b_half = .5 * b
            if gammaL != 0:
                C1 = twist / (1 - np.exp(gammaL * b_half))
                C2 = C1 * np.exp(gammaL * b_half)
                return C1 - C2 * np.exp(-gammaL * abs(y))
            else:
                return twist * (1 - abs(y) / b)

        def _calculateZeroLiftAngle(theta, b, twist, a0_r, a0_t, gammaL):
            alpha_geometric = _calculateTwistAngle(theta, b, twist, gammaL)
            y = self.transformTheta(theta, b)
            return 2. / b * (a0_t - a0_r) * abs(y) + a0_r - alpha_geometric

        def _calculateFuselageContribution():

            return 0  # Rdu - 1

        matrix = np.ndarray((N, N))  # Create sample matrix
        column2 = np.zeros((N, 1))  # Create column for twist and fuselage contributions

        samplepoints = np.linspace((self.b / 2 * tipCutoff - np.pi / 2) / N, np.pi / 2, N)

        for i in range(N):
            theta_sample = samplepoints[i]  # Use sample point i

            a0 = _calcLiftSlope(theta_sample, self.b, self.Cla_r,
                                self.Cla_t)  # Calculate the lift slope of sample point i
            c = self.calculateChord(theta_sample, self.taper, self.S, self.b)  # Calculate the chord of sample point i

            zeroliftangle = _calculateZeroLiftAngle(theta_sample, self.b, self.twist, self.a0_r, self.a0_t,
                                                    self.gamma)  # Calculate the zero lift angle
            fuselageangle = FuselageIncluded * _calculateFuselageContribution()  # Calculate the fuselage contribution to the angle of attack.
            column2[i] = - zeroliftangle - fuselageangle  # Calculate element (i,1) of the coefficient matrix

            for j in range(N):
                element = np.sin(theta_sample * (2 * j + 1)) * (4 * self.b / (a0 * c) + (2 * j + 1) / np.sin(
                    theta_sample))  # Calculate element (i,j) of the matrix

                # Add element to matrix; if element is close to zero, add zero.
                if abs(element) > 1e-4:
                    matrix[i, j] = element
                else:
                    matrix[i, j] = 0.

        matrix_inverse = la.inv(matrix)  # Calculate inverse of the matrix
        column1 = np.matmul(matrix_inverse, np.ones((N, 1)))
        column2 = np.matmul(matrix_inverse, column2)

        coefficientmatrix = np.concatenate((column1, column2), axis=1)  # Merge columns into one matrix

        self.coeff = coefficientmatrix

    def calcCL(self, alpha):

        A1 = self.coeff[0][0] * alpha + self.coeff[0][1]
        return np.pi * self.A * A1

    def calcLiftDistribution(self, alpha, N, showCircComponents=False):

        A = np.array([A[0] * alpha + A[1] for A in self.coeff])

        def _showCircComponents():
            testfig = plt.figure(figsize=(10, 4.5))
            testax = testfig.add_subplot(111)
            testy = np.linspace(-self.b / 2, self.b / 2, 1000)
            print(testy)
            for n in range(A.size):
                i = n
                n = 2 * n + 1

                term = lambda theta: 2 * self.b * A[i] * np.sin(n * theta)
                testax.plot(testy, term(self.transformSpan(testy, self.b)))

            plt.show()

        if showCircComponents: _showCircComponents()

        def _circ(theta):
            nsum = 0
            for n in range(A.size):
                i = n
                n = 2 * n + 1
                nsum += A[i] * np.sin(n * theta)
            return 2 * self.b * nsum

        def _alphai(theta):
            nsum = 0
            for n in range(A.size):
                i = n
                n = 2 * n + 1
                nsum += n * A[i] * np.sin(n * theta) / np.sin(theta)
            return nsum

        Cl_distr = []
        yPnts = np.linspace(-self.b / 2, self.b / 2, N)
        for y in yPnts:
            theta = -self.transformSpan(y, self.b)
            Cl = (2 * _circ(theta)) / (self.calculateChord(theta, self.taper, self.S, self.b))
            Cl_distr.append(Cl)

        return Cl_distr, yPnts

    def calcCLa(self):
        self.CL_alpha_wing = np.pi * self.A * self.coeff[0][0]

    def calcCLmax(self, plotProgression=False, printMaxLoc=False):

        alphaStep = 0.2
        alphaRange = np.radians(np.arange(8, 20, alphaStep))
        ClmaxDistr = lambda y: (self.Clmax_t - self.Clmax_r) / (self.b / 2) * abs(y) + self.Clmax_r

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

            ax1.plot(stallProgression[:, 0], np.degrees(stallProgression[:, 1]), linewidth=2, color='red', marker='',
                     fillstyle='none', markevery=4, label='stall onset')
            ax1.plot(-1 * stallProgression[:, 0], np.degrees(stallProgression[:, 1]), linewidth=2, color='red',
                     marker='', fillstyle='none', markevery=4)

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
        
        self.wing_CL_max = CLmax*(self.b-self.fuselagewidth)/self.b
        self.wing_alpha_max = alphaMax
        self.stallloc = alphaMaxLoc
        #return CLmax, alphaMax, Cl_distrMax, yPntsMax, ClmaxDistr, alphaMaxLoc

    def calcDelta(self, alpha):
        A = np.array([A[0] * alpha + A[1] for A in self.coeff])
        deltasum = 0
        for n in range(1, A.size):
            i = n
            n = 2 * n + 1
            deltasum += n * (A[i] / A[0]) ** 2
        return deltasum

    def calcCDi(self, alpha):
        ###!!! ONLY RUN FOR TIPCUTOFF < 0.8 !!!###

        A = np.array([A[0] * alpha + A[1] for A in self.coeff])

        def _nsum():
            nsum = 0
            for n in range(A.size):
                i = n
                n = 2 * n + 1
                nsum += n * A[i] ** 2
            return nsum

        CDi = np.pi * self.A * _nsum()

        return CDi

    def calcespan(self):
        elist = []
        for alpha in np.radians(np.linspace(2, 10, 3)):
            elist.append(1 / (1 + self.calcDelta(alpha)))

        return sum(elist) / len(elist)

    def calcAlphai(self, alpha, N):

        def _alphai(A, theta):
            nsum = 0
            for n in range(A.size):
                i = n
                n = 2 * n + 1
                nsum += n * A[i] * np.sin(n * theta) / np.sin(theta)
            return nsum

        A = np.array([A[0] * alpha + A[1] for A in self.coeff])

        alphai_distr = []
        tipCutoff = 0.9
        yPnts = np.linspace(-self.b / 2 * tipCutoff, self.b / 2 * tipCutoff, N)
        for y in yPnts:
            theta = self.transformSpan(y, self.b)
            alphai_distr.append(-_alphai(A, theta))

        return alphai_distr, yPnts

    def calcCD0(self,S_wet_fus,l_fus,fus_A_max,w_fuselage,S_h,S_v,MAC_emp,BLturbratio_fus, BLturbratio_wing, BLturbratio_emp,l_gear,w_gear,dCD_gear,MAC,flap_area_ratio,tc_airfoil=0.15,xc_airfoil=0.3,tc_emp=0.12,xc_emp=0.3,V_stall=23.15,rho_cruise=1.04,clean_config=True,visc=1.8e-5):
        
        def _CfLaminar(rho,V,L,visc):
            Re = rho*V*L/visc
            return 1.328/np.sqrt(Re)

        def _CfTurbulent(rho,V,L,visc):
            Re = rho*V*L/visc
            return 0.445/(np.log10(Re)**2.58)
        
        Cf_fus = (1-BLturbratio_fus)*_CfLaminar(rho_cruise,V_stall,l_fus,visc) + BLturbratio_fus*_CfTurbulent(rho_cruise,V_stall,l_fus,visc)
        ld_fus = l_fus/np.sqrt(4*fus_A_max/np.pi)
        FF_fus = 1 + 60./ld_fus**3 + ld_fus/400.
        IF_fus = 1.
        
        S_wet_wing = (self.S - self.calculateChord(self.transformSpan(.25*w_fuselage,self.b),self.taper,self.S,self.b)*w_fuselage)*2.06
        Cf_wing = (1-BLturbratio_wing)*_CfLaminar(rho_cruise,V_stall,MAC,visc) + BLturbratio_wing*_CfTurbulent(rho_cruise,V_stall,MAC,visc)
        FF_wing = 1. + 0.6*tc_airfoil/xc_airfoil + 100*tc_airfoil**4
        IF_wing = 1.25

        S_wet_emp = (S_h + S_v)*2.04
        Cf_emp = (1-BLturbratio_emp)*_CfLaminar(rho_cruise,V_stall,MAC_emp,visc) + BLturbratio_emp*_CfTurbulent(rho_cruise,V_stall,MAC_emp,visc)
        FF_emp = 1. + 0.6*tc_emp/xc_emp + 100*tc_emp**4
        IF_emp = 1.05

        CDS_wet_fus =  Cf_fus*  FF_fus*  IF_fus*  S_wet_fus
        CDS_wet_wing = Cf_wing* FF_wing* IF_wing* S_wet_wing
        CDS_wet_emp =  Cf_emp*  FF_emp*  IF_emp*  S_wet_emp
        CDS_ref_gear = dCD_gear*l_gear*w_gear
        
        dCD_flap = 0.0144*0.2*flap_area_ratio*(40-10)

        if clean_config:
            return 1.05*(CDS_wet_fus + CDS_wet_wing + CDS_wet_emp + CDS_ref_gear)/self.S

        if not clean_config:
            return 1.05*(CDS_wet_fus + CDS_wet_wing + CDS_wet_emp + CDS_ref_gear)/self.S + dCD_flap

    def calcOswald(self,S_wet_fus,l_fus,fus_A_max,w_fuselage,S_h,S_v,MAC_emp,BLturbratio_fus, BLturbratio_wing, BLturbratio_emp, l_gear,w_gear,dCD_gear,MAC,flap_area_ratio,tc_airfoil=0.15,xc_airfoil=0.3,tc_emp=0.12,xc_emp=0.3,V_stall=23.15,rho_cruise=1.04,clean_config=True,visc=1.8e-5,hasWinglets=False):
        k_fuselage = 1-2*(w_fuselage/self.b)**2
        Q = 1/(self.calcespan()*k_fuselage)
        P = 0.38*self.calcCD0(S_wet_fus,l_fus,fus_A_max,w_fuselage,S_h,S_v,MAC_emp,BLturbratio_fus, BLturbratio_wing, BLturbratio_emp,l_gear,w_gear,dCD_gear,flap_area_ratio,tc_airfoil,xc_airfoil,MAC,tc_emp,xc_emp,V_stall,rho_cruise,clean_config,visc)
        k_winglet = (1+2*self.hwl/(self.kwl*self.b))**2
        
        if not self.hasWinglets:
            self.e = 1/(Q+P*np.pi*self.A)

        else:
            self.e = k_winglet/(Q+P*np.pi*self.A)

    def flap_sizing(self):
        # Inputs
        S = self.S                      # [m2]      Wing surface area
        b = self.b                      # [m]       Wing span
        sweepc4 = 0                     # [rad]     Wing quarter chord sweep angle
        taper = self.taper              # [rad]     Wing taper ratio
        cr = self.c_r                   # [m]       Wing root chord

        CLmax_req = self.CLmax_req      # [-]       Required maximum lift coefficient
        CLmax_wing = self.calcCLmax()   # [-]       Wing maximum lift coefficient
        CLa = self.calcCLa()            # [/rad]    Wing lift curve slope

        dClmax = self.dClmax            # [-]
        da0l_airfoil = self.da0l_airfoil# [rad]

        cfc = self.cfc                  # [-]       Start of the flap as percentage of the chord

        sm = 0.1                        # [-]       Safety margin

        # Parameter calculations
        # Chord at flap start/end location
        bf = self.fuselagewidth         # [m]       Fuselage width
        d_ff = self.d_ff                # [m]       Spacing between fuselage and flap

        # Flap start location
        f1 = bf / 2 + d_ff
        cf1 = self.calcchord(f1, b, sweepc4, taper, cr)

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
        # cf1*bfl - 0.5*bf^2*tan(sweepLE) + 0.5*bf^2*tan(sweepTE) = Swf/2
        # A*bfl^2 + B*bfl + C = 0
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

            # Calculate flap end
            f2 = f1 + bfl

            # Update variables
            self.flapstart = f1
            self.flapend = f2
            self.flapspan = bfl
            self.flapaffectedarea = Swf

    def aileron_sizing(self):
        # Iteration parameters
        sizing = True
        i = 0                       # [-]       Number of iterations
        max_i = 1000                # [-]       Maximum number of iterations
        step = 0.01                 # [m]       Increase/Decrease of aileron length at every iteration

        # Inputs
        dphi = self.dphi            # [rad]     Bank angle
        dtTO = self.dtTO            # [s]       Bank time take-off
        dtL = self.dtL              # [s]       Bank time landing
        VTO = self.VTO              # [m/s]     Take-off speed
        VL = self.VL                # [m/s]     Landing speed

        bf = self.fuselagewidth     # [m]       Fuselage width

        S = self.S                  # [m2]      Wing surface area
        b = self.b                  # [m]       Wing span
        taper = self.taper          # [-]       Wing taper ratio
        cr = self.c_r               # [m]       Wing root chord
        cla = self.calcCLa()        # [/rad]    Take-off configuration lift curve slope
        cd0_TO = self.CD0clean      # [-]       Take-off configuration zero lift drag coefficient
        cd0_L = self.CD0flap        # [-]       Landing configuration zero lift drag coefficient

        b1 = self.b1                # [m]       Aileron start
        clear_tip = self.clear_tip  # [m]       Distance from the tip that should be clear of control surfaces
        da = self.da                # [rad]     Aileron deflection angle
        clda = self.clda            # [/rad]    Take-off configuration change in the airfoil’s lift coefficient with aileron deflection

        sm = 0.1                    # [-]       Safety margin

        # Parameter calculations
        # Required roll rate
        p_reqTO = dphi / dtTO * (1 + sm)
        p_reqL = dphi / dtL * (1 + sm)

        # Aileron end
        b2 = b / 2 - clear_tip

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

        # Transform to numpy arrays
        # b1lst = np.array(b1lst)
        # b2lst = np.array(b2lst)

        # Update values in variables class
        self.b1 = b1
        self.b2 = b2


def sys_Aerodynamics_wing(v,resolution):
    v.setAirfoils(v.Clmax_r, v.Clmax_t, v.Cla_r, v.Cla_t, v.a0_r, v.a0_t)
    v.calcCoefficients(resolution,0.7)
    v.calcCLa()
    v.calcCLmax()
    return v

def sys_Aerodynamics_total(v):    
    v.CD0clean = calcCD0(v.fuselagewetted,v.fuselagelength,v.fuselagefrontalarea,v.fuselagewidth,v.Sh,v.Sv,v.MAC_emp,v.BLturbratio_fus,v.BLturbratio_wing,v.BLturbratio_emp,v.h_landinggear,v.w_landinggear,v.dCD_landinggear,v.MAC,v.flapaffectedarea)
    v.CD0flaps = calcCD0(v.fuselagewetted,v.fuselagelength,v.fuselagefrontalarea,v.fuselagewidth,v.Sh,v.Sv,v.MAC_emp,v.BLturbratio_fus,v.BLturbratio_wing,v.BLturbratio_emp,v.h_landinggear,v.w_landinggear,v.dCD_landinggear,v.MAC,v.flapaffectedarea,clean_config = False)
    v.e_clean = calcOswald(v.fuselagewetted,v.fuselagelength,v.fuselagefrontalarea,v.fuselagewidth,v.Sh,v.Sv,v.MAC_emp,v.BLturbratio_fus,v.BLturbratio_wing,v.BLturbratio_emp,v.h_landinggear,v.w_landinggear,v.dCD_landinggear,v.MAC,v.flapaffectedarea,hasWinglets=v.haswinglets)
    v.e_flaps = calcOswald(v.fuselagewetted,v.fuselagelength,v.fuselagefrontalarea,v.fuselagewidth,v.Sh,v.Sv,v.MAC_emp,v.BLturbratio_fus,v.BLturbratio_wing,v.BLturbratio_emp,v.h_landinggear,v.w_landinggear,v.dCD_landinggear,v.MAC,v.flapaffectedarea,hasWinglets=v.haswinglets,clean_config=True)
    return v

