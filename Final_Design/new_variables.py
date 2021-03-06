import numpy as np
from math import sqrt,cos,pi, floor
import matplotlib.pyplot as plt
from numpy import linalg as la
try:
    from materials import materials
except ModuleNotFoundError:
    from Structures.materials import materials

class NewVariables:
    def __init__(self,haswinglets=True,wingletheight=0.43):
        self.init_general()
        self.init_aerodynamics(haswinglets,wingletheight)
        self.init_wing()
        self.init_fuselage()
        self.init_propulsion()
        self.init_sc()
        self.init_weight()

    def init_general(self):
        self.WTO = 750*9.81
        self.rho0 = 1.225
        self.g0 = 9.80665
        self.R = 287.05
        self.gamma = 1.4 
        self.lmbda = -0.0065
        self.T0 = 288.15

        self.Vcruise = 50           # [m/s]         Cruise speed
        self.V = self.Vcruise
        self.Vs = 23.15    
        self.hcruise = 914.4        # [m]           Cruise altitude
        self.rhocruise = 1.12
        self.rhotakeoff = 1.04
        self.sigma = self.rhotakeoff/self.rhocruise

        self._sto = 500
        self.k = np.sqrt(5647.9 + 17.331 * self.sto) - 75.153
        self.sland = 500

        self.designpointfactor = 1


        self.c = 2
        self.n_ult = 4.81 # TO BE OVERWRITTEN. ASK CASPER HOW TO IMPLEMENT THIS

        # Landing gear
        self.xng = None             # [m]   Nose gear location
        self.xmg = None             # [m]   Main gear location
        self.zmg = None             # [m]   Main gear vertical location
        self.propclear = 0.23
        self.h_landinggear = 0.8
        self.w_landinggear = 0.175

        self._xtail = 8.82

        self.fuselagelength = self.xtail + 0.42
        self.fuselagewidth = 1.05
        self.fuselageheight = 1.211
        self.fuselagefrontalarea = 1.218
        self.fuselagewettedarea = 17.507

        self.proplength = 0.3       # [m] Length of the propeller part, in front of the fuselage

        self.la = self.fuselagelength + self.proplength

        self.h_htail = self.h_landinggear + .5*self.fuselageheight

        self.cg_wing = 0.8         # [m]           Distance LE root chord - wing cg
        self.xwing = None

        self.xcg_fgroup = 2.3

        self.xcgPL = 1.5
        self.xcg_min = None
        self.xcg_max = None

        self.zcg = self.h_landinggear+0.5*self.fuselageheight    # [m]           Centre of gravity height

        self._WP = 0.121
        self._WS = 592

        self.printing = False # Whether to print the structures calculations results or not

    @property
    def xtail(self):
        return self._xtail

    @xtail.setter
    def xtail(self, val):
        self._xtail = val
        self.fuselagelength = val+0.42
        self.framelocs = np.array([self._cockpitbulkhead + (self.fuselagelength-self.cockpitbulkhead)*n/(self.framesamount+1) for n in range(self.framesamount+1)])[::-1]
        self.la = self.fuselagelength + self.proplength

    @property
    def sto(self):
        return self._sto

    @sto.setter
    def sto(self, val):
        self._sto = val
        self.k = np.sqrt(5647.9 + 17.331 * self._sto) - 75.153

    @property
    def WP(self):
        return self._WP

    @WP.setter
    def WP(self, val):
        self._WP = val
        # P1 = 65 * 1000      # Maximum power produced by the engine in W
        # P2 = self.WTO/val
        # self.P_max = max(P1, P2)

    @property
    def WS(self):
        return self._WS

    @WS.setter
    def WS(self, val):
        self._WS = val
        self.S = self.WTO/val

    def init_aerodynamics(self,haswinglets,wingletheight):
        # Conditions
        self.clean_config = None
        self.hasWinglets = haswinglets
        
        # General wing geometry
        self._S = 16.5
        self._A = 10.1
        self._taper = 0.45
        self._twist = np.radians(5)
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

        # Aileron geometry
        self.aileronstart = 4.9
        self.aileronend = None
        self.aileronspan = None
        
        # Horizontal tail geometry
        self._Sh = 0.4*self.S
        self.A_h = 3
        self.taper_h = 0.7
        self.sweeph = np.radians(10)
        self._b_h = np.sqrt(self.A_h * self.Sh)
        self._c_r_h = (2 * self.Sh) / (self.b_h * (1 + self.taper_h))
        self.c_t_h = self.taper_h * self.c_r_h
        self.MAC_h = (2 / 3) * self.c_r_h * ((1 + self.taper_h + self.taper_h ** 2) / (1 + self.taper_h))
        self._sweepLE_h = np.arctan( 4/self.A_h * (1-self.taper_h)/(1+self.taper_h) )
        self._YMAC_h = self.b / 6 * (1 + 2 * self.taper_h) / (1 + self.taper_h)
        self.XMAC_h = self.YMAC_h * np.tan(self.sweepLE_h)        


        # Exterior design inputs
        self.BLturbratio_wing = 0.65
        self.BLturbratio_fus = 1
        self.BLturbratio_emp = 0.65
        
        # Requirements
        self.CL_landing = 2.0
        self.CL_takeoff = self.CL_landing/(1.1**2)
        self.CL_climb = 1.8

        # Constants & statistical variables
        self.kwl = 2.1
        self.dCD_landinggear = 0.15
        self.going_mad = True

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
        self.clear_tip = 0.2   # [m]       Distance from the tip that should be clear of control surfaces
        self.da = 20*pi/180                 # [rad]     Aileron deflection angle
        self.clda = 0.046825*180/pi         # [/rad]    Take-off configuration change in the airfoil’s lift coefficient with aileron deflection


        # Output variables
        self.wing_CL_alpha = None
        self.wing_CL_max = None
        self.wing_alpha_max = None
        self._CD0clean = 0.01912
        self._CD0flap = 0.04
        self.CD0to = self.CD0flap/(1.1**2)
        if haswinglets:
            self._eclean = 0.84
        if not haswinglets:
            self._eclean = 0.79
        self.eflaps = None

        self.CD_climb = self.CD0clean + self.CL_climb**2/(np.pi*self.A*self.eclean)

        # Working variables
        self.coeff = None

    @property
    def CD0clean(self):
        return self._CD0clean

    @CD0clean.setter
    def CD0clean(self, val):
        self._CD0clean = val
        self.CD0clean + self.CL_climb ** 2 / (np.pi * self.A * self.eclean)

    @property
    def eclean(self):
        return self._eclean

    @eclean.setter
    def eclean(self, val):
        self._eclean = val
        self.CD_climb = self.CD0clean + self.CL_climb ** 2 / (np.pi * self.A * self.eclean)

    def init_wing(self):
        self.w_n_stiff_u = 3
        self.w_n_stiff_l = 3
        self.w_n_ribs = 3
        self.spar_le_locs = 0.25
        self.spar_te_locs = 0.75
        self.n_discretizations = 5 # Don't make more than 10, contact Yann if you want to change this.

    def init_fuselage(self):
        self._max_fuselage_iterations = None # Don't change without contacting Max
        self._Wfus_aft = 50 # weight of the aft section of the fuselage
        self.Wfus_aft_xbar = None # centroid of the aft section of the fuselage
        self._Wfus_aft_ybar = None # centroid of the aft section of the fuselage
        self.Wfus_aft_zbar = None # centroid of the aft section of the fuselage
        self._cockpitbulkhead = 2.2+self.proplength # m, back of the cockpit
        self._framesamount = 8 # [-] amount of frames behind the cockpit bulkhead in the fuselage
        self.framelocs = np.array([self._cockpitbulkhead + (self.fuselagelength-self.cockpitbulkhead)*n/(self.framesamount+1) for n in range(self.framesamount+1)])[::-1]
        self.skin_t = 2.5 # 2.138574136 # 2.0 # [mm] thickness of the skin
        self.skin_t_func = None # Only change if a custom skin thickness function is required. Two arguments: First is skin_t, second is y position where the thickness should be taken.
        self.mats = materials()
        self.stringermod = 1.0 # 1.5 # 2.23937285795561 # 1.0 # one-dimensional scaling parameter for stringers, geometry stays the same.
        self.circstringermod = 1.0 # 1.5 # 2.62916210608877 # 1.0 # one-dimensional scaling parameter for stringers in circular area, geometry stays the same.
        self.longeronmod = 4. # 3.36724386380556 # 4.0 # one-dimensional scaling parameter for longerons, geometry stays the same.
        self.stringermat = self.mats['alu2024']
        self.circstringermat = self.mats['alu2024']
        self.longeronmat = self.mats['carbonfibre']
        self._n_stiff = 6 # Amount of stiffeners per panel
        self.n_stiff_circ = 12 # Amount of stringers total for circular section
        self.n_longs = 4 # Can't change, only to 0 if longerons should be disabled
        self.n_stringers = np.ones(8, dtype=int) * self.n_stiff
        self.n_stringers[0]=self.n_stiff//2
        self.n_stringers[4]=self.n_stiff//2
        self.material_regular = self.mats['alu7075']
        self.material_circular = self.mats['carbonfibre']
        self.batteryinfuselage = True
        self._batteryoffset = 0.1 # [m] battery distance behind the cockpit aft bulkhead
        self._batterywidth = 0.2 # [m] distance of second attachment point of the battery from the first
        self.yout = None
        self._fuselage = None # Fuselage object, don't change without contacting Max
        self.loads = None # Loads acting on the fuselage
        if self.batteryinfuselage:
            self.xcgbat = self._cockpitbulkhead + self._batteryoffset + self._batterywidth*0.5
        else:
            self.xcgbat = 0.7

        self.emergencyparachute = self.cockpitbulkhead + self.batteryoffset + self.batterywidth + 0.1

        self.nosegearfraction = 0.08
    @property
    def cockpitbulkhead(self):
        return self._cockpitbulkhead

    @cockpitbulkhead.setter
    def cockpitbulkhead(self, val):
        self._cockpitbulkhead = val
        if self.batteryinfuselage:
            self.xcgbat = val + self._batteryoffset + self._batterywidth*0.5


    @property
    def batteryoffset(self):
        return self._batteryoffset

    @property
    def batterywidth(self):
        return self._batterywidth

    @batteryoffset.setter
    def batteryoffset(self, val):
        self._batteryoffset = val
        if self.batteryinfuselage:
            self.xcgbat = self._cockpitbulkhead + val + self._batterywidth*0.5

    @batterywidth.setter
    def batterywidth(self, val):
        self._batterywidth = val
        if self.batteryinfuselage:
            self.xcgbat = self._cockpitbulkhead + self._batteryoffset + val*0.5


    def init_weight(self):        
        self.WPL = 1961.33
        
        self.W_wing    = None           # Wing weight

        self.W_batt    = None           # Battery weight in Newtons
        self.W_motor   = 30 * 9.81      # Motor weight in Newtons
        self.W_shaft   = 4.48 * 9.81    # Engine shaft weight in Newtons
        self.W_prop    = 12 * 9.81      # Propeller weight in Newtons
        
        self.W_syscomp = 69.2*9.81      # System component weight (TE package + avionics + electronics)

        # self.Wfus_fwd = max(1228.755-self._Wfus_aft, 500)
        self.Wfus_fwd = self._Wfus_aft*5/7
        self.W_fgroup = None

        self.W_htail = None
        self.W_vtail = None
        
        self.W_OEW = None

        self.xcg_fus_fwd = 1.3 # 1.811
        self.xcg_fus_aft = 4.1
        self.xprop = 0.15
        self.xmotor = 0.6
        self.xshaft = 0.525

    @property
    def Wfus_aft_ybar(self):
        return self._Wfus_aft_ybar

    @Wfus_aft_ybar.setter
    def Wfus_aft_ybar(self, val):
        self._Wfus_aft_ybar = val
        self.xcg_fus_aft = val

    @property
    def Wfus_aft(self):
        return self._Wfus_aft

    @Wfus_aft.setter
    def Wfus_aft(self, val):
        self._Wfus_aft = val
        # self.Wfus_fwd = max(1228.755-val, 500)
        self.Wfus_fwd = self._Wfus_aft*5/7

    @property
    def fuselage(self):
        return self._fuselage

    @fuselage.setter
    def fuselage(self, object):
        self._fuselage = object
        self.Wfus_aft = object.mass * 9.81
        self.Wfus_aft_xbar = object.xbar
        self.Wfus_aft_ybar = object.ybar
        self.Wfus_aft_zbar = object.zbar

    @property
    def n_stiff(self):
        return self._n_stiff

    @n_stiff.setter
    def n_stiff(self, val):
        self._n_stiff = val
        self.n_stringers = np.ones(8)*val
        self.n_stringers[0]=val//2
        self.n_stringers[4]=val//2

    @property
    def framesamount(self):
        return self._framesamount

    @framesamount.setter
    def framesamount(self, val):
        self._framessamount = val
        self.framelocs = np.array([2.2 + self.proplength + (self.fuselagelength-self.cockpitbulkhead)*n/(self.framesamount+1) for n in range(self.framesamount+1)])[::-1]

    def init_propulsion(self):
        # Sizing
        self.v_batt       = None           # Battery volume in liters
        self.P_max        = 65 * 1000      # Maximum power produced by the engine in W

        # Battery characteristics
        self.batt_E_prod  = None           # Energy produced by the battery in Wh
        self.batt_Ns      = None           # Number of battery cells in series
        self.batt_Np      = None           # Number of battery cells in parallel
        self.batt_N       = None           # Number of battery cells
        self.eff_batt     = 0.95           # Battery efficiency
        self._eff_prop    = 0.85           # Propeller efficiency
        self.eff_tot_prop = 0.95 * self.eff_prop    # Total propulsive efficiency
        self.V_req_batt   = 400            # Required voltage in Volts
        self.I_req_batt   = 0.0 #189.75         # Required current in Amps UPDATE THIS
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
        self._batt_cell_C_Ah       = 5.0                                                             # Capacity [Ah]
        self.batt_cell_C_Wh       = self.batt_cell_C_Ah * self.batt_cell_V_nom                      # Capacity [Wh]

        # WARNING!: do not change batt_cell_E_spec! -> Change batt_cell_Ah or batt_cell_V_nom
        self.batt_cell_E_spec     = self.batt_cell_C_Wh / self.batt_cell_mass                       # Capacity per kg of weight [Wh/kg]
        self.batt_cell_E_vol_spec = self.batt_cell_C_Wh / self.batt_cell_volume                     # Capacity per unit volume [Wh/m^3]
        self.batt_cell_P          = self.batt_cell_I_max * self.batt_cell_V_nom                     # Maximum power [W]

    @property
    def batt_cell_C_Ah(self):
        return self._batt_cell_C_Ah

    @batt_cell_C_Ah.setter
    def batt_cell_C_Ah(self, val):
        self._batt_cell_C_Ah = val
        self.batt_cell_C_Wh = val * self.batt_cell_V_nom


    @property
    def eff_prop(self):
        return self._eff_prop

    @eff_prop.setter
    def eff_prop(self, val):
        self._eff_prop = val
        self.eff_tot_prop = self.eff_batt * val

    def init_sc(self):
        # Fuselage variables
        self.Sfs = 5.258            # [m2]          Fuselage lateral surface area
        self.hf1 = 1.146            # [m]           Fuselage nose height
        self.hf2 = 0.306            # [m]           Fuselage tail height
        self.bf1 = 0.960            # [m]           Fuselage nose width
        self.bf2 = 0.243            # [m]           Fuselage tail width


        # Propeller geometry parameters
        self.Bp = 2.5               # [-]           Number of blades per propeller
        self.Dp1 = 1.8                # [m]           1st propeller disk diameter
        self.Dp2 = 1.8                # [m]           2nd propeller disk diameter


        # Wing variables
        self.xlemac = None          # [m]           Distance nose - leading edge mean aerodynamic chord

        self.VhV = sqrt(0.85)       # [-]           Tail/wing speed ratio
        self.eta = 0.95             # [-]           Airfoil efficiency coefficient
        self.Cm0af = -0.1           # [-]           Airfoil zero lift pitching moment coefficient
        self.mu1 = 0.3              # [-]           Flap coefficient 1
        self.deda = None            # [-]           Downwash gradient
        self.cc = 1                 # [-]           Chord ratio (extended flap/clean)
        self.CnBi = 0.024           # [-]           Wing configuration lateral stability component


        # Horizontal tail variables
        self.lh = 6.9               # [m]           Horizontal tail arm
        self.ih = 0                 # [rad]         Horizontal tail incidence angle
        
        self.CLh_L = -0.8           # [-]           Horizontal tail landing configuration lift coefficient
        self.CLh_TO = None          # [-]           Horizontal tail take-off configuration lift coefficient
        

        # Vertical tail variables
        self._Sv = 3               # [m2]          Vertical tail surface

        self.CnB = None             # [-]           Directional stability coefficient


        # Flight performance parameters
        self.a_pitch = 12*pi/180    # [rad/s]       Take-off pitch angular velocity
        self.VTO = 1.05 * 25.2      # [m/s]         Take-off velocity
        self.rhoTO = 1.04           # [kg/m3]       Take-off density
        self.VL = 1.1 * 25.2
        self.mu = 0.05              # [-]           Take-off friction factor

        # Component locations
        self.zT = self.h_landinggear + 0.61
                                    # [m]           Thrust vector height
        self.zD = self.h_landinggear + .5*self.fuselageheight
                                    # [m]           Drag vector height

        # Elevator geometry parameters
        self.bebh = 1               # [-]           Elevator span
        self.de_max = 20*pi/180     # [rad]         Maximum elevator deflection
        self.de_min = -25*pi/180    # [rad]         Minimum elevator deflection (maximum negative)
        
        # Horizontal tail aerodynamic parameters
        self.CLh_TO = None          # [-]           Horizontal tail take-off lift coefficient


    @property
    def Sh(self):
        return self._Sh

    @Sh.setter
    def Sh(self,val):
        self._Sh = val
        self.b_h = np.sqrt(self.A_h * self.Sh)
        self.c_r_h = (2 * self.Sh) / (self.b_h * (1 + self.taper_h))

    @property
    def b_h(self):
        return self._b_h

    @b_h.setter
    def b_h(self, val):
        self._b_h = val
        self.c_r_h = (2 * self.Sh) / (self.b_h * (1 + self.taper_h))
        self.YMAC_h = self.b_h / 6 * (1 + 2 * self.taper_h) / (1 + self.taper_h)
    
    @property
    def c_r_h(self):
        return self._c_r_h

    @c_r_h.setter
    def c_r_h(self, val):
        self._c_r_h = val
        self.c_t_h = self.taper_h * self.c_r_h
        self.sweepLE_h = np.arctan(-self.c_r_h / (2 * self.b_h) * (self.taper_h - 1))

    @property
    def YMAC_h(self):
        return self._YMAC_h

    @YMAC_h.setter
    def YMAC_h(self,val):
        self._YMAC_h = val
        self.XMAC_h = self.YMAC_h * np.tan(self.sweepLE_h) 

    @property
    def sweepLE_h(self):
        return self._sweepLE_h

    @sweepLE_h.setter
    def sweepLE_h(self,val):
        self._sweepLE_h = val
        self.XMAC_h = self.YMAC_h * np.tan(self.sweepLE_h) 
    
    @property
    def Sv(self):
        return self._Sv

    @Sv.setter
    def Sv(self,val):
        self._Sv = val

    @property
    def S(self):
        return self._S

    @S.setter
    def S(self, val):
        self._S = val
        self.b = np.sqrt(self.A * self.S)
        self.Snet = self.S - self.calculateChord(self.transformSpan(.25*self.fuselagewidth,self.b),self.taper,self.S,self.b)*self.fuselagewidth
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
        self.Snet = self.S - self.calculateChord(self.transformSpan(.25*self.fuselagewidth,self.b),self.taper,self.S,self.b)*self.fuselagewidth
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
        self.Snet = self.S - self.calculateChord(self.transformSpan(.25*self.fuselagewidth,self.b),self.taper,self.S,self.b)*self.fuselagewidth
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

    @property
    def CD0flap(self):
        return self._CD0flap

    @CD0flap.setter
    def CD0flap(self,val):
        self._CD0flap = val
        self.CD0to = self.CD0flap/(1.1**2)


    ###################################

    def update_WTO(self):
        self.WTO = self.W_OEW + self.WPL + self.W_batt
        self.S = self.WTO/self.WS
        # P1 = 65 * 1000      # Maximum power produced by the engine in W
        # P2 = self.WTO/self.WP
        # self.P_max = max(P1, P2)

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

    def calcsweep(self,pc,b,sweepc4,taper,cr):
        sweepLE = np.arctan(np.tan(sweepc4) - cr/(2*b)*(taper-1))
        return np.arctan(np.tan(sweepLE) + 2*pc*cr*(taper-1)/b)

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
        self.wing_CL_alpha= np.pi * self.A * self.coeff[0][0]

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

    def calcCD0(self,S_wet_fus,l_fus,fus_A_max,w_fuselage,S_h,S_v,MAC_emp,BLturbratio_fus, BLturbratio_wing, BLturbratio_emp,l_gear,w_gear,dCD_gear,MAC,flap_area_ratio,tc_airfoil=0.15,xc_airfoil=0.3,tc_emp=0.12,xc_emp=0.3,V_stall=23.15,rho_cruise=1.12,clean_config=True,visc=1.8e-5):
        
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

    def calcOswald(self,S_wet_fus,l_fus,fus_A_max,w_fuselage,S_h,S_v,MAC_emp,BLturbratio_fus, BLturbratio_wing, BLturbratio_emp, l_gear,w_gear,dCD_gear,MAC,flap_area_ratio,tc_airfoil=0.15,xc_airfoil=0.3,tc_emp=0.12,xc_emp=0.3,V_stall=23.15,rho_cruise=1.12,clean_config=True,visc=1.8e-5,hasWinglets=False):
        k_fuselage = 1-2*(w_fuselage/self.b)**2
        Q = 1/(self.calcespan()*k_fuselage)
        P = 0.38*self.calcCD0(S_wet_fus,l_fus,fus_A_max,w_fuselage,S_h,S_v,MAC_emp,BLturbratio_fus, BLturbratio_wing, BLturbratio_emp,l_gear,w_gear,dCD_gear,flap_area_ratio,tc_airfoil,xc_airfoil,MAC,tc_emp,xc_emp,V_stall,rho_cruise,clean_config,visc)
        k_winglet = (1+2*self.hwl/(self.kwl*self.b))**2
        
        if not self.hasWinglets:
            return 1/(Q+P*np.pi*self.A)

        else:
            return k_winglet/(Q+P*np.pi*self.A)

    def flap_sizing(self):
        # Inputs
        S = self.S                      # [m2]      Wing surface area
        b = self.b                      # [m]       Wing span
        sweepc4 = 0                     # [rad]     Wing quarter chord sweep angle
        taper = self.taper              # [rad]     Wing taper ratio
        cr = self.c_r                   # [m]       Wing root chord

        #CLmax_req = self.CLmax_req      # [-]       Required maximum lift coefficient
        #CLmax_wing = self.calcCLmax()   # [-]       Wing maximum lift coefficient
        #CLa = self.calcCLa()            # [/rad]    Wing lift curve slope

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
        cf1 = self.calculateChord(self.transformSpan(f1,self.b), self.taper, self.S, self.b)

        # Leading edge sweep angle
        sweepLE = self.calcsweep(0, b, sweepc4, taper, cr)

        # Hinge line sweep angle
        sweep_hinge = self.calcsweep(cfc, b, sweepc4, taper, cr)

        # Trailing edge sweep angles
        sweepTE = self.calcsweep(1, b, sweepc4, taper, cr)

        # Increase in lift coefficient
        dCLmax = (2.0 - self.wing_CL_max) * (1 + sm)

        # Required flapped surface
        SwfS = dCLmax / (0.9 * dClmax * cos(sweep_hinge))
        Swf = SwfS * S

        # Shift in zero lift angle of attack
        da0L = da0l_airfoil * SwfS * cos(sweep_hinge)

        # Change in lift curve slope
        CLa_flapped = self.wing_CL_alpha

        # Flap span calculation
        # Solving the equation:
        # cf1*bfl - 0.5*bf^2*tan(sweepLE) + 0.5*bf^2*tan(sweepTE) = Swf/2
        # A*bfl^2 + B*bfl + C = 0
        A = 0.5 * (-np.tan(sweepLE) + np.tan(sweepTE))
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
        cla = self.wing_CL_alpha    # [/rad]    Take-off configuration lift curve slope
        cd0_TO = self.CD0clean      # [-]       Take-off configuration zero lift drag coefficient
        cd0_L = self.CD0flap        # [-]       Landing configuration zero lift drag coefficient

        b1 = self.aileronstart      # [m]       Aileron start
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
        self.aileronstart = b1
        self.aileronend = b2
        self.aileronspan = b2 - b1


def sys_Aerodynamics_wing(v,resolution):
    v.setAirfoils(v.Clmax_r, v.Clmax_t, v.Cla_r, v.Cla_t, v.a0_r, v.a0_t)
    v.calcCoefficients(resolution,0.7)
    v.calcCLa()
    v.calcCLmax()
    v.CL0clean= v.calcCL(0)
    v.flap_sizing()
    v.aileron_sizing()
    v.CL0flap = v.CL0clean + (2.0 - v.wing_CL_max)
    return v

def sys_Aerodynamics_total(v):    
    v.CD0clean = v.calcCD0(v.fuselagewettedarea,v.fuselagelength,v.fuselagefrontalarea,v.fuselagewidth,v.Sh,v.Sv,v.MAC_h,v.BLturbratio_fus,v.BLturbratio_wing,v.BLturbratio_emp,v.h_landinggear,v.w_landinggear,v.dCD_landinggear,v.MAC,v.flapaffectedarea)
    v.CD0flaps = v.calcCD0(v.fuselagewettedarea,v.fuselagelength,v.fuselagefrontalarea,v.fuselagewidth,v.Sh,v.Sv,v.MAC_h,v.BLturbratio_fus,v.BLturbratio_wing,v.BLturbratio_emp,v.h_landinggear,v.w_landinggear,v.dCD_landinggear,v.MAC,v.flapaffectedarea,clean_config = False)
    v.eclean = v.calcOswald(v.fuselagewettedarea,v.fuselagelength,v.fuselagefrontalarea,v.fuselagewidth,v.Sh,v.Sv,v.MAC_h,v.BLturbratio_fus,v.BLturbratio_wing,v.BLturbratio_emp,v.h_landinggear,v.w_landinggear,v.dCD_landinggear,v.MAC,v.flapaffectedarea,hasWinglets=v.hasWinglets)
    v.eflaps = v.calcOswald(v.fuselagewettedarea,v.fuselagelength,v.fuselagefrontalarea,v.fuselagewidth,v.Sh,v.Sv,v.MAC_h,v.BLturbratio_fus,v.BLturbratio_wing,v.BLturbratio_emp,v.h_landinggear,v.w_landinggear,v.dCD_landinggear,v.MAC,v.flapaffectedarea,hasWinglets=v.hasWinglets,clean_config=False)
    return v



def calcFusgroup(v):
    Wlist = np.array([v.Wfus_fwd,v.Wfus_aft,v.W_htail,v.W_vtail,v.W_prop,v.W_shaft,v.W_motor,v.W_syscomp])
    Xlist = np.array([v.xcg_fus_fwd,v.xcg_fus_aft,v.xtail,v.xtail,v.xprop,v.xshaft,v.xmotor,v.xcg_fus_fwd])

    Wlistreduced = np.array([v.Wfus_fwd,v.Wfus_aft,v.W_shaft,v.W_motor,v.W_syscomp])
    Xlistreduced = np.array([v.xcg_fus_fwd,v.xcg_fus_aft,v.xshaft,v.xmotor,v.xcg_fus_fwd])

    v.xcg_fuselage = np.sum(Wlistreduced*Xlistreduced)/np.sum(Wlistreduced)

    v.xcg_fgroup = np.sum(Wlist*Xlist)/np.sum(Wlist)
    v.W_fgroup = np.sum(Wlist)
    return v


def CalcOEW(v):
    v.W_OEW = v.W_wing + v.W_fgroup #v.W_fuselage + v.W_motor + v.W_shaft + v.W_prop + v.W_syscomp + v.W_htail + v.W_vtail
    return v

def CalcMTOWnew(v):
    v.WTO = v.W_OEW + v.WPL + v.W_batt
    return v

def CalcTTO(v):
    s           = 250                           #[m] assumed take-off roll (NOT 500!!!! Need to climb to 50 ft!)
    V_lof       = 45*1.05*1.851/3.6      #[m/s] lift of speed
    rho_sea     = v.rho0                        #[kg/m3] air density at sea level
    S_wing      = v.S                           #[m2] wing surface
    W           = v.WTO
    CD_to       = 0.0414                        #[-] DO NOT CHANGE THIS VARIABLE!!!
    CL_to       = 0.628                         #[-] DO NOT CHANGE THIS VARIABLE!!!
    mu          = 0.05                          #[-] DO NOT CHANGE THIS VARIABLE!!!

    a_avg = V_lof*2/(2*s)  #[m/s*2] average acceleration
    
    V_avg = V_lof/np.sqrt(2) #[m/s] average velocity
    
    D_avg = 0.5*rho_sea*V_avg**2*S_wing*CD_to #[N] average Drag
    
    L_avg = 0.5*rho_sea*V_avg**2*S_wing*CL_to #[N] average Lift
    
    
    v.TTO = a_avg*W/9.80665 + mu*(W-L_avg) + D_avg #[N] average Take off thrust
    
    return v
     #[N]

if __name__ == "__main__":
    v = NewVariables(False,0)
    sys_Aerodynamics_wing(v,10)
    sys_Aerodynamics_total(v)
    print(v.eclean)
