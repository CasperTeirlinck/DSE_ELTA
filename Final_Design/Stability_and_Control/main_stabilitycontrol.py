from math import sqrt,pi
from horizontaltail_design import sizing_htail_wingpos

class stabilitycontrol_variables():
    def __init__(self):
        # Constants
        self.R = 287.05         # [J/kg K]  Gas constant
        self.gamma = 1.4        # [-]       Heat capacity ratio
        self.T0 = 288.15        # [K]       Base temperature
        self.lmbda = -0.0065    # [degC/m]  Lapse rate

        # Flight performance parameters
        self.Vcruise = 50       # [m/s]     Cruise speed
        self.hcruise = 914.4    # [m]       Cruise altitude

        # Safty margins
        #sm = 0.05              # [-]       Safety Margin
        #sm_free = 0.05         # [-]       Fraction neutral point shift for stick-free stability
        #sm = 0.05              # [-]       Stability margin
        #sm = 0.1               # [-]       Safety margin

        # Weight parameters
        self.m_fgroup = 100     # [kg]      Fuselage group mass
        self.m_wing = 150       # [kg]      Wing mass
        self.m_bat = 300        # [kg]      Battery mass
        self.m_repbat = 100     # [kg]      Replaceable battery mass
        self.m_pax = 100        # [kg]      Single Pilot/Passenger mass

        # Centre of gravity parameters
        self.xcg_min = None     # [m]       Minimum center of gravity location
        self.xcg_max = None     # [m]       Maximum center of gravity location
        self.xcg_fgroup = 3     # [m]       Fuselage group center of gravity
        self.xcg_f = 2.3        # [m]       Fuselage centre of gravity location
        self.xcg_wing = 1.7     # [m]       Wing center of gravity
        self.xcg_bat = 3        # [m]       Battery center of gravity
        self.xcg_repbat = 3     # [m]       Replaceable battery center of gravity
        self.xcg_pax = 1        # [m]       Pilot/Passenger center of gravity

        # Fuselage geometry parameters
        self.lf = 9.420         # [m]       Fuselage length
        self.bf = 1.6           # [m]       Fuselage width
        self.hf = 1.213         # [m]       Fuselage height
        self.hfmax = 1.213      # [m]       Maximum fuselage height
        self.Sfs = 5.258        # [m2]      Fuselage lateral surface area
        self.hf1 = 1.146        # [m]       Fuselage nose height
        self.hf2 = 0.306        # [m]       Fuselage tail height
        self.bf1 = 0.960        # [m]       Fuselage nose width
        self.bf2 = 0.243        # [m]       Fuselage tail width

        # Wing geometry parameters
        self.lfn = 1.5          # [m]       Distance nose - wing
        self.hw = 0.5           # [m]       Height of the wing, from ground
        self.Sw = 15.6          # [m2]      Horizontal tail surface area
        self.Snet = 10          # [m2]      Net wing surface area
        self.bw = 12.6          # [m]       Wing span
        self.Aw = 10.1          # [-]       Wing aspect ratio
        self.sweepw = 0         # [rad]     Wing quarter chord sweep angle
        self.taperw = 0.45      # [-]       Wing taper ratio
        self.twistwr = 5*pi/180 # [rad]     Wing twist at the root
        self.MAC = 1.30         # [m]       Mean Aerodynamic Chord
        self.crw = 1.98         # [m]       Wing root chord

        # Horizontal tail geometry parameters
        self.lh = 6.2           # [m]       Horizontal tail arm
        self.Sh_min = None      # [m2]      Minimum required horizontal tail surface
        self.hh = 1.5           # [m]       Height horizontal tail from ground
        self.Ah = 3             # [-]       Horizontal tail aspect ratio
        self.sweeph = 0         # [rad]     Horizontal tail half chord sweep

        # Vertical tail geometry parameters
        self.lv = 9             # [m]       Vertical tail arm

        # Propeller geometry parameters
        self.Bp = 3             # [-]       Number of blades per porpeller
        self.lp1 = 2            # [m]       Distance 1st propeller plane - aircraft centre of gravity
        self.lp2 = 2            # [m]       Distance 1st propeller plane - aircraft centre of gravity
        self.Dp1 = 2            # [m2]      1st propeller disk diameter
        self.Dp2 = 2            # [m2]      1st propeller disk diameter

        # Wing aerodynamic parameters
        self.eta = 0.95         # [-]       Airfoil efficiency coefficient TODO Check if aerodynamics guys determine this value
        self.CLaw = 5.4         # [/rad]    Wing lift rate coefficient TODO Check if aerodynamics guys determine this value
        self.Cm0af = 0.05       # [-]       Airfoil zero lift pitching moment coefficient TODO Check this value
        self.mu1 = 1            # [-]       Flap coefficient 1 TODO Check this value
        self.mu2 = 1            # [-]       Flap coefficient 2 TODO Check this value
        self.mu3 = 1            # [-]       Flap coefficient 3 TODO Check this value
        self.dClmax = 1         # [-]       Airfoil lift coefficient increase at landing TODO Check this value
        self.cc = 1             # [-]       Chord ratio (extended flap/clean) TODO Check this value
        self.CL_landing = 2     # [-]       Wing lift coefficient at landing (all flaps deployed) TODO Check this value
        self.Swf = 10           # [m2]      Reference wing flapped surface area TODO Check this value
        self.CL0 = 1            # [-]       Flapped wing lift coefficient at zero angle of attack TODO Check this value
        self.CLA_h = 3          # [-]       Aircraft less tail lift coefficient TODO Check this value
        self.VhV = sqrt(0.85)   # [-]       Tail/wing speed ratio
        self.CnBi = 0.024       # [-]       Wing configuration stability component

        # Horizontal tail aerodynamic parameters
        self.CLh = -0.8         # [-]       Horizontal tail lift coefficient TODO Check this value

        # Vertical tail aerodynamic parameters


sc_v = stabilitycontrol_variables()

# Horizontal tail sizing and wing position
sc_v = sizing_htail_wingpos(sc_v,plot=True)