from math import sqrt,pi
from horizontaltail_design import sizing_htail_wingpos
from verticaltail_design import verticaltail_sizing
from empennagecontrol_design import elevator_sizing

class stabilitycontrol_variables:
    def __init__(self):
        # Constants
        self.R = 287.05             # [J/kg K]  Gas constant
        self.gamma = 1.4            # [-]       Heat capacity ratio
        self.T0 = 288.15            # [K]       Base temperature
        self.lmbda = -0.0065        # [degC/m]  Lapse rate
        self.g = 9.80665            # [m/s2]    Gravitational acceleration

        # Flight performance parameters
        self.Vcruise = 50           # [m/s]         Cruise speed
        self.hcruise = 914.4        # [m]           Cruise altitude
        self.a_pitch = 12*pi/180    # [rad/s]       Take-off pitch angular velocity
        self.VTO = 1.05 * 25.2      # [m/s]         Take-off velocity
        self.rhoTO = 1.225          # [kg/m3]       Take-off density
        self.mu = 0.05              # [-]           Take-off friction factor

        # Weight parameters
        self.WTO = 7354.9875        # [N]       Take-off weight
        self.m_fgroup = 100         # [kg]      Fuselage group mass
        self.m_wing = 150           # [kg]      Wing mass
        self.m_bat = 300            # [kg]      Battery mass
        self.m_repbat = 100         # [kg]      Replaceable battery mass
        self.m_pax = 100            # [kg]      Single Pilot/Passenger mass

        # Centre of gravity parameters
        self.xcg_min = None         # [m]       Minimum center of gravity location
        self.xcg_max = None         # [m]       Maximum center of gravity location
        self.xcg_fgroup = 3         # [m]       Fuselage group center of gravity
        self.xcg_f = 2.3            # [m]       Fuselage centre of gravity location
        self.xcg_wing = 1.7         # [m]       Wing center of gravity
        self.xcg_bat = 3            # [m]       Battery center of gravity
        self.xcg_repbat = 3         # [m]       Replaceable battery center of gravity
        self.xcg_pax = 1            # [m]       Pilot/Passenger center of gravity
        self.xmg = 2                # [m]       Main gear location
        self.xacw = 1               # [m]       Wing/fuselage aerodynamic centre location
        self.xach = 6               # [m]       Horizontal tail aerodynamic centre location

        self.zcg = 1                # [m]       Centre of gravity height
        self.zmg = 0                # [m]       Main gear height
        self.zT = 1                 # [m]       Thrust vector height
        self.zD = 0.5               # [m]       Drag vector height

        # Fuselage geometry parameters
        self.la = 9                 # [m]       Aircraft length
        self.lf = 9.420             # [m]       Fuselage length
        self.bf = 1.6               # [m]       Fuselage width
        self.hf = 1.213             # [m]       Fuselage height
        self.hfmax = 1.213          # [m]       Maximum fuselage height
        self.Sfs = 5.258            # [m2]      Fuselage lateral surface area
        self.hf1 = 1.146            # [m]       Fuselage nose height
        self.hf2 = 0.306            # [m]       Fuselage tail height
        self.bf1 = 0.960            # [m]       Fuselage nose width
        self.bf2 = 0.243            # [m]       Fuselage tail width

        # Wing geometry parameters
        self.lfn = 1.5              # [m]       Distance nose - wing
        self.xlemacw = None         # [m]       Distance nose - leading edge of the MAC
        self.hw = 0.5               # [m]       Height of the wing, from ground
        self.Sw = 15.6              # [m2]      Wing surface area
        self.Snet = 10              # [m2]      Net wing surface area
        self.bw = 12.6              # [m]       Wing span
        self.Aw = 10.1              # [-]       Wing aspect ratio
        self.sweepw = 0             # [rad]     Wing quarter chord sweep angle
        self.taperw = 0.45          # [-]       Wing taper ratio
        self.twistwr = 5*pi/180     # [rad]     Wing twist at the root
        self.MAC = 1.30             # [m]       Mean Aerodynamic Chord
        self.crw = 1.98             # [m]       Wing root chord
        self.f1 = 0.8               # [m]       Flap start
        self.f2 = 4.42              # [m]       Flap end

        # Horizontal tail geometry parameters
        self.lh = 6.2               # [m]       Horizontal tail arm
        self.Sh = None              # [m2]      Minimum required horizontal tail surface
        self.bh = 5                 # [m]       Horizontal tail span
        self.hh = 1.5               # [m]       Height horizontal tail from ground
        self.Ah = 3                 # [-]       Horizontal tail aspect ratio
        self.sweeph = 0             # [rad]     Horizontal tail half chord sweep
        self.chr = 1                # [m]       Horizontal tail root chord
        self.ih = 0                 # [rad]     Horizontal tail incidence angle

        # Vertical tail geometry parameters
        self.lv = 6.2               # [m]       Vertical tail arm
        self.Sv = None              # [m2]      Vertical tail surface

        # Elevator geometry parameters
        self.bebh = 1               # [-]       Elevator span
        self.de_max = 20*pi/180     # [rad]     Maximum elevator deflection
        self.de_min = -25*pi/180    # [rad]     Minimum elevator deflection (maximum negative)

        # Propeller geometry parameters
        self.Bp = 2.5               # [-]       Number of blades per porpeller
        self.lp1 = 2                # [m]       Distance 1st propeller plane - aircraft centre of gravity
        self.lp2 = 2                # [m]       Distance 1st propeller plane - aircraft centre of gravity
        self.Dp1 = 3.14             # [m2]      1st propeller disk diameter
        self.Dp2 = 3.14             # [m2]      1st propeller disk diameter

        # Wing aerodynamic parameters
        self.eta = 0.95             # [-]       Airfoil efficiency coefficient
        self.CLaw = 5.4             # [/rad]    Wing lift rate coefficient
        self.Cm0af = 0.05           # [-]       Airfoil zero lift pitching moment coefficient TODO Check this value
        self.mu1 = 1                # [-]       Flap coefficient 1 TODO Check this value
        self.mu2 = 0.25             # [-]       Flap coefficient 2
        self.mu3 = 1                # [-]       Flap coefficient 3 TODO Check this value
        self.dClmax = 1             # [-]       Airfoil lift coefficient increase at landing TODO Check this value
        self.cc = 1                 # [-]       Chord ratio (extended flap/clean)
        self.CL_landing = 2         # [-]       Wing lift coefficient at landing (all flaps deployed)
        self.Swf = 10               # [m2]      Reference wing flapped surface area TODO Check this value
        self.CL0 = 1                # [-]       Flapped wing lift coefficient at zero angle of attack TODO Check this value
        self.CLA_h = 3              # [-]       Aircraft less tail lift coefficient TODO Check this value
        self.VhV = sqrt(0.85)       # [-]       Tail/wing speed ratio
        self.CnBi = 0.024           # [-]       Wing configuration stability component
        self.CLTO = 0.628           # [-]       Take-off lift coefficient
        self.CDTO = 0.0414          # [-]       Take-off drag coefficient
        self.Cmacwf = -0.565        # [-]       Wing-fuselage pitching moment coefficient around the aerodynamic centre
        self.deda = None            # [-]       Downwash gradient
        self.a0 = -7.255*pi/180     # [rad]     Zero lift angle of attack

        # Horizontal tail aerodynamic parameters
        self.CLh_L = -0.8           # [-]       Horizontal tail landing lift coefficient
        self.CLh_TO = None          # [-]       Horizontal tail take-off lift coefficient
        self.CLah = 4               # [/rad]    Horizontal lift curve slope

        # Vertical tail aerodynamic parameters
        self.CnB = None             # [-]       Directional stability coefficient

        # Power parameters
        self.TTO = 1500             # [N]       Take-off thrust


def stability_control(variables):
    variables = sizing_htail_wingpos(variables)
    variables = verticaltail_sizing(variables)
    variables = elevator_sizing(variables)
    return variables

# Test
if __name__ ==  "__main__":
    sc_v = stabilitycontrol_variables()

    sc_v = stability_control(sc_v)