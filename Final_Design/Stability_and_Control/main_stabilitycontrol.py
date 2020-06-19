from math import sqrt,pi
from horizontaltail_design import sizing_htail_wingpos
from verticaltail_design import verticaltail_sizing
from empennagecontrol_design import elevator_sizing


class stabilitycontrol_variables:
    def __init__(self):
        # Constants
        self.R = 287.05                 # [J/kg K]  Gas constant
        self.gamma = 1.4                # [-]       Heat capacity ratio
        self.T0 = 288.15                # [K]       Base temperature
        self.lmbda = -0.0065            # [degC/m]  Lapse rate
        self.g0 = 9.80665               # [m/s2]    Gravitational acceleration

        # Flight performance parameters
        self.Vcruise = 50               # [m/s]         Cruise speed
        self.hcruise = 914.4            # [m]           Cruise altitude
        self.a_pitch = 12*pi/180        # [rad/s]       Take-off pitch angular velocity
        self.VTO = 1.05 * 25.2          # [m/s]         Take-off velocity
        self.rhoTO = 1.225              # [kg/m3]       Take-off density
        self.mu = 0.05                  # [-]           Take-off friction factor

        # Weight parameters
        self.WTO = 7354.9875            # [N]       Take-off weight
        self.W_fgroup = 100*9.80665     # [kg]      Fuselage group weight
        self.W_wing = 150*9.80665       # [kg]      Wing weight
        self.W_batt = 300*9.80665       # [kg]      Battery weight
        self.WPL = 1961.33              # [kg]      Payload weight

        # Centre of gravity parameters
        self.xcg_min = None             # [m]       Minimum center of gravity location
        self.xcg_max = None             # [m]       Maximum center of gravity location
        self.xcg_fgroup = 3             # [m]       Fuselage group center of gravity
        self.xcg_f = 2.3                # [m]       Fuselage centre of gravity location
        self.cg_wing = 1                # [m]       Wing center of gravity
        self.xwing = None
        self.xcgbat = 4                 # [m]       Battery center of gravity
        self.xcgPL = 1                  # [m]       Payload center of gravity
        self.xmg = 2                    # [m]       Main gear location
        self.xacw = 1                   # [m]       Wing/fuselage aerodynamic centre location
        self.xach = 6                   # [m]       Horizontal tail aerodynamic centre location

        self.zcg = 1                    # [m]       Centre of gravity height
        self.zmg = 0                    # [m]       Main gear height
        self.zT = 1                     # [m]       Thrust vector height
        self.zD = 0.5                   # [m]       Drag vector height

        # Fuselage geometry parameters
        self.la = 9                     # [m]       Aircraft length
        self.fuselagelength = 9.420     # [m]       Fuselage length
        self.fuselagewidth = 1.6        # [m]       Fuselage width
        self.fuselageheight = 1.213     # [m]       Fuselage height
        self.hfmax = 1.213              # [m]       Maximum fuselage height
        self.Sfs = 5.258                # [m2]      Fuselage lateral surface area
        self.hf1 = 1.146                # [m]       Fuselage nose height
        self.hf2 = 0.306                # [m]       Fuselage tail height
        self.bf1 = 0.960                # [m]       Fuselage nose width
        self.bf2 = 0.243                # [m]       Fuselage tail width

        # Wing geometry parameters
        self.xlemac = None              # [m]       Distance nose - leading edge of the MAC
        self.h_landinggear = 0.5        # [m]       Height of the wing, from ground
        self.S = 15.6                   # [m2]      Wing surface area
        self.Snet = 10                  # [m2]      Net wing surface area
        self.b = 12.6                   # [m]       Wing span
        self.A = 10.1                   # [-]       Wing aspect ratio
        self.sweep = 0                  # [rad]     Wing quarter chord sweep angle
        self.taper = 0.45               # [-]       Wing taper ratio
        self.twist = 5*pi/180           # [rad]     Wing twist at the root
        self.MAC = 1.30                 # [m]       Mean Aerodynamic Chord
        self.c_r = 1.98                 # [m]       Wing root chord
        self.flapspan = 4.42-0.8        # [m]       Flap span

        # Horizontal tail geometry parameters
        self.lh = 6.2                   # [m]       Horizontal tail arm
        self.Sh = None                  # [m2]      Minimum required horizontal tail surface
        self.bh = 5                     # [m]       Horizontal tail span
        self.h_htail = 1.5              # [m]       Height horizontal tail from ground
        self.Ah = 3                     # [-]       Horizontal tail aspect ratio
        self.sweeph = 0                 # [rad]     Horizontal tail half chord sweep
        self.cr_h = 1                   # [m]       Horizontal tail root chord
        self.ih = 0                     # [rad]     Horizontal tail incidence angle
        self.taperh = 0.4               # [-]       Horizontal tail taper ratio

        # Vertical tail geometry parameters
        self.lv = 6.2                   # [m]       Vertical tail arm
        self.Sv = None                  # [m2]      Vertical tail surface

        # Elevator geometry parameters
        self.bebh = 1                   # [-]       Elevator span
        self.de_max = 20*pi/180         # [rad]     Maximum elevator deflection
        self.de_min = -25*pi/180        # [rad]     Minimum elevator deflection (maximum negative)

        # Propeller geometry parameters
        self.Bp = 2.5                   # [-]       Number of blades per propeller
        self.lp1 = 2                    # [m]       Distance 1st propeller plane - aircraft centre of gravity
        self.lp2 = 2                    # [m]       Distance 1st propeller plane - aircraft centre of gravity
        self.Dp1 = 2                    # [m2]      1st propeller disk diameter
        self.Dp2 = 2                    # [m2]      1st propeller disk diameter

        # Wing aerodynamic parameters
        self.eta = 0.95                 # [-]       Airfoil efficiency coefficient
        self.Cm0af = -0.1               # [-]       Airfoil zero lift pitching moment coefficient TODO Check this value
        self.mu1 = 0.3                  # [-]       Flap coefficient 1
        self.dClmax = 1.25 * 0.96       # [-]       Airfoil lift coefficient increase at landing TODO Check this value
        self.cc = 1                     # [-]       Chord ratio (extended flap/clean)
        self.CL_landing = 2             # [-]       Wing lift coefficient at landing (all flaps deployed)
        self.flapaffectarea = 7         # [m2]      Reference wing flapped surface area TODO Check this value
        self.CL0 = 1                    # [-]       Flapped wing lift coefficient at zero angle of attack TODO Check this value
        self.CLA_h = 2                  # [-]       Aircraft less tail lift coefficient TODO Check this value
        self.VhV = sqrt(0.85)           # [-]       Tail/wing speed ratio
        self.CnBi = 0.024               # [-]       Wing configuration stability component
        self.CLTO = 0.628               # [-]       Take-off lift coefficient
        self.CDTO = 0.0414              # [-]       Take-off drag coefficient
        self.Cmacwf = -0.565            # [-]       Wing-fuselage pitching moment coefficient around the aerodynamic centre
        self.deda = None                # [-]       Downwash gradient
        self.a0 = -7.255*pi/180         # [rad]     Zero lift angle of attack

        # Horizontal tail aerodynamic parameters
        self.CLh_L = -0.8               # [-]       Horizontal tail landing lift coefficient
        self.CLh_TO = None              # [-]       Horizontal tail take-off lift coefficient
        self.CLah = 4                   # [/rad]    Horizontal lift curve slope

        # Vertical tail aerodynamic parameters
        self.CnB = None                 # [-]       Directional stability coefficient

        # Power parameters
        self.TTO = 1500                 # [N]       Take-off thrust

    def calcCLa(self):
        return 5.4

'''
def stability_control(variables):
    variables = sizing_htail_wingpos(variables)
    variables = verticaltail_sizing(variables)
    variables = elevator_sizing(variables)
    return variables
'''

# Test
if __name__ ==  "__main__":
    test_v = stabilitycontrol_variables()

    test_v = sizing_htail_wingpos(test_v,plot=True)
    print(test_v.xcg_wing)