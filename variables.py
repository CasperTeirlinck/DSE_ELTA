import numpy as np
from math import sqrt, pi, atan, sin, cos, tan, log

def IAS_TAS(h, V_IAS):
    """
    Calculates the true airspeed at for a certain indicated airspeed at a altitude
    :param h: alititude [m]
    :param V_IAS: indicated airspeed [m/s]
    :return: V_TAS: true airspeed [m/s]
    """
    p0 = 101325
    rho0 = 1.225
    T0 = 288.15
    g0 = 9.80665
    lmbda = -0.0065
    R = 287.05
    gamma = 1.4

    p = p0*(1 + lmbda*(h/T0))**(-g0/(lmbda*R))
    T = T0 + lmbda*h
    a = np.sqrt(gamma*R*T)

    M = np.sqrt((2/(gamma - 1))*((1 + (p0/p)*((1 + (gamma - 1)/(2*gamma)*(rho0/p0)*(V_IAS**2))**(gamma/(gamma-1)) - 1))**((gamma - 1)/gamma) - 1))
    V_TAS = M*a

    return V_TAS


def glauert_function(phi):
    return (2+5*tan(phi)**2)/(8*cos(phi)) - 3 / 16 * tan(phi)**4*log((1-cos(phi))/(1+cos(phi)), 10)


class CurrentVariables():
    def __init__(self, n_engines=1, wing_mounted=False, T_tail=False, x_cg_pass=0.97, x_cg_batt=1.96,
                 x_cg_f = 2.33, ducted=False):
        ## Requirements
        self.sto         = 500                                          # [m] take-off distance
        self.WPL         = 200*9.80665                                  # [N] Payload weight
        self.c           = 2                                            # [m/s] climb rate
        self.phi         = 60                                           # [deg] Maximum bank angle
        self.climbrate   = 2                                            # [m/s] climb rate
        self.Vmax_kts    = 120                                          # [KIAS] Max speed in knots
        self.Vmax        = self.Vmax_kts*0.514444444                    # IAS [m/s] max speed
        self.Vs          = IAS_TAS(1700, 23.15)                         # [m/s] stall speed (45 kts calibrated)
        self.rho         = 1.04                                         # [kg/m3] airdensity take-off and landing
        self.M_tip       = 0.7                                          # [-] Tip mach (non-helical, only rotation)
        self.helicaltipmach = 0.75                                      # [-] Helical tip mach number, or critical mach of propeller tip airfoil
        self.rho0        = 1.225                                        # [kg/m3] density at sealvl
        self.sland       = 500                                          # [m] landing distance
        self.f           = 1                                            # WL/WTO

        ## Design concept parameters (Concept number: )
        self.n_engines   = n_engines                                    # [-] number of engines
        self.wing_mounted_engine = wing_mounted                         # Condition whether engines are mounted on the wing (otherwise, fuselage mounted)
        self.T_tail = T_tail                                            # Condition whether a t-tail configuration is used (otherwise, conventional tail)
        self.x_cg_passenger = x_cg_pass                                 # [m] Passenger CG location as measured from aircraft nose
        self.x_cg_battery = x_cg_batt                                   # [m] Battery CG location as measured from aircraft nose
        self.x_cg_fuselage = x_cg_f                                     # [m] Empty fuselage CG location as measured from the aircraft nose
        self.ducted = ducted
        if self.wing_mounted_engine:
            self.x_cg_engine = 0.755                                    # [m] Engine CG location as measured from aircraft nose
        if not self.wing_mounted_engine:
            self.chordwise_cg_engine = 000                              # [-] Engine CG location as measured from LEMAC%



        ## Statistical values
        self.htail_volume = 0.8                                         # [-] Horizontal tail volume
        self.vtail_volume = 0.4                                         # [-] Vertical tail volume
        self.Especif_bat = 900000                                       # J/kg Li-ion from Maarten
        self.rho_bat     = 500*3600                                     # J/L  Li-ion from Maarten
        self.CD0to       = 0.0380                                       # drag constant
        self.CD0clean    = 0.0280                                       # drag constant


        # Free design choices (None means TBD)
        self.init_single_engine() if n_engines == 1 else self.init_multi_engine()
        self.sweep_h     = 0                                            # [deg] (0- 0) Quarter chord sweep angle of the horizontal tailplane
        self.sweep_v     = None                                         # [deg] (0-50) Quarter chord sweep angle of the vertical tailplane
        self.A_h         = 3                                            # [-] (3-5)aspect ratio of the horizontal tailplane
        self.A_v         = 12                                           # [-] (1-2) aspect ratio of the vertical tailplane
        self.x_htail     = 4.7                                          # [m] location of the ac of the horizontal tail
        self.x_vtail     = 4.7                                          # [m] location of the ac of the vertical tail
        self.taper_h     = None                                         # [-] (0.3-1.0) Taper ratio of the horizontal tailplane
        self.taper_v     = None                                         # [-] (0.3-0.7) Taper ratio of the vertical tailplane
        self.A           = 12                                           # aspect ratio of the main wing
        self.e           = 0.83                                         # oswald efficiency factor
        self.n_blades    = 4                                            # [-] Number of propeller blades single prop
        self.rps_TO      = 2400 / 60                                    # [rps] revolution speed of prop at TO
        self.rps_cruise  = 2700 / 60                                    # [rps] revolution speed of prop at cruise
        self.coverR      = 0.5                                          # [-] Chord of duct to radius of fan ratio
        self.angleduct   = 2                                            # [deg] angle of the chordline of the duct
        self.initial_etap        = 0.7                                  # propeller efficiency
        self.V           = IAS_TAS(914.4, 48.87)                        # [m/s] cruise velocity
        self.rhocruise   = 1.12                                         # [kg/m3] airdensity cruise

        # Miscellaneous
        self.eff_propeller = None
        self.T_propeller = None

        ## Calculated values
        self.CLto        = np.array(self.CLmaxto / (1.1 ** 2)).round(1) # take-off CL, = CLmax,TO/1.1^2
        self.CDto = self.CD0to + self.CLto**2/(pi*self.A*self.e)        # drag during takeoff
        self.CD0 = self.CD0clean
        self.CLclimb     = self.CLmaxto - 0.2                           # CLmax - safety margin fo crimb gradient
        self.CDclimb     = self.CD0to + (self.CLclimb**2)/(np.pi*self.A*self.e) # CD for climb gradient

        self.sigma       = self.rho/self.rho0                           # density ratio rho/rho0
        if n_engines == 1:
            self.WP      = 0.1218                                       # [N/W] power loading
            self.WS      = 465                                          # [N/m2] wing loading
        else:
            self.WP      = 0.149                                        # [N/W] power loading
            self.WS      = 592                                          # [N/m2] wing loading
        self.WTO         = 750*9.81                                     # [N] take-off weight
        self.Wbat        = 0                                            # [N] Battery weight
        self.Woew        = 0                                            # [N] Operational empty weight
        self.Woew_classII= 1000 * 9.81                                  # [N] Operational empty weight
        self.Wprop       = 9.81*48.2                                    # [N] Propeller weight
        self.Wmotor      = 9.81*19.75                                   # [N] Motor weight

        self.sweep       = 0                                            # [deg] Quarter chord sweep angle of the main wing
        self.taper       = 0.4                                          # [-] Taper ratio of the main wing
        self.b           = 12.2                                         # [m] Wing span
        self.cr          = 1.5                                          # [m] Root chord
        self.ct          = 0.6                                          # [m] Tip chord
        self.MAC         = 1.1                                          # [m] Mean aerodynamic chord
        self.chordwise_cg_oew = 0.25                                    # Position of the OEW aircraft CG as measured from the chordwise OEW

        self.Weng        = self.Wprop + self.Wmotor                     # [N] Total engine weight
        self.S           = self.WTO/self.WS                             # [m] Wing surface area
        self.do_engine_sizing()
        self.eff_tot_prop= 0.95*self.eff_propeller                      # Total propulsion efficiency (motor and bat)
        self.x_maingear = None                                          # [m] from nose, more negative is further from nose
        self.y_maingear = None                                          # [m] from nose, right wing positive
        self.z_maingear = None                                          # [m] from top of fuselage, down positive
        self.x_nosegear = None                                          # [m] from nose, more negative is further from nose
        self.maingeardiameter = None                                    # [m] diameter of main gear wheel
        self.maingearwidth = None                                       # [m] width of main gear wheel
        self.nosegeardiameter = None                                    # [m] diameter of nose gear wheel
        self.nosegearwidth = None                                       # [m] width of nose gear wheel






    def init_single_engine(self):
        # self.CLmaxto = np.array([1.7, 1.8, 1.9])                      # CLmax take-off
        self.CLmaxto = 1.8                                              # CLmax take-off
        self.CLmaxland = np.array([1.7, 2.0, 2.3])                      # CLmax landing
        self.CLmaxclean = np.array([1.7, 1.8, 1.9])                     # CLmax clean
        self.k = np.sqrt(5647.9 + 17.331 * self.sto) - 75.153           # [N2/m2W] take-off parameter

    def init_multi_engine(self):
        # self.CLmaxto     = np.array([ 1.8, 1.9, 2.0 ])                # CLmax take-off
        self.CLmaxto     = 1.9                                          # CLmax take-off
        self.CLmaxland   = np.array([ 1.8, 2.1, 2.4 ])                  # CLmax landing
        self.CLmaxclean  = np.array([ 1.8, 1.8, 1.8 ])                  # CLmax clean
        self.k           = 93                                           # [N2/m2W] take-off parameter


    def do_engine_sizing(self):
        # Using methods from Rik's book
        self.P_total           = self.WTO/self.WP                       # [W] Engine power
        self.P           = self.P_total/self.n_engines                  # Single engine thrust

        if self.n_blades == 2:
            factor = 0.56
        elif self.n_blades == 3:
            factor = 0.52
        else:
            factor = 0.49

        self.prop_d      = factor*(self.P*0.001)**0.25
        self.estimate_eff_T()

    def estimate_eff_T(self):
        if self.eff_propeller == None:
            self.eff_propeller = 0.8
            self.T_propeller = 800
        else:
            if self.ducted:
                print("Use Javaprop to find efficiency. Parameters to use are P={} W, blades={} [-], V={} m/s, rpm={}, c/R={} [-], angle={} deg".format(self.P, self.n_blades, self.Vmax, self.rps_TO*60, self.coverR, self.angleduct))
            else:
                print("Use Javaprop to find efficiency. Parameters to use are P={} W, blades={} [-], V={} m/s, rpm={}.".format(self.P, self.n_blades, self.Vmax, self.rps_TO*60))
            self.eff_propeller = float(input("Efficiency [%]: "))/100       # [-] Efficiency of a propeller
            self.T_propeller = input("Thrust [N]: ")                        # [N] Thrust of a single engine
        self.T_total = self.T_propeller * self.n_engines                # [N] Total thrust of aircraft


if __name__ == "__main__":
    variables = CurrentVariables(n_engines=1)