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


class CurrentVariables:
    def __init__(self, conceptnumber=2, n_engines=1, wing_mounted=False, T_tail=False, x_cg_pass=0.97, x_cg_batt=1.96,
                 x_cg_f = 1.965, ducted=False, lowwing=True, strutted_wing=False, l_fuselage=6., d_fuselage=1.2,
                 spinner_length=0.365, bulkhead_loc=0.875, cabin_length=1.3, tailheight=1.2, n_blades=2,):
        ## Requirements
        self.sto         = 500                                          # [m] take-off distance
        self.WPL         = 200*9.80665                                  # [N] Payload weight
        self.c           = 2                                            # [m/s] climb rate
        self.phi         = 60                                           # [deg] Maximum bank angle
        self.climbrate   = 2                                            # [m/s] climb rate
        self.Vmax_kts    = 120                                          # [KIAS] Max speed in knots
        self.Vmax        = self.Vmax_kts*0.514444444                    # IAS [m/s] max speed
        # self.Vmax_lvl_kts= IAS_TAS(914.4, 110*0.514444444)/0.514444444
        self.Vs_eas      = 23.15                                        # [m/s] stall speed (45 kts calibrated)
        self.Vs          = IAS_TAS(1700, self.Vs_eas)                         # [m/s] stall speed (45 kts calibrated)
        self.rho         = 1.04                                         # [kg/m3] airdensity take-off and landing
        self.M_tip       = 0.7                                          # [-] Tip mach (non-helical, only rotation)
        self.helicaltipmach = 0.75                                      # [-] Helical tip mach number, or critical mach of propeller tip airfoil
        self.rho0        = 1.225                                        # [kg/m3] density at sealvl
        self.sland       = 500                                          # [m] landing distance
        self.f           = 1                                            # [-] WL/WTO
        self.n_ult       = 5.7                                          # [-] Ultimate load factor in g
        self.rollrate    = 15*np.pi/180                                 # [rad/s] roll rate requirement
        self.propclear   = 0.23                                         # [m] Minimum clearance of the propeller or fuselage
        self.endurance_s  = 3600*2.5                                     # [s] Endurance in seconds
        self.range_m     = 250000                                       # [m] Range in meters

        ## Design concept parameters (Concept number: )
        self.concept_number = conceptnumber                             # [-] Number of the concept design
        self.n_engines   = n_engines                                    # [-] number of engines
        self.wing_mounted_engine = wing_mounted                         # Condition whether engines are mounted on the wing (otherwise, fuselage mounted)
        self.T_tail = T_tail                                            # Condition whether a t-tail configuration is used (otherwise, conventional tail)
        self.x_cg_passenger = x_cg_pass                                 # [m] Passenger CG location as measured from aircraft nose
        self.x_cg_battery = x_cg_batt                                   # [m] Battery CG location as measured from aircraft nose
        self.x_cg_fuselage = x_cg_f                                     # [m] Empty fuselage CG location as measured from the aircraft nose
        self.ducted = ducted
        if self.wing_mounted_engine:
            self.chordwise_cg_engine = 000                              # [-] Engine CG location as measured from LEMAC%
        if not self.wing_mounted_engine:
            self.x_cg_engine = 0.755                                    # [m] Engine CG location as measured from aircraft nose
        self.strutted_wing = strutted_wing
        self.l_fus = l_fuselage                                         # [m] Length of fuselage without tail
        self.d_fus = d_fuselage                                         # Depth of fuselage (?)

        self.tailtiplength= 3.325                                       # [m] from LE of main wing to most aft point of aircraft
        self.prop_spin   = 0.3650                                       # [m] length of the propeller spinner
        self.fuselage_len= 4.300                                        # [m] Length of the fuselage, nose to tail.
        self.bulkhead    = bulkhead_loc                                 # [m] Location of the bulkhead in the nose.
        self.cabinlength = cabin_length                                 # [m] Length of the cabin of the pilot from bulkhead
        self.enginelength= 0.6                                          # [m] Length of the propellers on the wing
        self.eng_perc_rootchord = 20                                    # [%] engine cg in percentage of root chord
        self.eng_height_above_w = 0.3                                   # [m] engine cg in m above the main wing
        self.tail_height = tailheight                                   # [m] Tail height above ground
        if self.ducted:
            self.AR_duct = 5.0
            self.D_fan = 1.09
            self.rho_ductmaterial = 1.5

        ## Statistical values
        self.htail_volume = 0.8                                         # [-] Horizontal tail volume
        self.vtail_volume = 0.05                                        # [-] Vertical tail volume
        self.Especif_bat = 900000                                       # J/kg Li-ion from Maarten
        self.rho_bat     = 500*3600                                     # J/L  Li-ion from Maarten
        self.motor_spec_mass = 2.5                                      # kW/kg from Maarten
        self.motor_spec_volume = 7                                      # kW/L from Maarten
        self.CD0to       = 0.0380                                       # drag constant
        self.CD0clean    = 0.0280                                       # drag constant
        self.chordwise_wing_cg = 0.30                                   # [-] 100*%MAC of wing centre of gravity
        self.aileroncutoff = 0.5                                        # [m] fixed position of outer aileron location from wingtip
        self.max_controlsurface_deflection = 25                         # [deg] maximum deflection of control surfaces
        self.c_l_a       = 5.729578                                     # [1/rad] lift slope
        self.c_l_a_flaps = 6.231                                        # [1/rad] lift slope of airfoil with flaps deployed
        self.c_l_delta_a = 0.046825*180/(np.pi)                         # [1/rad] change in the airfoil’s lift coefficient with aileron deflection
        self.fus_height  = 1.725                                        # [m] height from lowes point floor until heighest point of the fuselage
        self.prop_spin   = spinner_length                               # [m] Length of the propeller spinner

        self.W_wsn   = 5                                                #[kg]   Nose wheel weight + its strut assembly
        self.W_wsm   = 7                                                #[kg]   Main wheel weight + its strut assembly
        self.l_sn    = 0.2                                              #[m]    Shock strut length nose wheel
        self.l_sm    = 0.2                                              #[m]    Shock strut length main wheel

        self.W_avion = 15                                               #[kg] avionic weight

        # Free design choices (None means TBD)
        self.init_single_engine() if n_engines == 1 else self.init_multi_engine()
        self.lowwing     = lowwing                                      # [-] Boolean low wing or high wing.
        self.sweep_h     = 0                                            # [deg] (0- 0) Quarter chord sweep angle of the horizontal tailplane
        self.sweep_v     = 10                                         # [deg] (0-50) Quarter chord sweep angle of the vertical tailplane
        self.A_h         = 3                                            # [-] (3-5)aspect ratio of the horizontal tailplane
        self.A_v         = 2                                           # [-] (1-2) aspect ratio of the vertical tailplane
        self.x_htail     = 6.2                                          # [m] location of the ac of the horizontal tail
        self.x_vtail     = 6.2                                          # [m] location of the ac of the vertical tail
        self.taper_h     = None                                         # [-] (0.3-1.0) Taper ratio of the horizontal tailplane
        self.taper_v     = None                                         # [-] (0.3-0.7) Taper ratio of the vertical tailplane
        self.A           = 12                                            # aspect ratio of the main wing
        self.e           = 0.83                                         # oswald efficiency factor
        self.dihedral    = 0                                            # [deg] dihedral, positive upwards.
        self.n_blades    = n_blades                                     # [-] Number of propeller blades single prop
        self.rps_TO      = 2400 / 60                                    # [rps] revolution speed of prop at TO
        self.rps_cruise  = 2700 / 60                                    # [rps] revolution speed of prop at cruise
        self.coverR      = 0.5                                          # [-] Chord of duct to radius of fan ratio
        self.angleduct   = 2                                            # [deg] angle of the chordline of the duct
        self.V_eas       = 48.87                                        # [m/s] cruise velocity EAS
        self.V           = IAS_TAS(914.4, self.V_eas)                        # [m/s] cruise velocity
        self.rhocruise   = 1.12                                         # [kg/m3] airdensity cruise
        self.tcr_h        = 0.14                                        # [-] Thickness-to-rootchord ratio of the horizontal stabiliser
        self.tcr_v        = 0.14                                        # [-] Thickness-to-rootchord ratio of the vertical stabiliser    
        self.chordwise_cg_oew = 0.30                                    # Position of the OEW aircraft CG as measured from the chordwise OEW


        # Miscellaneous
        self.eff_propeller = 0.8
        self.T_propeller = None
        self.tail_ready  = False
        self.tcwing = 0.12

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
        self.bmf         = None                                         # [-] Battery mass fraction
        self.Wbat        = None                                         # [N] Battery weight
        self.Woew        = 0                                            # [N] Operational empty weight
        self.Woew_classII= None                                         # [N] Operational empty weight
        self.Wprop       = 9.81*48.2                                    # [N] Propeller weight
        self.Wmotor      = 9.81*19.75                                   # [N] Motor weight
        self.Wwing       = None                                         # [N] Wing weight
        self.Wgear_front = None                                         # [N] Front landing gear weight       
        self.Wgear_main  = None                                         # [N] Main landing gear weight
        self.Wgear       = None                                         # [N] Total landing gear weight
        self.W_htail     = None                                         # [N] Horizontal tailplane weight
        self.W_vtail     = None                                         # [N] Vertical tailplane weight
        self.xcg_aft     = None                                         # [m] Aft-most cg location
        self.xcg_fwr     = None                                         # [m] Forward-most cg location

        self.sweep       = 0                                            # [deg] Quarter chord sweep angle of the main wing
        self.taper       = 0.4                                          # [-] Taper ratio of the main wing
        self.b           = 12.2                                         # [m] Wing span
        self.cr          = 1.5                                          # [m] Root chord
        self.ct          = 0.6                                          # [m] Tip chord
        self.MAC         = 1.1                                          # [m] Mean aerodynamic chord

        self.wingpos     = 1.27                                         # [m] Location of root chord LE from nose
        self.sweep_LE    = atan((self.cr-self.ct)*0.25/self.b*2)        # [rad] Sweep of wing leading edge
        self.y_mac       = -(self.MAC-self.cr)/(self.cr-self.ct)*self.b*0.5 # [m] spanwise position of MAC
        self.x_lemac     = self.wingpos-self.y_mac*tan(self.sweep_LE)   # [m] x position from nose of LEMAC.


        self.aileronend  = None                                         # [m] End of the aileron from the aircraft center
        self.aileronstart= None                                         # [m] Start of the aileron from the aircraft center
        self.flapstart   = None                                         # [m] Start of the flap from the aircraft center
        self.Swf         = None                                         # [m^2] Wetted surface by the flaps
        self.deltaC_l_max= 1.13                                         # [-] change in C_l due to flaps at max deflection (0.2 chord @ 40 deg flap)

        self.Weng        = self.Wprop + self.Wmotor                     # [N] Total engine weight
        self.S           = self.WTO/self.WS                             # [m2] Wing surface area
        self.Sh          = None                                         # [m2] Horizontal tail surface area
        self.Sv          = None                                         # [m2] Vertical tail surface area
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

        self.propcg_x    = None
        self.propcg_z    = None
        self.enginecg_x  = None
        self.enginecg_z  = None
        self.batterycg_x = None
        self.batterycg_z = None
        self.baggagecg_x = None
        self.baggagecg_z = None
        self.payloadcg_x = None
        self.payloadcg_z = None
        self.fuselagecg_x= None
        self.fuselagecg_z= None

        self.do_engine_sizing()


    def init_single_engine(self):
        # self.CLmaxto = np.array([1.7, 1.8, 1.9])                      # CLmax take-off
        self.CLmaxto = 1.8                                              # CLmax take-off
        # self.CLmaxland = np.array([1.7, 2.0, 2.3])                      # CLmax landing
        self.CLmaxland = 2                    # CLmax landing
        # self.CLmaxclean = np.array([1.7, 1.8, 1.9])                     # CLmax clean
        self.CLmaxclean = 1.5                     # CLmax clean
        self.k = np.sqrt(5647.9 + 17.331 * self.sto) - 75.153           # [N2/m2W] take-off parameter

    def init_multi_engine(self):
        # self.CLmaxto     = np.array([ 1.8, 1.9, 2.0 ])                # CLmax take-off
        self.CLmaxto     = 1.9                                          # CLmax take-off
        # self.CLmaxland   = np.array([ 1.8, 2.1, 2.4 ])                  # CLmax landing
        self.CLmaxland   = 2.1                  # CLmax landing
        # self.CLmaxclean  = np.array([ 1.8, 1.8, 1.8 ])                  # CLmax clean
        self.CLmaxclean = 1.4
        self.k           = 93                                           # [N2/m2W] take-off parameter


    def update_WTO(self):
        self.WTO = self.Woew_classII + self.Wbat + self.WPL

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
        # self.estimate_eff_T()

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
    dictionary = vars(variables)
    print(dictionary['rhocruise'])