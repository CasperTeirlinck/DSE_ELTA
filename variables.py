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


# Single engine
# A           = np.array([8, 10, 12])                     # aspect ratio
# e           = 0.83                                      # oswald efficiency factor
# CLmaxto     = np.array([ 1.7, 1.8, 1.9 ])               # CLmax take-off
# CLmaxland   = np.array([ 1.7, 2.0, 2.3 ])               # CLmax landing
# CLmaxclean  = np.array([ 1.7, 1.8, 1.9 ])               # CLmax clean
# CLto        = np.array(CLmaxto/(1.1**2)).round(1)       # take-off CL, = CLmax,TO/1.1^2
# CD0to       = 0.0380                                    # drag constant
# CD0clean    = 0.0280                                    # drag constant
# CLclimb     = CLmaxto[1] - 0.2                          # CLmax - safety margin fo crimb gradient
# CDclimb     = CD0to + (CLclimb**2)/(np.pi*A*e)          # CD for climb gradient
# Vs          = IAS_TAS(1700, 23.15)                      # [m/s] stall speed (45 kts calibrated)
# V           = IAS_TAS(914.4, 48.87)                     # [m/s] velocity
# rho         = 1.04                                      # [kg/m3] airdensity take-off and landing
# rhocruise   = 1.12                                      # [kg/m3] airdensity cruise
# rho0        = 1.225                                     # [kg/m3] density at sealvl
# sland       = 500                                       # [m] landing distance
# sto         = 500                                       # [m] take-off distance
# k           = np.sqrt(5647.9 + 17.331 * sto) - 75.153   # [N2/m2W] take-off parameter
# sigma       = rho/rho0                                  # density ratio rho/rho0
# f           = 1                                         # WL/WTO
# etap        = 0.7                                       # propeller efficiency
# c           = 2                                         # [m/s] climb rate
# phi         = 60                                        # [deg] bank angle

# Twin engine
# A           = np.array([8, 10, 12])                     # aspect ratio
# e           = 0.83                                      # oswald efficiency factor
# CLmaxto     = np.array([ 1.8, 1.9, 2.0 ])               # CLmax take-off
# CLmaxland   = np.array([ 1.8, 2.1, 2.4 ])               # CLmax landing
# CLmaxclean  = np.array([ 1.8, 1.8, 1.8 ])               # CLmax clean
# CLto        = np.array(CLmaxto/(1.1**2)).round(1)       # take-off CL, = CLmax,TO/1.1^2
# CD0to       = 0.0380                                    # drag constant
# CD0clean    = 0.0280                                    # drag constant
# CLclimb     = CLmaxto[1] - 0.2                          # CLmax - safety margin fo crimb gradient
# CDclimb     = CD0to + (CLclimb**2)/(np.pi*A*e)          # CD for climb gradient
# Vs          = IAS_TAS(1700, 23.15)                      # [m/s] stall speed (45 kts calibrated)
# V           = IAS_TAS(914.4, 48.87)                     # [m/s] velocity
# rho         = 1.04                                      # [kg/m3] airdensity take-off and landing
# rhocruise   = 1.12                                      # [kg/m3] airdensity cruise
# rho0        = 1.225                                     # [kg/m3] density at sealvl
# sland       = 500                                       # [m] landing distance
# sto         = 500                                       # [m] take-off distance
# k           = 93                                        # [N2/m2W] take-off parameter
# sigma       = rho/rho0                                  # density ratio rho/rho0
# f           = 1                                         # WL/WTO
# etap        = 0.7                                       # propeller efficiency
# c           = 2                                         # [m/s] climb rate
# phi         = 60                                        # [deg] bank angle


class CurrentVariables():
    def __init__(self, n_engines=1):
        self.n_engines   = n_engines                                    # [-] amount of engines
        self.sto         = 500                                          # [m] take-off distance
        self.init_single_engine() if n_engines == 1 else self.init_multi_engine()
        self.A           = 12                                           # aspect ratio
        self.e           = 0.83                                         # oswald efficiency factor
        self.CLto        = np.array(self.CLmaxto / (1.1 ** 2)).round(1) # take-off CL, = CLmax,TO/1.1^2
        self.CD0to       = 0.0380                                       # drag constant
        self.CDto = self.CD0to + self.CLto**2/(pi*self.A*self.e)        # drag during takeoff
        self.CD0clean    = 0.0280                                       # drag constant
        self.CD0 = self.CD0clean
        self.CLclimb     = self.CLmaxto - 0.2                           # CLmax - safety margin fo crimb gradient
        self.CDclimb     = self.CD0to + (self.CLclimb**2)/(np.pi*self.A*self.e) # CD for climb gradient
        self.climbrate   = 2                                            # [m/s] climb rate
        self.Vs          = IAS_TAS(1700, 23.15)                         # [m/s] stall speed (45 kts calibrated)
        self.Vmax_kts    = 120                                          # [KIAS] Max speed in knots
        self.Vmax        = self.Vmax_kts*0.514444444                    # IAS [m/s] max speed
        self.V           = IAS_TAS(914.4, 48.87)                        # [m/s] cruise velocity
        self.rho         = 1.04                                         # [kg/m3] airdensity take-off and landing
        self.rhocruise   = 1.12                                         # [kg/m3] airdensity cruise
        self.rho0        = 1.225                                        # [kg/m3] density at sealvl
        self.sland       = 500                                          # [m] landing distance
        self.sigma       = self.rho/self.rho0                           # density ratio rho/rho0
        self.f           = 1                                            # WL/WTO
        self.etap        = 0.7                                          # propeller efficiency
        self.c           = 2                                            # [m/s] climb rate
        self.phi         = 60                                           # [deg] bank angle
        self.WP          = 0.1218                                       # [N/W] power loading
        self.WS          = 465                                          # [N/m2] wing loading
        self.WTO         = 750*9.80665                                  # [N] take-off weight
        self.Wbat        = 0                                            # [N] Battery weight
        self.Woew        = 0                                            # [N] Operational empty weight
        self.WPL         = 200*9.80665                                  # [N] Payload weight
        self.S           = self.WTO/self.WS                             # [m] Wing surface area
        self.n_blades    = 4                                            # [-] Number of propeller blades single prop
        self.M_tip       = 0.7                                          # [-] Tip mach (non-helical, only rotation)
        self.helicaltipmach = 0.75                                      # [-] Helical tip mach number, or critical mach of propeller tip airfoil
        self.rps_TO      = 2400 / 60                                    # [rps] revolution speed of prop at TO
        self.rps_cruise  = 2700 / 60                                    # [rps] revolution speed of prop at cruise
        self.C_d_prop    = 0.012                                        # [-] 2D airfoil drag of propeller
        self.prop_A_ratio= (81836-12446)/372230                         # [-] propeller area ratio: Lockheed YO-3 estimation
        self.contrarotate= True
        self.duct_t_over_c= 0.0001
        self.do_engine_sizing(self.contrarotate, self.duct_t_over_c)



    def init_single_engine(self):
        # self.CLmaxto = np.array([1.7, 1.8, 1.9])                # CLmax take-off
        self.CLmaxto = 1.8                                      # CLmax take-off
        self.CLmaxland = np.array([1.7, 2.0, 2.3])              # CLmax landing
        self.CLmaxclean = np.array([1.7, 1.8, 1.9])             # CLmax clean
        self.k = np.sqrt(5647.9 + 17.331 * self.sto) - 75.153   # [N2/m2W] take-off parameter

    def init_multi_engine(self):
        # self.CLmaxto     = np.array([ 1.8, 1.9, 2.0 ])               # CLmax take-off
        self.CLmaxto     = 1.9                                       # CLmax take-off
        self.CLmaxland   = np.array([ 1.8, 2.1, 2.4 ])               # CLmax landing
        self.CLmaxclean  = np.array([ 1.8, 1.8, 1.8 ])               # CLmax clean
        self.k           = 93                                        # [N2/m2W] take-off parameter

    def calc_engine_efficiency(self, tolerance=1e-3, maxiterations=1000):
        for _ in range(maxiterations):
            eta2 = 1 - 4 / pi**3 * self.eta1 * self.C_P / self.J
            print("eta2", eta2)
            phi = atan(self.J/pi/self.eta1/eta2)
            eta3 = 1 - pi*pi*pi*pi*eta2*eta2*self.prop_A_ratio*self.C_d_prop*glauert_function(phi)/8/self.C_P
            print("eta3", eta3)
            eta1 = 1 - 2/pi * self.C_P * eta2 * eta3 * (self.eta1/self.J)**3
            print("eta1", eta1)
            if abs(eta2-self.eta2) < tolerance and abs(eta1-self.eta1) < tolerance and abs(eta3-self.eta3) < tolerance:
                self.eta1, self.eta2, self.eta3 = eta1, eta2, eta3
                return eta1*eta2*eta3
            self.eta1, self.eta2, self.eta3 = eta1, eta2, eta3

        else:
            print("Process did not converge.")
            return 0.0


    def do_engine_sizing(self, contra, t_over_c_duct, s=1.2, z=0.075):
        # Using methods from Rik's book
        self.P_total           = self.WTO/self.WP                           # [W] Engine power
        self.P           = self.P_total/self.n_engines                            # Single engine thrust
        PBHP             = self.P * 0.0013410220888           # [BPH] Engine power
        CL, CD = sqrt(pi * self.A * self.e * self.CD0clean), 2 * self.CD0clean

        if self.n_blades == 2:
            factor = 0.56
        elif self.n_blades == 3:
            factor = 0.52
        else:
            factor = 0.49
        self.prop_d      = factor*(self.P*0.001)**0.25
        self.T_static = 78.034*PBHP**(2/3)*(self.rho*self.prop_d*self.prop_d)**(1/3)
        self.T = 0.5*self.rho*self.V**2*self.WTO/self.WS*CD
        self.J           = self.Vmax / self.rps_TO / self.prop_d        # [-] Advance ratio
        self.CP          = self.P / (self.rho * self.rps_cruise**3 * self.prop_d**5)
        self.CT_static   = self.T_static / (self.rho * self.rps_cruise**2 * self.prop_d**4)
        self.CT          = self.T / (self.rho * self.rps_cruise**2 * self.prop_d**4)
        self.CQ          = self.CP*0.5/pi
        self.prop_eff    = self.J * self.CT / self.CP
        print(self.J, "This")

        self.speedfactor = 1
        if t_over_c_duct != 0:
            c = 0.5*self.prop_d*(1/(1/s-t_over_c_duct))
            Re = 0.5*self.prop_d+0.5*c*t_over_c_duct
            delta_d = 1 - sqrt(2*Re/self.prop_d)*((0.458+4.431*s)/(1+1.089*s)*z+(2.033+4.88*s)/(1+0.893*s)*s*z*z)
            K = 1
            CT_noduct = self.T/(0.5*self.rho*self.V**2*pi*self.prop_d**2*0.25)
            delta_i = K * (sqrt(1+CT_noduct)-1)
            print(delta_i)
            print(delta_d)
            self.speedfactor = 1+(delta_d+delta_i)



        self.propefficiency = self.WP**(sqrt(2*self.WS/self.rho/self.CLto)+self.climbrate) # [-] propeller efficiency (Pa to thrust)
        self.Vclimb      = sqrt(self.WS/(0.5*self.rho*self.CLto))
        self.M_to        = self.Vclimb/331                              # [-] Take-off mach number (critical for noise)
        # self.prop_d      = 331/pi/self.engine_rps*sqrt(self.helicaltipmach**2-self.M_to**2) # [m] Propeller diameter
        self.prop_d     = 0.49 * (self.P/1000)**0.25
        self.C_P         = self.P / self.rho / self.prop_d ** 3 / self.rps_TO ** 5
        # Using 'Maximum Propeller Efficiency Estimation' - Ohad Gur
        self.eta1, self.eta2, self.eta3 = 0.5, 0.5, 0.5
        # self.propefficiency = self.calc_engine_efficiency()


if __name__ == "__main__":
    variables = CurrentVariables()
    # print(variables.prop_d3)
    print(variables.speedfactor)
    variables.do_engine_sizing(variables.contrarotate, variables.duct_t_over_c, s=2., z=0.1)
    print(variables.speedfactor)
    print()
