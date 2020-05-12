import numpy as np
from math import sqrt, pi

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
        self.A           = np.array([8, 10, 12])                        # aspect ratio
        self.e           = 0.83                                         # oswald efficiency factor
        self.CLto        = np.array(self.CLmaxto / (1.1 ** 2)).round(1) # take-off CL, = CLmax,TO/1.1^2
        self.CD0to       = 0.0380                                       # drag constant
        self.CDto = self.CD0to + self.CLto**2/(pi*self.A*self.e)        # drag during takeoff
        self.CD0clean    = 0.0280                                       # drag constant
        self.CD0 = self.CD0clean
        self.CLclimb     = self.CLmaxto[1] - 0.2                        # CLmax - safety margin fo crimb gradient
        self.CDclimb     = self.CD0to + (self.CLclimb**2)/(np.pi*self.A*self.e) # CD for climb gradient
        self.climbrate   = 2                                            # [m/s] climb rate
        self.Vs          = IAS_TAS(1700, 23.15)                         # [m/s] stall speed (45 kts calibrated)
        self.Vmax        = 120*0.514444444                              # IAS [m/s] max speed
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
        self.WTO         = 750*9.81                                     # [N] take-off weight
        self.Wbat        = 0                                            # [N] Battery weight
        self.Woew        = 0                                            # [N] Operational empty weight
        self.WPL         = 200*9.80665                                  # [N] Payload weight
        self.P           = self.WTO/self.WP                             # [W] Engine power
        #self.propefficiency = self.WP*(sqrt(2*self.WS/self.rho/self.CLto)+self.climbrate) # [-] propeller efficiency (Pa to thrust)
        self.update_engine_power()
        self.n_blades    = 6                                            # [-] Number of propeller blades single prop
        self.prop_d      = 20 * (self.P*0.00134102209)**0.25 * 0.3048   # [m] Propeller diameter
        self.sweep       = 0                                            # [deg] Quarter chord sweep angle
        self.taper       = 0.4                                          # [-] Taper ratio
        self.b           = 13.78                                        # [m] Wing span
        self.cr          = 1.64                                         # [m] Root chord
        self.ct          = 0.66                                         # [m] Tip chord
        self.MAC         = 1.22                                         # [m] Mean aerodynamic chord


    def init_single_engine(self):
        self.CLmaxto = np.array([1.7, 1.8, 1.9])                # CLmax take-off
        self.CLmaxland = np.array([1.7, 2.0, 2.3])              # CLmax landing
        self.CLmaxclean = np.array([1.7, 1.8, 1.9])             # CLmax clean
        self.k = np.sqrt(5647.9 + 17.331 * self.sto) - 75.153   # [N2/m2W] take-off parameter

    def init_multi_engine(self):
        self.CLmaxto     = np.array([ 1.8, 1.9, 2.0 ])               # CLmax take-off
        self.CLmaxland   = np.array([ 1.8, 2.1, 2.4 ])               # CLmax landing
        self.CLmaxclean  = np.array([ 1.8, 1.8, 1.8 ])               # CLmax clean
        self.k           = 93                                        # [N2/m2W] take-off parameter

    def update_engine_power(self):
        self.P_single_engine = self.P/self.n_engines


if __name__ == "__main__":
    variables = CurrentVariables()
    print(variables.WS)