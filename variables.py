import numpy as np

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


CLmaxto     = np.array([ 1.3, 1.6, 1.9 ])   # CLmax take-off
CLmaxland   = np.array([ 1.6, 2.0, 2.3 ])   # CLmax landing
CLto        = CLmaxland/(1.1**2)            # take-off CL, = CLmax,TO/1.1^2
Vs          = IAS_TAS(0, 23.15)             # stall speed (45 kts calibrated)
rho         = 1.225                         # airdensity
rho0        = 1.225                         # density at sealvl
k           = 0                             # take-off parameter
sigma       = 0                             # density ratio rho/rho0
WS          = 0                             # wing loading
sland       = 0                             # landing distance
f           = 0                             # WL/WTO
etap        = 0                             # propeller efficiency
V           = 0                             # velocity
CD0         = 0                             # drag constant
A           = 0                             # aspect ratio
e           = 0                             # oswald efficiency factor
c           = 0                             # climb gradient







