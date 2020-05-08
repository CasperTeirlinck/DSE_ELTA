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


CLmax =                     # CLmax of the aircraft
Vs    = IAS_TAS(0, 23.15)   # stall speed (45 kts calibrated)
rho   =                     # airdensity
rho0  = 1.225               # density at sealvl
k     =                     # take-off parameter
CLTO  =                     # take-off CL, = CLmax,TO/1.1^2
sigma =                     # density ratio rho/rho0
WS    =                     # wing loading
sland =                     # landing distance
f     =                     # WL/WTO
etap  =                     # propeller efficiency
V     =                     # velocity
CD0   =                     # drag constant
A     =                     # aspect ratio
e     =                     # oswald efficiency factor
c     =                     # climb gradient







