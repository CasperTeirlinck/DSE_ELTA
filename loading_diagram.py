import numpy as np

def stallspeed(CLmax, Vs, rho):
    """
    :param CLmax: CLmax of the aircraft, scalar
    :param Vs: Stall speed, scalar
    :param rho: Airdensity, scalar
    :return: Stall speed requirement, W/S scalar
    """
    return 0.5*rho*(Vs**2)*CLmax


def takeoff(k, CLTO, sigma, WS):
    """
    :param k: take-off parameter, scalar
    :param CLTO: take-off CL, scalar
    :param sigma: density ratio rho/rho0, scalar
    :param WS: plot x values, array
    :return: WP: plot y values, array
    """
    return (k/WS)*CLTO*sigma


def landing(CLmax, rho, sland, f):
    """
    :param CLmax: CLmax of the aircraft, scalar
    :param rho: air density, scalar
    :param sland: landing distance, scalar
    :param f: WL/WTO
    :return: landing speed requirement, W/S scalar
    """
    return (CLmax*rho*(sland/0.5915))/(2*f)


def cruisspeed(etap, rho, rho0, CD0, V, A, e, WS):
    """
    :param etap: NO IDEA WHAT THIS IS, scalar
    :param rho: density, scalar
    :param rho0: density at sealvl, scalar
    :param CD0: drag constant, scalar
    :param V: Velocity, scalar
    :param A: Aspect ratio, scalar
    :param e: oswald efficiency factor
    :param WS: plot x values, array
    :return: WP: plot y values, array
    """
    return (0.9/0.8)*etap*(rho/rho0)**(3/4)*((CD0*0.5*rho*V**3)/(0.8*WS) + (0.8*WS)*(1/(np.pi*A*e*0.5*rho*V)))


def climbrate(etap, rho, A, e, CD0, c, WS):
    """
    :param etap: NO IDEA WHAT THIS IS, scalar
    :param rho: density, scalar
    :param A: Aspect ratio, scalar
    :param e: oswald efficiency factor
    :param CD0: drag constant, scalar
    :param c: NO IDEA WHAT THIS IS, scalar
    :param WS: plot x values, array
    :return: plot y values, array
    """
    return etap/(c+ (np.sqrt(WS*(2/rho)))/(1.345*((A*e)**(3/4))/(CD0**(1/4))))


def climbgradient(etap, cV, CD, CL, rho, WS):
    """
    :param etap: NO IDEA WHAT THIS IS, scalar
    :param cV: NO IDEA WHAT THIS IS, scalar
    :param CD: drag coefficient, scalar
    :param CL: lift coefficient, scalar
    :param rho: density, scalar
    :param WS: plot x values, array
    :return: plot y values, array
    """
    return etap/(np.sqrt(WS)*(cV + CD/CL)*np.sqrt((2/rho)*(1/CL))