import numpy as np
import matplotlib.pyplot as plt
import variables

def stallspeed(CLmax, Vs, rho):
    """
    :param CLmax: CLmax of the aircraft, scalar
    :param Vs: Stall speed, scalar
    :param rho: Airdensity, scalar; sea level can be assumed
    :return: Stall speed requirement, W/S scalar
    """
    return 0.5*rho*(Vs**2)*CLmax


def takeoff(k, CLTO, sigma, WS):
    """
    :param k: take-off parameter, scalar
    :param CLTO: take-off CL, = CLmax,TO/1.1^2 scalar
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
    :param etap: propeller efficiency, scalar
    :param rho: density, scalar
    :param rho0: density at sealvl, scalar
    :param CD0: drag constant, scalar
    :param V: Velocity, scalar
    :param A: Aspect ratio, scalar
    :param e: oswald efficiency factor
    :param WS: plot x values, array
    :return: WP: plot y values, array
    """
    return (0.9/0.8)*etap*(rho/rho0)**(3/4)*((((CD0*0.5*rho*V**3)/(0.8*WS)) + (0.8*WS)*(1/(np.pi*A*e*0.5*rho*V)))**(-1))


def climbrate(etap, rho, A, e, CD0, c, WS):
    """
    :param etap: propeller efficiency, scalar
    :param rho: density, scalar
    :param A: Aspect ratio, scalar
    :param e: oswald efficiency factor
    :param CD0: drag constant, scalar
    :param c: NO IDEA WHAT THIS IS, scalar
    :param WS: plot x values, array
    :return: plot y values, array
    """
    return etap/(c + (np.sqrt(WS*(2/rho)))/(1.345*((A*e)**(3/4))/(CD0**(1/4))))


def climbgradient(etap, cV, CD, CL, rho, WS):
    """
    :param etap: propeller efficiency, scalar
    :param cV: NO IDEA WHAT THIS IS, scalar
    :param CD: drag coefficient, scalar
    :param CL: lift coefficient, scalar
    :param rho: density, scalar
    :param WS: plot x values, array
    :return: plot y values, array
    """
    return etap/(np.sqrt(WS)*(cV + CD/CL)*np.sqrt((2/rho)*(1/CL)))


## Plotting
WS_plot = np.arange(0.1, 2000, 0.1)

# stallspeed
WP_s = stallspeed(CLmax = 1.3, Vs = 23.15, rho = 1.225)
plt.plot([WP_s, WP_s], [0, 0.5], label = "Stall, CLmax = 1.3" )

# take off
WP_to = takeoff(k = 100, CLTO= 1.3, sigma = 1, WS= WS_plot)
plt.plot(WS_plot, WP_to, label = "Take-off = 1.3")

# landing
WP_landing = landing(CLmax = 1.3, rho = 1.225, sland=500, f=1)
plt.plot([WP_landing, WP_landing], [0, 0.5], label = "Landing = 1.3")

# cruise
WP_cruise = cruisspeed(etap=0.7, rho=1.225, rho0=1.225, CD0=0.05, V=70, A=10, e=0.7, WS=WS_plot)
plt.plot(WS_plot, WP_cruise, label = "Cruise")

# climbrate
WP_climbrate = climbrate(etap=0.7, rho=1.225, A=10, e=0.7, CD0=0.05, c=2, WS=WS_plot)
plt.plot(WS_plot, WP_climbrate, label = "Climb rate")

# climbgrad
WP_climbgrad = climbgradient(etap=0.7, cV=2/70, CD=0.7, CL=1.2, rho=1.225, WS=WS_plot)
plt.plot(WS_plot, WP_climbgrad, label = "Climb gradient")

plt.ylim(0, 0.4)
plt.xlim(0, 2000)
plt.legend()
plt.show()

