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


def takeoff(k, CLto, sigma, WS):
    """
    :param k: take-off parameter, scalar
    :param CLto: take-off CL, = CLmax,TO/1.1^2 scalar
    :param sigma: density ratio rho/rho0, scalar
    :param WS: plot x values, array
    :return: WP: plot y values, array
    """
    return (k/WS)*CLto*sigma


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

## Determining design point
# Calculate W/S value of vertical limits
WS_s = stallspeed(CLmax=variables.CLmaxclean[1], Vs=variables.Vs, rho=variables.rho) # Stallspeed
WS_landing1 = landing(CLmax=variables.CLmaxland[0], rho=variables.rho, sland=variables.sland, f=variables.f) # Pessimistic landing 
WS_landing2 = landing(CLmax=variables.CLmaxland[1], rho=variables.rho, sland=variables.sland, f=variables.f) # Neutral landing 
WS_landing3 = landing(CLmax=variables.CLmaxland[2], rho=variables.rho, sland=variables.sland, f=variables.f) # Optimistic landing 

# Calculate limiting W/S value
WS_limit = min([WS_s,WS_landing1,WS_landing2,WS_landing3])

# Calculate W/P value of curved limits
WP_to1 = takeoff(k=variables.k, CLto=variables.CLto[0], sigma=variables.sigma, WS= WS_limit) # Pessimistic take-off
WP_to2 = takeoff(k=variables.k, CLto=variables.CLto[1], sigma=variables.sigma, WS= WS_limit) # Neutral take-off
WP_to3 = takeoff(k=variables.k, CLto=variables.CLto[2], sigma=variables.sigma, WS= WS_limit) # Optimistic take-off

WP_cruise1 = cruisspeed(etap=variables.etap, rho=variables.rho, rho0=variables.rho, CD0=variables.CD0clean, V=variables.V, A=variables.A[0], e=variables.e, WS=WS_limit) # Pessimistic cruise
WP_cruise2 = cruisspeed(etap=variables.etap, rho=variables.rho, rho0=variables.rho, CD0=variables.CD0clean, V=variables.V, A=variables.A[1], e=variables.e, WS=WS_limit) # Neutral cruise
WP_cruise3 = cruisspeed(etap=variables.etap, rho=variables.rho, rho0=variables.rho, CD0=variables.CD0clean, V=variables.V, A=variables.A[2], e=variables.e, WS=WS_limit) # Optimistic cruise

WP_climbrate1 = climbrate(etap=variables.etap, rho=variables.rho, A=variables.A[0], e=variables.e, CD0=variables.CD0to, c=variables.c, WS=WS_limit) # Pessimistic climb rate
WP_climbrate2 = climbrate(etap=variables.etap, rho=variables.rho, A=variables.A[1], e=variables.e, CD0=variables.CD0to, c=variables.c, WS=WS_limit) # Neutral climb rate
WP_climbrate3 = climbrate(etap=variables.etap, rho=variables.rho, A=variables.A[2], e=variables.e, CD0=variables.CD0to, c=variables.c, WS=WS_limit) # Optimistic climb rate

WP_climbgrad1 = climbgradient(etap=variables.etap, cV=variables.c/variables.V, CD=variables.CDclimb[0], CL=variables.CLclimb, rho=variables.rho, WS=WS_limit) # Pessimistic climb gradient
WP_climbgrad2 = climbgradient(etap=variables.etap, cV=variables.c/variables.V, CD=variables.CDclimb[1], CL=variables.CLclimb, rho=variables.rho, WS=WS_limit) # Neutral climb rate
WP_climbgrad3 = climbgradient(etap=variables.etap, cV=variables.c/variables.V, CD=variables.CDclimb[2], CL=variables.CLclimb, rho=variables.rho, WS=WS_limit) # Optimistic climb rate

# Calculate limiting W/P values
WP_limit = min([WP_to1,WP_to2,WP_to3,WP_cruise1,WP_cruise2,WP_cruise3,WP_climbrate1,WP_climbrate2,WP_climbrate3,WP_climbgrad1,WP_climbgrad2,WP_climbgrad3])


## Plotting
WS_plot = np.arange(100, 1500, 0.1)

# stallspeed
WP_s = stallspeed(CLmax=variables.CLmaxclean[1], Vs=variables.Vs, rho=variables.rho)
plt.plot([WP_s, WP_s], [0, 0.5], label = f"Stall, CLmax = {variables.CLmaxclean[1]}", color='c' )

# take off
WP_to = takeoff(k=variables.k, CLto=variables.CLto[0], sigma=variables.sigma, WS= WS_plot)
plt.plot(WS_plot, WP_to, label = f"Take-off, CLto = {variables.CLto[0]}", color='firebrick')

WP_to = takeoff(k=variables.k, CLto=variables.CLto[1], sigma=variables.sigma, WS= WS_plot)
plt.plot(WS_plot, WP_to, label = f"Take-off, CLto = {variables.CLto[1]}", color='indianred')

WP_to = takeoff(k=variables.k, CLto=variables.CLto[2], sigma=variables.sigma, WS= WS_plot)
plt.plot(WS_plot, WP_to, label = f"Take-off, CLto = {variables.CLto[2]}", color='lightcoral')

# landing
WP_landing = landing(CLmax=variables.CLmaxland[0], rho=variables.rho, sland=variables.sland, f=variables.f)
plt.plot([WP_landing, WP_landing], [0, 0.5], label = f"Landing, CL = {variables.CLmaxland[0]}", color='forestgreen')

WP_landing = landing(CLmax=variables.CLmaxland[1], rho=variables.rho, sland=variables.sland, f=variables.f)
plt.plot([WP_landing, WP_landing], [0, 0.5], label = f"Landing, CL = {variables.CLmaxland[1]}", color='limegreen')

WP_landing = landing(CLmax=variables.CLmaxland[2], rho=variables.rho, sland=variables.sland, f=variables.f)
plt.plot([WP_landing, WP_landing], [0, 0.5], label = f"Landing, CL = {variables.CLmaxland[2]}", color='darkgreen')

# cruise
WP_cruise = cruisspeed(etap=variables.etap, rho=variables.rho, rho0=variables.rho, CD0=variables.CD0clean, V=variables.V, A=variables.A[0], e=variables.e, WS=WS_plot)
plt.plot(WS_plot, WP_cruise, label = f"Cruise, A = {variables.A[0]}", color='darkorchid')

WP_cruise = cruisspeed(etap=variables.etap, rho=variables.rho, rho0=variables.rho, CD0=variables.CD0clean, V=variables.V, A=variables.A[1], e=variables.e, WS=WS_plot)
plt.plot(WS_plot, WP_cruise, label = f"Cruise, A = {variables.A[1]}", color='mediumorchid')

WP_cruise = cruisspeed(etap=variables.etap, rho=variables.rho, rho0=variables.rho, CD0=variables.CD0clean, V=variables.V, A=variables.A[2], e=variables.e, WS=WS_plot)
plt.plot(WS_plot, WP_cruise, label = f"Cruise, A = {variables.A[2]}", color='plum')

# climbrate
WP_climbrate = climbrate(etap=variables.etap, rho=variables.rho, A=variables.A[0], e=variables.e, CD0=variables.CD0to, c=variables.c, WS=WS_plot)
plt.plot(WS_plot, WP_climbrate, label = f"Climb rate, A = {variables.A[0]}", color='mediumblue')

WP_climbrate = climbrate(etap=variables.etap, rho=variables.rho, A=variables.A[1], e=variables.e, CD0=variables.CD0to, c=variables.c, WS=WS_plot)
plt.plot(WS_plot, WP_climbrate, label = f"Climb rate, A = {variables.A[1]}", color='royalblue')

WP_climbrate = climbrate(etap=variables.etap, rho=variables.rho, A=variables.A[2], e=variables.e, CD0=variables.CD0to, c=variables.c, WS=WS_plot)
plt.plot(WS_plot, WP_climbrate, label = f"Climb rate, A = {variables.A[2]}", color='cornflowerblue')

# climbgrad
WP_climbgrad = climbgradient(etap=variables.etap, cV=variables.c/variables.V, CD=variables.CDclimb[0], CL=variables.CLclimb, rho=variables.rho, WS=WS_plot)
plt.plot(WS_plot, WP_climbgrad, label = f"Climb gradient, A = {variables.A[0]}", color='darkorange')

WP_climbgrad = climbgradient(etap=variables.etap, cV=variables.c/variables.V, CD=variables.CDclimb[1], CL=variables.CLclimb, rho=variables.rho, WS=WS_plot)
plt.plot(WS_plot, WP_climbgrad, label = f"Climb gradient, A = {variables.A[1]}", color='orange')

WP_climbgrad = climbgradient(etap=variables.etap, cV=variables.c/variables.V, CD=variables.CDclimb[2], CL=variables.CLclimb, rho=variables.rho, WS=WS_plot)
plt.plot(WS_plot, WP_climbgrad, label = f"Climb gradient, A = {variables.A[2]}", color='gold')

plt.ylim(0, 0.4)
plt.xlim(0, 1500)
plt.legend(loc=1)
plt.show()