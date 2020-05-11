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

def turnrate(rho, V, CD0, phi, A, e, WS):
    """
    :param rho: air density cruise, scalar
    :param V: cruise speed, scalar
    :param CD0: zero drag coefficient, scalar
    :param phi: bank angle, scalar
    :param A: aspect ratio, scalar
    :param e: oswald efficiency factor
    :param WS: plot x values, array
    :return: plot y values, array
    """
    
    return WS/( 0.5*rho*V**3*( CD0 + WS**2/( (0.5*rho*V**2*np.cos( np.radians(phi)))**2 * np.pi*A*e ) ) )

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

WP_cruise1 = cruisspeed(etap=variables.etap, rho=variables.rhocruise, rho0=variables.rho0, CD0=variables.CD0clean, V=variables.V, A=variables.A[0], e=variables.e, WS=WS_limit) # Pessimistic cruise
WP_cruise2 = cruisspeed(etap=variables.etap, rho=variables.rhocruise, rho0=variables.rho0, CD0=variables.CD0clean, V=variables.V, A=variables.A[1], e=variables.e, WS=WS_limit) # Neutral cruise
WP_cruise3 = cruisspeed(etap=variables.etap, rho=variables.rhocruise, rho0=variables.rho0, CD0=variables.CD0clean, V=variables.V, A=variables.A[2], e=variables.e, WS=WS_limit) # Optimistic cruise

WP_climbrate1 = climbrate(etap=variables.etap, rho=variables.rho, A=variables.A[0], e=variables.e, CD0=variables.CD0to, c=variables.c, WS=WS_limit) # Pessimistic climb rate
WP_climbrate2 = climbrate(etap=variables.etap, rho=variables.rho, A=variables.A[1], e=variables.e, CD0=variables.CD0to, c=variables.c, WS=WS_limit) # Neutral climb rate
WP_climbrate3 = climbrate(etap=variables.etap, rho=variables.rho, A=variables.A[2], e=variables.e, CD0=variables.CD0to, c=variables.c, WS=WS_limit) # Optimistic climb rate

WP_climbgrad1 = climbgradient(etap=variables.etap, cV=variables.c/variables.V, CD=variables.CDclimb[0], CL=variables.CLclimb, rho=variables.rho, WS=WS_limit) # Pessimistic climb gradient
WP_climbgrad2 = climbgradient(etap=variables.etap, cV=variables.c/variables.V, CD=variables.CDclimb[1], CL=variables.CLclimb, rho=variables.rho, WS=WS_limit) # Neutral climb rate
WP_climbgrad3 = climbgradient(etap=variables.etap, cV=variables.c/variables.V, CD=variables.CDclimb[2], CL=variables.CLclimb, rho=variables.rho, WS=WS_limit) # Optimistic climb rate

WP_turnrate1 = turnrate(rho=variables.rhocruise, V=variables.V, CD0=variables.CD0clean, phi=variables.phi, A=variables.A[0], e=variables.e, WS=WS_limit)
WP_turnrate2 = turnrate(rho=variables.rhocruise, V=variables.V, CD0=variables.CD0clean, phi=variables.phi, A=variables.A[1], e=variables.e, WS=WS_limit)
WP_turnrate3 = turnrate(rho=variables.rhocruise, V=variables.V, CD0=variables.CD0clean, phi=variables.phi, A=variables.A[2], e=variables.e, WS=WS_limit)

# Calculate limiting W/P values
WP_limit = min([WP_to1,WP_to2,WP_to3,WP_cruise1,WP_cruise2,WP_cruise3,WP_climbrate1,WP_climbrate2,WP_climbrate3,WP_climbgrad1,WP_climbgrad2,WP_climbgrad3,WP_turnrate1,WP_turnrate2,WP_turnrate3])


## Plotting
plt.figure(figsize=(9,4.5))

WS_plot = np.arange(100, 1500, 0.1)

# stallspeed
WS_s = stallspeed(CLmax=variables.CLmaxclean[1], Vs=variables.Vs, rho=variables.rho)
plt.plot([WS_s, WS_s], [0, 0.5], label = f"Stall, CLmax = {variables.CLmaxclean[1]}", color='c' )

# take off
WP_to_1 = takeoff(k=variables.k, CLto=variables.CLto[0], sigma=variables.sigma, WS= WS_plot)
plt.plot(WS_plot, WP_to_1, label = f"Take-off, CLto = {variables.CLto[0]}", color='firebrick')

WP_to_2 = takeoff(k=variables.k, CLto=variables.CLto[1], sigma=variables.sigma, WS= WS_plot)
plt.plot(WS_plot, WP_to_2, label = f"Take-off, CLto = {variables.CLto[1]}", color='indianred')

WP_to_3 = takeoff(k=variables.k, CLto=variables.CLto[2], sigma=variables.sigma, WS= WS_plot)
plt.plot(WS_plot, WP_to_3, label = f"Take-off, CLto = {variables.CLto[2]}", color='lightcoral')

# landing
WP_landing_1 = landing(CLmax=variables.CLmaxland[0], rho=variables.rho, sland=variables.sland, f=variables.f)
plt.plot([WP_landing_1, WP_landing_1], [0, 0.5], label = f"Landing, CL = {variables.CLmaxland[0]}", color='forestgreen')

WP_landing_2 = landing(CLmax=variables.CLmaxland[1], rho=variables.rho, sland=variables.sland, f=variables.f)
plt.plot([WP_landing_2, WP_landing_2], [0, 0.5], label = f"Landing, CL = {variables.CLmaxland[1]}", color='limegreen')

WP_landing_3 = landing(CLmax=variables.CLmaxland[2], rho=variables.rho, sland=variables.sland, f=variables.f)
plt.plot([WP_landing_3, WP_landing_3], [0, 0.5], label = f"Landing, CL = {variables.CLmaxland[2]}", color='darkgreen')

# cruise
WP_cruise_1 = cruisspeed(etap=variables.etap, rho=variables.rhocruise, rho0=variables.rho0, CD0=variables.CD0clean, V=variables.V, A=variables.A[0], e=variables.e, WS=WS_plot)
plt.plot(WS_plot, WP_cruise_1, label = f"Cruise, A = {variables.A[0]}", color='darkorchid')

WP_cruise_2 = cruisspeed(etap=variables.etap, rho=variables.rhocruise, rho0=variables.rho0, CD0=variables.CD0clean, V=variables.V, A=variables.A[1], e=variables.e, WS=WS_plot)
plt.plot(WS_plot, WP_cruise_2, label = f"Cruise, A = {variables.A[1]}", color='mediumorchid')

WP_cruise_3 = cruisspeed(etap=variables.etap, rho=variables.rhocruise, rho0=variables.rho0, CD0=variables.CD0clean, V=variables.V, A=variables.A[2], e=variables.e, WS=WS_plot)
plt.plot(WS_plot, WP_cruise_3, label = f"Cruise, A = {variables.A[2]}", color='plum')

# climbrate
WP_climbrate_1 = climbrate(etap=variables.etap, rho=variables.rho, A=variables.A[0], e=variables.e, CD0=variables.CD0to, c=variables.c, WS=WS_plot)
plt.plot(WS_plot, WP_climbrate_1, label = f"Climb rate, A = {variables.A[0]}", color='mediumblue')

WP_climbrate_2 = climbrate(etap=variables.etap, rho=variables.rho, A=variables.A[1], e=variables.e, CD0=variables.CD0to, c=variables.c, WS=WS_plot)
plt.plot(WS_plot, WP_climbrate_2, label = f"Climb rate, A = {variables.A[1]}", color='royalblue')

WP_climbrate_3 = climbrate(etap=variables.etap, rho=variables.rho, A=variables.A[2], e=variables.e, CD0=variables.CD0to, c=variables.c, WS=WS_plot)
plt.plot(WS_plot, WP_climbrate_3, label = f"Climb rate, A = {variables.A[2]}", color='cornflowerblue')

# climbgrad
WP_climbgrad_1 = climbgradient(etap=variables.etap, cV=variables.c/variables.V, CD=variables.CDclimb[0], CL=variables.CLclimb, rho=variables.rho, WS=WS_plot)
plt.plot(WS_plot, WP_climbgrad_1, label = f"Climb gradient, A = {variables.A[0]}", color='darkorange')

WP_climbgrad_2 = climbgradient(etap=variables.etap, cV=variables.c/variables.V, CD=variables.CDclimb[1], CL=variables.CLclimb, rho=variables.rho, WS=WS_plot)
plt.plot(WS_plot, WP_climbgrad_2, label = f"Climb gradient, A = {variables.A[1]}", color='orange')

WP_climbgrad_3 = climbgradient(etap=variables.etap, cV=variables.c/variables.V, CD=variables.CDclimb[2], CL=variables.CLclimb, rho=variables.rho, WS=WS_plot)
plt.plot(WS_plot, WP_climbgrad_3, label = f"Climb gradient, A = {variables.A[2]}", color='gold')

# turnrate
WP_turnrate_1 = turnrate(rho=variables.rhocruise, V=variables.V, CD0=variables.CD0clean, phi=variables.phi, A=variables.A[0], e=variables.e, WS=WS_plot)
plt.plot(WS_plot, WP_turnrate_1, label = f"Turn rate, A = {variables.A[0]}", color='red')

WP_turnrate_2 = turnrate(rho=variables.rhocruise, V=variables.V, CD0=variables.CD0clean, phi=variables.phi, A=variables.A[1], e=variables.e, WS=WS_plot)
plt.plot(WS_plot, WP_turnrate_2, label = f"Turn rate, A = {variables.A[1]}", color='orangered')

WP_turnrate_3 = turnrate(rho=variables.rhocruise, V=variables.V, CD0=variables.CD0clean, phi=variables.phi, A=variables.A[2], e=variables.e, WS=WS_plot)
plt.plot(WS_plot, WP_turnrate_3, label = f"Turn rate, A = {variables.A[2]}", color='tomato')


# plot DP
WP_1 = np.minimum(WP_to_1, WP_cruise_1)
i = (np.abs(WS_plot - WS_s)).argmin()
plt.fill_between(WS_plot[:i], WP_1[:i], 1000*(WS_plot[:i] - WS_s), color='green', alpha=0.5)

print(f'\nDesign point: W/S = {WS_limit} \t W/P = {WP_limit}\n')

plt.plot(WS_limit, WP_limit, marker='o', markersize='10', markerfacecolor='limegreen', markeredgewidth=1, markeredgecolor='black')

plt.ylim(0, 0.4)
plt.xlim(0, 1500)
plt.ylabel('W/P [N/W]')
plt.xlabel('W/S [N/m2]')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fancybox=False, shadow=False, edgecolor='black')
plt.tight_layout()
plt.show()