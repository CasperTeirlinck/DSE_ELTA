import numpy as np
from math import pi,sqrt,tan
import matplotlib.pyplot as plt

# Input parameters
lh = 6.2        # [m]       Tail arm
bf = 1.6        # [m]       Fuselage width
hf = 2          # [m]       Fuselage height
lfn = 2         # [m]       Distance nose - wing
MAC = 1.47      # [m]       Mean Aerodynamic Chord
Sw = 23         # [m2]      Horizontal tail surface area
Snet = 20       # [m2]      Net wing surface area
bw = 16.6       # [m]       Wing span
sweepw = 0      # [rad]     Wing quarter chord sweep angle
taperw = 0.4    # [-]       Wing taper ratio
Ah = 3          # [-]       Horizontal tail aspect ratio
sweeph = 0      # [rad]     Horizontal tail half chord sweep

deda = 0.3      # [-]       Downwash effect
VhV = 0.85      # [-]       Tail/wing speed ratio

Vcruise = 50    # [m/s]     Cruise speed
a = 300         # [m/s]     Speed of sound at cruise altitude

xacw = 0.6      # [-]
eta = 0.95      # [-]       Airfoil efficiency coefficient
CLaw = 0.6      # [/rad]    Wing lift rate coefficient
Cmac = -0.05    # [-]       Aircraft less tail pitching moment coefficient
CLA_h = 0.5     # [-]       Aircraft less tail lift coefficient
CLh = 1         # [-]       Horizontal tail lift coefficient

sm = 0.1        # [-]       Safety margin

# Parameter calculations

# Tail lift rate coefficient
Vh = VhV*Vcruise
Mh = Vh/a
beta = sqrt(1 - Mh**2)
CLah = 2*pi*Ah/(2 + sqrt(4 + (Ah*beta/eta)**2 * (1 + tan(sweeph)**2/(beta**2))))

# Lift rate coefficient of the aircraft less tail
CLaA_h = CLaw * (1 + 2.15*bf/bw) * Snet/Sw + pi/2*bf**2/Sw

# Aerodynamic center of the aircraft less tail
cg = Sw/bw
xacf2 = 0.273/(1 + taperw) * bf*cg*(bw-bf)/(MAC**2*(bw + 2.15*bf)) * tan(sweepw)
xacf1 = -1.8/CLaA_h * bf*hf*lfn/(Sw*MAC)
xacwf = xacw + xacf1 + xacf2
xacn = 1
xac = xacwf + xacn

# Horizontal tail surface area range
ShS = np.arange(0,1.001,0.001)

# Stability Analysis
def stability_curve(ShS):
    return xac/MAC + CLah/CLaA_h * (1-deda) * ShS*lh/MAC * VhV**2 - sm

xcg_stability = stability_curve(ShS)

# Controllability Analysis
def control_curve(ShS):
    return xac/MAC - Cmac/CLA_h + CLh/CLA_h * ShS*lh/MAC * VhV**2

xcg_control = control_curve(ShS)

# Create Scissor Plot
plt.plot(xcg_stability,ShS)
plt.plot(xcg_control,ShS)
plt.show()