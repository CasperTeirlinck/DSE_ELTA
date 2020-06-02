import numpy as np
from math import pi,sqrt,tan,cos
import matplotlib.pyplot as plt

# Input parameters TODO Check all input values
R = 287.05      # [J/kg K]  Gas constant
gamma = 1.4     # [-]       Heat capacity ratio
T0 = 288.15     # [K]       Base temperature
lmbda = -0.0065 # [degC/m]  Lapse rate

bf = 1.6        # [m]       Fuselage width
hf = 2          # [m]       Fuselage height

lfn = 2         # [m]       Distance nose - wing
hw = 0.5        # [m]       Height of the wing, from ground
MAC = 1.47      # [m]       Mean Aerodynamic Chord
Sw = 23         # [m2]      Horizontal tail surface area
Snet = 20       # [m2]      Net wing surface area
bw = 16.6       # [m]       Wing span
Aw = 12         # [-]       Wing aspect ratio
sweepw = 0      # [rad]     Wing quarter chord sweep angle
taperw = 0.4    # [-]       Wing taper ratio
twistwr = 0     # [deg]     Wing twist at the root

lh = 6.2        # [m]       Tail arm
hh = 1.5        # [m]       Height horizontal tail from ground
Ah = 3          # [-]       Horizontal tail aspect ratio
sweeph = 0      # [rad]     Horizontal tail half chord sweep

VhV = 0.85      # [-]       Tail/wing speed ratio

Vcruise = 50    # [m/s]     Cruise speed
hcruise = 914.4 # [m]       Cruise altitude

xacw = 0.5      # [%MAC]    Wing location aerodynamic center TODO Check if aerodynamics guys determine this value
eta = 0.95      # [-]       Airfoil efficiency coefficient TODO Check if aerodynamics guys determine this value
CLaw = 5.4      # [/rad]    Wing lift rate coefficient TODO Check if aerodynamics guys determine this value
Cmac = 1        # [-]       Aircraft less tail pitching moment coefficient TODO Check this value
CLA_h = 3       # [-]       Aircraft less tail lift coefficient TODO Check this value
CLh = -0.5      # [-]       Horizontal tail lift coefficient TODO Check this value

sm = 0.1        # [-]       Safety margin

# Parameter calculations

# Tail lift rate coefficient
Vh = VhV*Vcruise
Tcruise = T0 + lmbda*hcruise
a = sqrt(gamma*R*Tcruise)
Mh = Vh/a
betah = sqrt(1 - Mh**2)
CLah = 2*pi*Ah/(2 + sqrt(4 + (Ah*betah/eta)**2 * (1 + tan(sweeph)**2/(betah**2))))

# Lift rate coefficient of the aircraft less tail
CLaA_h = CLaw * (1 + 2.15*bf/bw) * Snet/Sw + pi/2*bf**2/Sw

# Aerodynamic center of the aircraft less tail
cg = Sw/bw
xacf1 = -1.8/CLaA_h * bf*hf*lfn/(Sw*MAC)
xacf2 = 0.273/(1 + taperw) * bf*cg*(bw-bf)/(MAC**2*(bw + 2.15*bf)) * tan(sweepw)
xacwf = xacw + xacf1 + xacf2
xacn = 0 # TODO Check this value
xac = xacwf + xacn

# Wing downwash gradient
r = lh*2/bw
mtv = 2/bw * ((hh-hw)+lh*tan(twistwr)) * cos(twistwr)
KeLambda = (0.1124 + 0.1265*sweepw + sweepw**2)/(r**2) + 0.1025/r + 2
KeLambda0 = 0.1124/(r**2) + 0.1024/r + 2
deda = KeLambda/KeLambda0 * (r/(r**2 + mtv**2)*0.4876/sqrt(r**2+0.6319+mtv**2)+
                             (1+(r**2/(r**2+0.7915+5.0734*mtv**2))**0.3113)
                             *(1-sqrt(mtv**2/(1+mtv**2)))) * CLaw/(pi*Aw)

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
plt.plot(xcg_stability*100,ShS,label='Stability curve')
plt.plot((xcg_stability+sm)*100,ShS,'k--')
plt.plot(xcg_control*100,ShS,label='Control curve')
plt.plot((xcg_control-sm)*100,ShS,'k--')
plt.title('Scissor Plot')
plt.xlabel('$x_{cg}$/MAC (%)')
plt.ylabel('$S_h/S$ [-]')
plt.legend()
plt.grid()
plt.show()