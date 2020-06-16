from math import pi

g = 9.80665

a_pitch = 12.5*pi/180       # [deg/s]       Take-off pitch angular velocity
bebh = 1                    # [-]           Elevator span
de_max = 20*pi/180          # [rad]         Maximum elevator deflection
de_min = -25*pi/180         # [rad]         Minimum elevator deflection (maximum negative)
Sw = 21.09                  # [m2]          Wing surface area
VTO = 1.05*25.2             # [m/s]         Take-off velocity
rhoTO = 1.225               # [kg/m3]       Take-off density
MAC = 1.47                  # [m]           Mean aerodynamic chord
CLTO = 1.5                  # [-]           Take-off lift coefficient
CDTO = 0.5                  # [-]           Take-off drag coefficient
Cmacwf = 1                  # [-]           Wing-fuselage pitching moment coefficient around the aerodynamic centre
TTO = 100                   # [N]           Take-off thrust
WTO = 750*g                 # [N]           Take-off weight
xcg = 1.5                   # [m]           Most forward centre of gravity
zcg = 1                     # [m]           Centre of gravity height
xmg = 2                     # [m]           Main gear location
xacwf = 1                   # [m]           Wing/fuselage aerodynamic centre location
zmg = 0                     # [m]           Main gear height
zD = 0.5                    # [m]           Drag vector height
zT = 1                      # [m]           Thrust vector height
zwac = 0.5                  # [m]           Aerodynamic centre height
mu = 0.3                    # [-]           Friction factor
la = 9                      # [m]           Aircraft length
Sh = 10                     # [m2]          Horizontal tail surface area
xach = 6                    # [m]           Horizontal tail aerodynamic centre location

# Wing/fuselage lift, drag and pitching moment
Lwf = CLTO*0.5*rhoTO*VTO**2*Sw
DTO = CDTO*0.5*rhoTO*VTO**2*Sw
Macwf = Cmacwf*0.5*rhoTO*VTO**2*Sw*MAC

# Aircraft linear acceleration
LTO = Lwf
N = WTO-LTO
ma = (TTO-DTO-mu*N)

# Contributing pitching moments during take-off rotation
Mw = WTO*(xmg-xcg)
MD = DTO*(zD-zmg)
MT = TTO*(zT-zmg)
MLwf = Lwf*(xmg-xacwf)
Ma = ma*(zcg-zmg)

# Aircraft pitch moment of inertia
C = 0.176
ky = C*la
Iyy = (WTO/g)*ky**2
Iyymg = Iyy + (zcg-zmg)**2*(WTO/g)

# Desired horizontal tail lift
Lh = (MLwf+Macwf+Ma-Mw+MD-MT-Iyymg*a_pitch)/(xach-xmg)

# Desired horizontal tail lift coefficient
CLh = 2*Lh/(rhoTO*VTO**2*Sh)

# Elevator angle of attack effectiveness
tau_e = (CLh/CLah-ah)/de_max
