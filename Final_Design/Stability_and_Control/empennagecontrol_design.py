from math import pi

def elevator_sizing(variables):
    g = 9.80665

    a_pitch = 12*pi/180         # [rad/s]       Take-off pitch angular velocity
    VTO = 1.05*25.2             # [m/s]         Take-off velocity
    rhoTO = 1.225               # [kg/m3]       Take-off density
    mu = 0.05                   # [-]           Friction factor

    la = 9                      # [m]           Aircraft length

    WTO = 750*g                 # [N]           Take-off weight
    TTO = 1500                  # [N]           Take-off thrust

    xcg = 1.5                   # [m]           Most forward centre of gravity
    xmg = 2                     # [m]           Main gear location
    xacw = 1                    # [m]           Wing/fuselage aerodynamic centre location
    xach = 6                    # [m]           Horizontal tail aerodynamic centre location

    zcg = 1                     # [m]           Centre of gravity height
    zmg = 0                     # [m]           Main gear height
    zT = 1                      # [m]           Thrust vector height
    zD = 0.5                    # [m]           Drag vector height

    Sw = 21.09                  # [m2]          Wing surface area
    MAC = 1.47                  # [m]           Mean aerodynamic chord
    CLTO = 0.628                # [-]           Take-off lift coefficient
    CDTO = 0.0414               # [-]           Take-off drag coefficient
    Cmacwf = -0.565             # [-]           Wing-fuselage pitching moment coefficient around the aerodynamic centre
    deda = 0.0689               # [-]           Downwash gradient
    a0 = -7.255*pi/180          # [rad]         Zero lift angle of attack

    Sh = 10                     # [m2]          Horizontal tail surface area
    ih = 0                      # [rad]         Horizontal tail incidence angle
    bh = 5                      # [m]           Horizontal tail surface area
    chr = 1                     # [m]           Horizontal tail root chord
    CLah = 4                    # [/rad]        Horizontal lift curve slope

    bebh = 1                    # [-]           Elevator span
    de_max = 20*pi/180          # [rad]         Maximum elevator deflection
    de_min = -25*pi/180         # [rad]         Minimum elevator deflection (maximum negative)

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
    MLwf = Lwf*(xmg-xacw)
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

    # Horizontal tail angle of attack
    ah = ih - a0*deda

    # Elevator angle of attack effectiveness
    tau_e = (CLh/CLah-ah)/de_min
    print('Elevator angle of attack effectiveness =',round(tau_e,2))

    # Elevator-to-tail chord ratio
    cech = float(input('\nCE/Ch = '))

    # Elevator size
    be = bebh*bh
    ce = cech*chr
    Se = be*ce

    # Update in variables class
    variables.be = be
    variables.ce = ce
    variables.Se = Se

    return variables

# Test
class Test_variables_sc:
    def __init__(self):
        self.Se = None      # [m2]      Elevator surface area
        self.be = None      # [m]       Elevator span
        self.ce = None      # [m]       Elevator chord

if __name__ ==  "__main__":
    test_v = Test_variables_sc()

    test_v = elevator_sizing(test_v)
    print('\nSe =',round(test_v.Se,2),'m2')
    print('be =',round(test_v.be,2),'m')
    print('ce =',round(test_v.ce,2),'m')