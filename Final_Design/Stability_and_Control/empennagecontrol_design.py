from math import pi,cos,tan

def elevator_sizing(variables):
    g = variables.g0

    a_pitch = variables.a_pitch     # [rad/s]       Take-off pitch angular velocity
    VTO = variables.VTO             # [m/s]         Take-off velocity
    rhoTO = variables.rhoTO         # [kg/m3]       Take-off density
    mu = variables.mu               # [-]           Friction factor

    la = variables.la               # [m]           Aircraft length

    WTO = variables.WTO             # [N]           Take-off weight
    TTO = variables.TTO             # [N]           Take-off thrust

    xcg = variables.xcg_min         # [m]           Most forward centre of gravity
    xmg = variables.xmg             # [m]           Main gear location
    xacw = variables.xacw           # [m]           Wing/fuselage aerodynamic centre location
    xach = variables.xach           # [m]           Horizontal tail aerodynamic centre location

    zcg = variables.zcg             # [m]           Centre of gravity height
    zmg = variables.zmg             # [m]           Main gear height
    zT = variables.zT               # [m]           Thrust vector height
    zD = variables.zD               # [m]           Drag vector height

    Sw = variables.S                # [m2]          Wing surface area
    MAC = variables.MAC             # [m]           Mean aerodynamic chord
    CLTO = variables.CL_takeoff     # [-]           Take-off lift coefficient
    CDTO = variables.CD0to + CLTO/(pi*variables.A*variables.eflaps)
                                    # [-]           Take-off drag coefficient
    Cm0af = variables.Cm0af
    Aw = variables.A
    sweepw = 0
    bf = variables.fuselagewidth
    lf = variables.fuselagelength
    hf = variables.fuselageheigth
    mu1 = variables.mu1
    dClmax = 0.9*variables.deClmax
    CL0 = variables.CL0flap
    CLaw = variables.wing_CL_alpha
    bw = variables.b
    Snet = variables.Snet
    cc = variables.cc
    Swf = variables.flapaffectedarea
    bfl = variables.flapspan            # [m]       Flap span
    mu1 = variables.mu1                 # [-]       Flap coefficient 1
    mu2 = 1.2*(bfl/bw)+0.13             # [-]       Flap coefficient 2
    mu3 = 0.06*(bfl/bw)+0.0335          # [-]       Flap coefficient 3
    CLaA_h = CLaw * (1 + 2.15*bf/bw) * Snet/Sw + pi/2*bf**2/Sw
    Cmacw = Cm0af*(Aw*cos(sweepw)**2/(Aw + 2*cos(sweepw)))
    dfusCmac = -1.8*(1-2.5*bf/lf) * pi*bf*hf*lf/(4*Sw*MAC) * CL0/CLaA_h
    dfCmac = mu2*(-mu1*dClmax*cc-(CLTO+dClmax*(1-Swf/Sw))*(1/8)*cc*(cc-1))\
             + 0.7*Aw/(1+2/Aw)*mu3*dClmax*tan(sweepw)
    dnacCmac = 0
    Cmacwf = Cmacw + dfCmac + dfusCmac + dnacCmac

    deda = variables.deda           # [-]           Downwash gradient
    a0 = variables.a0               # [rad]         Zero lift angle of attack

    Sh = variables.Sh               # [m2]          Horizontal tail surface area
    ih = variables.ih               # [rad]         Horizontal tail incidence angle
    bh = variables.bh               # [m]           Horizontal tail surface area
    chr = variables.chr             # [m]           Horizontal tail root chord
    CLah = variables.CLah           # [/rad]        Horizontal lift curve slope

    bebh = variables.bebh           # [-]           Elevator span
    de_max = variables.de_max       # [rad]         Maximum elevator deflection
    de_min = variables.de_min       # [rad]         Minimum elevator deflection (maximum negative)

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
    print('\n----- Elevator design -----')
    print('Elevator angle of attack effectiveness =',round(tau_e,2))

    # Elevator-to-tail chord ratio
    cech = float(input('CE/Ch = '))

    # Elevator size
    be = bebh*bh
    ce = cech*chr
    Se = be*ce

    # Update in variables class
    variables.be = be
    variables.ce = ce
    variables.Se = Se
    variables.CLh_TO = CLh

    return variables
'''
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

'''