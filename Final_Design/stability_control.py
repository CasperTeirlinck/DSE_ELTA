'''
This script calculates the minimum horizontal tail size, based on the center of gravity range
Author: Bob
'''

import numpy as np
from math import pi,sqrt,tan,cos
import matplotlib.pyplot as plt

'''
scissor_plot()

Inputs:
    xcg_min [float]:    Minimum x cg location [%MAC]
    xcg_max [float]:    Maximum x center of gravity location [%MAC]

Outputs:
    ShS_min [float]:    Minimum required horizontal tail surface for stability and controllability [m2]
'''

def scissor_plot(xcg_min,xcg_max,plot=True):
    # Input parameters TODO Check all input values
    R = 287.05      # [J/kg K]  Gas constant
    gamma = 1.4     # [-]       Heat capacity ratio
    T0 = 288.15     # [K]       Base temperature
    lmbda = -0.0065 # [degC/m]  Lapse rate

    bf = 1.6        # [m]       Fuselage width
    hf = 2          # [m]       Fuselage height
    lf = 6          # [m]       Fuselage length

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
    Cm0af = 0.05    # [-]       Airfoil zero lift pitching moment coefficient TODO Check this value
    CL0 = 1         # [-]       Flapped wing lift coefficient at zero angle of attack TODO Check this value
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
    xacn = 0 # TODO Check this
    xac = xacwf + xacn

    # Wing downwash gradient
    r = lh*2/bw
    mtv = 2/bw * ((hh-hw)+lh*tan(twistwr)) * cos(twistwr)
    KeLambda = (0.1124 + 0.1265*sweepw + sweepw**2)/(r**2) + 0.1025/r + 2
    KeLambda0 = 0.1124/(r**2) + 0.1024/r + 2
    deda = KeLambda/KeLambda0 * (r/(r**2 + mtv**2)*0.4876/sqrt(r**2+0.6319+mtv**2)+
                                 (1+(r**2/(r**2+0.7915+5.0734*mtv**2))**0.3113)
                                 *(1-sqrt(mtv**2/(1+mtv**2)))) * CLaw/(pi*Aw)

    # Pitching moment coefficient
    Cmacw = Cm0af*(Aw*cos(sweepw)**2/(Aw + 2*cos(sweepw)))
    dfusCmac = -1.8*(1-2.5*bf/lf) * pi*bf*hf*lf/(4*Sw*MAC) * CL0/CLaA_h
    dfCmac = 0 # TODO Implement this
    dnacCmac = 0 # TODO Check this
    Cmac = Cmacw + dfCmac + dfusCmac + dnacCmac

    # Stability Analysis
    def stability(xcg):
        a = 1/(CLah/CLaA_h * (1-deda) * lh/MAC * VhV**2)
        b = -(xac-sm)/(CLah/CLaA_h * (1-deda) * lh/MAC * VhV**2)
        return a*xcg + b

    ShS_stability = stability(xcg_max)

    # Controllability Analysis
    def control(xcg):
        a = 1/(CLh/CLA_h * lh/MAC * VhV**2)
        b = (Cmac/CLA_h - xac - sm)/(CLh/CLA_h * lh/MAC * VhV**2)
        return a*xcg + b

    ShS_Control = control(xcg_min)

    # Horizontal tail surface
    ShS_min = max(ShS_stability,ShS_Control)
    print(ShS_min)

    if plot:
        # Create Scissor Plot
        def stability_curve(ShS):
            return xac + CLah/CLaA_h * (1-deda) * ShS*lh/MAC * VhV**2 - sm

        def control_curve(Shs):
            return xac - Cmac/CLA_h + CLh/CLA_h * ShS*lh/MAC * VhV**2 + sm

        ShS = np.arange(0,1.001,0.001)

        xcg_stability = stability_curve(ShS)
        xcg_control = control_curve(ShS)

        plt.plot(xcg_stability*100,ShS,label='Stability curve')
        plt.plot((xcg_stability+sm)*100,ShS,'k--')
        plt.plot(xcg_control*100,ShS,label='Control curve')
        plt.plot((xcg_control-sm)*100,ShS,'k--')
        plt.plot([xcg_min*100,xcg_max*100],[ShS_min,ShS_min],'r',label='Center of Gravity range')
        plt.title('Scissor Plot')
        plt.xlabel('$x_{cg}$/MAC (%)')
        plt.ylabel('$S_h/S$ [-]')
        plt.legend()
        plt.grid()
        plt.show()

    else:
        pass

    return ShS_min

# Test
if __name__ ==  "__main__":
    xcg_min = 0.3
    xcg_max = 1.5

    ShS_min = scissor_plot(xcg_min,xcg_max,plot=True)