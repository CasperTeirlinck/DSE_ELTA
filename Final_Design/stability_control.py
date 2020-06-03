'''
This script calculates the minimum horizontal tail size, based on the center of gravity range
Author: Bob
'''

import numpy as np
from math import pi,sqrt,tan,cos
import matplotlib.pyplot as plt

def percMAC(xcg,MAC,xlemac,percentage=False):
    if percentage:
        return (xcg-xlemac)/MAC * 100
    else:
        return (xcg-xlemac)/MAC


payload = 'Battery'
payload_lst = ['Pilot','Battery']
plot = True

if payload in payload_lst:
    # Inputs TODO Get values out of variables class
    MAC = 1.47          # [m]       Mean Aerodynamic Chord
    xlemac = 1.5        # [m]       Distance nose, leading edge mean aerodynamic chord

    m_oe = 400          # [kg]      Operational Empty mass
    xcg_oew = 2         # [m]       OEW center of gravity
    m_bat = 150         # [kg]      Battery mass
    xcg_bat = 3         # [m]       Battery center of gravity
    m_repbat = 100      # [kg]      Replaceable battery mass
    xcg_repbat = 3      # [m]       Replaceable battery center of gravity
    m_pax = 100         # [kg]      Pilot/Passenger mass
    xcg_pax = 1         # [m]       Pilot/Passenger center of gravity

    # Payload mass and center of gravity
    if payload == 'Pilot':
        m_repbat = 0    # [kg]      Replaceable battery mass
        m_pax = 2*m_pax # [kg]      Passenger/Pilot mass

    elif payload == 'Battery':
        m_repbat = 100    # [kg]    Replaceable battery mass

    # Create mass and center of gravity lists
    m_lst = [m_oe+m_bat]
    xcg_lst = [(xcg_oew*m_oe + xcg_bat*m_bat)/m_lst[0]]

    # Create loading diagram
    m_repbat_lst = []
    xcg_repbat_lst = []
    m_pax_lst = []
    xcg_pax_lst = []

    # Battery loading
    def battery_loading(m,xcg):
        m_new = m + m_repbat
        xm_new = xcg*m + xcg_repbat*m_repbat
        xcg_new = xm_new/m_new

        return m_new,xcg_new

    # Passenger/Pilot loading
    def pax_loading(m,xcg):
        m_new = m + m_pax
        xm_new = xcg*m + xcg_pax*m_pax
        xcg_new = xm_new/m_new

        return m_new,xcg_new

    # First Battery, then Pilot
    m_new,xcg_new = battery_loading(m_lst[-1],xcg_lst[-1])

    bp_m_repbat = [m_lst[-1]]
    bp_xcg_repbat = [xcg_lst[-1]]

    m_lst.append(m_new)
    bp_m_repbat.append(m_new)
    xcg_lst.append(xcg_new)
    bp_xcg_repbat.append(xcg_new)

    m_new,xcg_new = pax_loading(m_lst[-1],xcg_lst[-1])

    bp_m_pax = [m_lst[-1]]
    bp_xcg_pax = [xcg_lst[-1]]

    m_lst.append(m_new)
    bp_m_pax.append(m_new)
    xcg_lst.append(xcg_new)
    bp_xcg_pax.append((xcg_new))

    # First Pilot, then Battery
    m_new,xcg_new = pax_loading(m_lst[0],xcg_lst[0])

    pb_m_pax = [m_lst[0]]
    pb_xcg_pax = [xcg_lst[0]]

    m_lst.append(m_new)
    pb_m_pax.append(m_new)
    xcg_lst.append(xcg_new)
    pb_xcg_pax.append(xcg_new)

    m_new,xcg_new = battery_loading(m_lst[-1],xcg_lst[-1])

    pb_m_repbat = [m_lst[-1]]
    pb_xcg_repbat = [xcg_lst[-1]]

    m_lst.append(m_new)
    pb_m_repbat.append(m_new)
    xcg_lst.append(xcg_new)
    pb_xcg_repbat.append(xcg_new)

    # Transform to numpy arrays
    m_lst = np.array(m_lst)
    bp_m_repbat = np.array(bp_m_repbat)
    bp_m_pax = np.array(bp_m_pax)
    pb_m_repbat = np.array(pb_m_repbat)
    pb_m_pax = np.array(pb_m_pax)

    xcg_lst = percMAC(np.array(xcg_lst),MAC,xlemac)
    bp_xcg_repbat = percMAC(np.array(bp_xcg_repbat),MAC,xlemac)
    bp_xcg_pax = percMAC(np.array(bp_xcg_pax),MAC,xlemac)
    pb_xcg_repbat = percMAC(np.array(pb_xcg_repbat),MAC,xlemac)
    pb_xcg_pax = percMAC(np.array(pb_xcg_pax),MAC,xlemac)

    # Get Maximum and minimum center of gravity location
    xcg_min = min(xcg_lst)
    xcg_max = max(xcg_lst)

    # Create plot
    if plot:
        plt.plot(bp_xcg_repbat,bp_m_repbat,color='#1f77b4',label='Replaceable Battery Loading')
        plt.plot(pb_xcg_repbat,pb_m_repbat,color='#1f77b4')
        plt.plot(bp_xcg_pax,bp_m_pax,color='#ff7f0e',label='Passenger/Pilot Loading')
        plt.plot(pb_xcg_pax,pb_m_pax,color='#ff7f0e')
        plt.title('Loading Diagram')
        plt.xlabel('$x_{cg}$/MAC (%)')
        plt.ylabel('Weight [kg]')
        plt.legend()
        plt.grid()
        plt.show()

    else:
        pass

else:
    print("Incorrect payload input ('"+payload+"'). Choose 'Pilot' or 'Battery'")

'''
scissor_plot()

Inputs:
    xcg_min [float]:    Minimum x cg location [%MAC]
    xcg_max [float]:    Maximum x center of gravity location [%MAC]

Outputs:
    ShS_min [float]:    Minimum required horizontal tail surface for stability and controllability [m2]
'''

def scissor_plot(xcg_min,xcg_max,plot=True):
    # Input parameters TODO Get values out of variables class
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
'''
if __name__ ==  "__main__":
    xcg_min = 0.3
    xcg_max = 1.5

    ShS_min = scissor_plot(xcg_min,xcg_max,plot=True)
'''