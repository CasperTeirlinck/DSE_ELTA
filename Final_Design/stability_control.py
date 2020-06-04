'''
This script calculates the minimum horizontal tail size, based on the center of gravity range
Author: Bob
'''

import numpy as np
from math import pi,sqrt,tan,cos,atan
import matplotlib.pyplot as plt


'''
XLEMAC() :              Calculates the distance nose - leading edge of the MAC

Inputs:
    lfn [float]:        Distance nose - wing leading edge root chord [m]
    sweep [float]:      Quarter chord sweep angle [rad]
    cr [float]:         Root chord [m]
    b [float]:          Span [m]
    taper [float]:      Taper ratio [-]

Outputs:
    xmac [float]:       Distance LE root chord - LE of the MAC [m]
    xlemac [float]:     Distance nose - leading edge of the MAC [m]
'''

def XLEMAC(lfn,sweep,cr,b,taper):
    sweepLE = atan(tan(sweep) - cr/(2*b)*(taper-1))
    ymac = (b/6)*((1+2*taper)/(1+taper))
    xmac = ymac*tan(sweepLE)
    xlemac = lfn + xmac

    return xlemac,xmac


'''
percMAC() :             Calculates the center of gravity as a fraction of the MAC

Inputs:
    xcg [float/array]:  Center of gravity location [m]
    MAC [float]:        Mean aerodynamic chord [m]
    xlemac [float]:     Distance nose - leading edge mean aerodynamic chord [m]
    percentage [bool]:  Return percentage (True) or fraction (False)

Outputs:
    [float]:            Center of gravity location [%MAC]
'''

def percMAC(xcg,MAC,xlemac,percentage=False):
    if percentage:
        return (xcg-xlemac)/MAC * 100
    else:
        return (xcg-xlemac)/MAC


'''
loading_diagram() :     Creates the loading diagram

Inputs:
    variables [class]:  
    xcg_wing [float]:   Wing center of gravity location [m]
    xlemac [float]:     Distance nose - leading edge Mean Aerodynamic Chord [m]
    plot [bool]:        Create a plot (True or False), default is False

Outputs:
    xcg_min [float]:    Minimum cg location [%MAC]
    xcg_max [float]:    Maximum cg location [%MAC]
'''

def loading_diagram(variables,xcg_wing,xlemac,plot=False):
    # Inputs TODO Get values out of variables class
    MAC = 1.47                          # [m]       Mean Aerodynamic Chord

    m_wing = 100                        # [kg]      Wing mass
    m_fgroup = 300                      # [kg]      Fuselage group mass
    xcg_fgroup = 3                      # [m]       Fuselage group center of gravity

    m_oe = m_wing + m_fgroup            # [kg]      Operational Empty mass
    xcg_oew = (xcg_wing*m_wing + xcg_fgroup*m_fgroup)/m_oe
                                        # [m]       OEW center of gravity

    m_bat = 150                         # [kg]      Battery mass
    xcg_bat = 3                         # [m]       Battery center of gravity
    m_repbat = 100                      # [kg]      Replaceable battery mass
    xcg_repbat = 3                      # [m]       Replaceable battery center of gravity
    m_pax = 100                         # [kg]      Pilot/Passenger mass
    xcg_pax = 1                         # [m]       Pilot/Passenger center of gravity

    sm = 0.05                           # [-]       Safety Margin

    # Create mass and center of gravity lists
    m_lst = [m_oe+m_bat]
    xcg_lst = [(xcg_oew*m_oe + xcg_bat*m_bat)/m_lst[0]]

    # Battery loading
    def battery_loading(initial_m,initial_xcg,m_repbat,xcg_repbat):
        m_new = initial_m + m_repbat
        xm_new = initial_xcg*initial_m + xcg_repbat*m_repbat
        xcg_new = xm_new/m_new

        return m_new,xcg_new

    # Pilot loading
    def pax_loading(initial_m,initial_xcg,m_pax,xcg_pax):
        m_new = initial_m + m_pax
        xm_new = initial_xcg*initial_m + xcg_pax*m_pax
        xcg_new = xm_new/m_new

        return m_new,xcg_new


    # 1 Pilot + Replaceable Battery
    # First Battery, then Pilot
    # Battery loading
    m_new,xcg_new = battery_loading(m_lst[-1],xcg_lst[-1],m_repbat,xcg_repbat)

    # Create loading lists
    bp_m_repbat = [m_lst[-1]]
    bp_xcg_repbat = [xcg_lst[-1]]

    # Append new mass and center of gravity
    m_lst.append(m_new)
    bp_m_repbat.append(m_new)
    xcg_lst.append(xcg_new)
    bp_xcg_repbat.append(xcg_new)

    # Pilot loading
    m_new,xcg_new = pax_loading(m_lst[-1],xcg_lst[-1],m_pax,xcg_pax)

    # Create loading lists
    bp_m_pax = [m_lst[-1]]
    bp_xcg_pax = [xcg_lst[-1]]

    # Append new mass and center of gravity
    m_lst.append(m_new)
    bp_m_pax.append(m_new)
    xcg_lst.append(xcg_new)
    bp_xcg_pax.append((xcg_new))

    # First Pilot, then Battery
    # Pilot loading
    m_new,xcg_new = pax_loading(m_lst[0],xcg_lst[0],m_pax,xcg_pax)

    # Create loading lists
    pb_m_pax = [m_lst[0]]
    pb_xcg_pax = [xcg_lst[0]]

    # Append new mass and center of gravity
    m_lst.append(m_new)
    pb_m_pax.append(m_new)
    xcg_lst.append(xcg_new)
    pb_xcg_pax.append(xcg_new)

    # Battery loading
    m_new,xcg_new = battery_loading(m_lst[-1],xcg_lst[-1],m_repbat,xcg_repbat)

    # Create loading lists
    pb_m_repbat = [m_lst[-1]]
    pb_xcg_repbat = [xcg_lst[-1]]

    # Append new mass and center of gravity
    m_lst.append(m_new)
    pb_m_repbat.append(m_new)
    xcg_lst.append(xcg_new)
    pb_xcg_repbat.append(xcg_new)


    # 2 Pilots
    # Calculate mass of 2 pilots
    m_pax = 2*m_pax

    # Pilot loading
    m_new,xcg_new = pax_loading(m_lst[0],xcg_lst[0],m_pax,xcg_pax)

    # Create loading lists
    p_m = [m_lst[0]]
    p_xcg = [xcg_lst[0]]

    # Append new mass and center of gravity
    m_lst.append(m_new)
    p_m.append(m_new)
    xcg_lst.append(xcg_new)
    p_xcg.append(xcg_new)


    # Transform to numpy arrays and calculate fraction of MAC
    m_lst = np.array(m_lst)
    bp_m_repbat = np.array(bp_m_repbat)
    bp_m_pax = np.array(bp_m_pax)
    pb_m_repbat = np.array(pb_m_repbat)
    pb_m_pax = np.array(pb_m_pax)
    p_m = np.array(p_m)

    xcg_lst = percMAC(np.array(xcg_lst),MAC,xlemac)
    bp_xcg_repbat = percMAC(np.array(bp_xcg_repbat),MAC,xlemac)
    bp_xcg_pax = percMAC(np.array(bp_xcg_pax),MAC,xlemac)
    pb_xcg_repbat = percMAC(np.array(pb_xcg_repbat),MAC,xlemac)
    pb_xcg_pax = percMAC(np.array(pb_xcg_pax),MAC,xlemac)
    p_xcg = percMAC(np.array(p_xcg),MAC,xlemac)

    # Get Maximum and minimum center of gravity location and apply safety margin
    xcg_min = min(xcg_lst) - sm
    xcg_max = max(xcg_lst) + sm

    # Create loading diagram
    if plot:
        plt.plot(bp_xcg_repbat*100,bp_m_repbat,color='#1f77b4',label='Pilot + Battery: Replaceable Battery Loading')
        plt.plot(pb_xcg_repbat*100,pb_m_repbat,color='#1f77b4')
        plt.plot(bp_xcg_pax*100,bp_m_pax,color='#ff7f0e',label='Pilot + Battery: Pilot Loading')
        plt.plot(pb_xcg_pax*100,pb_m_pax,color='#ff7f0e')
        plt.plot(p_xcg*100,p_m,color='#2ca02c',label='2 Pilots: Pilot Loading')
        plt.plot([xcg_min*100,xcg_min*100],[min(m_lst),max(m_lst)],'r',label='Minimum/Maximum Center of Gravity')
        plt.plot([xcg_max*100,xcg_max*100],[min(m_lst),max(m_lst)],'r')
        plt.title('Loading Diagram')
        plt.xlabel('$x_{cg}$/MAC (%)')
        plt.ylabel('Weight [kg]')
        plt.legend()
        plt.grid()
        plt.show()

    else:
        pass

    return xcg_min,xcg_max


'''
scissor_plot()          Creates the scissor plot

Inputs:
    variables[class]:   
    lfn [float]         Distance nose - wing leading edge root chord [m]
    xcg_min [float]:    Minimum center of gravity location [%MAC]
    xcg_max [float]:    Maximum center of gravity location [%MAC]
    plot [bool]:        Create a plot (True or False), default is False

Outputs:
    ShS_min [float]:    Minimum required horizontal tail surface [-]
'''

def scissor_plot(variables,lfn,xcg_min,xcg_max,plot=False):
    # Input parameters TODO Get values out of variables class
    R = 287.05                  # [J/kg K]  Gas constant
    gamma = 1.4                 # [-]       Heat capacity ratio
    T0 = 288.15                 # [K]       Base temperature
    lmbda = -0.0065             # [degC/m]  Lapse rate

    bf = 1.6                    # [m]       Fuselage width
    hf = 2                      # [m]       Fuselage height
    lf = 6                      # [m]       Fuselage length

    hw = 0.5                    # [m]       Height of the wing, from ground
    MAC = 1.47                  # [m]       Mean Aerodynamic Chord
    Sw = 23                     # [m2]      Horizontal tail surface area
    Snet = 20                   # [m2]      Net wing surface area
    bw = 16.6                   # [m]       Wing span
    Aw = 12                     # [-]       Wing aspect ratio
    sweepw = 0                  # [rad]     Wing quarter chord sweep angle
    taperw = 0.4                # [-]       Wing taper ratio
    twistwr = 0                 # [deg]     Wing twist at the root

    lh = 6.2                    # [m]       Tail arm
    hh = 1.5                    # [m]       Height horizontal tail from ground
    Ah = 3                      # [-]       Horizontal tail aspect ratio
    sweeph = 0                  # [rad]     Horizontal tail half chord sweep

    VhV = 0.85                  # [-]       Tail/wing speed ratio

    Vcruise = 50                # [m/s]     Cruise speed
    hcruise = 914.4             # [m]       Cruise altitude

    xacw = 0.25                 # [%MAC]    Wing location aerodynamic center TODO Check if aerodynamics guys determine this value
    eta = 0.95                  # [-]       Airfoil efficiency coefficient TODO Check if aerodynamics guys determine this value
    CLaw = 5.4                  # [/rad]    Wing lift rate coefficient TODO Check if aerodynamics guys determine this value
    Cm0af = 0.05                # [-]       Airfoil zero lift pitching moment coefficient TODO Check this value
    CL0 = 1                     # [-]       Flapped wing lift coefficient at zero angle of attack TODO Check this value
    CLA_h = 3                   # [-]       Aircraft less tail lift coefficient TODO Check this value
    CLh = -0.5                  # [-]       Horizontal tail lift coefficient TODO Check this value

    sm = 0.1                    # [-]       Safety margin

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

    # Minimum required horizontal tail surface
    ShS_min = max(ShS_stability,ShS_Control)

    # Create Scissor Plot
    if plot:
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


'''
sizing_htail_wingpos()  Sizing of the horizontal tail and wing position

Inputs:
    variables[class]:   
    plot [bool]:        Create a plot (True or False), default is False

Outputs:
    variables [class]:  The class contains updated values of:
                         - lfn [float]:         Distance nose - wing leading edge root chord [m]
                         - xcg_wing [float]:    Wing center of gravity location [m]
                         - xlemacw [float]:     Distance nose - leading edge Mean Aerodynamic Chord [m]
                         - xcg_min [float]:     Minimum center of gravity location [%MAC]
                         - xcg_max [float]:     Maximum center of gravity location [%MAC]
                         - ShS_min [float]:     Minimum required horizontal tail surface [-]
'''

def sizing_htail_wingpos(variables,plot=False):
    # Inputs
    lf = 6                          # [m]   Fuselage length
    lfn = 1.5                       # [m]   Distance nose - wing
    sweepw = 0                      # [rad] Wing quarter chord sweep angle
    taperw = 0.4                    # [-]   Wing taper ratio
    bw = 16.6                       # [m]   Wing span
    crw = 1.98                      # [m]   Wing root chord
    MAC = 1.47                      # [m]   Mean Aerodynamic Chord
    xcg_wing = 1.5                  # [m]   Wing center of gravity

    # Parameter calculations
    xlemacw,xmacw = XLEMAC(lfn,sweepw,crw,bw,taperw)
    pxcg_wing = percMAC(xcg_wing,MAC,xlemacw)

    # Create lists
    xlemaclf_lst = np.arange(0,1.001,0.001)
    xcgmin_lst = []
    xcgmax_lst = []
    ShSmin_lst = []

    # Perform wing shift
    for xlemaclf in xlemaclf_lst:
        xlemacw_i = xlemaclf*lf
        lfn_i = xlemacw_i - xmacw
        xcg_wing_i = lfn_i + pxcg_wing*MAC

        xcg_min_i,xcg_max_i = loading_diagram(variables,xcg_wing_i,xlemacw_i)
        ShS_min_i = scissor_plot(variables,lfn_i,xcg_min_i,xcg_max_i)

        xcgmin_lst.append(xcg_min_i)
        xcgmax_lst.append(xcg_max_i)
        ShSmin_lst.append(ShS_min_i)

    # Determine minimum horizontal tail surface
    ShS_min = min(ShSmin_lst)
    # Get index of minimum horizontal tail surface
    i = ShSmin_lst.index(ShS_min)

    # Determine center of gravity range for minimum horizontal tail surface
    xcg_min = xcgmin_lst[i]
    xcg_max = xcgmax_lst[i]

    # Determine wing location parameters for minimum horizontal tail surface
    xlemaclf = xlemaclf_lst[i]
    xlemacw = xlemaclf*lf
    lfn = xlemacw - xmacw
    xcg_wing = lfn + pxcg_wing*MAC

    # Update values in variables class
    #variables.lfn = lfn
    #variables.xcg_wing = xcg_wing
    #variables.xlemacw = xlemacw
    variables.xcg_min = xcg_min
    variables.xcg_max = xcg_max
    variables.ShS_min = ShS_min

    # Create plots
    if plot:
        # Transform lists to numpy arrays
        xlemaclf_lst = np.array(xlemaclf_lst)
        xcgmin_lst = np.array(xcgmin_lst)
        xcgmax_lst = np.array(xcgmax_lst)

        # Loading diagram
        loading_diagram(variables,xcg_wing,xlemacw,plot)

        # Scissor plot
        scissor_plot(variables,lfn,xcg_min,xcg_max,plot)

        # Center of gravity range plot
        plt.plot(xcgmin_lst*100, xlemaclf_lst*100,label='Most forward center of gravity')
        plt.plot(xcgmax_lst*100, xlemaclf_lst*100,label='Most aft center of gravity')
        plt.plot([xcg_min*100,xcg_max*100],[xlemaclf*100,xlemaclf*100],'r',label='Center of Gravity range')
        plt.title('Center of gravity range')
        plt.xlabel('$x_{cg}$/MAC (%)')
        plt.ylabel('$x_{lemac}/l_{fus}$ (%)')
        plt.legend()
        plt.grid()
        plt.show()
    else:
        pass

    return variables


# Test
class Test_variables_sc:
    def __init__(self):
        self.xcg_min = None     # [%MAC]    Minimum center of gravity location
        self.xcg_max = None     # [%MAC]    Maximum center of gravity location
        self.ShS_min = None     # [m2]      Minimum required horizontal tail surface

if __name__ ==  "__main__":
    payload = 'Battery'
    test_v = Test_variables_sc()

    test_v = sizing_htail_wingpos(test_v,plot=True)
    print('xcg_min =',test_v.xcg_min)
    print('xcg_max =',test_v.xcg_max)
    print('Sh/S =',test_v.ShS_min)