'''
This script sizes the empennage.
Author: Bob
'''

import numpy as np
from math import pi,sqrt,tan,cos
import matplotlib.pyplot as plt
from Stability_and_Control.wing_properties import XMAC,XLEMAC,sweep,downwash,percMAC
#

'''
loading_diagram() :     Creates the loading diagram

Inputs:
    variables [class]:  
    xcg_wing [float]:   Wing center of gravity location [m]
    plot [bool]:        Create a plot (True or False), default is False

Outputs:
    xcg_min [float]:    Minimum cg location [%MAC]
    xcg_max [float]:    Maximum cg location [%MAC]

V&V:    Verified
'''

def loading_diagram(variables,xcg_wing,plot=False):
    # Inputs
    g = variables.g0                                        # [m/s2]    Gravitational acceleration

    cg_wing = variables.cg_wing                             # [m]       Distance LE root chord - cg wing
    b = variables.b                                         # [m]       Wing sweep
    sweepc4 = 0                                             # [rad]     Wing quarter chord sweep angle
    taper = variables.taper                                 # [-]       Wing taper ratio
    cr = variables.c_r                                      # [m]       Wing root chord
    xmac = XMAC(b,sweepc4,taper,cr)                         # [m]       Distance LE root chord - LE MAC
    xlemac = xcg_wing-cg_wing + xmac                        # [m]       Distance nose - LE MAC
    MAC = variables.MAC                                     # [m]       Mean Aerodynamic Chord

    m_wing = variables.W_wing/g                             # [kg]      Wing mass
    m_fgroup = variables.W_fgroup/g                         # [kg]      Fuselage group mass
    xcg_fgroup = variables.xcg_fgroup                       # [m]       Fuselage group center of gravity

    m_oe = m_wing + m_fgroup                                # [kg]      Operational Empty mass
    xcg_oew = (xcg_wing*m_wing + xcg_fgroup*m_fgroup)/m_oe  # [m]       OEW center of gravity

    m_bat = variables.W_batt/g                              # [kg]      Battery mass
    xcg_bat = variables.xcgbat                              # [m]       Battery center of gravity
    m_PL = variables.WPL/g                                  # [kg]      Payload mass
    xcg_PL = variables.xcgPL                                # [m]       Payload center of gravity

    sm = 0.05                                               # [-]       Safety Margin

    # Create mass and center of gravity lists
    m_lst = [m_oe]
    xcg_lst = [xcg_oew]


    # First Battery, then Payload

    # Battery loading
    bp_b_m = [m_lst[0]]
    bp_b_xcg =[xcg_lst[0]]

    m_new = bp_b_m[-1] + m_bat
    xm_new = bp_b_xcg[-1]*bp_b_m[-1] + xcg_bat*m_bat
    xcg_new = xm_new/m_new

    m_lst.append(m_new)
    bp_b_m.append(m_new)
    xcg_lst.append(xcg_new)
    bp_b_xcg.append(xcg_new)

    # Payload loading
    bp_p_m = [bp_b_m[-1]]
    bp_p_xcg = [bp_b_xcg[-1]]

    m_new = bp_p_m[-1] + m_PL
    xm_new = bp_p_xcg[-1]*bp_p_m[-1] + xcg_PL*m_PL
    xcg_new = xm_new/m_new

    m_lst.append(m_new)
    bp_p_m.append(m_new)
    xcg_lst.append(xcg_new)
    bp_p_xcg.append(xcg_new)


    # First payload, then Battery

    # Payload loading
    pb_p_m = [m_lst[0]]
    pb_p_xcg = [xcg_lst[0]]

    m_new = pb_p_m[-1] + m_PL
    xm_new = pb_p_xcg[-1]*pb_p_m[-1] + xcg_PL*m_PL
    xcg_new = xm_new/m_new

    m_lst.append(m_new)
    pb_p_m.append(m_new)
    xcg_lst.append(xcg_new)
    pb_p_xcg.append(xcg_new)

    # Battery loading
    pb_b_m = [pb_p_m[-1]]
    pb_b_xcg =[pb_p_xcg[-1]]

    m_new = pb_b_m[-1] + m_bat
    xm_new = pb_b_xcg[-1]*pb_b_m[-1] + xcg_bat*m_bat
    xcg_new = xm_new/m_new

    m_lst.append(m_new)
    pb_b_m.append(m_new)
    xcg_lst.append(xcg_new)
    pb_b_xcg.append(xcg_new)

    # Transform to numpy arrays and calculate fraction of MAC
    m_lst = np.array(m_lst)
    xcg_lst = np.array(xcg_lst)/MAC

    # Get Maximum and minimum center of gravity location and apply safety margin
    xcg_min = xcg_lst[-1] - sm/2
    xcg_min_ground = min(xcg_lst) - sm/2
    xcg_max = max(xcg_lst) + sm/2

    # Create loading diagram
    if plot:
        # Transform lists to numpy arrays
        bp_b_m = np.array(bp_b_m)
        bp_p_m = np.array(bp_p_m)
        pb_p_m = np.array(pb_p_m)
        pb_b_m = np.array(pb_b_m)

        # Calculate fraction of MAC
        bp_b_xcg = percMAC(np.array(bp_b_xcg), xlemac, MAC)
        bp_p_xcg = percMAC(np.array(bp_p_xcg), xlemac, MAC)
        pb_p_xcg = percMAC(np.array(pb_p_xcg), xlemac, MAC)
        pb_b_xcg = percMAC(np.array(pb_b_xcg), xlemac, MAC)

        xcg_min_plot = percMAC(xcg_min*MAC,xlemac,MAC)
        xcg_min_ground_plot = percMAC(xcg_min_ground*MAC,xlemac,MAC)
        xcg_max_plot = percMAC(xcg_max * MAC, xlemac, MAC)

        # Create the plot
        plt.plot(bp_b_xcg*100,bp_b_m,color='#1f77b4',label='Battery loading')
        plt.plot(pb_b_xcg*100,pb_b_m,'--',color='#1f77b4')
        plt.plot(bp_p_xcg*100,bp_p_m,color='#ff7f0e',label='Payload loading')
        plt.plot(pb_p_xcg*100,pb_p_m,'--',color='#ff7f0e')
        plt.plot([xcg_min_plot*100,xcg_min_plot*100],[min(m_lst),max(m_lst)],'r',label='Minimum/Maximum Center of Gravity')
        plt.plot([xcg_max_plot*100,xcg_max_plot*100],[min(m_lst),max(m_lst)],'r')
        plt.plot([xcg_min_ground_plot*100,xcg_min_ground_plot*100],[min(m_lst),max(m_lst)],'r--')
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
    variables [class]:  The class contains updated values of:
                         - deda [float]:    Wing downwash gradient [-]
    ShS_min [float]:    Minimum required horizontal tail surface [-]

V&V:    Verified
'''

def scissor_plot(variables,lh,xlemac,xcg_min,xcg_max,plot=False):
    # Input parameters
    R = variables.R                     # [J/kg K]  Gas constant
    gamma = variables.gamma             # [-]       Heat capacity ratio
    T0 = variables.T0                   # [K]       Base temperature
    lmbda = variables.lmbda             # [degC/m]  Lapse rate

    bf = variables.fuselagewidth        # [m]       Fuselage width
    hf = variables.fuselageheight       # [m]       Fuselage height
    lf = variables.fuselagelength       # [m]       Fuselage length

    zw = variables.h_landinggear        # [m]       Height of the wing, from ground
    MAC = variables.MAC                 # [m]       Mean Aerodynamic Chord
    Sw = variables.S                    # [m2]      Horizontal tail surface area
    Snet = variables.Snet               # [m2]      Net wing surface area
    bw = variables.b                    # [m]       Wing span
    Aw = variables.A                    # [-]       Wing aspect ratio
    sweepw = 0                          # [rad]     Wing quarter chord sweep angle
    taperw = variables.taper            # [-]       Wing taper ratio
    twistwr = variables.twist           # [rad]     Wing twist at the root
    crw = variables.c_r                 # [m]       Wing root chord
    bfl = variables.flapspan            # [m]       Flap span

    zh = variables.h_htail              # [m]       Height horizontal tail from ground
    Ah = variables.A_h                  # [-]       Horizontal tail aspect ratio
    bh = variables.b_h                  # [m]       Horizontal tail span
    taperh = variables.taper_h          # [-]       Horizontal tail taper ratio
    sweeph = variables.sweeph           # [rad]     Horizontal tail quarter chord sweep angle
    crh = variables.c_r_h               # [m]       Horizontal tail root chord

    VhV = variables.VhV                 # [-]       Tail/wing speed ratio

    Vcruise = variables.Vcruise         # [m/s]     Cruise speed
    hcruise = variables.hcruise         # [m]       Cruise altitude

    eta = variables.eta                 # [-]       Airfoil efficiency coefficient
    CLaw = variables.wing_CL_alpha      # [/rad]    Wing lift rate coefficient
    Cm0af = variables.Cm0af             # [-]       Airfoil zero lift pitching moment coefficient
    mu1 = variables.mu1                 # [-]       Flap coefficient 1
    mu2 = 1.2*(bfl/bw)+0.13             # [-]       Flap coefficient 2
    mu3 = 0.06*(bfl/bw)+0.0335          # [-]       Flap coefficient 3
    dClmax = variables.dClmax           # [-]       Airfoil lift coefficient increase at landing
    cc = variables.cc                   # [-]       Chord ratio (extended flap/clean)
    CL_landing = variables.CL_landing   # [-]       Wing lift coefficient at landing (all flaps deployed)
    Swf = variables.flapaffectedarea    # [m2]      Reference wing flapped surface area
    CL0 = variables.CL0flap             # [-]       Flapped wing lift coefficient at zero angle of attack
    CLA_h = CL_landing                  # [-]       Aircraft less tail lift coefficient
    CLh = variables.CLh_L               # [-]       Horizontal tail landing configuration lift coefficient

    sm_free = 0.05                      # [-]       Fraction neutral point shift for stick-free stability
    sm = 0.05                           # [-]       Stability margin

    # Parameter calculations
    # lfn
    xmac = XMAC(bw,sweepw,taperw,crw)
    xwing = xlemac-xmac
    sweepLEw = sweep(0,bw,sweepw,taperw,crw)
    lfn = xwing + bf/2*tan(sweepLEw)

    # Tail lift rate coefficient
    Vh = VhV*Vcruise
    Tcruise = T0 + lmbda*hcruise
    a = sqrt(gamma*R*Tcruise)
    Mh = Vh/a
    betah = sqrt(1 - Mh**2)
    sweepc2h = sweep(0.5,bh,sweeph,taperh,crh)
    CLah = 2*pi*Ah/(2 + sqrt(4 + (Ah*betah/eta)**2 * (1 + tan(sweepc2h)**2/(betah**2))))
    # Lift rate coefficient of the aircraft less tail
    CLaA_h = CLaw * (1 + 2.15*bf/bw) * Snet/Sw + pi/2*bf**2/Sw

    # Aerodynamic center
    xacw = (0.25*MAC+xlemac)/MAC
    #cg = Sw/bw
    xacf1 = -1.8/CLaA_h * bf*hf*lfn/(Sw*MAC)
    xacf2 = 0#0.273/(1 + taperw) * bf*cg*(bw-bf)/(MAC**2*(bw + 2.15*bf)) * tan(sweepw)
    xacwf = xacw + xacf1 + xacf2
    xacn = 0
    xac = xacwf + xacn

    # Wing downwash gradient
    deda = downwash(lh,bw,zh,zw,twistwr,sweepw,CLaw,Aw)

    # Pitching moment coefficient
    Cmacw = Cm0af*(Aw*cos(sweepw)**2/(Aw + 2*cos(sweepw)))
    dfusCmac = -1.8*(1-2.5*bf/lf) * pi*bf*hf*lf/(4*Sw*MAC) * CL0/CLaA_h
    dfCmac = mu2*(-mu1*dClmax*cc-(CL_landing+dClmax*(1-Swf/Sw))*(1/8)*cc*(cc-1))\
             + 0.7*Aw/(1+2/Aw)*mu3*dClmax*tan(sweepw)
    dnacCmac = 0
    Cmac = Cmacw + dfCmac + dfusCmac + dnacCmac

    # Stability Analysis
    def stability(xcg):
        a = 1/(CLah/CLaA_h * (1-deda) * lh/MAC * VhV**2)
        b = -(xac-sm)/(CLah/CLaA_h * (1-deda) * lh/MAC * VhV**2)
        return (a*xcg/(1-sm_free) + b)

    def stability_curve(ShS):
        return (1-sm_free)*(xac+CLah/CLaA_h*(1-deda)*ShS*lh/MAC*VhV**2)-sm

    ShS_stability = stability(xcg_max)

    # Controllability Analysis
    def control(xcg):
        a = 1/(CLh/CLA_h * lh/MAC * VhV**2)
        b = (Cmac/CLA_h - xac)/(CLh/CLA_h * lh/MAC * VhV**2)
        return a*xcg + b

    def control_curve(ShS):
        return xac-Cmac/CLA_h+CLh/CLA_h*ShS*lh/MAC*VhV**2

    ShS_Control = control(xcg_min)

    # Minimum required horizontal tail surface
    ShS_min = max(ShS_stability,ShS_Control)

    # Neutral point
    xnp = stability_curve(ShS_min) + sm

    # Create Scissor Plot
    if plot:

        ShS = np.arange(0,0.401,0.001)

        xcg_stability = stability_curve(ShS)*MAC
        xcg_control = control_curve(ShS)*MAC

        xcg_stability = percMAC(xcg_stability,xlemac,MAC)
        xcg_control = percMAC(xcg_control,xlemac,MAC)

        xcg_min_plot = percMAC(xcg_min*MAC,xlemac,MAC)
        xcg_max_plot = percMAC(xcg_max*MAC,xlemac,MAC)

        plt.plot(xcg_stability*100,ShS,label='Stability curve')
        plt.plot((xcg_stability+sm)*100,ShS,'k--')
        plt.plot(xcg_control*100,ShS,label='Control curve')
        plt.plot([xcg_min_plot*100,xcg_max_plot*100],[ShS_min,ShS_min],'r',label='Center of Gravity range')
        plt.title('Scissor Plot')
        plt.xlabel('$x_{cg}$/MAC (%)')
        plt.ylabel('$S_h/S$ [-]')
        plt.legend()
        plt.grid()
        plt.show()

    else:
        pass

    return ShS_min,xnp


'''
sizing_htail_wingpos()  Sizing of the horizontal tail and wing position

Inputs:
    variables[class]:   
    plot [bool]:        Create a plot (True or False), default is False

Outputs:
    variables [class]:  The class contains updated values of:
                         - xcg_min [float]:     Minimum center of gravity location [m]
                         - xcg_max [float]:     Maximum center of gravity location [m]
                         - xwing [float]:       Wing location [m]
                         - xlemac [float]:      Distance nose - leading edge mean aerodynamic chord [m]
                         - Sh_min [float]:      Minimum required horizontal tail surface [m2]
                         - deda [float]:        Wing downwash gradient [-]

V&V:    Verified
'''

def sizing_htail_wingpos(variables,plot=False):
    # Inputs
    lf = variables.fuselagelength       # [m]       Fuselage length

    cg_wing = variables.cg_wing         # [m]       Wing center of gravity
    zw = variables.h_landinggear        # [m]       Wing height
    Sw = variables.S                    # [m2]      Wing surface area
    sweepw = 0                          # [rad]     Wing quarter chord sweep angle
    taperw = variables.taper            # [-]       Wing taper ratio
    Aw = variables.A                    # [-]       Wing aspect ratio
    bw = variables.b                    # [m]       Wing span
    twistwr = variables.twist           # [rad]     Wing twist
    MAC = variables.MAC                 # [m]       Mean aerodynamic chord
    crw = variables.c_r                 # [m]       Wing root chord
    CLaw = variables.wing_CL_alpha      # [/rad]    Wing lift rate coefficient

    xtail = variables.xtail             # [m]       Tail location
    lh = variables.lh                   # [m]       Tail arm
    zh = variables.h_htail              # [m]       Tail height

    sm = 0.1                            # [-]       Safety margin

    # Parameter calculations
    # xlemac and xmac
    xmac = XMAC(bw,sweepw,taperw,crw)

    # Create lists
    xlemaclf_lst = np.arange(0,1.001,0.001)
    xcgmin_lst = []
    xcgmax_lst = []
    ShSmin_lst = []
    lh_lst = []

    # Perform wing shift
    for xlemaclf in xlemaclf_lst:
        xlemacw_i = xlemaclf*lf + variables.proplength
        xwing_i = xlemacw_i-xmac
        xcg_wing_i = xwing_i + cg_wing

        xcg_min_i,xcg_max_i = loading_diagram(variables,xcg_wing_i)
        ShS_min_i,xnp = scissor_plot(variables,lh,xlemacw_i,xcg_min_i,xcg_max_i)
        lh = xtail - xnp*MAC

        xcgmin_lst.append(xcg_min_i)
        xcgmax_lst.append(xcg_max_i)
        ShSmin_lst.append(ShS_min_i)
        lh_lst.append(lh)

    # Determine minimum horizontal tail surface
    ShS_min = min(ShSmin_lst)
    Sh = ShS_min*Sw*(1+sm)
    # Get index of minimum horizontal tail surface
    i = ShSmin_lst.index(ShS_min)

    # Determine center of gravity range for minimum horizontal tail surface
    xcg_min = xcgmin_lst[i]
    xcg_max = xcgmax_lst[i]

    # Determine tail length
    lh = lh_lst[i-1]

    # Determine wing location parameters for minimum horizontal tail surface
    xlemaclf = xlemaclf_lst[i]
    xlemacw = xlemaclf*lf + variables.proplength
    xwing = xlemacw - xmac
    xcg_wing = xwing + cg_wing

    # Determine downwash
    deda = downwash(lh,bw,zh,zw,twistwr,sweepw,CLaw,Aw)

    # Update values in variables class
    variables.xwing = xwing
    variables.xlemac = xlemacw
    variables.deda = deda
    variables.lh = lh
    variables.xcg_min = xcg_min*MAC
    variables.xcg_max = xcg_max*MAC
    variables.Sh = Sh

    # Create plots
    if plot:
        # Transform lists to numpy arrays
        xlemac_lst = xlemaclf_lst*lf
        xcgmin_lst = percMAC(np.array(xcgmin_lst)*MAC-variables.proplength,xlemac_lst,MAC)
        xcgmax_lst = percMAC(np.array(xcgmax_lst)*MAC-variables.proplength,xlemac_lst,MAC)

        xcg_min_plot = percMAC(xcg_min*MAC,xlemacw,MAC)
        xcg_max_plot = percMAC(xcg_max*MAC,xlemacw,MAC)

        # Loading diagram
        loading_diagram(variables,xcg_wing,plot)

        # Scissor plot
        scissor_plot(variables,lh,xlemacw,xcg_min,xcg_max,plot)

        # Center of gravity range plot
        plt.plot(xcgmin_lst*100, xlemaclf_lst*100,label='Most forward center of gravity')
        plt.plot(xcgmax_lst*100, xlemaclf_lst*100,label='Most aft center of gravity')
        plt.plot([xcg_min_plot*100,xcg_max_plot*100],[xlemaclf*100,xlemaclf*100],'r',label='Center of Gravity range')
        plt.title('Center of gravity range')
        plt.xlabel('$x_{cg}$/MAC (%)')
        plt.ylabel('$x_{lemac}/l_{fus}$ (%)')
        plt.legend()
        plt.grid()
        plt.show()
    else:
        pass

    return variables