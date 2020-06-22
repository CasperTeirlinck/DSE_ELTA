import numpy as np
import variables as v
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.lines import Line2D
from collections import OrderedDict
import os
from utils import plotLiftDistribution, readWinglist, plotDesignParams, readAeroLoads, plotPlanform, plotHtail

from WingPlanform import WingPlanform

def optimize():
    WingList = []
    for taper in np.linspace(0.2,1.0,10):
        for twist in np.linspace(0, np.radians(6), 7):
            if twist == 0:
                print(taper*100,"% completed")

            wing = WingPlanform(v.S, v.A, taper, twist, v.gamma, v.CD0)
            wing.setAirfoils(v.Clmax_r, v.Clmax_t, v.Cla_r, v.Cla_t, v.a0_r, v.a0_t, v.Cd0_r, v.Cd0_t, v.deltaAlphaStall_r, v.deltaAlphaStall_t)
            wing.setWinglets(v.hwl, v.kwl)
            wing.calcCoefficients(200, tipCutoff=0.6)
            
            CLmax, alphaMax, Cl_distrMax, yPntsMax, ClmaxDistr, stallpos = wing.calcCLmax()
            espan = wing.calcespan()
            e_oswald = wing.calcOswald(v.w_fuselage, hasWinglets=True)
            if True: 
            # if stallpos <= 0.4*wing.b/2 and CLmax >= 1.40:
                WingList.append([taper, twist, CLmax, stallpos, alphaMax, espan, e_oswald])

    with open('Final_Design/aerodynamics/winglist.csv', 'w') as file:
        file.write('taper, twist, CLmax, stallpos, alphaMax, espan, e_oswald\n')
        for wingOption in WingList:
            for prop in wingOption:
                file.write(f'{round(prop, 3)}, ')
            file.write('\n')

    WingList.sort(key= lambda x: x[5], reverse = True)

    with open('Final_Design/aerodynamics/winglist_sorted.csv', 'w') as file:
        file.write('taper, twist, CLmax, stallpos, alphaMax, espan, e_oswald\n')
        for wingOption in WingList:
            for prop in wingOption:
                file.write(f'{round(prop, 3)}, ')
            file.write('\n')

if __name__ == "__main__":
    
    """ === OPTIMIZE === """

    # optimize()
    taper, CLmax, espan = readWinglist()
    # plotDesignParams(taper, espan, CLmax, 'Taper', 'e span', 'CLmax')

    # y_list, cl_list, cd_list, xcp_list = readAeroLoads(5)

    taper = 0.45
    twist = np.radians(5)

    """ === SHOW === """

    if False:
        S_h = 20
        A_h = 3
        taper_h = 0.7
        b_h = np.sqrt(A_h*S_h)
        c_r_h = (2*S_h)/(b_h*(1 + taper_h))
        c_t_h = taper_h*c_r_h
        plotHtail(c_r_h, c_t_h, b_h)

    wing = WingPlanform(v.S, v.A, taper, twist, v.gamma)
    wing.setAirfoils(v.Clmax_r, v.Clmax_t, v.Cla_r, v.Cla_t, v.a0_r, v.a0_t, v.tc_airfoil, v.xc_airfoil, v.Cd0_r, v.Cd0_t)
    wing.setWinglets(v.hwl, v.kwl)
    wing.calcCoefficients(200, tipCutoff=0.6)

    print(f'\n=== WINGPLANFORM ===\n')
    print(f'A = {wing.A}')
    print(f'S = {wing.S} m')
    print(f'b = {round(wing.b, 2)} m')
    print(f'taper = {wing.taper}')
    print(f'twist = {round(np.degrees(wing.twist), 1)} deg')
    print(f'cr = {round(wing.c_r, 2)} m')
    print(f'ct = {round(wing.c_t, 2)} m')
    print(f'mac = {round(wing.MAC, 2)} m')
    print(f'Xmac = {round(wing.XMAC, 3)} m')
    print(f'Ymac = {round(wing.YMAC, 3)} m')
    # print(f'e = {round(wing.calcOswald(v.w_fuselage, hasWinglets=True), 2)}')
    CLa = wing.calcCLa()
    print(f'CLa = {round(CLa, 2)} 1/rad or {round(CLa*np.pi/180, 2)} /deg')
    alpha0 = wing.calcAlpha0()
    print(f'alpha_L=0 = {round(np.degrees(alpha0), 2)} degrees')

    if False:
        CLmax, alphaMax, Cl_distrMax, yPntsMax, ClmaxDistr, stallpos = wing.calcCLmax(plotProgression=True)
        print(f'CLmax = {round(CLmax, 2)} @ a = {round(np.degrees(alphaMax), 2)} deg')
        plotLiftDistribution(yPntsMax, [Cl_distrMax], ClmaxDistr=ClmaxDistr, legend=True)

    print(f'\n=== ============ ===\n')

    if False:
        plotPlanform(wing.c_r, wing.c_t, wing.b, wing.MAC, wing.XMAC, wing.YMAC)

    if False:
        alpha = np.radians(10.2)
        Cl_distr, yPnts = wing.calcLiftDistribution(alpha, 100)
        plotLiftDistribution(yPnts, [Cl_distr])

    CD0 = wing.calcCD0(v.S_wet_fus, v.l_fus, v.fus_A_max, v.w_fuselage, v.S_h, v.S_v, v.MAC_emp, v.BLturbratio_fus, v.BLturbratio_wing, v.BLturbratio_emp, v.l_gear, v.w_gear, v.dCD_gear, v.flap_area_ratio, v.tc_emp, v.xc_emp, v.V_stall, v.rho_cruise, v.visc,clean_config=False)

    CD0wing = wing.calcCD0wing(1.,0,0)

    e = wing.calcOswald(CD0, v.w_fuselage, hasWinglets=False)

    print(CD0)
    print(CD0wing)
    print(e)