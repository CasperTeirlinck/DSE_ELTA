import numpy as np
import variables as v
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.lines import Line2D
from collections import OrderedDict
import os
from utils import plotLiftDistribution, readWinglist, plotDesignParams, readAeroLoads, plotPlanform

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

    y_list, cl_list, cd_list = readAeroLoads()

    taper = 0.55
    twist = np.radians(5)

    """ === SHOW === """

    wing = WingPlanform(v.S, v.A, taper, twist, v.gamma, v.CD0)
    wing.setAirfoils(v.Clmax_r, v.Clmax_t, v.Cla_r, v.Cla_t, v.a0_r, v.a0_t, v.Cd0_r, v.Cd0_t, v.deltaAlphaStall_r, v.deltaAlphaStall_t)
    wing.setWinglets(v.hwl, v.kwl)
    wing.calcCoefficients(200, tipCutoff=0.6)

    print(f'\n=== WINGPLANFORM ===\n')
    print(f'A = {wing.A}')
    print(f'S = {wing.S} m')
    print(f'b = {round(wing.b, 2)}')
    print(f'taper = {wing.taper}')
    print(f'twist = {round(np.degrees(wing.twist), 1)} deg')
    cr = wing.calculateChord(np.pi/2, wing.taper, wing.S, wing.b)
    ct = wing.calculateChord(0, wing.taper, wing.S, wing.b)
    print(f'cr = {round(cr, 2)} m')
    print(f'ct = {round(ct, 2)} m')
    CLa = wing.calcCLa()
    print(f'CLa = {round(CLa, 2)} 1/rad or {round(CLa*np.pi/180, 2)} /deg')

    if False:
        CLmax, alphaMax, Cl_distrMax, yPntsMax, ClmaxDistr, stallpos = wing.calcCLmax(plotProgression=True)
        print(f'CLmax = {round(CLmax, 2)} @ a = {round(np.degrees(alphaMax), 2)} deg')
        # plotLiftDistribution(yPntsMax, [Cl_distrMax], ClmaxDistr=ClmaxDistr, legend=True)

    print(f'\n=== ============ ===\n')

    plotPlanform(cr, ct, wing.b)

    if False:
        alpha = np.radians(10.2)
        Cl_distr, yPnts = wing.calcLiftDistribution(alpha, 100)
        plotLiftDistribution(yPnts, [Cl_distr])