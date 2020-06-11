import numpy as np
import variables as v
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.lines import Line2D
from collections import OrderedDict
import os
from utils import plotLiftDistribution, readWinglist, plotDesignParams

from WingPlanform import WingPlanform

def optimize():
    WingList = []
    for taper in np.linspace(0.1,1.0,11):
        for twist in np.linspace(0,np.radians(6),7):
            if twist == 0:
                print(taper*100,"% completed")

            wing = WingPlanform(v.S, v.A, taper, twist, v.gamma)
            wing.setAirfoils(v.Clmax_r, v.Clmax_t, v.Cla_r, v.Cla_t, v.a0_r, v.a0_t, v.Cd0_r, v.Cd0_t, v.deltaAlphaStall_r, v.deltaAlphaStall_t)
            wing.calcCoefficients(200, tipCutoff=0.6)
            
            CLmax, alphaMax, Cl_distrMax, yPntsMax, ClmaxDistr, stallpos = wing.calcCLmax()
            espan = wing.calcespan()

            # if stallpos <= 0.5*wing.b and CLmax >= 1.40:
            if True:
                WingList.append([taper, twist, CLmax, stallpos, alphaMax, espan])

    with open('Final_Design/aerodynamics/winglist.csv', 'w') as file:
        file.write('taper, twist, CLmax, stallpos, alphaMax, espan\n')
        for wingOption in WingList:
            for prop in wingOption:
                file.write(f'{round(prop, 3)}, ')
            file.write('\n')

if __name__ == "__main__":
    
    # optimize()
    bestWing, taper, CLmax, espan = readWinglist()
    print(f'\nOptimized wing:\t taper={bestWing[0]} \t twist={bestWing[1]} \t CLmax={bestWing[2]} \t espan={bestWing[3]} \n')
    
    plotDesignParams(taper, espan, CLmax, 'Taper', 'e span', 'CLmax')



    # wing = WingPlanform(v.S, v.A, taper, twist, v.gamma)
    # wing.setAirfoils(v.Clmax_r, v.Clmax_t, v.Cla_r, v.Cla_t, v.a0_r, v.a0_t, v.Cd0_r, v.Cd0_t, v.deltaAlphaStall_r, v.deltaAlphaStall_t)
    # wing.calcCoefficients(200, tipCutoff=0.6)

    # alpha = np.radians(5)
    # Cl_distr, yPnts = wing.calcLiftDistribution(alpha, 100)
    # plotLiftDistribution(yPnts, [Cl_distr])

    """ CLa """
    # if False:
    #     CLa = wing.calcCLa()
    #     print(f'CLa = {CLa} [/rad] or {CLa*np.pi/180} [deg]')

    # """ CLmax """
    # if False:
    #     CLmax, alphaMax, Cl_distrMax, yPntsMax, ClmaxDistr, stallpos = wing.calcCLmax(plotProgression=False, printMaxLoc=True)
    #     print(f'CLmax = {round(CLmax, 2)} @ a = {round(np.degrees(alphaMax), 2)}')
    #     plotLiftDistribution(yPntsMax, [Cl_distrMax], ClmaxDistr=ClmaxDistr, legend=True)



# wing.calcCoefficients(500, tipCutoff=0.9)
    
# alpha = np.radians(6)
# Cl_distr, yPnts, CDi_distr, yPntsCDi = wing.calcLiftDistribution(alpha, 100)
# plotLiftDistribution(yPnts, [Cl_distr])
# plotLiftDistribution(yPntsCDi, [CDi_distr])

# # alphai_distr, yPnts = wing.calcAlphai(alpha, 100)
# # plotLiftDistribution(yPnts, [np.degrees(alphai_distr)])

# """ CLa """
# if False:
#     CLa = wing.calcCLa()
#     print(f'CLa = {CLa} [/rad] or {CLa*np.pi/180} [deg]')

# """ CLmax """
# if False:
#     CLmax, alphaMax, Cl_distrMax, yPntsMax, ClmaxDistr = wing.calcCLmax(plotProgression=False, printMaxLoc=True)
#     print(f'CLmax = {round(CLmax, 2)} @ a = {round(np.degrees(alphaMax), 2)}')
#     plotLiftDistribution(yPntsMax, [Cl_distrMax], ClmaxDistr=ClmaxDistr, legend=True)