import numpy as np
import variables as v
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os

from WingPlanform import WingPlanform

def plotLiftDistribution(y, Cl_range, ClmaxDistr=None, legend=False):
    fig = plt.figure(figsize=(10, 4.5))
    ax1 = fig.add_subplot(111)

    i = 0
    for distr in Cl_range:
        marker = 'o' if i == 0 else 'v' if i == 1 else 's'
        alpha = 2 if i == 0 else 6 if i == 1 else 10
        ax1.plot(y, distr, linewidth=2, color='green', marker=marker, fillstyle='full', markevery=5, label=f'Lift ditribution @ alpha={alpha}')
        if ClmaxDistr: ax1.plot(y, ClmaxDistr(y), linewidth=1.5, linestyle='-.', color='black', label='Local stall limit')
        i += 1

    ax1.axvline(x=0, linewidth=2, color='black')
    ax1.axhline(y=0, linewidth=2, color='black')
    ax1.set_xlabel('Wingspan [m]')
    ax1.set_ylabel('Cl [-]')
    ax1.xaxis.grid(color='black', linestyle='--')
    ax1.yaxis.grid(color='black', linestyle='--')
    if legend: plt.legend(loc='lower center')
    fig.suptitle('Model Lift Distribution', fontsize=16, y=0.97)
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.show()

if __name__ == "__main__":
    # OPTIMISATION LOOP
    # WingList = np.array([])
    WingList = []
    for taper in np.linspace(0.1,1.0,11):
        for twist in np.linspace(0,np.radians(6),7):
            if twist == 0:
                print(taper*100,"% completed")

            wing = WingPlanform(v.S, v.A, taper, twist, v.gamma)
            wing.setAirfoils(v.Clmax_r, v.Clmax_t, v.Cla_r, v.Cla_t, v.a0_r, v.a0_t, v.Cd0_r, v.Cd0_t, v.deltaAlphaStall_r, v.deltaAlphaStall_t)
    # wing.calcCoefficients(1000, tipCutoff=0.5) # use to get sensible CDi
            wing.calcCoefficients(200, tipCutoff=0.6)
            
            # alpha = np.radians(5)
            # Cl_distr, yPnts = wing.calcLiftDistribution(alpha, 100)
            # plotLiftDistribution(yPnts, [Cl_distr])


            CLmax, alphaMax, Cl_distrMax, yPntsMax, ClmaxDistr, stallpos = wing.calcCLmax()
            # plotLiftDistribution(yPntsMax, [Cl_distrMax], ClmaxDistr=ClmaxDistr, legend=True)

            CDi1,e1 = wing.calcCDi(np.radians(2))
            CDi2,e2 = wing.calcCDi(np.radians(6))
            CDi3,e3 = wing.calcCDi(np.radians(10))
            if stallpos <= 0.5*wing.b and CLmax >= 1.40:
                # WingList = np.append(WingList, [[taper,twist,CLmax,stallpos,alphaMax,e1,e2,e3]])
                WingList.append([taper,twist,CLmax,stallpos,alphaMax,e1,e2,e3])

    with open('Final_Design/aerodynamics/winglist.csv', 'w') as file:
        file.write('taper, twist, CLmax, stallpos, alphaMax, e1, e2, e3\n')
        for wingOption in WingList:
            for prop in wingOption:
                file.write(f'{round(prop, 3)}, ')
            file.write('\n')

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