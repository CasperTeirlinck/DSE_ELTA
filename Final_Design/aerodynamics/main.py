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

    WingList = np.array([])
    for taper in np.linspace(0,1.0,11):
        for twist in linspace(0,6,7):

            wing = WingPlanform(v.S, v.A, taper, twist, v.gamma)
            wing.setAirfoils(v.Clmax_r, v.Clmax_t, v.Cla_r, v.Cla_t, v.a0_r, v.a0_t, v.Cd0_r, v.Cd0_t, v.deltaAlphaStall_r, v.deltaAlphaStall_t)
    # wing.calcCoefficients(1000, tipCutoff=0.5) # use to get sensible CDi
            wing.calcCoefficients(100, tipCutoff=0.9)

            CLmax, alphaMax, Cl_distrMax, yPntsMax
    # alpha = np.radians(5)
    # Cl_distr, yPnts = wing.calcLiftDistribution(alpha, 100)
    # plotLiftDistribution(yPnts, [Cl_distr])

    """ CLa """
    if False:
        CLa = wing.calcCLa()
        print(f'CLa = {CLa} [/rad] or {CLa*np.pi/180} [deg]')

    """ CLmax """
    if False:
        CLmax, alphaMax, Cl_distrMax, yPntsMax, ClmaxDistr = wing.calcCLmax(plotProgression=False)
        print(f'CLmax = {round(CLmax, 2)} @ a = {round(np.degrees(alphaMax), 2)}')
        plotLiftDistribution(yPntsMax, [Cl_distrMax], ClmaxDistr=ClmaxDistr, legend=True)

    