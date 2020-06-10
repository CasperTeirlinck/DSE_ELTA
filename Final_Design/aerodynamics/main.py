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

def get_xfoil_data(alpha):
    y = []
    Cl = []
    with open(os.path.join(os.path.dirname(__file__), f'wingValidationData/wing{alpha}.dat'), 'r') as f:
        for line in f.readlines()[21:58]:
            columns = line.strip().split()
            if len(columns) > 0:
                y.append(float(columns[0]))
                Cl.append(float(columns[3]))
    return y, Cl

def calc_err(x1, y1, x2, y2, errrange, rangestep):
    ipx = np.arange(*errrange, rangestep)
    ipy1 = interp1d(x1, y1, kind='linear', fill_value='extrapolate')(ipx)
    ipy2 = interp1d(x2, y2, kind='linear', fill_value='extrapolate')(ipx)
    err = np.abs(ipy2 - ipy1)
    rmsd = np.sqrt(np.sum(err**2 / err.size)) / (np.max(y2) - np.min(y2)) * 100 # https://en.wikipedia.org/wiki/Root-mean-square_deviation
    return rmsd

if __name__ == "__main__":

    wing = WingPlanform(v.S, v.A, v.taper, v.twist, v.gamma)
    wing.setAirfoils(v.Clmax_r, v.Clmax_t, v.Cla_r, v.Cla_t, v.a0_r, v.a0_t, v.Cd0_r, v.Cd0_t, v.deltaAlphaStall_r, v.deltaAlphaStall_t)
    # wing.calcCoefficients(1000, tipCutoff=0.5) # use to get sensible CDi
    wing.calcCoefficients(500, tipCutoff=0.9)
    
    # alpha = np.radians(5)
    # Cl_distr, yPnts = wing.calcLiftDistribution(alpha, 100)
    # plotLiftDistribution(yPnts, [Cl_distr])

    """ CLa """
    if False:
        CLa = wing.calcCLa()
        print(f'CLa = {CLa} [/rad] or {CLa*np.pi/180} [deg]')

    """ CLmax """
    if True:
        CLmax, alphaMax, Cl_distrMax, yPntsMax, ClmaxDistr = wing.calcCLmax(plotProgression=True)
        print(f'CLmax = {round(CLmax, 2)} @ a = {round(np.degrees(alphaMax), 2)}')
        plotLiftDistribution(yPntsMax, [Cl_distrMax], ClmaxDistr=ClmaxDistr, legend=True)