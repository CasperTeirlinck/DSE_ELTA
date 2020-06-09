import numpy as np
import variables_aero as v
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os

from WingPlanform import WingPlanform

def plotLiftDistribution(y, Cl_range, y2=None, Cl_range2=None, ClmaxDistr=None, legend=False):
    fig = plt.figure(figsize=(10, 4.5))
    ax1 = fig.add_subplot(111)

    if (y2 and Cl_range2):
        i = 0
        for y2, distr2 in zip(y2, Cl_range2):
            marker = 'o' if i == 0 else 'v' if i == 1 else 's'
            alpha = 2 if i == 0 else 6 if i == 1 else 10
            ax1.plot(y2, distr2, linewidth=2, color='black', marker=marker, fillstyle='none', markevery=2, label=f'Xflr5 @ alpha={alpha}')
            i += 1

    i = 0
    for distr in Cl_range:
        marker = 'o' if i == 0 else 'v' if i == 1 else 's'
        alpha = 2 if i == 0 else 6 if i == 1 else 10
        ax1.plot(y, distr, linewidth=2, color='green', marker=marker, fillstyle='full', markevery=5, label=f'Model @ alpha={alpha}')
        if ClmaxDistr: ax1.plot(y, ClmaxDistr(y), linewidth=1.5, linestyle='-.', color='black')
        i += 1

    ax1.axvline(x=0, linewidth=2, color='black')
    ax1.axhline(y=0, linewidth=2, color='black')
    ax1.set_xlabel('Wingspan [m]')
    ax1.set_ylabel('Cl [-]')
    ax1.xaxis.grid(color='black', linestyle='--')
    ax1.yaxis.grid(color='black', linestyle='--')
    if legend: plt.legend(loc='lower center')
    fig.suptitle('Model Lift Distribution vs Xflr5', fontsize=16, y=0.97)
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.show()

def get_xfoil_data(alpha):
    y = []
    Cl = []
    with open(os.path.join(os.path.dirname(__file__), f'wing{alpha}.dat'), 'r') as f:
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
    wing.setAirfoils(v.Clmax_r, v.Clmax_t, v.Cla_r, v.Cla_t, v.a0_r, v.a0_t, v.Cd0_r, v.Cd0_t)
    wing.calcCoefficients(1000, tipCutoff=0.5) # use for CDi
    # wing.calcCoefficients(1000, tipCutoff=0.9)

    """ Test convergence """
    if False:
        CLlist = []
        alpha = np.radians(13)
        Nmax = 100
        for N in range(1, Nmax):
            coeff = wing.calcCoefficients(N)
            A1 = coeff[0][0]*alpha + coeff[0][1]
            CL = np.pi * v.A * A1
            print(f'N: {N}')
            CLlist.append(CL)
        plt.plot(range(1, Nmax), CLlist)
        plt.show()

    """ Lift Distribution Validation """
    if True:
        Cl_distr_range = []
        Cl_distr2_range = []
        yPnts2_range = []
        for a in [2, 6, 10]:
            alpha = np.radians(a)
            Cl_distr, yPnts = wing.calcLiftDistribution(alpha, 100)
            yPnts2, Cl_distr2 = get_xfoil_data(a)
            Cl_distr_range.append(Cl_distr)
            Cl_distr2_range.append(Cl_distr2)
            yPnts2_range.append(yPnts2)

            CL = wing.calcCL(alpha)
            print(f'CL: {round(CL, 2)} @ a = {a}')

            CDi1, e = wing.calcCDi(alpha)
            # wing.calcCoefficients(1000, tipCutoff=0.5)
            print(f'CDi: {round(CDi1, 3)} @ a = {a}')
            print(f'e: {round(e, 2)} @ a = {a}')

            err_range = (yPnts[0], yPnts[-1])
            err = calc_err(yPnts, Cl_distr, yPnts2, Cl_distr2, err_range, 0.1)
            print(f'error: {round(err, 1)}% @ a = {a}\n')

        plotLiftDistribution(yPnts, Cl_distr_range, y2=yPnts2_range, Cl_range2=Cl_distr2_range, legend=True)

    """ CLa """
    CLa = wing.calcCLa()
    print(f'CLa = {CLa} [/rad] or {CLa*np.pi/180} [deg]')

    """ CLmax """
    if False:
        CLmax, alphaMax, Cl_distrMax, yPntsMax, ClmaxDistr = wing.calcCLmax()
        print(f'CLmax = {round(CLmax, 2)} @ a = {round(np.degrees(alphaMax), 2)}')
        plotLiftDistribution(yPntsMax, [Cl_distrMax], ClmaxDistr=ClmaxDistr)