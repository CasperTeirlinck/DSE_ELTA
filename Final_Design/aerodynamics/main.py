import numpy as np
import variables_aero as v
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os

from WingPlanform import WingPlanform

def plotLiftDistribution(y, Cl_range, y2=None, Cl_range2=None, ClmaxDistr=None):
    fig = plt.figure(figsize=(10, 4.5))
    ax1 = fig.add_subplot(111)

    if (y2 and Cl_range2):
        for distr2 in Cl_range2:
            ax1.plot(y2, distr2, linewidth=2, color='black', marker='o', fillstyle='none', markevery=2)

    for distr in Cl_range:
        ax1.plot(y, distr, linewidth=2, color='green', marker='o', fillstyle='full', markevery=5)
        if ClmaxDistr: ax1.plot(y, ClmaxDistr(y), linewidth=1.5, linestyle='-.', color='black')

    ax1.axvline(x=0, linewidth=2, color='black')
    ax1.axhline(y=0, linewidth=2, color='black')
    ax1.set_xlabel('Wingspan [m]')
    ax1.set_ylabel('')
    ax1.xaxis.grid(color='black', linestyle='--')
    ax1.yaxis.grid(color='black', linestyle='--')
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
    wing.calcCoefficients(1000)

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
        alpha = np.radians(5)
        CL = wing.calcCL(alpha)
        print(f'CL = {round(CL, 2)} @ a = {round(np.degrees(alpha), 2)}')

        CDi, e = wing.calcCDi(alpha)
        print(f'CDi = {CDi} @ a = {round(np.degrees(alpha), 2)}')
        print(f'e = {e} @ a = {round(np.degrees(alpha), 2)}')

        Cl_distr, yPnts = wing.calcLiftDistribution(alpha, 100)
        yPnts2, Cl_distr2 = get_xfoil_data(5)
        plotLiftDistribution(yPnts, [Cl_distr], y2=yPnts2, Cl_range2=[Cl_distr2])

        # err_range = (yPnts[0], yPnts[-1])
        # err = ( calc_err(Cla_x, Cla_y, naca_Cla_x, naca_Cla_y, err_Cla_range1, 1) )

        # print(f'error Cla:  {round(err_Cla[0], 1)}% full range  {round(err_Cla[1], 1)}% limited range')

    # Cl_distr_range = []
    # # alphaRange = [5]
    # alphaRange = np.arange(-15, 15, 1)
    # for alpha in alphaRange:
    #     alpha = np.radians(alpha)
    #     Cl_distr, yPnts = wing.calcLiftDistribution(alpha)
    #     Cl_distr_range.append(Cl_distr)
    # plotLiftDistribution(yPnts, Cl_distr_range)

    """ CLa """
    CLa = wing.calcCLa()
    print(f'CLa = {CLa}')

    """ CLmax """
    if False:
        CLmax, alphaMax, Cl_distrMax, yPntsMax, ClmaxDistr = wing.calcCLmax()
        print(f'CLmax = {round(CLmax, 2)} @ a = {round(np.degrees(alphaMax), 2)}')
        plotLiftDistribution(yPntsMax, [Cl_distrMax], ClmaxDistr=ClmaxDistr)