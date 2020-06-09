import os
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

kts_to_ms = 0.514444444
ft_to_m = 0.3048

""" Variables """
h_cruise = 3000 * ft_to_m # [m] ref midterm
rho_cruise = 1.12 # [kg/m3] ref midterm

V_s = 45 * kts_to_ms # [m/s] ref midterm
V_cruise = 100 * kts_to_ms # [m/s] ref midterm
V_max = 120 * kts_to_ms # [m/s] ref midterm

mac_initial = 1.45 # [m] ref midterm

def Re_min(c):
    kin_visc_min = 13.28e-6 # [m2/s] @ 0C
    return (V_s*c)/kin_visc_min

def Re_max(c):
    kin_visc_max = 15.52e-6 # [m2/s] @ 25C
    return (V_max*c)/kin_visc_max

# c_mean = (1.16 + 1.289)/2
# print('test:', (Re_max(c_mean) + Re_min(c_mean))/2)

""" Get Xfoil Data """
def get_data_xfoil(foil, Re):

    a = []
    Cl = []
    Cd = []
    with open(os.path.join(os.path.dirname(__file__), f'NACA{foil}_T1_Re{Re}.000_M0.00_N9.0.txt'), 'r') as f:
        for line in f.readlines()[11:]:
            columns = line.strip().split()
            if len(columns) > 0:
                a.append(float(columns[0]))
                Cl.append(float(columns[1]))
                Cd.append(float(columns[2]))
        
    return a, Cl, Cd    

""" Get NACA Data """
def get_data_naca(foil, Re):
    a = [-16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16]
    if foil == '4415': # pdf p 499
        if Re == 3:
            Cl = [-0.85, -0.9, -0.85, -0.65,  -0.4,  -0.2,      0,    0.2,   0.45,    0.6,   0.85,  1.05,   1.2,  1.31,  1.41, 1.4, 1.31]
            Cd = [    0,    0, 0.013, 0.011, 0.009, 0.008, 0.0077, 0.0075, 0.0072, 0.0073, 0.0075, 0.008, 0.011, 0.015, 0.021,   0,    0]
        if Re == 6:
            Cl = [-0.85, -0.9, -0.85, -0.65,  -0.4,  -0.2,      0,    0.2,   0.45,    0.6,   0.85,  1.05,   1.2,   1.4,  1.51, 1.55, 1.49]
            Cd = [    0,    0, 0.0105, 0.0095, 0.008, 0.0075, 0.007, 0.0068, 0.0065, 0.0067, 0.0072, 0.0082, 0.0097, 0.014, 0.017, 0.018,    0]
    if foil == '652415': # pdf p 639
        a = [-16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
        if Re == 3:
            Cl = [-0.6, -1.0, -0.95,    -0.8,   -0.6,   -0.4,  -0.18,   0.05,    0.3,    0.5,   0.75,   0.9,  1.05,  1.15, 1.25, 1.34, 1.38, 1.37, 1.3]
            Cd = [   0,    0,     0,  0.0134, 0.0105,  0.009, 0.0078, 0.0075,  0.005, 0.0056, 0.0057, 0.0105, 0.016, 0.024,    0,    0,    0,    0,  0]
        if Re == 6:
            Cl = [-0.6, -1.0, -0.95,    -0.8,   -0.6,   -0.4,  -0.18,   0.05,    0.3,    0.5,   0.75,   0.9,   1.1,  1.25,  1.35, 1.43, 1.5, 1.53, 1.53]
            Cd = [   0,    0,     0,   0.011, 0.0092, 0.0079,  0.007, 0.0067, 0.0044, 0.0045, 0.0048, 0.010, 0.013, 0.018, 0.028,   0,    0,    0,    0]

    return a, Cl, Cd

if __name__ == "__main__":
    
    """ Data Setup """
    # foil = '4415' # why?: midterm choice
    foil = '652415' # why?: 6-series: good lift and drag (ref midterm), potential choice in airfoil selection
    Re = 6

    a, Cl, Cd = get_data_xfoil(foil, Re)
    naca_a, naca_Cl, naca_Cd = get_data_naca(foil, Re)

    # Cl alpha
    Cla_x = a
    Cla_y = Cl
    naca_Cla_x = naca_a
    naca_Cla_y = naca_Cl
    # Cl Cd
    if foil == '4415':
        if Re == 3:
            ClCd_x = Cl[:]
            ClCd_y = Cd[:]
            naca_ClCd_x = naca_Cl[2:-2]
            naca_ClCd_y = naca_Cd[2:-2]
        if Re == 6:
            ClCd_x = Cl[2:-1]
            ClCd_y = Cd[2:-1]
            naca_ClCd_x = naca_Cl[2:-1]
            naca_ClCd_y = naca_Cd[2:-1]
    if foil == '652415':
        if Re == 3:
            ClCd_x = Cl[:]
            ClCd_y = Cd[:]
            naca_ClCd_x = naca_Cl[3:-5]
            naca_ClCd_y = naca_Cd[3:-5]
        if Re == 6:
            ClCd_x = Cl[:]
            ClCd_y = Cd[:]
            naca_ClCd_x = naca_Cl[3:-4]
            naca_ClCd_y = naca_Cd[3:-4]

    """ Errors """
    def calc_err(x1, y1, x2, y2, errrange, rangestep):
        ipx = np.arange(*errrange, rangestep)
        ipy1 = interp1d(x1, y1, kind='linear', fill_value='extrapolate')(ipx)
        ipy2 = interp1d(x2, y2, kind='linear', fill_value='extrapolate')(ipx)
        err = np.abs(ipy2 - ipy1)
        rmsd = np.sqrt(np.sum(err**2 / err.size)) / (np.max(y2) - np.min(y2)) * 100 # https://en.wikipedia.org/wiki/Root-mean-square_deviation
        return rmsd

    # Cl alpha
    err_Cla = []
    err_Cla_range1 = (max(Cla_x[0], naca_Cla_x[0]), min(Cla_x[-1], naca_Cla_x[-1]))
    err_Cla_range2 = (-10, 8)

    err_Cla.append( calc_err(Cla_x, Cla_y, naca_Cla_x, naca_Cla_y, err_Cla_range1, 1) )
    err_Cla.append( calc_err(Cla_x, Cla_y, naca_Cla_x, naca_Cla_y, err_Cla_range2, 1) )

    print(f'error Cla:  {round(err_Cla[0], 1)}% full range  {round(err_Cla[1], 1)}% limited range')

    # Cl Cd
    err_ClCd = []
    err_ClCd_range1 = (max(ClCd_x[0], naca_ClCd_x[0]), min(ClCd_x[-1], naca_ClCd_x[-1]))
    err_ClCd_range2 = (-0.8, 1.2)

    err_ClCd.append( calc_err(ClCd_x, ClCd_y, naca_ClCd_x, naca_ClCd_y, err_ClCd_range1, 0.1) )
    err_ClCd.append( calc_err(ClCd_x, ClCd_y, naca_ClCd_x, naca_ClCd_y, err_ClCd_range2, 0.1) )

    print(f'error ClCd:  {round(err_ClCd[0], 1)}% full range  {round(err_ClCd[1], 1)}% limited range')

    """ Plotting """
    fig = plt.figure(figsize=(10, 4.5))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    # Cl alpha
    ax1.plot(Cla_x, Cla_y, linewidth=1, color='black', marker='o', fillstyle='full', markevery=1, label=f'Xfoil - Re={Re}.0e6')
    ax1.plot(naca_Cla_x, naca_Cla_y, linewidth=1, color='black', marker='o', fillstyle='none', markevery=1, label=f'NACA - Re={Re}.0e6')

    ax1.axvline(x=0, linewidth=2, color='black')
    ax1.axhline(y=0, linewidth=2, color='black')
    ax1.set_xlabel('alpha [deg]')
    ax1.set_ylabel('Cl [-]')
    ax1.xaxis.grid(color='black', linestyle='--')
    ax1.yaxis.grid(color='black', linestyle='--')

    ax1.legend()

    # Cl Cd
    ax2.plot(ClCd_x, ClCd_y, linewidth=1, color='black', marker='o', fillstyle='full', markevery=1, label=f'Xfoil - Re={Re}.0e6')
    ax2.plot(naca_ClCd_x, naca_ClCd_y, linewidth=1, color='black', marker='o', fillstyle='none', markevery=1, label=f'NACA - Re={Re}.0e6')

    ax2.axvline(x=0, linewidth=2, color='black')
    ax2.axhline(y=0, linewidth=2, color='black')
    ax2.set_xlabel('Cl [-]')
    ax2.set_ylabel('Cd [-]')
    ax2.xaxis.grid(color='black', linestyle='--')
    ax2.yaxis.grid(color='black', linestyle='--')

    fig.suptitle(f'NACA {foil}: Experimental vs Xfoil Polars', fontsize=16, y=0.97)
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.show()