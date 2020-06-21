import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os

# Geometry derived with 0 TE sweep constraint
def sweepLE(A, taper):
    return np.arctan( 4/A * (1-taper)/(1+taper) )

def C1(taper):
    x_raymerC1 = np.arange(0, 1.1, 0.1)
    y_raymerC1 = [0, 0.25, 0.49, 0.49, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0]
    return interp1d(x_raymerC1, y_raymerC1, kind='cubic')(taper)
def C2(taper):
    x_raymerC2 = np.arange(0, 1.1, 0.1)
    y_raymerC2 = [0, 0.25, 0.5, 0.8, 1.1, 1.05, 1, 0.95, 0.9, 0.86, 0.85]
    return interp1d(x_raymerC2, y_raymerC2, kind='cubic')(taper)

# Raymer Fig 12.15
def x_raymer1215(A, sweepLE, C2):
    return (C2 + 1) * A * np.tan(sweepLE)
def line_raymer1215(A, sweepLE, taper):
    return A * np.cos(sweepLE)*(1 + (2*taper)**2)
# Raymer Fig 12.14
def x_amaxBase(A, sweepLE, C1):
    beta = 1
    return (C1 + 1) * A/beta * np.cos(sweepLE)

def amaxBase(x):
    x_raymeramaxBase = np.arange(0, 4.4+0.4, 0.4)
    y_raymeramaxBase = [35, 35, 35, 32, 28, 25, 23, 22, 21, 21, 20.5, 20.5, 20.5]
    return float(interp1d(x_raymeramaxBase, y_raymeramaxBase, kind='cubic')(x))

if __name__ == "__main__":

    taper = 0.7
    A = 3

    sweepLE = sweepLE(A, taper)
    print(np.degrees(sweepLE))

    x_amaxBase = x_amaxBase(A, sweepLE, C1(taper))
    if x_amaxBase > 4.4:
        print('too high aspect ratio')
    else:
        amaxBase = amaxBase(x_amaxBase)
        print(f'amax base = {round(amaxBase, 1)} deg')

    x_raymer1215 = x_raymer1215(A, sweepLE, C2(taper))
    line_raymer1215 = line_raymer1215(A, sweepLE, taper)
    print( f'amax raymer Fig 12.15: x: {round(x_raymer1215, 1)}, line: {round(line_raymer1215, 1)}' )


    print(f'CLmax base raymer Fig 12.12: x: {round(x_amaxBase, 2)}')
    print(f'CLmax delta raymer Fig: x: {round(x_raymer1215, 1)}')



















    


    # fig = plt.figure(figsize=(4.5, 2))
    # ax1 = fig.add_subplot(111)
    # x = np.arange(0, 3.6+0.1, 0.01)
    # ax1.plot(x, test(x))
    # ax1.set_ylim([0, 50])
    # ax1.set_yticks(np.arange(0, 50+10, 10))
    # ax1.set_xticks(np.arange(0, 3.6+0.4, 0.4))
    # ax1.xaxis.grid(color='black', linestyle='--')
    # ax1.yaxis.grid(color='black', linestyle='--')
    # plt.show()

    
