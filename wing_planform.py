import numpy as np
import matplotlib.pyplot as plt
from variables import *

def wing_planform(variables,A,drawing=False):
    S = variables.WTO/variables.WS

    # Sweep
    sweep = np.rad2deg(np.arccos(1)) # [deg]

    # Taper
    taper = 0.2*(2-sweep*np.pi/180)

    # Span
    b = np.sqrt(S*A)

    # Chord
    cr = 2*S/((1+taper)*b)
    ct = taper*cr
    MAC = (2/3)*cr*((1+taper+taper**2)/(1+taper))

    # Wing drawing
    if drawing:

        xc4 = [-b/2, b/2]
        yc4 = [0.75*cr, 0.75*cr]

        xctl = [-b/2, -b/2]
        yctl = [yc4[0]-0.75*ct, yc4[0]+0.25*ct]
        xctr = [b/2, b/2]
        yctr = [yc4[0]-0.75*ct, yc4[0]+0.25*ct]

        xcr = [0, 0]
        ycr = [0, cr]

        xle = [-b/2, 0, b/2]
        yle = [yctl[1], cr, yctr[1]]

        xte = [-b/2, 0, b/2]
        yte = [yctl[0], 0, yctr[0]]

        plt.plot(xctl,yctl,color='k')
        plt.plot(xctr,yctr,color='k')
        plt.plot(xcr,ycr,color='k')
        plt.plot(xle,yle,color='k')
        plt.plot(xte,yte,color='k')
        plt.plot(xc4,yc4,'k:')
        plt.xlim(-b/2-2, b/2+2)
        plt.ylim(-b/2-2, b/2+2)
        plt.show()

    else:
        pass

    return sweep,taper,b,cr,ct,MAC

if __name__ == "__main__":
    n_engines = 1
    variables = CurrentVariables(n_engines=n_engines)
    A = 12

    sweep, taper, b, cr, ct, MAC = wing_planform(variables,A,drawing=False)

    variables.sweep = sweep
    variables.taper = taper
    variables.b = b
    variables.cr = cr
    variables.ct = ct
    variables.MAC = MAC