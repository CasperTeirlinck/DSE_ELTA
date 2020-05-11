import numpy as np

def wing_planform(S,A):
    # Sweep
    sweep = np.rad2deg(np.arccos(1)) # [deg]

    # Taper
    taper = 0.2*(2-sweep*np.pi/180)

    # Span
    b = np.sqrt(S*A)

    # Chord
    cr = 2*S/((1+taper)*b)
    ct = taper*cr
    mac = (2/3)*cr*((1+taper+taper**2)/(1+taper))

    return sweep,taper,b,cr,ct,mac