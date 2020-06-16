import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os

def optimize():

    CLa_h_range = []
    A_h_range = np.linspace(1, 4, 100)
    sweepHalfc_h_range = np.radians(np.linspace(0, 10, 100))
    # for A_h in A_h_range:
    for sweepHalfc_h in sweepHalfc_h_range:
        A_h = 2
        # sweepHalfc_h = np.radians(0)
        empng = Empennage(A_h, sweepHalfc_h)

        CLa_h_range.append( empng.CLa_h() )

    fig = plt.figure(figsize=(10, 4.5))
    ax1 = fig.add_subplot(111)

    ax1.plot(A_h_range, CLa_h_range)

    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.show()

class Empennage:
    def __init__(self, A_h, sweepHalfc_h):
        self.A_h = A_h
        self.sweepHalfc_h = sweepHalfc_h

    def setAirfoils(self, Cla_h):
        self.Cla_h = Cla_h

    def CLa_h(self):
        beta = 1
        eta_h = 0.95
        return 2*np.pi*self.A_h / ( 2 + np.sqrt( 4 + (self.A_h*beta/eta_h)**2 * ( 1 + (np.tan(self.sweepHalfc_h)/beta)**2 ) ) )


if __name__ == "__main__":
    optimize()
    
    
