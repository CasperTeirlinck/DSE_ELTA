import numpy as np
from math import pi, sqrt
from stringerclass import *


def make_stringers(positions, stringerinstance):
    stringers = []
    for pos in list(positions):
        stringers.append(Stringer(*pos, stringerinstance))
    return stringers



def make_panel(panel, *args):
    for arg in args:
        panel.addStringer(arg)
    return panel


class Sheet:
    def __init__(self, ts):
        self.ts = ts
        self.stringers = []

    def addStringer(self, stringer: Stringer):
        self.stringers.append(stringer)

    # def calc_we(self, C=4.0):
    #     for stringer in Stringers:
    #         self.ts * sqrt(C*pi*pi/)

if __name__ == "__main__":
    print("Hello world!")
    aluminium = MatProps(sigma_y=450000000, E=72000000000, poisson=0.3, alpha=0.8, n=0.6)

    j = J_Stringer(Le=0.4, material=aluminium)
    stringers = make_stringers([(1,2),(2,3),(3,4)], j)

    panel = Sheet(ts = 0.001)

    panel = make_panel(panel, *stringers)
    print(panel.stringers[0].properties.)
