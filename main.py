from variables import *
import numpy as np
import math as m

from loading_diagram import *
from class_I import *
from wing_planform import *
from AileronflapsizingUSETHISONE import *
from cg_determination import *
from gear import *
from ClassIIWeightEstimation import *
from propulsion import *


def loop(v: CurrentVariables):
    v.WS, v.WP = get_design_point(v)
    v = classIestimation(v)
    v = wing_planform(v)
    # TODO: add Wwing, Wfus, etc. to variables class. Also fix Class I bug.
    # Control surface sizing here
    v = cg_calculations(v)
    v = size_gear(v)
    print("Here")
    v = function2(v)
    v = function3(v)
    v = classII(v)
    return v


def do_loop(v: CurrentVariables, difference=35, maxiterations=10):
    for iteration in range(maxiterations):
        if v.Woew_classII != None and abs(v.Woew - v.Woew_classII) < difference*9.81:
            print(v.Woew)
            print(v.Woew_classII)
            print("This")
            return v
        else:
            v = loop(v)
    else:
        print("Did not converge within {} iterations".format(maxiterations))
        return v


if __name__ == "__main__":

    n_engines = 1
    wing_mounted = False
    T_tail = False
    x_cg_pass = 0.97
    x_cg_batt = 1.96
    x_cg_f = 2.33
    ducted=False


    v = CurrentVariables(n_engines=n_engines, wing_mounted=wing_mounted, T_tail=T_tail, x_cg_pass=x_cg_pass,
                         x_cg_batt=x_cg_batt, x_cg_f=x_cg_f, ducted=ducted)


    v = do_loop(v)
    print("Done!")