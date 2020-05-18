from variables import *
import numpy as np
import math as m

from loading_diagram import *
from class_I import *
from wing_planform import *
from FuselageCGs import *
from AileronflapsizingUSETHISONE import *
from cg_determination import *
from gear import *
from ClassIIWeightEstimation import *
from propulsion import *
from Empennage import *


def subloop(v):
    v = cg_calculations_total(v)
    v = size_gear(v)
    v = empennage_sizing(v)
    v = classIIestimation(v)
    return v


def dosubloop(v: CurrentVariables, subdifference=0.05, maxsubiterations=20):
    XTAIL = [0]
    for subiteration in range(maxsubiterations):
        XTAIL.append(v.x_htail)
        if abs(XTAIL[subiteration+1] - XTAIL[subiteration]) < subdifference:
            return v            
        else:
            v = subloop(v)
    else: 
        print(print("Subloop did not converge within {} iterations".format(maxsubiterations)))
        return v


def loop(v: CurrentVariables):
    v.WS, v.WP = get_design_point(v)
    v = classIestimation_alt(v)
    v = wing_planform(v)
    v = calculate_cg_groups(v)
    # TODO: Discuss TE weight addition in code.
    print("before")
    v = size_control_surfaces(v)
    print("After")
    v = subloop(v)
    v = dosubloop(v)    
    print("here2")
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

    conceptnumber = 1
    n_engines = 1
    wing_mounted = False
    T_tail = False
    x_cg_pass = 0.97
    x_cg_batt = 1.96
    x_cg_f = 2.33
    ducted=False


    v = CurrentVariables(conceptnumber=conceptnumber, wing_mounted=wing_mounted, T_tail=T_tail, x_cg_pass=x_cg_pass,
                         x_cg_batt=x_cg_batt, x_cg_f=x_cg_f, ducted=ducted)


    v = do_loop(v)
    print("Done!")