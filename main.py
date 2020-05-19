from variables import *
import numpy as np
import math as m
import matplotlib.pyplot as plt
import csv

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
from conceptparameters import conceptparameters


def subloop(v):

    v = cg_calculations_total(v)
    v = size_gear(v)
    v = empennage_sizing(v)
    v = CalculateClassII(v)

    v.update_WTO()

    return v


def dosubloop(v: CurrentVariables, subdifference=0.1, maxsubiterations=40):
    XTAIL = [0]
    for subiteration in range(maxsubiterations):
        XTAIL.append(v.x_htail)
        if abs(XTAIL[subiteration+1] - XTAIL[subiteration]) < subdifference:
            return v            
        else:
            v = subloop(v)
    else: 
        print("Subloop did not converge within {} iterations".format(maxsubiterations))
        return v


def loop(v: CurrentVariables):
    v = get_design_point(v)
    v = classIestimation_alt(v)
    v = motor_mass_volume(v)
    v = wing_planform(v)
    v = calculate_cg_groups(v)
    v = size_control_surfaces(v)

    # Do the c.g. related positioning (wing/tail/gear) subloop
    v = dosubloop(v)

    return v


def do_loop(v: CurrentVariables, difference=0.1, maxiterations=500):
    OEWS = []
    for iteration in range(maxiterations):
        if v.Woew_classII != None and abs(v.Woew - v.Woew_classII) < difference*9.81:
            print("Converged! Woew_classII={} N and Woew={} N".format(v.Woew_classII, v.Woew))
            return v, OEWS
        else:
            v = loop(v)
            OEWS.append([v.Woew_classII/9.81, v.Woew/9.81])
    else:
        print("Did not converge within {} iterations".format(maxiterations))
        return v, OEWS


if __name__ == "__main__":

    conceptnumberlist = [1]

    for conceptnumber in conceptnumberlist:

        v = CurrentVariables(*conceptparameters(conceptnumber))
        v.range_m *= 0.001
        v.endurance_s *= 0.46
        v, oews = do_loop(v)
        oews = np.array(oews)
        v_dict = vars(v)
        print(v_dict)
        with open('concept_{}.csv'.format(conceptnumber), 'w') as file:
            for key in v_dict.keys():
                file.write("%s, %s\n" % (key, v_dict[key]))
        print("Final weight = ",v.WTO/9.81," kg")
        plt.plot(oews[:,0])
        plt.plot(oews[:,1])
        plt.show()
