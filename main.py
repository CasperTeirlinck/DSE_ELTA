from variables import *
import numpy as np
import math as m
import matplotlib.pyplot as plt
import csv
from collections import OrderedDict

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
from engine import *


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
    v = propeller_weight(v)
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
            # print("Converged! Woew_classII={} N and Woew={} N".format(v.Woew_classII, v.Woew))
            print(f'concept {v.concept_number} converged!')
            return v, OEWS
        else:
            v = loop(v)
            OEWS.append([v.Woew_classII/9.81, v.Woew/9.81])
    else:
        print("Did not converge within {} iterations".format(maxiterations))
        return v, OEWS


if __name__ == "__main__":

    sensAnalysis = True

    if not sensAnalysis:
        conceptnumberlist = [1,2,3,4,5]
        for conceptnumber in conceptnumberlist:

            v = CurrentVariables(*conceptparameters(conceptnumber))
            # v.range_m *= 0.001
            # v.endurance_s *= 0.46
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
            
    else:
        conceptnumber = 3

        v = CurrentVariables(*conceptparameters(conceptnumber))
        v, _ = do_loop(v)
        v_dict = vars(v)

        v1 = CurrentVariables(*conceptparameters(conceptnumber))

        """ Test single variable: """
        """ ===================== """
        change = -10

        ## Reqs
        # v1.WPL = v1.WPL*(1 + change/100)                              # <--!!!
        # v1.Vmax_kts = v1.Vmax_kts*(1 + change/100)
        # v1.rho = v1.rho*(1 + change/100)
        # v1.n_ult = v1.n_ult*(1 + change/100)                          # <--!!!
        # v1.rollrate = v1.rollrate*(1 + change/100)
        # v1.propclear = v1.propclear*(1 + change/100)
        # v1.endurance_s = v1.endurance_s*(1 + change/100)                # <--!!! non linear?
        # v1.range_m = v1.range_m*(1 + change/100)                      # <--!!!


        ## Free variables
        v1.A = v1.A*(1 + change/100)
        # v1.e = v1.e*(1 + change/100)
        # v1.A_h = v1.A_h*(1 + change/100)
        # v1.A_v = v1.A_v*(1 + change/100)
        # v1.x_htail = v1.x_htail*(1 + change/100)
        # v1.x_vtail = v1.x_vtail*(1 + change/100)
        # v1.n_blades = v1.n_blades*(1 + change/100)
        # v1.eff_propeller = v1.eff_propeller*(1 + change/100)
        # v1.rhocruise = v1.rhocruise*(1 + change/100)
        # v1.tcr_h = v1.tcr_h*(1 + change/100)
        # v1.tcr_v = v1.tcr_v*(1 + change/100)
        
        ## Statistical values
        # v1.htail_volume = v1.htail_volume*(1 + change/100)                                      # <--!!! expected
        # v1.vtail_volume = v1.vtail_volume*(1 + change/100)                                          # <--!!! expected
        # v1.Especif_bat = v1.Especif_bat*(1 + change/100)                                        # <--!!! expected ELABORATE
        # v1.motor_spec_mass = v1.motor_spec_mass*(1 + change/100)                                # <--!!! expected
        # v1.CD0clean = v1.CD0clean*(1 + change/100)                                                  # <--!!! expected ELABORATE
        # v1.max_controlsurface_deflection = v1.max_controlsurface_deflection*(1 + change/100)
        # v1.c_l_a_flaps = v1.c_l_a_flaps*(1 + change/100)                                            # <--!!! unexpected ELABORATE
        # v1.c_l_delta_a = v1.c_l_delta_a*(1 + change/100)                                            # <--!!! unexpected ELABORATE
        # v1.fus_height = v1.fus_height*(1 + change/100)
        # v1.prop_spin = v1.prop_spin*(1 + change/100)
        # v1.W_wsn = v1.W_wsn*(1 + change/100)
        # v1.W_wsm = v1.W_wsm*(1 + change/100)
        # v1.W_avion = v1.W_avion*(1 + change/100)
        """ ===================== """

        v1, _ = do_loop(v1)
        v1_dict = vars(v1)
        with open(f'concept_{conceptnumber}_sens.csv', 'w') as file:
            for (key, value), (key1, value1) in zip(v_dict.items(), v1_dict.items()):
                
                if not (isinstance(value, list) or isinstance(value, np.ndarray)):
                    if value:
                        if value != 0:
                            diff = round((value1-value)/value*100, 1)
                            warning = ', <--!!!' if np.abs(diff) >= np.abs(change) else ''
                            if np.abs(diff) >= np.abs(change)/10:
                                file.write(f'{key}, {diff}%{warning}\n')
                        else:
                            file.write(f'{key}, {value1}\n')

        """ Tornado Chart: """
        freeVars = [
            ['A', 'A'],
            ['e', 'e'],
            ['A_h', 'A_htail'],
            ['A_v', 'A_vtail'],
            ['x_htail', 'x_htail'],
            ['x_vtail', 'x_vtail'],
            ['n_blades', 'n_blades'],
            ['eff_propeller', 'eff_propeller'],
            ['rhocruise', 'rhocruise'],
            ['tcr_h', 'tcr_h'],
            ['tcr_v', 'tcr_v']
        ]

        variation = [5, 10] # [%]
        sensDict = {}

        # for var in np.array(freeVars)[:,0]:
        for var, label in freeVars:
            v2 = CurrentVariables(*conceptparameters(conceptnumber))            
            setattr(v2, var, getattr(v2, var)*(1 + variation[0]/100))
            v2, _ = do_loop(v2)
            WTOchangePos1 =  (v2.WTO-v.WTO)/v.WTO*100

            v2 = CurrentVariables(*conceptparameters(conceptnumber))            
            setattr(v2, var, getattr(v2, var)*(1 + variation[1]/100))
            v2, _ = do_loop(v2)
            WTOchangePos2 =  (v2.WTO-v.WTO)/v.WTO*100

            v2 = CurrentVariables(*conceptparameters(conceptnumber))
            setattr(v2, var, getattr(v2, var)*(1 - variation[0]/100))
            v2, _ = do_loop(v2)
            WTOchangeNeg1 =  (v2.WTO-v.WTO)/v.WTO*100

            v2 = CurrentVariables(*conceptparameters(conceptnumber))
            setattr(v2, var, getattr(v2, var)*(1 - variation[1]/100))
            v2, _ = do_loop(v2)
            WTOchangeNeg2 =  (v2.WTO-v.WTO)/v.WTO*100

            sensDict[label] = [[WTOchangePos1, WTOchangeNeg1], [WTOchangePos2, WTOchangeNeg2]]

        # Plotting
        orderedSensDict = OrderedDict(sorted(sensDict.items(), key=lambda item: item[1][-1][0], reverse=False))

        fig = plt.figure()
        ax = fig.add_subplot(111)
        width = 0.4
        xAx = np.arange(len(list(orderedSensDict.keys())))

        ax.barh(xAx, [value[0][0] for value in orderedSensDict.values()], width, color='lightcoral', label='+5%')
        ax.barh(xAx, [value[0][1] for value in orderedSensDict.values()], width, color='cornflowerblue', label='-5%')
        ax.barh(xAx+width, [value[1][0] for value in orderedSensDict.values()], width, color='indianred', label='+10%')
        ax.barh(xAx+width, [value[1][1] for value in orderedSensDict.values()], width, color='royalblue', label='-10%')

        ax.axvline(linewidth=1, color='black')

        ax.set_xlabel('WTO [%]')
        ax.set_yticks(xAx + width / 2)
        ax.set_yticklabels( list(orderedSensDict.keys()) )

        plt.gca().invert_yaxis()
        plt.legend(loc='lower right')
        plt.tight_layout()
        plt.show()


                