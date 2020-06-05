import Midterm_Design
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
from main import *

if __name__ == "__main__":
    conceptnumber = 3

    v = CurrentVariables(*conceptparameters(conceptnumber))
    v1 = CurrentVariables(*conceptparameters(conceptnumber))

    basefactor = 0.765*0.7
    variation = [1,2]  # [%]

    v.designpointfactor = basefactor
    v.designpointfactor = basefactor

    v, _ = do_loop(v)
    v_dict = vars(v)

    freeVars = [
        ['A', 'aspect ratio'],
        ['e', 'oswald efficiency factor'],
        ['A_h', 'aspect ratio htail'],
        ['A_v', 'aspect ratio vtail'],
        ['x_htail', 'aerodynamic center htail'],
        ['x_vtail', 'aerodynamic center htail'],
        ['n_blades', 'nr. of propeller blades'],
        ['eff_propeller', 'propeller eff.'],
        ['tcr_h', 'Thickness-to-rootchord ratio htail'],
        ['tcr_v', 'Thickness-to-rootchord ratio vtail'],
        ['V', 'cruise speed'],
        ['designpointfactor', 'Design point choice, WS vs WP']
    ]

    statVars = [
        ['htail_volume', 'htail volume'],
        ['vtail_volume', 'vtail volume'],
        ['Especif_bat', 'specific E batteries'],
        ['motor_spec_mass', 'battery specific mass'],
        ['CD0clean', 'drag constant clean'],
        ['W_wsn', 'nose wheel weight'],
        ['W_wsm', 'main wheel weight']
    ]

    reqVars = [
        ['WPL', 'payload weight'],
        ['Vmax_kts', 'max speed'],
        ['n_ult', 'ultimate load factor'],
        ['endurance_s', 'endurance'],
        ['range_m', 'range']
    ]

    sensDict = {}

    for var, label in freeVars:
    # for var, label in statVars:
    # for var, label in reqVars:
    #     if var == 'Vmax_kts':
    #         variation[0] = 4.7
    #     else:
    #         variation[0] = 5

        v2 = CurrentVariables(*conceptparameters(conceptnumber))
        setattr(v2, var, getattr(v2, var) * (1 + variation[0] / 100))
        v2, _ = do_loop(v2)
        WTOchangePos1 = (v2.WTO - v.WTO) / v.WTO * 100

        v2 = CurrentVariables(*conceptparameters(conceptnumber))
        setattr(v2, var, getattr(v2, var) * (1 + variation[1] / 100))
        v2, _ = do_loop(v2)
        WTOchangePos2 = (v2.WTO - v.WTO) / v.WTO * 100

        v2 = CurrentVariables(*conceptparameters(conceptnumber))
        setattr(v2, var, getattr(v2, var) * (1 - variation[0] / 100))
        v2, _ = do_loop(v2)
        WTOchangeNeg1 = (v2.WTO - v.WTO) / v.WTO * 100

        v2 = CurrentVariables(*conceptparameters(conceptnumber))
        setattr(v2, var, getattr(v2, var) * (1 - variation[1] / 100))
        v2, _ = do_loop(v2)
        WTOchangeNeg2 = (v2.WTO - v.WTO) / v.WTO * 100

        sensDict[label] = [[WTOchangePos1, WTOchangeNeg1], [WTOchangePos2, WTOchangeNeg2]]

    # Plotting
    orderedSensDict = OrderedDict(sorted(sensDict.items(), key=lambda item: item[1][-1][0], reverse=False))

    fig = plt.figure(figsize=(6, 4.5))
    ax = fig.add_subplot(111)
    width = 0.2
    xAx = np.arange(len(list(orderedSensDict.keys())))

    ax.barh(xAx, [value[0][0] for value in orderedSensDict.values()], width, color='lightcoral', label='+{}%'.format(variation[0]))
    ax.barh(xAx + width, [value[0][1] for value in orderedSensDict.values()], width, color='cornflowerblue', label='-{}%'.format(variation[0]))
    ax.barh(xAx + width*2, [value[1][0] for value in orderedSensDict.values()], width, color='indianred', label='+{}%'.format(variation[1]))
    ax.barh(xAx + width*3, [value[1][1] for value in orderedSensDict.values()], width, color='royalblue', label='-{}%'.format(variation[1]))

    ax.axvline(linewidth=1, color='black')

    ax.set_xlabel('WTO [%]')
    ax.set_yticks(xAx + width / 2)
    ax.set_yticklabels(list(orderedSensDict.keys()))
    ax.xaxis.grid(color='black', linestyle='--')
    fig.suptitle('Original WTO={} kg'.format(v.WTO/9.81), fontsize=16)

    plt.gca().invert_yaxis()
    # plt.legend(loc='lower right')
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.show()