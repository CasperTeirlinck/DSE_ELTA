from collections import OrderedDict
import matplotlib.pyplot as plt
import numpy as np

from new_main import do_loop
from new_variables import NewVariables

def sensAnalysis(refVar='WTO'):
    # sensVars = [
    #     ['A', 'aspect ratio'],
    #     ['e', 'oswald efficiency factor'],
    #     ['A_h', 'aspect ratio htail'],
    #     ['A_v', 'aspect ratio vtail'],
    #     ['x_htail', 'aerodynamic center htail'],
    #     ['x_vtail', 'aerodynamic center htail'],
    #     ['n_blades', 'nr. of propeller blades'],
    #     ['eff_propeller', 'propeller eff.'],
    #     ['tcr_h', 'Thickness-to-rootchord ratio htail'],
    #     ['tcr_v', 'Thickness-to-rootchord ratio vtail'],
    #     ['V', 'cruise speed']
    # ]
    sensVars = [
        # ['designpointfactor', 'Design point choice, WS vs WP'],
        # ['V', 'cruise speed'],
        # ['sto', 'Take-off field length'],
        # ['eff_prop', 'Propeller efficiency'],
        # ['endurance_s', 'endurance'],
        # ['range_m', 'range'],
        # ['batt_cell_C_Ah', 'specific E batteries'],
        # ['DoD', 'Depth of discharge'],
        ['hwl', 'Winglet height'],
        # ['xtail', 'Tail position from nose']
    ]

    variation = [5, 10] # [%]
    sensDict = {}

    v = NewVariables(True, 0.3)
    v, _, _ = do_loop(v)
    refValue = getattr(v, refVar)

    def _do_loop_sens(var: str, change: int):
        if var != 'hwl':
            v2 = NewVariables(True, 0.3)
            setattr(v2, var, getattr(v2, var)*(1 + change/100))
            v2, _, _ = do_loop(v2)
            return (getattr(v2, refVar) - refValue) / refValue *100
        else:
            v2 = NewVariables(True, 0.3*(1+change/100))
            print(v2.winglet)
            v2, _, _ = do_loop(v2)
            return (getattr(v2, refVar) - refValue) / refValue *100

    for var, label in sensVars:
        print("Now doing var={}, label={}".format(var, label))
        refValueChange_pos_1 = _do_loop_sens(var, variation[0])
        refValueChange_pos_2 = _do_loop_sens(var, variation[1])
        refValueChange_neg_1 = _do_loop_sens(var, -variation[0])
        refValueChange_neg_2 = _do_loop_sens(var, -variation[1])

        sensDict[label] = [[refValueChange_pos_1, refValueChange_pos_2], [refValueChange_neg_1, refValueChange_neg_2]]

    orderedSensDict = OrderedDict(sorted(sensDict.items(), key=lambda item: item[1][-1][0], reverse=False))

    """ PLOTTING """
    fig = plt.figure(figsize=(6, 4.5))
    ax = fig.add_subplot(111)
    width = 0.2
    xAx = np.arange(len(list(orderedSensDict.keys())))

    ax.barh(xAx, [value[0][0] for value in orderedSensDict.values()], width, color='lightcoral', label=f'+{variation[0]}%')
    ax.barh(xAx+width, [value[0][1] for value in orderedSensDict.values()], width, color='indianred', label=f'+{variation[1]}%')
    ax.barh(xAx+width*2, [value[1][0] for value in orderedSensDict.values()], width, color='cornflowerblue', label=f'-{variation[0]}%')
    ax.barh(xAx+width*3, [value[1][1] for value in orderedSensDict.values()], width, color='royalblue', label=f'-{variation[1]}%')

    ax.axvline(linewidth=1, color='black')

    ax.set_xlabel(f'{refVar} [%]')
    ax.set_yticks(xAx + width / 2)
    ax.set_yticklabels( list(orderedSensDict.keys()) )
    ax.xaxis.grid(color='black', linestyle='--')

    plt.gca().invert_yaxis()
    # plt.legend(loc='lower right')
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.show()
    print("End")

if __name__ == "__main__":
    sensAnalysis("WTO")

