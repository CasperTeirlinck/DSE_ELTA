from new_main import *

from copy import *
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    voriginal = NewVariables(True,0.3)
    variation = [5,10]  # [%]

    v = copy.deepcopy(voriginal)

    tweakVars = [
        ['designpointfactor', 'Design point choice, WS vs WP'],
        ['A_v', 'aspect ratio vtail'],
        ['V', 'cruise speed'],
        ['sto', 'Take-off field length'],
        ['endurance_s', 'endurance'],
        ['range_m', 'range'],
        ['batt_cell_E_spec', 'specific E batteries'],
    ]

    sensDict = {}

    v, _ = do_loop(v)

    print("WTO={} kg, b={} m".format(v.WTO/9.81, v.b))
    v_dict = vars(v)

    def changefunction(key, x1, x2):
        return (getattr(x2, key) - getattr(x1, key)) / getattr(x1, key) * 100

    changevariable = "WTO"
    # changevariable = "b"
    func = partial(changefunction, changevariable)

    for var, label in tweakVars:

        v2 = copy.deepcopy(voriginal)
        setattr(v2, var, getattr(v2, var) * (1 + variation[0] / 100))
        v2, _ = do_loop(v2)
        WTOchangePos1 = func(v, v2)

        v2 = copy.deepcopy(voriginal)
        setattr(v2, var, getattr(v2, var) * (1 + variation[1] / 100))
        v2, _ = do_loop(v2)
        WTOchangePos2 = func(v, v2)

        v2 = copy.deepcopy(voriginal)
        setattr(v2, var, getattr(v2, var) * (1 - variation[0] / 100))
        v2, _ = do_loop(v2)
        WTOchangeNeg1 = func(v, v2)

        v2 = copy.deepcopy(voriginal)
        setattr(v2, var, getattr(v2, var) * (1 - variation[1] / 100))
        v2, _ = do_loop(v2)
        WTOchangeNeg2 = func(v, v2)

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

    ax.set_xlabel('{} [%]'.format(changevariable))
    ax.set_yticks(xAx + width / 2)
    ax.set_yticklabels(list(orderedSensDict.keys()))
    ax.xaxis.grid(color='black', linestyle='--')
    # fig.suptitle('Original WTO={} kg'.format(v.WTO/9.81), fontsize=16)

    plt.gca().invert_yaxis()
    # plt.legend(loc='lower right')
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.show()