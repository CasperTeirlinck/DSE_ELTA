from new_main import *
from copy import deepcopy
import numpy as np

if __name__ == "__main__":

    voriginal = NewVariables(True,0.3)
    variation = [5,10]  # [%]

    v = copy.deepcopy(voriginal)

    tweakVars = [
        ['designpointfactor', 'Design point choice, WS vs WP'],
        ['A', 'aspect ratio'],
        ['e', 'oswald efficiency factor'],
        ['A_h', 'aspect ratio htail'],
        ['A_v', 'aspect ratio vtail'],
        ['x_htail', 'aerodynamic center htail'],
        ['x_vtail', 'aerodynamic center htail'],
        ['n_blades', 'nr. of propeller blades'],
        ['eff_propeller', 'propeller eff.'],
        # ['tcr_h', 'Thickness-to-rootchord ratio htail'],
        # ['tcr_v', 'Thickness-to-rootchord ratio vtail'],
        ['CD0', 'Zero-lift drag coefficient'],
        ['V', 'cruise speed'],
        ['CLmaxclean', 'CLmax clean configuration'],
        ['CLmaxto', 'CLmax take-off configuration'],
        ['CLmaxland', 'CLmax landing configuration'],
        ['sto', 'Take-off field length'],
        # ['endurance_s', 'endurance'],
        # ['range_m', 'range']
    ]

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
        ['CD0', 'Zero-lift drag coefficient'],
        ['V', 'cruise speed']
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

    v, _ = do_loop(v)

    print("WTO={} kg, b={} m".format(v.WTO/9.81, v.b))
    v_dict = vars(v)
    print(v_dict)

    def changefunction(key, x1, x2):
        return (getattr(x2, key) - getattr(x1, key)) / getattr(x1, key) * 100

    changevariable = "WTO"
    # changevariable = "b"
    func = partial(changefunction, changevariable)

    for var, label in tweakVars:
        # for var, label in freeVars:
        # for var, label in statVars:
        # for var, label in reqVars:
        #     if var == 'Vmax_kts':
        #         variation[0] = 4.7
        #     else:
        #         variation[0] = 5

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