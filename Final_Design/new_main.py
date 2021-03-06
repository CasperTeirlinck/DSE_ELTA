import numpy as np

from loading_diagram2 import *
from new_variables import *
from Stability_and_Control.horizontaltail_design import *
from Stability_and_Control.verticaltail_design import *
from Stability_and_Control.empennagecontrol_design import *
from Power_and_Propulsion.battery_main import *
from ClassIIWeightEstimationnew import *
from Midterm_Design.cg_determination import *
from Structures.fuselagedesign import *
from Structures.landing_gear_pos import *
from Structures.stress_analysis import *

wingresolution = 100
iterating_designpoint = False


def subloop(v):
    v.W_htail,v.W_vtail = EmpennageEstimation(v)
    v = calcFusgroup(v)

    v = sizing_htail_wingpos(v)
    v = verticaltail_sizing(v)
    return v

def dosubloop(v):
    for i in range(5):
        v = subloop(v)
    return v

def loop(v):
    if iterating_designpoint:
        get_design_point(v)

    sys_Aerodynamics_wing(v ,wingresolution)
    sys_Aerodynamics_total(v)
    v.W_wing = MainWingEstimationNew(v)

    v = power_calculation(v)
    v = main_bat(v)

    v = dosubloop(v)
    v = CalcTTO(v)
    v = size_gear(v)
    v = elevator_sizing(v)

    v = design_fuselage(v)
    v = size_wing(v)

    v = CalcOEW(v)
    v.update_WTO()
    return v

def do_loop(v, difference=1, maxiterations=30):
    WTOS = []
    for iteration in range(maxiterations):
        if len(WTOS) > 2 and abs(WTOS[-1] - WTOS[-2]) < difference*9.81:
            return v, WTOS, iteration
        else:
            v = loop(v)
            print("Iteration {} done, mtom = {} kg".format(iteration, v.WTO/9.81))
            # print("skin_t={}, n_stiff={}, n_stiff_circ={}, stringermod={}, circstringermod={}, longeronmod={}".format(v.skin_t, v.n_stiff, v.n_stiff_circ, v.stringermod, v.circstringermod, v.longeronmod))
            WTOS.append(v.W_OEW)
    else:
        print("Did not converge within {} iterations".format(maxiterations))
        return v, WTOS, iteration

if __name__ == "__main__":
    v = NewVariables(True,0.51)
    var = 'endurance_s'
    change = -35
    # change = 0
    setattr(v, var, getattr(v, var)*(1 + change/100))
    # print(v.endurance_s)
    # print("Original mass = 875.3897121968254 kg")
    # for i in range(5):
    v, wtos, i = do_loop(v)
    wtos = np.array(wtos)
    v_dict = vars(v)
    fus_dict = vars(v.fuselage)
    print(v_dict)
    print("Final mass = {}, b={}, S={}, WS={}, WP={}".format(v.WTO/9.81, v.b, v.S, v.WS, v.WP))

    print("Iteration ", i + 1, " results:")
    print("Oswald:", v.eclean)
    print("Take-off: ", v.WTO / 9.81)
    print("Battery: ", v.W_batt / 9.81)
    print("Horizontal tail: ", v.W_htail / 9.81)
    print("Vertical tail: ", v.W_vtail / 9.81)
    print("Wing: ", v.W_wing / 9.81)
    print("Fuselage aft: ", v.Wfus_aft / 9.81)
    print("Take-off weight history: ", wtos)
    print("_______________")

    with open('finaldesign.csv', 'w') as file:
        for key in v_dict.keys():
            file.write("%s, %s\n" % (key, v_dict[key]))
    print("Done!")