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

wingresolution = 100
iterating_designpoint = True


def subloop(v):
    v.W_htail,v.W_vtail = EmpennageEstimation(v)
    v = calcFusgroup(v)

    v = sizing_htail_wingpos(v)
    v = verticaltail_sizing(v)
    return v

def dosubloop(v):
    for i in range(5):
        v = subloop(v)

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
    v = elevator_sizing(v)
    
    v = design_fuselage(v)

    v = CalcOEW(v)
    v.update_WTO()
    return v

def do_loop(v, difference=0.1, maxiterations=500):
    WTOS = []
    for iteration in range(maxiterations):
        if len(WTOS) > 2 and abs(WTOS[-1] - WTOS[-2]) < difference*9.81:
            return v, WTOS, iteration
        else:
            v = loop(v)
            print("\nIteration done, mtom = {} kg\n".format(v.WTO/9.81))
            WTOS.append(v.W_OEW)
    else:
        print("Did not converge within {} iterations".format(maxiterations))
        return v, WTOS, iteration

if __name__ == "__main__":
    v = NewVariables(True,0.3)
    # for i in range(5):
    v, wtos, i = do_loop(v)
    wtos = np.array(wtos)
    v_dict = vars(v)
    print(v_dict)
    print("Final mass = {}, b={}, S={}, WS={}, WP={}".format(v.WTO/9.81, v.b, v.S, v.WS, v.WP))

    print("Iteration ", i + 1, " results:")
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