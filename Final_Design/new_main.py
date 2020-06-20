from loading_diagram2 import *
from new_variables import *
from Stability_and_Control.horizontaltail_design import *
from Stability_and_Control.verticaltail_design import *
from Stability_and_Control.empennagecontrol_design import *
from Power_and_Propulsion.battery_main import *
from ClassIIWeightEstimationnew import *

wingresolution = 100
iterating_designpoint = True

def do_subloop1(airplane):
    sizing_htail_wingpos(v)
    sys_Aerodynamics_total(v)
    



def do_loop(v):
    if iterating_designpoint:
        get_design_point(v)

    sys_Aerodynamics_wing(v ,wingresolution)
    sys_Aerodynamics_total(v)
    v.W_wing = MainWingEstimationNew(v)

    v = power_calculation(v)
    v = main_bat(v)


    v.W_htail,v.W_vtail = EmpennageEstimation(v)


    v = sizing_htail_wingpos(v)
    v = verticaltail_sizing(v)
    v = elevator_sizing(v)
    v.W_ = EmpennageEstimation(v)


if __name__ == "__main__":
    v = NewVariables(False,0.)
    print(v.Snet)
    do_loop(v)
    print("IT WORKS!!")