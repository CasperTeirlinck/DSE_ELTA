from loading_diagram2 import *
from new_variables import *
from Stability_and_Control.horizontaltail_design import *
from Power_and_Propulsion.battery_main import *
from Midterm_Design.ClassIIWeightEstimation import *

wingresolution = 100
iterating_designpoint = False

def do_subloop1(airplane):
    sizing_htail_wingpos(v)
    sys_Aerodynamics_total(v)
    



def do_loop(v):
    if iterating_designpoint:
        get_design_point(v)

    sys_Aerodynamics_wing(v,wingresolution)
    sys_Aerodynamics_total(v)
    v.Wwing = MainWingEstimation(v)
    
    v = main_bat(v)
    v = power_calculation(v)
    
    v = sizing_htail_wingpos(v)
    v = verticaltail_sizing(v)
    v = elevator_sizing(v)


if __name__ == "__main__":
    v = NewVariables(False,0.)
    do_loop(v)