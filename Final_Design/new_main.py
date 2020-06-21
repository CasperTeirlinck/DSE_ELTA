from loading_diagram2 import *
from new_variables import *
from Stability_and_Control.horizontaltail_design import *
from Stability_and_Control.verticaltail_design import *
from Stability_and_Control.empennagecontrol_design import *
from Power_and_Propulsion.battery_main import *
from ClassIIWeightEstimationnew import *
from Midterm_Design.cg_determination import *
 
wingresolution = 100
iterating_designpoint = True

    



def do_loop(v):
    if iterating_designpoint:
        get_design_point(v)

    sys_Aerodynamics_wing(v ,wingresolution)
    sys_Aerodynamics_total(v)
    v.W_wing = MainWingEstimationNew(v)

    v = power_calculation(v)
    v = main_bat(v)

    v.W_htail,v.W_vtail = EmpennageEstimation(v)
    v = calcFusgroup(v)


    v = sizing_htail_wingpos(v)
    v = verticaltail_sizing(v)
    v = CalcTTO(v)
    v = elevator_sizing(v)    
    



    v = CalcOEW(v)
    v = CalcMTOWnew(v)

if __name__ == "__main__":
    v = NewVariables(False,0.)
    print(v.Snet)
    for i in range(5):
        do_loop(v)
        print(v.WTO/9.81)

    print("IT WORKS!!")