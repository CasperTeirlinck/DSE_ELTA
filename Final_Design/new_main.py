from loading_diagram2 import *
from new_variables import *


wingresolution = 100
iterating_designpoint = True

def do_loop(airplane):
    if iterating_designpoint:
        get_design_point(airplane)

    sys_Aerodynamics_wing(airplane,wingresolution)
    
