from variables import *
from cg_determination import *

def wing_surface_htail(variables):
    return variables.htail_volume * variables.S * variables.MAC/(variables.x_htail - x_aftCG(variables))

def wing_surface_vtail(variables):
    return variables.vtail_volume * variables.S * variables.b/(variables.x_vtail - x_aftCG(variables))

def empennage_sizing(variables):
    variables.Sh = wing_surface_htail(variables)
    variables.Sv = wing_surface_vtail(variables)
    return variables