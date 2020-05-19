from variables import *

def wing_surface_htail(variables):
    return variables.htail_volume * variables.S * variables.MAC/(variables.x_htail - variables.xcg_aft)

def wing_surface_vtail(variables):
    return variables.vtail_volume * variables.S * variables.b/(variables.x_vtail - variables.xcg_aft)

def empennage_sizing(variables):
    variables.Sh = wing_surface_htail(variables)
    print(wing_surface_htail(variables))
    variables.Sv = wing_surface_vtail(variables)
    return variables