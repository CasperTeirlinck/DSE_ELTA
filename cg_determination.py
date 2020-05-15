import numpy as np 
from variables import *


# Calculate x-coordinate of the fuselage group CG
def X_fuselagegroupCG(variables):
    if variables.wing_mounted_engine:
        return (variables.Wfus * variables.x_cg_fuselage + variables.Wbat * x_cg_battery)/(variables.Wfus + variables.Wbat) 
    else:
        return (variables.Weng * variables.x_cg_engine + variables.Wfus * variables.x_cg_fuselage + variables.Wbat * x_cg_battery)/(Weng + variables.Wfus + variables.Wbat) 

# Calculate the weight of the fuselage group
def fuselagegroup_weight(variables):
    if variables.wing_mounted_engine:
        return variables.Wfus + variables.Wbat
    else:
        return Weng + variables.Wfus + variables.Wbat

# Calculate chordwise location of the wing group CG
def chordwise_winggroupCG(variables):
    if variables.engine_on_wings:
        return (variables.Wwing*variables.chordwise_wing_cg + variables.Weng*variables.chordwise_cg_engine)/(variables.Wwing + variables.Weng)
    else:
        return variables.chordwise_wing_cg

# Calculate the weight of the wing group
def winggroup_weight(variables):
    if variables.wing_mounted_engine:
        return variables.Wwing + variables.Weng
    else:
        return variables.Wwing


# Calculate exact wing position for a chosen chordwise OEW CG
def wing_position(variables):
    return fuselagegroupCG(variables) + variables.MAC * (chordwise_winggroupCG(variables) * winggroup_weight(variables)/fuselagegroup_weight(variables) - variables.chordwise_cg_OEW * (1 + winggroup_weight(variables) / fuselagegroup_weight(variables)))


# Calculate x-coordinate of total aircraft CG OEW
def x_OEWCG(variables):
    return (winggroup_weight(variables)*wing_position(variables) + fuselagegroup_weight(variables)*X_fuselagegroupCG(variables)) / (winggroup_weight(variables) + fuselagegroup_weight(variables))



def x_forwardCG(variables):
    return min((variables.WPL*variables.x_cg_passenger + (winggroup_weight(variables) + fuselagegroup_weight(variables)) * x_OEWCG(variables) )/ (winggroup_weight(variables) + fuselagegroup_weight(variables) + variables.WPL), x_OEWCG)

def x_aftCG(variables):
    return max((variables.WPL*variables.x_cg_passenger + (winggroup_weight(variables) + fuselagegroup_weight(variables) * x_OEWCG(variables) / (winggroup_weight(variables) + fuselagegroup_weight(variables) + variables.WPL), x_OEWCG)

