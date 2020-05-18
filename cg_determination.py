from variables import *

## AIRCRAFT GROUP CONTRIBUTIONS
# Calculate x-coordinate of the fuselage group CG
def X_fuselagegroupCG(variables):
    if variables.wing_mounted_engine:
        return (variables.Wfus * variables.x_cg_fuselage + variables.Wbat * variables.x_cg_battery)/(variables.Wfus + variables.Wbat) 
    else:
        return (variables.Weng * variables.x_cg_engine + variables.Wfus * variables.x_cg_fuselage + variables.Wbat * variables.x_cg_battery)/(variables.Weng + variables.Wfus + variables.Wbat) 

# Calculate the weight of the fuselage group
def fuselagegroup_weight(variables):
    if variables.wing_mounted_engine:
        return variables.Wfus + variables.Wbat
    else:
        return variables.Weng + variables.Wfus + variables.Wbat

# Calculate chordwise location of the wing group CG
def chordwise_winggroupCG(variables):
    if variables.wing_mounted_engine:
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
    return X_fuselagegroupCG(variables) + variables.MAC * (chordwise_winggroupCG(variables) * winggroup_weight(variables)/fuselagegroup_weight(variables) - variables.chordwise_cg_OEW * (1 + winggroup_weight(variables) / fuselagegroup_weight(variables)))

## TAILLESS AIRCRAFT CG
# Calculate x-coordinate of tailless aircraft CG OEW
def x_OEWCG_tailless(variables):
    return (winggroup_weight(variables)*wing_position(variables) + fuselagegroup_weight(variables)*X_fuselagegroupCG(variables)) / (winggroup_weight(variables) + fuselagegroup_weight(variables))

# Calculate most forward CG of tailless aircraft
def x_forwardCG_tailless(variables):
    return min((variables.WPL*variables.x_cg_passenger + (winggroup_weight(variables) + fuselagegroup_weight(variables)) * x_OEWCG_tailless(variables))/ (winggroup_weight(variables) + fuselagegroup_weight(variables) + variables.WPL), x_OEWCG_tailless(variables))

# Calculate most aft CG of tailless aircraft
def x_aftCG_tailless(variables):
    return max((variables.WPL*variables.x_cg_passenger + (winggroup_weight(variables) + fuselagegroup_weight(variables)) * x_OEWCG_tailless(variables)) / (winggroup_weight(variables) + fuselagegroup_weight(variables) + variables.WPL), x_OEWCG_tailless(variables))

def cg_calculations_tailless(variables):
    variables.xcg_frw = x_forwardCG_tailless(variables)
    variables.xcg_aft = x_aftCG_tailless(variables)
    return variables

## TOTAL AIRCRAFT CG
# Calculate x-coordinate of total aircraft CG OEW
def x_OEW_total(variables):
    W_htail = Empennage_Estimation(variables)[0]
    return (winggroup_weight(variables)*wing_position(variables) + fuselagegroup_weight(variables)*X_fuselagegroupCG(variables) + W_htail*variables.x_htail) / (winggroup_weight(variables) + fuselagegroup_weight(variables) + W_htail)

# Calculate most forward CG of total aircraft
def x_forwardCG_total(variables):
    return min((variables.WPL*variables.x_cg_passenger + (winggroup_weight(variables) + fuselagegroup_weight(variables)) * x_OEWCG_total(variables)) / (winggroup_weight(variables) + fuselagegroup_weight(variables) + variables.WPL), x_OEWCG_total(variables))

# Calculate most aft CG of total aircraft
def x_aftCG_total(variables):
    return max((variables.WPL*variables.x_cg_passenger + (winggroup_weight(variables) + fuselagegroup_weight(variables)) * x_OEWCG_total(variables)) / (winggroup_weight(variables) + fuselagegroup_weight(variables) + variables.WPL), x_OEWCG_total(variables))

def cg_calculations_total(variables):
    variables.xcg_frw = x_forwardCG_total(variables)
    variables.xcg_aft = x_aftCG_total(variables)
    return variables    