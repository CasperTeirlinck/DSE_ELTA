from Midterm_Design.variables import *
from Midterm_Design.ClassIIWeightEstimation import FuselageEstimation, MainWingEstimation, EngineEstimation

## AIRCRAFT GROUP CONTRIBUTIONS
# Calculate x-coordinate of the fuselage group CG
def X_fuselagegroupCG(variables):
    if variables.wing_mounted_engine:
        return (FuselageEstimation(variables) * variables.fuselagecg_x + variables.Wbat * variables.batterycg_x)/(FuselageEstimation(variables) + variables.Wbat)
    else:
        return (EngineEstimation(variables) * variables.enginecg_x + FuselageEstimation(variables) * variables.fuselagecg_x + variables.Wbat * variables.batterycg_x)/(EngineEstimation(variables) + FuselageEstimation(variables) + variables.Wbat)

# Calculate the weight of the fuselage group
def fuselagegroup_weight(variables):
    if variables.wing_mounted_engine:
        return FuselageEstimation(variables) + variables.Wbat
    else:
        return EngineEstimation(variables) + FuselageEstimation(variables) + variables.Wbat

# Calculate chordwise location of the wing group CG
def chordwise_winggroupCG(variables):
    if variables.wing_mounted_engine:
        return (MainWingEstimation(variables)*variables.chordwise_wing_cg + EngineEstimation(variables)*variables.chordwise_cg_engine)/(MainWingEstimation(variables) + EngineEstimation(variables))
    else:
        return variables.chordwise_wing_cg

# Calculate the weight of the wing group
def winggroup_weight(variables):
    if variables.wing_mounted_engine:
        return MainWingEstimation(variables) + EngineEstimation(variables)
    else:
        return MainWingEstimation(variables)

# Calculate exact wing position for a chosen chordwise OEW CG
def wing_position(variables):
    return X_fuselagegroupCG(variables) + variables.MAC * (chordwise_winggroupCG(variables) * winggroup_weight(variables)/fuselagegroup_weight(variables) - variables.chordwise_cg_oew * (1 + winggroup_weight(variables) / fuselagegroup_weight(variables)))

## TAILLESS AIRCRAFT CG
# Calculate x-coordinate of tailless aircraft CG OEW
def x_OEWCG_tailless(variables):
    return (winggroup_weight(variables)*wing_position(variables) + fuselagegroup_weight(variables)*X_fuselagegroupCG(variables)) / (winggroup_weight(variables) + fuselagegroup_weight(variables))

# Calculate most forward CG of tailless aircraft
def x_forwardCG_tailless(variables):
    return min((variables.WPL*variables.payloadcg_x + (winggroup_weight(variables) + fuselagegroup_weight(variables)) * x_OEWCG_tailless(variables))/ (winggroup_weight(variables) + fuselagegroup_weight(variables) + variables.WPL), x_OEWCG_tailless(variables))

# Calculate most aft CG of tailless aircraft
def x_aftCG_tailless(variables):
    return max((variables.WPL*variables.payloadcg_x + (winggroup_weight(variables) + fuselagegroup_weight(variables)) * x_OEWCG_tailless(variables)) / (winggroup_weight(variables) + fuselagegroup_weight(variables) + variables.WPL), x_OEWCG_tailless(variables))

def cg_calculations_tailless(variables):
    variables.xcg_frw = x_forwardCG_tailless(variables)
    variables.xcg_aft = x_aftCG_tailless(variables)
    variables.tail_ready = True
    return variables

## TOTAL AIRCRAFT CG
# Calculate x-coordinate of total aircraft CG OEW
def x_OEWCG_total(variables):   
    return (winggroup_weight(variables)*wing_position(variables) + fuselagegroup_weight(variables)*X_fuselagegroupCG(variables) + variables.W_htail*variables.x_htail + variables.W_vtail*variables.x_vtail + variables.Wgear_main*variables.x_maingear + variables.Wgear_front*variables.x_nosegear) / (winggroup_weight(variables) + fuselagegroup_weight(variables) + variables.W_htail + variables.W_vtail + variables.Wgear_front + variables.Wgear_main)

# Calculate most forward CG of total aircraft
def x_forwardCG_total(variables):
    return min((variables.WPL*variables.payloadcg_x + (winggroup_weight(variables) + fuselagegroup_weight(variables) + variables.W_htail + variables.W_vtail + variables.Wgear_main + variables.Wgear_front) * x_OEWCG_total(variables)) / (winggroup_weight(variables) + fuselagegroup_weight(variables) + variables.WPL + variables.W_htail + variables.W_vtail + variables.Wgear_main + variables.Wgear_front), x_OEWCG_total(variables))

# Calculate most aft CG of total aircraft
def x_aftCG_total(variables):
    return max((variables.WPL*variables.payloadcg_x + (winggroup_weight(variables) + fuselagegroup_weight(variables) + variables.W_htail + variables.W_vtail + variables.Wgear_main + variables.Wgear_front) * x_OEWCG_total(variables)) / (winggroup_weight(variables) + fuselagegroup_weight(variables) + variables.WPL + variables.W_htail + variables.W_vtail + variables.Wgear_main + variables.Wgear_front), x_OEWCG_total(variables))

def cg_calculations_total(variables):
    if variables.tail_ready:
        variables.xcg_frw = x_forwardCG_total(variables)
        variables.xcg_aft = x_aftCG_total(variables)
        # print(variables.xcg_aft)
        return variables
    else:
        return cg_calculations_tailless(variables)


## TESTING
if __name__ == "__main__":
    # Test variables
    variables = CurrentVariables()
    variables.Wbat = 200
    variables.Wwing = 100 
    variables.Wfus = 300
    variables.x_cg_passenger = 1.97
    variables.x_cg_battery = 1.96
    variables.x_cg_fuselage = 2.33
    variables.chordwise_cg_oew = 0.25

    variables.x_htail = 6
    variables.x_vtail = variables.x_htail
    W_htail = 35
    W_vtail = 20

    # Test:                                 # Correct behaviour indicating verification:
    print(wing_position(variables))         # Yields realistic, feasible and consistent values
    print(x_OEWCG_tailless(variables))      # Yields realistic, feasible and consistent values, placed behind wing, moves correctly with wing&battery position, does not move with passenger CG
    print(x_forwardCG_tailless(variables))  # Yields realistic, feasible and consistent values, placed behind wing, moves correctly with wing, battery and passenger position,  either forward or aft cg is always equal to OEW CG, always behind aft cg
    print(x_aftCG_tailless(variables))      # Yields realistic, feasible and consistent values, placed behind wing, moves correctly with wing, battery and passenger position,  either forward or aft cg is always equal to OEW CG, always further than forward cg
    print()                                 # -
    print(x_OEWCG_total(variables))         # Yields realistic, feasible and consistent values, placed behind wing, moves correctly with wing and battery position, always further aft than tailless CG
    print(x_forwardCG_total(variables))     # Yields realistic, feasible and consistent values, placed behind wing, moves correctly with wing, battery and passenger position, either forward or aft cg is always equal to OEW CG, always behind aft cg, always further aft than tailless CG
    print(x_aftCG_total(variables))         # Yields realistic, feasible and consistent values, placed behind wing, moves correctly with wing, battery and passenger position, either forward or aft cg is always equal to OEW CG, always further than forward cg, always further aft than tailless CG

    # Untested: Tail weight integration