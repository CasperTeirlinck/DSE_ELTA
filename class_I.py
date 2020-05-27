# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
from propulsion import *
from variables import *
from math import exp

################################################################################
################################# CLASS I TOOL #################################
################################################################################


g0 = 9.80665


def battery_ratio(R, A, e, E, eff, n=1, SwetS=3.7, Cfe_dict = [0.0055, 0.0045]):
    #R: Range [m], A: Aspect ratio, e: Oswald efficiency, n: number of eng (1 or 2)
    # E: Battery specific energy [Wh/kg], eff: Efficiency of propulsion system, W_PL: Payload [kg]
    #Returns: Take-off, OE, PL and bat weights
    CD_0 = Cfe_dict[n-1]*SwetS
    CD_optr = 2*CD_0
    CL_optr = np.sqrt(np.pi*A*e*CD_0)
    LD_optr = CL_optr/CD_optr
    E = E*3600
    return R/(LD_optr*1/g0*E*eff)

def calc_W_TO(W_PL, WbatWTO, a_classI=0.548, b_classI=55.899):
    # Returns: Take-off, OE, and bat weights in N, W_PL input in kg and WbatWTO dimensionless
    WTO = (W_PL+b_classI)/(1-a_classI-WbatWTO)
    return g0*WTO, g0*(a_classI*WTO+b_classI), g0*WbatWTO*WTO


def iterate_wo(wpl, bmf, WeWo):
    return wpl/(1-bmf-WeWo)


def iterate_wo2(wpl, bmf, WeWo):
    return wpl/(1-bmf-WeWo)


def classIestimation(variables, a=0.656, b=-109, maxiterations=100):
    bat_range = flight_profile_energy_per_WTO(variables, range_m=range_m, endurance_s=0.0)
    bat_endurance = flight_profile_energy_per_WTO(variables, range_m=0.0, endurance_s=endurance_s)
    print(bat_range, bat_endurance)
    variables.bmf = max(bat_range, bat_endurance)
    # variables.bmf = flight_profile_energy_per_WTO(variables, range_m=0.0, endurance_s=0.6*endurance_s)
    # print(variables.bmf)
    if variables.Woew_classII == 0 or variables.Woew_classII == None:
        WeWo = lambda Wo: (a*Wo + b)/Wo
        # WeWo = lambda Wo: (68.3*exp(0.00227*Wo))/Wo
    else:
        WeWo = lambda Wo: variables.Woew_classII/Wo
    Wo_previous = variables.WTO
    Wo_new = iterate_wo(variables.WPL, variables.bmf, WeWo(Wo_previous))
    print(variables.bmf)
    for i in range(maxiterations):
        if abs(Wo_new-Wo_previous)<0.001*Wo_previous:
            variables.WTO = Wo_new
            variables.Woew = a*Wo_new + b
            variables.Wbat = Wo_new*variables.bmf
            print("Class I converged, WTO={}, Woew={}, Wbat={}, BMF={}".format(Wo_new, variables.Woew,
                                                                               variables.Wbat, variables.bmf))
            return variables
        Wo_previous, Wo_new = Wo_new, iterate_wo(variables.WPL, variables.bmf, WeWo(Wo_new))
    else:
        print("Class I did not converge, please check the precision.")
        return variables

def classIestimation_alt(variables, a=0.656, b=-109, maxiterations=100):
    bat_endurance = flight_profile_energy_per_WTO(variables, range_m=0.0, endurance_s=variables.endurance_s)
    bat_range = flight_profile_energy_per_WTO(variables, range_m=variables.range_m, endurance_s=0.0)
    variables.bmf = max(bat_range, bat_endurance)
    # variables.bmf = bat_endurance
    if variables.Woew_classII != 0 and variables.Woew_classII != None:
        variables.Woew = variables.Woew_classII
    else:
        variables.Woew = a*variables.WTO + b
    variables.Wbat = variables.WTO*variables.bmf
    return variables



if __name__ == "__main__":
    v = CurrentVariables()
    print(classIestimation_alt(v))
    
    
    
    



