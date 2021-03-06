# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from Midterm_Design.propulsion import *
from Midterm_Design.variables import *

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

def classIestimation_alt(variables, a=0.656, b=-109):
    bat_endurance = flight_profile_energy_per_WTO(variables, range_m=0.0, endurance_s=variables.endurance_s)
    bat_range = flight_profile_energy_per_WTO(variables, range_m=variables.range_m, endurance_s=0.0)
    variables.bmf = max(bat_range, bat_endurance)
    # variables.bmf = bat_endurance
    if variables.Woew_classII != 0 and variables.Woew_classII != None:
        variables.Woew = variables.Woew_classII
        if bat_range > bat_endurance:
            variables.Wbat = (variables.WTO-variables.WPL*0.5)*variables.bmf
        else:
            variables.Wbat = variables.WTO*variables.bmf
    else:
        if bat_range > bat_endurance:
            variables.Woew = a*(variables.WTO-0.5*variables.WPL) + b
            variables.Wbat = (variables.WTO-variables.WPL*0.5)*variables.bmf
        else:
            variables.Woew = a*variables.WTO + b
            variables.Wbat = variables.WTO*variables.bmf
    return variables


def classIestimation_alt3(variables, a=0.656, b=-109):
    variables.bmf_e = flight_profile_energy_per_WTO(variables, range_m=0.0, endurance_s=variables.endurance_s)
    variables.bmf_r = flight_profile_energy_per_WTO(variables, range_m=variables.range_m, endurance_s=0.0)
    if variables.Woew_classII != 0 and variables.Woew_classII != None:
        variables.Woew = variables.Woew_classII
    else:
        variables.Woew = a*variables.WTO + b

    if variables.batterypilot:
        batteryrange = variables.WTO*variables.bmf_r-100*9.81
    else:
        batteryrange = (variables.WTO-variables.WPL*0.5)*variables.bmf_r
    batteryendurance = variables.WTO*variables.bmf_e
    variables.Wbat = max(batteryrange, batteryendurance)
    return variables


def classIestimation_alt2(variables, a=0.656, b=-109):
    if variables.rangecritical == True:
        variables.bmf = flight_profile_energy_per_WTO(variables, range_m=variables.range_m, endurance_s=0.0)
    else:
        variables.bmf = flight_profile_energy_per_WTO(variables, range_m=0.0, endurance_s=variables.endurance_s)

    if variables.Woew_classII != 0 and variables.Woew_classII != None:
        variables.Woew = variables.Woew_classII
        variables.Wbat = variables.WTO * variables.bmf
    else:
        variables.Woew = a * variables.WTO + b
        variables.Wbat = variables.WTO * variables.bmf
    return variables


# def check_range_endurance(variables):
#     bat_range = flight_profile_energy_per_WTO(variables, range_m=variables.range_m, endurance_s=0.0)
#     bat_endurance = flight_profile_energy_per_WTO(variables, range_m=0.0, endurance_s=variables.endurance_s)
#     if variables.rangecritical:
#         realbmfrange = variables.Wbat/variables.WTO
#         realbmfendurance = variables.Wbat/(variables.WTO+100*9.81)
#         if realbmf < bat_range:
#             print("Designed for range, but range requirement not met. Real bmf={}, target bmf={}".format(realbmf, bat_range))
#
#             print("Designed for range, but range requirement not met. Real bmf={}, target bmf={}".format(realbmf, bat_range))




    # bat_endurance = flight_profile_energy_per_WTO(variables, range_m=0.0, endurance_s=variables.endurance_s)
    # bat_range = flight_profile_energy_per_WTO(variables, range_m=variables.range_m, endurance_s=0.0)
    # # variables.bmf = bat_endurance
    # if variables.Woew_classII != 0 and variables.Woew_classII != None:
    #     variables.Woew = variables.Woew_classII
    # else:
    #     if bat_range > bat_endurance:
    #         variables.rangecritical = True
    #         variables.Woew = a * (variables.WTO - 0.5 * variables.WPL) + b
    #         variables.Wbat = (variables.WTO - variables.WPL * 0.5) * variables.bmf
    #     else:
    #         variables.rangecritical = False
    #         variables.Woew = a * variables.WTO + b
    #         variables.Wbat = variables.WTO * variables.bmf
    #
    # if variables.rangecritical == False:
    #     batcandirange = (variables.WTO - variables.WPL * 0.5) * bat_range
    #     batcandiendurance = variables.WTO * bat_endurance
    # else:
    #     batcandirange = variables.WTO * bat_range
    #     batcandiendurance = (variables.WTO + variables.WPL * 0.5) * bat_endurance
    # if batcandiendurance > batcandirange:
    #     variables.rangecritical = False
    # else:
    #     variables.rangecritical = True
    # variables.Wbat = max(batcandirange, batcandiendurance)
    # return variables


if __name__ == "__main__":
    v = CurrentVariables()
    print(classIestimation_alt(v))
    
    
    
    



