# -*- coding: utf-8 -*-
"""
Created on Fri May  8 09:00:06 2020

@author: HJ Hoogendoorn


Class II Weight Estimation Method based on the Cessna Method
Outputs only the empty weight, not yet the c.g. locations!
"""
import numpy as np
from variables import *
#----------------Input parameters----------------------------------
"""
#General aircraft parameters

W_to    = 750   #[kg]   take-off weight
W_l     = 750   #[kg]   landing weight
n_ult   = 5.7   #[-]    ultimate load factor
n_ult_l = 5.7   #[-]    ultimate load factor during landing, may be taken as 5.7 according to Cessna method from Roskam
pax     = 2     #[-]    number of pax
l_fn    = 5     #[m]    length fuselage excluding engine
d_fus   = 1.5   #[m]    depth of fuselage
V_max   = 120   #[kts]  max speed

#Powerplant
P_to    = 80    #[kW]   take-off power
W_b     = 30    #[kg]   battery weight

#Thrust
W_prop  = 5     #[kg]   propeller weight
K_p     = 2.5   #[kW/kg]Engine specific parameter
ductedfan = False

#Main Wing
S_w      = 14   #[m**2] surface
A_w      = 8    #[-]    Aspect ratio
Strut    = True #       Strut braced wing, if False: it is a cantilever wing

#Empennage
#Horizontal stabilizer
S_h      = 2    #[m**2] Surface
A_h      = 4    #[-]    Aspect ratio
t_r_h   = 0.1   #[m]    max root thickness

#Vertical tail
S_v      = 1.04 #[m**2] Surface
A_v      = 1.41 #[-]    Aspect ratio
t_r_v   = 0.15  #[m]    max root thickness
qsweep_v= 30    #[deg]  quarter chord sweep

#Landing Gear
Retract = True  #       retractable landing gear True or False
W_wsn   = 5     #[kg]   Nose wheel weight + its strut assembly
W_wsm   = 7     #[kg]   Main wheel weight + its strut assembly
l_sn    = 0.2   #[m]    Shock strut length nose wheel
l_sm    = 0.2   #[m]    Shock strut length main wheel

#Training Enhancement
W_avion = 15    #[kg] avionic weight

"""





#conversions
kglbs = 2.2046 #from kg to lbs
mft = 1/0.3048 #meter to foot




#---------------------functions for component weight estimation--------------------------

def MainWingEstimation(variables): #W_to,S_w,n_ult,A_w,Strut 
    if not variables.strutted_wing:
        W_w = (0.04674*(variables.WTO*kglbs)**0.397*(variables.S*mft**2)**0.360*variables.n_ult**0.397*variables.A**1.712)
    else:
        W_w = (0.002933*(variables.S*mft**2)**1.1018*variables.A**2.473*variables.n_ult**0.611)
    return W_w/kglbs

def EmpennageEstimation(variables):
    tr_h = variables.tcr_h*np.sqrt(variables.Sh/variables.A_h)
    tr_v = variables.tcr_v*np.sqrt(variables.Sv/variables.A_v)

    W_h = (3.183*(variables.WTO*kglbs)**0.887*(variables.Sh*mft**2)**0.101*variables.A_h**0.138)/(174.04*(tr_h*mft)**0.223)
    W_v = (1.68*(variables.WTO*kglbs)**0.567*(variables.Sv*mft**2)**1.249*variables.A_v**0.482)/(639.95*(tr_v*mft)**0.747*np.cos(variables.sweep_v/180*np.pi)**0.882)
    return W_h/kglbs, W_v/kglbs

def FuselageEstimation(variables):
    W_f = 200*((variables.n_ult*variables.WTO*kglbs/10**5)**0.286*(variables.l_fus*mft/10)**0.857*((variables.l_fus+variables.d_fus)*mft/10)*(variables.Vmax/100)**0.338)**1.1
    return W_f/kglbs
    
def EngineAndNacelleEstimation(W_to, ductedfan):
    #W_e = variables.P_total/K_p
    
    #if ductedfan == False:
        #W_n = 0.24*W_to*0.5 #0.5 is scaling factor due to electric engine! electric engine = 55 kg, piston = 144, more than 50% difference, but want to be conservative for the nacelle weight
    #else:
        #W_n = 0.24*W_to*0.6 #0.6 is engineering judgement
        
    #return W_e, W_n #already in kg
    return 100.

def LandingGearEstimation(variables):
    W_frontgear = 0.013*(variables.WTO*kglbs)+0.146*(variables.WTO*kglbs)**0.417*variables.n_ult**0.950*(variables.l_   sm*mft)**0.183 + variables.W_wsm*kglbs
    W_maingear = 6.2 + 0.0013*(variables.WTO*kglbs) + 0.000143*(variables.WTO*kglbs)**0.749*variables.n_ult_l*(variables.l_sn*mft)**0.788 + variables.W_wsn*kglbs
    W_g = W_frontgear + W_maingear
    if Retract == True:
        W_g += 0.014*(variables.WTO*kglbs)   
        
    return W_frontgear/kglbs, W_maingear/kglbs, W_g/kglbs
    
def FlightControlSystemEstimation(variables):
    W_fc = 0.0168*variables.WTO
    return W_fc #already in kg

def ElectricalSystemEstimation(variables):
    W_el = 0.0268*variables.WTO
    return W_el #already in kg
    

def EmptyWeight(variables):
    
    OEW = (MainWingEstimation(variables) + EmpennageEstimation(variables)[0] + EmpennageEstimation(variables)[1] +
           FuselageEstimation(variables) + EngineAndNacelleEstimation(variables)[0] + EngineAndNacelleEstimation(variables)[1] +
           LandingGearEstimation(variables)[3] + FlightControlSystemEstimation(variables) + ElectricalSystemEstimation(variables) +
           W_prop + W_avion)
    #W_b is not included as it is part of the payload
    return OEW

def CalculateClassII(variables):
    variables.Wwing = 9.81*MainWingEstimation(variables)
    variables.W_htail, variables.W_vtail = 9.81*EmpennageEstimation(variables)
    variables.Wfus = 9.81*FuselageEstimation(variables)
    variables.Weng = 9.81*Engine(variables)
    variables.Wgear_front,variables.Wgear_main,variables.Wgear = 9.81*LandingGearEstimation(variables) 
    variables.Wfcs = 9.81*FlightControlSystemEstimation(variables)
    variables.Wels = 9.81*ElectricalSystemEstimation(variables)
    variables.Woew_classII = 9.81*EmptyWeight(variables)

    return variables


## TEST
if __name__ == "__main__":
    variables = CurrentVariables()
    variables.Sh = 2
    variables.Sv = 2
    print(MainWingEstimation(variables))
    print(EmpennageEstimation(variables))
    print(FuselageEstimation(variables))
    print(LandingGearEstimation(variables))
    print(FlightControlSystemEstimation(variables))
    print(ElectricalSystemEstimation(variables))
    print(CalculateClassII(variables))

























