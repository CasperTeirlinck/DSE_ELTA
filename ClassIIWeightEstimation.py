# -*- coding: utf-8 -*-
"""
Created on Fri May  8 09:00:06 2020

@author: HJ Hoogendoorn


Class II Weight Estimation Method based on the Cessna Method
Outputs only the empty weight, not yet the c.g. locations!
"""
import numpy as np
#----------------Input parameters----------------------------------

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







#conversions
kglbs = 2.2046 #from kg to lbs
mft = 1/0.3048 #meter to foot




#---------------------functions for component weight estimation--------------------------

def MainWingEstimation(W_to,S_w,n_ult,A_w,Strut):
    if Strut == False:
        W_w = (0.04674*(W_to*kglbs)**0.397*(S_w*mft**2)**0.360*n_ult**0.397*A_w**1.712)
    else:
        W_w = (0.002933*(S_w*mft**2)**1.1018*A_w**2.473*n_ult**0.611)
    return W_w/kglbs


def EmpennageEstimation(W_to,S_h,S_v,A_h,A_v,t_r_h,t_r_v,qsweep_v):
    W_h = (3.183*(W_to*kglbs)**0.887*(S_h*mft**2)**0.101*A_h**0.138)/(174.04*(t_r_h*mft)**0.223)
    
    W_v = (1.68*(W_to*kglbs)**0.567*(S_v*mft**2)**1.249*A_v**0.482)/(639.95*(t_r_v*mft)**0.747*np.cos(qsweep_v/180*np.pi)**0.882)
    return W_h/kglbs, W_v/kglbs

def FuselageEstimation(W_to,pax,l_fn,Strut,n_ult, d_fus,V_max):
#    if Strut == False:
#        W_f = 0.04682*(W_to*kglbs)**0.692*pax**0.374*(l_fn*mft)**0.590
#    else:
#        W_f = 14.86*(W_to*kglbs)**0.144*((l_fn*mft)/pax)**0.778*(l_fn*mft)**0.383*pax**0.455
    W_f = 200*((n_ult*W_to*kglbs/10**5)**0.286*(l_fn*mft/10)**0.857*((l_fn+d_fus)*mft/10)*(V_max/100)**0.338)**1.1
    return W_f/kglbs
    
def EngineAndNacelleEstimation(W_to, ductedfan):
    W_e = P_to/K_p
    
    if ductedfan == False:
        W_n = 0.24*W_to*0.5 #0.5 is scaling factor due to electric engine! electric engine = 55 kg, piston = 144, more than 50% difference, but want to be conservative for the nacelle weight
    else:
        W_n = 0.24*W_to*0.6 #0.6 is engineering judgement
        
    return W_e, W_n #already in kg

def LandingGearEstimation(W_to,W_l, n_ult_l, l_sn,l_sm,W_wsn,W_wsm):
    W_frontgear = 0.013*(W_to*kglbs)+0.146*(W_l*kglbs)**0.417*n_ult_l**0.950*(l_sm*mft)**0.183+W_wsm*kglbs
    W_middlegear = 6.2 + 0.0013*(W_to*kglbs) + 0.000143*(W_l*kglbs)**0.749*n_ult_l*(l_sn*mft)**0.788 + W_wsn*kglbs
    W_g = W_frontgear + W_middlegear
    if Retract == True:
        W_g += 0.014*(W_to*kglbs)   
        
    return W_frontgear/kglbs, W_reargear/kglbs, W_g/kglbs
    
def FlightControlSystemEstimation(W_to):
    W_fc = 0.0168*W_to
    return W_fc #already in kg

def ElectricalSystemEstimation(W_to):
    W_el = 0.0268*W_to
    return W_el #already in kg
    

def EmptyWeight(W_to,W_l,n_ult,n_ult_l,pax,l_fn,P_to,W_b,W_prop,K_p,ductedfan,S_w,A_w,Strut,S_h,A_h,t_r_h,S_v,A_v,t_r_v,qsweep_v,Retract,W_wsn,W_wsm,l_sn,l_sm,V_max,d_fus):
    
    OEW = (MainWingEstimation(W_to,S_w,n_ult,A_w,Strut) + EmpennageEstimation(W_to,S_h,S_v,A_h,A_v,t_r_h,t_r_v,qsweep_v)[0] + EmpennageEstimation(W_to,S_h,S_v,A_h,A_v,t_r_h,t_r_v,qsweep_v)[1] +
           FuselageEstimation(W_to,pax,l_fn,Strut,n_ult, d_fus,V_max) + EngineAndNacelleEstimation(W_to, ductedfan)[0] + EngineAndNacelleEstimation(W_to, ductedfan)[1] +
           LandingGearEstimation(W_to,W_l, n_ult_l, l_sn,l_sm,W_wsn,W_wsm) + FlightControlSystemEstimation(W_to) + ElectricalSystemEstimation(W_to) +
           W_prop + W_avion)
    #W_b is not included as it is part of the payload
    return OEW





























