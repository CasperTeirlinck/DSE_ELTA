# -*- coding: utf-8 -*-
"""
Created on Fri May 15 16:46:57 2020

@author: HJ Hoogendoorn

CG POSITIONS OF FUSELAGE + ENGINE + PAX + BATTERY
"""
import numpy as np
##################INPUT PARAMETERS###############################
#All positions should be measured from the nose cone towards the back
LG_nose_pos     = 1.058 #[m] position of nose landing gear
LG_main_pos     = 2.395 #[m] position of main landing gear
prop_clearance  = 0.230 #[m] ground clearance of the propeller (concept 2 and 3) or fuselage (Concept 1, 4 and 5)
prop_diameter   = 1.090 #[m] diameter of propeller
W_main_pos      = 1.270 #[m] position of leading edge main wing at root
W_main_chord    = 1.6   #[m] root chord length of main wing
l_wt            = 3.325 #[m] leading edge main wing until trailing edge vertical tail measured at the root

#for concept 1:
engine_plac_concept_1 = 3.9 #[m] cg of engine measured from the nose

#for concept 4 and 5:
eng_perc_rootchord = 20 #[%] engine cg in percentage of root chord
eng_height_above_w = 0.3 #[m] engine cg in m above the main wing
l_eng              = 1.0 #[m] engine length

#Fixed dimensions
fus_height  = 1.725 #[m] height from lowes point floor until heighest point of the fuselage
cab_l       = 1.300 #[m] length of fire wall until rear of pilot's seat
prop_spin   = 0.3650 #[m] length of the propeller spinner



#Make a choice of concept:
Concept = 4



if Concept == 1:


    def PropellerCG(engine_plac_concept_1, prop_clearance, l_wt):
        CGy = engine_plac_concept_1 - 0.985/2
        CGz = prop_clearance + 0.5261904761904761 + np.tan(11.5/180*np.pi)*(0.402+l_wt)
        return CGy, CGz
    
    def EngineCG(engine_plac_concept_1, prop_clearance, l_wt):
        CGy = engine_plac_concept_1
        CGz = prop_clearance + 0.5261904761904761 + np.tan(11.5/180*np.pi)*(0.402+l_wt)
        return CGy, CGz
    
    def BatteryCG(LG_nose_pos, prop_clearance):
        CGy = LG_nose_pos/2
        CGz = prop_clearance + 0.402
        return CGy, CGz
    
    def Baggage(LG_nose_pos, prop_clearance, cab_l,l_wt):
        CGy = LG_nose_pos + cab_l + 0.155
        CGz = prop_clearance + 0.526 + np.tan(11.5/180*np.pi)*(0.402+l_wt)
        return CGy, CGz
    
    def Payload(prop_clearance, LG_nose_pos,l_wt):
        CGy = LG_nose_pos + 0.907
        CGz = prop_clearance + 0.526 + np.tan(11.5/180*np.pi)*(0.402+l_wt)
        return CGy, CGz
    
    def FuselageEmptyStructure(W_main_pos, l_wt, prop_clearance):
        CGy = (W_main_pos + l_wt )/2
        CGz = prop_clearance + 0.3095 + np.tan(10/180*np.pi)*(l_wt + W_main_pos)
        return CGy, CGz





if Concept == 2 or Concept == 3:

    def PropellerCG(prop_clearance,prop_diameter, prop_spin):
        CGy = prop_spin/2 #[m] fixed value
        CGz = prop_clearance + prop_diameter/2
        return CGy, CGz
    
    def EngineCG(prop_spin, LG_nose_pos):
        CGy = prop_spin + (LG_nose_pos - prop_spin)/2
        CGz = prop_clearance + prop_diameter/2 - 0.09
        return CGy, CGz
    
    def BatteryCG(W_main_pos, W_main_chord, prop_clearance,prop_diameter, prop_spin):
        CGy = W_main_pos + (W_main_chord/2)
        CGz = PropellerCG(prop_clearance,prop_diameter, prop_spin)[1] - 0.76 + 0.2 #+ 0.121 + prop_clearance
        return CGy, CGz
    
    def Baggage(LG_nose_pos, prop_clearance, cab_l,prop_diameter, prop_spin):
        CGy = LG_nose_pos + cab_l + 0.151
        CGz = PropellerCG(prop_clearance,prop_diameter, prop_spin)[1] - 0.76 + 0.302
        return CGy, CGz
    
    def Payload(LG_nose_pos, prop_clearance,prop_diameter, prop_spin):
        CGy = LG_nose_pos + 0.907
        CGz = PropellerCG(prop_clearance,prop_diameter, prop_spin)[1] - 0.76 + 0.7558
        return CGy, CGz
    
    def FuselageEmptyStructure(prop_spin, W_main_pos, l_wt, prop_clearance, prop_diameter):
        CGy = prop_spin + ((W_main_pos - prop_spin) + l_wt )/2
        CGz = prop_clearance + prop_diameter/2
        return CGy, CGz
    
    
    
    
    
if Concept == 4 or Concept == 5:


    def PropellerCG(l_eng,W_main_pos,eng_perc_rootchord,W_main_chord,fus_height,eng_height_above_w,  prop_clearance):
        CGy = W_main_pos + eng_perc_rootchord*W_main_chord/100 - l_eng
        CGz = prop_clearance + fus_height + eng_height_above_w
        return CGy, CGz
    
    def EngineCG(W_main_pos,eng_perc_rootchord,W_main_chord,fus_height,eng_height_above_w,  prop_clearance):
        CGy = W_main_pos + eng_perc_rootchord*W_main_chord/100
        CGz = prop_clearance + fus_height + eng_height_above_w
        return CGy, CGz
    
    def BatteryCG(prop_clearance):
        CGy = 1.8365384615384617
        CGz = prop_clearance + 0.2
        return CGy, CGz
    
    def Baggage(cab_l):
        CGy = 0.875 + cab_l + 0.385
        CGz = prop_clearance + 0.481
        return CGy, CGz
    
    def Payload(prop_clearance):
        CGy = 0.875 + 0.907
        CGz = prop_clearance + 0.526
        return CGy, CGz
    
    def FuselageEmptyStructure(W_main_pos , l_wt, prop_clearance):
        CGy = (W_main_pos + l_wt)/2
        CGz = prop_clearance + 0.526
        return CGy, CGz  

#def o(x):
#    return (1250/6.5*x)/1000
#
##def BatteryCG():
##    CGy = 
##    CGz = 
##    return CGy, CGz