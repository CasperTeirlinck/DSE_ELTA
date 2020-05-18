# -*- coding: utf-8 -*-
"""
Created on Fri May 15 16:46:57 2020

@author: HJ Hoogendoorn

CG POSITIONS OF FUSELAGE + ENGINE + PAX + BATTERY
"""
import numpy as np
from variables import *
##################INPUT PARAMETERS###############################
#All positions should be measured from the nose cone towards the back
LG_nose_pos     = 0.6 #[m] position of nose landing gear
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



def calculate_cg_groups(variables):
    Concept = variables.concept_number
    if Concept == 1:
        def PropellerCG(variables, engine_plac_concept_1=3.9):
            CGy = engine_plac_concept_1 - 0.985/2
            CGz = variables.propclear + 0.5261904761904761 + np.tan(11.5/180*np.pi)*(0.402+variables.tailtiplength)
            return CGy, CGz

        def EngineCG(variables, engine_plac_concept_1=3.9):
            CGy = engine_plac_concept_1
            CGz = variables.propclear + 0.5261904761904761 + np.tan(11.5/180*np.pi)*(0.402+variables.tailtiplength)
            return CGy, CGz

        def BatteryCG(variables):
            CGy = variables.bulkhead/2
            CGz = variables.propclear + 0.402
            return CGy, CGz

        def Baggage(variables):
            CGy = variables.bulkhead + variables.cabinlength + 0.155
            CGz = variables.propclear + 0.526 + np.tan(11.5/180*np.pi)*(0.402+variables.tailtiplength)
            return CGy, CGz

        def Payload(variables):
            CGy = variables.bulkhead + 0.907
            CGz = variables.propclear + 0.526 + np.tan(11.5/180*np.pi)*(0.402+variables.tailtiplength)
            return CGy, CGz

        def FuselageEmptyStructure(variables):
            CGy = variables.fuselage_len/2
            CGz = variables.propclear + 0.3095 + np.tan(10/180*np.pi)*(variables.fuselage_len)
            return CGy, CGz

    if Concept == 2 or Concept == 3:

        def PropellerCG(variables):
            CGy = variables.prop_spin/2 #[m] fixed value
            CGz = variables.propclear + variables.prop_d/2
            return CGy, CGz

        def EngineCG(variables):
            CGy = variables.prop_spin + (variables.bulkhead - variables.prop_spin)/2
            CGz = variables.propclear + variables.prop_d/2 - 0.09
            return CGy, CGz

        def BatteryCG(variables):
            CGy = variables.wingpos + (variables.cr/2)
            CGz = PropellerCG(variables)[1] - 0.76 + 0.2
            return CGy, CGz

        def Baggage(variables):
            CGy = variables.bulkhead + variables.cabinlength + 0.151
            CGz = PropellerCG(variables)[1] - 0.76 + 0.302
            return CGy, CGz

        def Payload(variables):
            CGy = variables.bulkhead + 0.907
            CGz = PropellerCG(variables)[1] - 0.76 + 0.7558
            return CGy, CGz

        def FuselageEmptyStructure(variables):
            CGy = variables.prop_spin + (variables.fuselage_len - variables.prop_spin)/2
            CGz = variables.propclear + variables.prop_d/2
            return CGy, CGz

    if Concept == 4 or Concept == 5:

        def PropellerCG(variables):
            CGy = variables.wingpos + variables.eng_perc_rootchord*variables.cr/100 - variables.enginelength
            CGz = variables.propclear + variables.fus_height + variables.eng_height_above_w
            return CGy, CGz

        def EngineCG(variables):
            CGy = variables.wingpos + variables.eng_perc_rootchord*variables.cr/100
            CGz = variables.propclear + variables.fus_height + variables.eng_height_above_w
            return CGy, CGz

        def BatteryCG(variables):
            CGy = 1.8365384615384617
            CGz = variables.propclear + 0.2
            return CGy, CGz

        def Baggage(variables):
            CGy = 0.875 + variables.cabinlength + 0.385
            CGz = variables.propclear + 0.481
            return CGy, CGz

        def Payload(variables):
            CGy = 0.875 + 0.907
            CGz = variables.propclear + 0.526
            return CGy, CGz

        def FuselageEmptyStructure(variables):
            CGy = variables.fuselage_len/2
            CGz = variables.propclear + 0.526
            return CGy, CGz

    variables.propcg_y, variables.propcg_z = PropellerCG(variables)
    variables.enginecg_y, variables.enginecg_z = EngineCG(variables)
    variables.batterycg_y, variables.batterycg_z = BatteryCG(variables)
    variables.baggagecg_y, variables.baggagecg_z = Baggage(variables)
    variables.payloadcg_y, variables.payloadcg_z = Payload(variables)
    variables.fuselagecg_y, variables.fuselagecg_z = FuselageEmptyStructure(variables)
    return variables