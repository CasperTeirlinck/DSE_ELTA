# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np

################################################################################
################################# CLASS I TOOL #################################
################################################################################

def classI(R, A, e, n, E, eff, W_PL): 
    #R: Range [m], A: Aspect ratio, e: Oswald efficiency, n: number of eng (1 or 2)
    # E: Battery specific energy [Wh/kg], eff: Efficiency of propulsion system, W_PL: Payload [kg]
    #Returns: Take-off, OE, PL and bat weights
    
    Cfe_dict = {'1': 0.0055,
                '2': 0.0045 }
    C_fe = Cfe_dict[str(n)]
    g = 9.80665
    SwetS = 2
    a, b = 0.548, 55.899
    CD_0 = C_fe*SwetS
    CD_optr = 2*CD_0
    CL_optr = np.sqrt(np.pi*A*e*CD_0)
    LD_optr = CL_optr/CD_optr
    E = E*3600
    WbatWTO = R/(LD_optr*1/g*E*eff)
    W_TO = (W_PL+b)/(1-a-WbatWTO)
    W_OE = a*W_TO+b
    W_bat = WbatWTO*W_TO
    
    return W_TO, W_OE, W_PL, W_bat
    
    
    
    



