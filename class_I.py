# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np

################################################################################
################################# CLASS I TOOL #################################
################################################################################


a_classI, b_classI = 0.548, 55.899
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

def calc_W_TO(W_PL, WbatWTO):
    # Returns: Take-off, OE, and bat weights in N, W_PL input in kg and WbatWTO dimensionless
    WTO = (W_PL+b_classI)/(1-a_classI-WbatWTO)
    return g0*WTO, g0*(a_classI*WTO+b_classI), g0*WbatWTO*WTO


if __name__ == "__main__":
    # TODO: Insert unit tests here
    print("Hello world!")
    
    
    
    



