# -*- coding: utf-8 -*-
"""
Created on Mon May 11 13:55:39 2020

@author: HJ Hoogendoorn and Yann Pejon

Aileron sizing
"""

import numpy as np



#----------------Input parameters----------------------------------
b           = 8.25                  #[m] wing span
S           = 10                    #[m2] wing surface
c_root      = 1.2                   #[m] root chord
taperrat    = 1.6                   #[-] taper ratio
b2          = b/2 - 0.15            #[m] fixed position of outer aileron location
delta_a     = 20*np.pi/180          #[rad] max aileron deflection both positive and negative
c_l_delta_a = 0.0145*180/(np.pi)    #[1/deg] change in the airfoil’s lift coefficient with aileron deflection
P           = 15*np.pi/180          #[deg/s] roll rate requirement
C_d0        = 0.01                  #[-] zero lift drag coefficient
c_l_a       = 2*np.pi               #[1/rad] lift slope
V           = 45  *1.1*1.852/3.6      #[m/s] 1.1*V_stall






    
  

def rolldampingcoef(c_l_a,C_d0,c_root,S,taperrat,b):
    C_l_p = -((c_l_a+C_d0)*c_root*b)/(24*S)*(1+3*taperrat)
    return C_l_p #[rad]

    
#def aileronauthorityderivative(c_l_a,C_d0,c_root,delta_a,V,b):
#    C_l_p = rolldampingcoef(c_l_a,C_d0,c_root)
#    C_l_delta_a = -(P*b)/(2*V*delta_a)*C_l_p
#
#    return C_l_delta_a #[1/rad]
    
    

def aileronstartinner(b,S,c_root,taperrat,b2,delta_a,c_l_delta_a,V,P):
    C_l_p = rolldampingcoef(c_l_a,C_d0,c_root,S,taperrat,b)
    
    
    
    if (P*b)/(2*V)*(C_l_p)/(c_l_delta_a*delta_a)*b**2+b2**2 < 0:
        return print('wing span should be decreased')
    else:
        b1 = np.sqrt((P*b)/(2*V)*(C_l_p)/(c_l_delta_a*delta_a)*b**2+b2**2)
        return b1
    
    
    
    
