# -*- coding: utf-8 -*-
"""
Created on Mon May 11 13:55:39 2020

@author: HJ Hoogendoorn and Yann Pejon

Aileron sizing and flap sizing
"""

import numpy as np



#----------------Input parameters aileron----------------------------------
b           = 12.2                   #[m] wing span
S           = 12.4                  #[m2] wing surface
c_root      = 1.5                    #[m] root chord
c_tip       = 0.6
taperrat    = 0.4                 #[-] taper ratio
b2          = b/2 - 0.5             #[m] fixed position of outer aileron location
delta_a     = 20*np.pi/180           #[rad] max aileron deflection both positive and negative
c_l_delta_a = 0.046825*180/(np.pi)   #[1/rad] change in the airfoil’s lift coefficient with aileron deflection
P           = 15*np.pi/180           #[rad/s] roll rate requirement
C_d0        = 0.007                  #[-] zero lift drag coefficient
c_l_a       = 5.729578               #[1/rad] lift slope
V           = 45  *1.1*1.852/3.6     #[m/s] 1.1*V_stall

#------------------Additional Input parameters Flap-------------------------
C_L_max_req = 1.8 #[-] max C_L needed
AR          = 12 #[-] Aspect ratio
a_max       = 14   /180*np.pi #[rad] Angle of attack at C_l_max
C_l_max     = 1.55 #[-] max C_l
deltaC_l_max= 1.13 #[-] change in C_l due to flaps at max deflection (0.2 chord @ 40 deg flap)
#deltaC_l_max= 1.43 #[-] change in C_l due to flaps at max deflection (0.3 chord @ 60 deg flap)



    
  

def rolldampingcoef(c_l_a,C_d0,c_root,S,taperrat,b):
    C_l_p = -((c_l_a+C_d0)*c_root*b)/(24*S)*(1+3*taperrat)
    return C_l_p #[rad]

    
    

def aileronstartinner(b,S,c_root,taperrat,b2,delta_a,c_l_delta_a,V,P): #outputs inner position of aileron
    C_l_p = rolldampingcoef(c_l_a,C_d0,c_root,S,taperrat,b)
    
    
    
    if (P*b)/(2*V)*(C_l_p)/(c_l_delta_a*delta_a)*b**2+b2**2 < 0:
        return print('wing span should be decreased')
    else:
        b1 = np.sqrt((P*b)/(2*V)*(C_l_p)/(c_l_delta_a*delta_a)*b**2+b2**2)
        return b1 #this value should be added or subtracted from b/2
    
    
    
    
    

def ClCLConverter(c_l_a,AR,a_max):
    C_L = c_l_a*AR/(AR+2)*a_max
    return C_L


def Flapsurface(c_l_a,AR,a_max,C_L_max_req,S,deltaC_l_max): #determines the wing surface that includes a flap, SO NOT ONLY THE FLAP SURFACE, BUT ENTIRE WING
    C_L_max = ClCLConverter(c_l_a,AR,a_max)
    C_L_max = C_L_max #- C_L_max/10
    deltaC_L_max = C_L_max_req- C_L_max
    Swf = (S*deltaC_L_max)/(deltaC_l_max*0.9)
    return Swf
    
def Flappos(b,S,c_root,taperrat,b2,delta_a,c_l_delta_a,V,P,c_tip,c_l_a,AR,a_max,C_L_max_req,deltaC_l_max):
    f2 = aileronstartinner(b,S,c_root,taperrat,b2,delta_a,c_l_delta_a,V,P) #- 1.5  #outer location of flap

    A = (c_root-c_tip)/b
    B = -c_root
    C = -(Flapsurface(c_l_a,AR,a_max,C_L_max_req,S,deltaC_l_max)/2 - c_root*f2 + A*f2**2)
    
    #f11 = (-B+np.sqrt(B**2-4*A*C))/(2*A)
    f12 = (-B-np.sqrt(B**2-4*A*C))/(2*A)
        
    return f12


        
    
    
    
    
    
