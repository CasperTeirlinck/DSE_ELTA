#import variables as v
#import cg_determination as cg
import numpy as np


class Test():
    def __init__(self):
        self.fus_height = 1.5
        self.tail_height = 1.15                         #m from ground
        self.l_tail = 3                                 #m from nose
        self.xcg_aft = 1                                #m from nose
        self.fuselagecg_z = self.fus_height*0.6
        self.x_lemac = -1
        self.z_lemac = 0.2
        self.ct = 0.6
        self.b = 12.2
        self.MAC = 1.1
        self.dihedral = 1
        self.WTO = 750*9.80665
        self.propclear = 0.23
        self.tailtiplength = 3.325
        self.lowwing = True
        
        
#variables = v.CurrentVariables()

variables = Test()

def size_gear(variables):
    h_cg   = (variables.fus_height + variables.propclear)-variables.fuselagecg_z
    z_tail = variables.tail_height- h_cg
    m1     = -1/np.tan(5.25*np.pi/180)
    m2     = np.tan(15*1.05*np.pi/180)
    c      = z_tail-m2*(variables.l_tail-variables.xcg_aft)
    x_int  = c/(m1-m2)
    move   = 0.1
    x_main = x_int+move
    z_main = m1*(x_main)
    factor = 0.12
    x_nose = -((1-factor)/factor)*x_main
    
    return x_main, z_main, x_nose

print(size_gear(variables))
    
    







