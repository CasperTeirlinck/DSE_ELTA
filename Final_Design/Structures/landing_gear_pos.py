#import variables as v
#import cg_determination as cg
import numpy as np

'''
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
'''

def size_gear(variables):
    sizing = True
    i = 0
    max_i = 1000
    move = 0.1
    nose_margin = 0.5 # [m]

    while sizing:
        h_cg   = (variables.fuselageheight + variables.propclear)-variables.zcg
        z_tail = variables.h_htail- h_cg
        m1     = -1/np.tan(5.25*np.pi/180)
        m2     = np.tan(15*1.05*np.pi/180)
        c      = z_tail-m2*(variables.xtail-variables.xcg_max)
        x_int  = c/(m1-m2)
        x_main = x_int+move
        z_main = m1*(x_main)
        factor = 0.12
        x_nose = -((1-factor)/factor)*x_main

        if (abs(x_nose)+nose_margin)>variables.xcg_max and i<max_i:
            move += 0.1
            i += 1
        else:
            sizing=False
            if i>=max_i:
                print('There is a problem with the landing gear sizing, ask Bob')
            else:
                pass

    variables.xmg = x_main
    #variables.zmg = z_main
    variables.xng = variables.xcg_max + x_nose

    return variables

#print(size_gear(variables))
    
    







