import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

class WingPlanform:
    def __init__(self):
        self.b          = 16.6                   #m
        self.S          = 20                     #m2
        self.qc_sweep   = 0                      #rad
        self.taper      = 1.4
        self.airfoil    = 'naca4415'
        self.cr         = 1.98                   #m
        self.ct         = 0.79

class StiffenedWing:
    def __init__(self, n, stringer, n_string_u, n_string_l):
        self.b               = 16.6                   #m
        self.S               = 20                     #m2
        self.qc_sweep        = 0                      #rad
        self.taper           = 1.4
        self.airfoil         = 'naca4415'
        self.cr              = 1.98                   #m
        self.ct              = 0.79
        self.stringer        = stringer
        self.n_string_u      = n_string_u
        self.n_string_l      = n_string_l
        self.cross_sections  = StiffenedWing.create_cross_sections(self,n)[0]
        self.ylst            = StiffenedWing.create_cross_sections(self,n)[1]
        self.clst            = StiffenedWing.create_cross_sections(self,n)[2]
    
    
    def create_cross_sections(self, n):
        wing = WingPlanform()
        theta = np.arctan((wing.cr-wing.ct)/2/wing.b)
        
        file_name = wing.airfoil+'.txt'
        
        file = open(file_name, 'r')
        lines = file.readlines()
        file.close()
    
        xlst, zlst = [], []
    
        for line in lines:
            xlst.append(line[:7])
            zlst.append(line[8:])
        
        for i in range(len(xlst)):
            xlst[i] = float(xlst[i])
            zlst[i] = float(zlst[i])
        
        clst = np.linspace(wing.cr, wing.ct, n)
        
        cross_sections = []
        
        for c in clst:
            xs = []
            zs = []
            for l in range(len(xlst)):
                xs.append(xlst[l]*c+0.25*(wing.cr-c))
                
                zs.append(zlst[l]*c)
    
            cross_sections.append([xs,zs])
        
        ylst = []
        
        theta = np.arctan((wing.cr-wing.ct)/2/wing.b)
        
        for k in range(len(clst)):
            y = 0.25*(wing.cr-clst[k])/np.tan(theta)
            ylst.append(y)
        
        return cross_sections, ylst, clst
    
    
        
    def position_stringers(self, n, stringer, n_string_u, n_string_l):
        
        wing = StiffenedWing(n, stringer, n_string_u, n_string_l)
        
        cross_sections = wing.cross_sections
        
        l_u = []
        l_l = []
        for cross_section in cross_sections:
            l_u_cs = []
            l_l_cs = []
            for i in range(len(cross_sections[0][0])-1):
                if cross_section[1][i] >=0 and cross_section[1][i+1] >=0:
                    l_u_cs.append(np.sqrt((cross_section[0][i+1]-cross_section[0][i])**2+(cross_section[0][i+1]-cross_section[0][i])**2))
                elif cross_section[1][i] <=0 and cross_section[1][i+1] <=0:
                    l_l_cs.append(np.sqrt((cross_section[0][i+1]-cross_section[0][i])**2+(cross_section[0][i+1]-cross_section[0][i])**2))
            
            l_u.append(sum(l_u_cs))
            l_l.append(sum(l_l_cs))
            
        
        
        return l_u, l_l
    
    
wing = StiffenedWing(5,None,None,None)

l_u, l_l = wing.position_stringers(5,None,None,None)
                    
                
        
        
    

#def calc_MoI(cross_sections):
#
#
#def create_half_wing(cross_sections, ylst):



#cross_sections, ylst, clst = create_cross_sections(5)
#
#fig = plt.figure()
#ax = plt.axes(projection='3d')
#ax.contour3D(cross_sections[0][0], cross_sections[0][1], ylst[0])
#plt.axes().set_aspect('equal')
#plt.show()


#plt.plot(cross_sections[-1][0], cross_sections[-1][1])
#plt.axes().set_aspect('equal')
#plt.show()


#def create_cross_sections(n):
#    wing = WingPlanform()
#    theta = np.arctan((wing.cr-wing.ct)/2/wing.b)
#    
#    file_name = wing.airfoil+'.txt'
#    
#    file = open(file_name, 'r')
#    lines = file.readlines()
#    file.close()
#
#    xlst, zlst = [], []
#
#    for line in lines:
#        xlst.append(line[:7])
#        zlst.append(line[8:])
#    
#    for i in range(len(xlst)):
#        xlst[i] = float(xlst[i])
#        zlst[i] = float(zlst[i])
#    
#    clst = np.linspace(wing.cr, wing.ct, n)
#    
#    cross_sections = []
#    
#    for c in clst:
#        xs = []
#        zs = []
#        for l in range(len(xlst)):
#            xs.append(xlst[l]*c+0.25*(wing.cr-c))
#            
#            zs.append(zlst[l]*c)
#
#        cross_sections.append([xs,zs])
#    
#    ylst = []
#    
#    theta = np.arctan((wing.cr-wing.ct)/2/wing.b)
#    
#    for k in range(len(clst)):
#        y = 0.25*(wing.cr-clst[k])/np.tan(theta)
#        ylst.append(y)
#    
#    return cross_sections, ylst, clst

    
    
    
