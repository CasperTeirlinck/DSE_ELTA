import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

#################################################################
### SOME LINES NEED TO BE REVIEWED IF QC SWEEP IS IMPLEMENTED ###
############## NEEDS REVIEW FOR SYMMETRIC AIRFOIL ###############
##### ENSURE DATA FOR AIRFOIL IS GIVEN FROM TE -> TOP -> BOT ####
#################################################################


class WingPlanform:
    def __init__(self):
        self.b          = 16.6                   #m
        self.S          = 20                     #m2
        self.qc_sweep   = 0                      #rad
        self.taper      = 1.4
        self.airfoil    = 'naca4415'
        self.cr         = 1.98                   #m
        self.ct         = 0.79
        

class StiffenedWing(WingPlanform):
    def __init__(self, n, stringer_u, stringer_l, n_string_u, n_string_l, spar_le_loc, spar_te_loc):
        super().__init__()

        self.stringer_u      = stringer_u                               #stringer type on upper sheet (class)
        self.stringer_l      = stringer_l                               #stringer type on lower sheet (class)
        self.n_string_u      = n_string_u                               #number of stringers on upper sheet
        self.n_string_l      = n_string_l                               #number of stringers on lower sheet
        self.spar_le_loc     = spar_le_loc                              #x/c location of leading edge spar
        self.spar_te_loc     = spar_te_loc                              #x/c location of trailing edge spar
        self.cross_sections  = self.create_cross_sections(n)[0]         #list of cross sections x,z coordinates (output: [[[x1,x2],[z1,z2]],...])
        self.ylst            = self.create_cross_sections(n)[1]         #list of all y coordinates of the cross sections
        self.clst            = self.create_cross_sections(n)[2]         #list of all chord lengths of cross sections
    
    
    def create_cross_sections(self, n):
        theta = np.arctan((self.cr-self.ct)/2/self.b)
        
        file_name = self.airfoil+'.txt'
        
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
        
        clst = np.linspace(self.cr, self.ct, n)
        
        cross_sections = []
        
        for c in clst:
            xs = []
            zs = []
            for l in range(len(xlst)):
                xs.append(xlst[l]*c+0.25*(self.cr-c))      #DOES NOT WORK IF THERES SWEEP
                
                zs.append(zlst[l]*c)
    
            cross_sections.append([xs,zs])
        
        ylst = []

        
        for k in range(len(clst)):
            y = 0.25*(self.cr-clst[k])/np.tan(theta)
            ylst.append(y)
        
        return cross_sections, ylst, clst
    
    
        
    def position_stringers(self, n_string_u, n_string_l):
        
        l_ucs  = []
        l_lcs  = []
        l_u    = []
        l_l    = []
        l_lst  = []
        
        for cross_section in self.cross_sections:
            l_u_cs    = []
            l_l_cs    = []
            l_lst_cs  = []
            
            for i in range(len(self.cross_sections[0][0])-1):
                x0 = cross_section[0][i]
                z0 = cross_section[1][i]
                x1 = cross_section[0][i+1]
                z1 = cross_section[1][i+1]
                l = np.sqrt((z1-z0)**2+(x1-x0)**2)
                l_lst_cs.append(l)
                
                if z0 >=0 and z1 >=0:
                    l_u_cs.append(l)
                    
                elif z0 <=0 and z0 <=0 or z0>=0 and z1<=0:
                    l_l_cs.append(l)
            
            l_ucs.append(l_u_cs)
            l_lcs.append(l_l_cs)
            l_lst.append(l_lst_cs)
            l_u.append(sum(l_u_cs))
            l_l.append(sum(l_l_cs))
        
        
        x_box_u = []
        z_box_u = []
        x_box_l = []
        z_box_l = []
        
        
        x_min_u_in = []
        x_max_u_in = []
        x_min_l_in = []
        x_max_l_in = []
        
        string_pos_u = []
        string_pos_l = []
        
        for j in range(len(self.cross_sections)):
            c     = self.clst[j]
            cross_section = self.cross_sections[j]
            x_min = self.spar_le_loc*c+0.25*(self.cr-c)             #DOES NOT WORK IF THERES QC SWEEP
            x_max = self.spar_te_loc*c+0.25*(self.cr-c)            #DOES NOT WORK IF THERES QC SWEEP
            x_box_ucs = []
            x_box_lcs = []
            z_box_ucs = []
            z_box_lcs = []
            
            for k in range(len(cross_section[0])):
                x_co = cross_section[0][k]
                z_co = cross_section[1][k]
                if x_co >= x_min and x_co <= x_max:
                    if z_co >= 0:
                        x_box_ucs.append(x_co)
                        z_box_ucs.append(z_co)
                    if z_co < 0:
                        x_box_lcs.append(x_co)
                        z_box_lcs.append(z_co)
            
            x_box_ucs.reverse()
            z_box_ucs.reverse()
            
            
            x_box_u.append(x_box_ucs)
            z_box_u.append(z_box_ucs)
            x_box_l.append(x_box_lcs)
            z_box_l.append(z_box_lcs)
            

            
            it = 0
            for x_coordinate in cross_section[0]:
                if x_coordinate == min(x_box_u[j]) and cross_section[1][it]>0:
                    x_min_u_in.append(it)
                if x_coordinate == max(x_box_l[j]) and cross_section[1][it]>0:
                    x_max_u_in.append(it)
                if x_coordinate == min(x_box_u[j]) and cross_section[1][it]<0:
                    x_min_l_in.append(it)
                if x_coordinate == max(x_box_l[j]) and cross_section[1][it]<0:
                    x_max_l_in.append(it)
                    
                it += 1
            
            l1 = sum(l_lst[j][:x_max_u_in[j]])
            l2 = sum(l_lst[j][:x_min_u_in[j]])
            string_pos_u_cs = list(np.linspace(l1,l2,n_string_u))
            
            l3 = sum(l_lst[j][:x_min_l_in[j]])
            l4 = sum(l_lst[j][:x_max_l_in[j]])
            string_pos_l_cs = list(np.linspace(l3,l4,n_string_l))
        
            string_pos_u.append(string_pos_u_cs)
            string_pos_l.append(string_pos_l_cs)
            
        #TODO: DEPENDING ON S-WISE LOCATION ALONG AIRFOIL GET X AND Z COORDINATE
        #TODO: MAKE CUMULATIVE LENGTH LIST, SEE WHERE STRING POINT LIES ON THAT LIST
        #TODO: INTERPOLATE BETWEEN TWO COORDINATES AND OBTAIN X AND Z
    
    
        
        return string_pos_u, string_pos_l
    
    
    
wing = StiffenedWing(5,None,None,None,None, 0.25, 0.75)

string_pos_u, string_pos_l = wing.position_stringers(5,5)

#plt.plot(x_box_u[0], z_box_u[0], x_box_l[0], z_box_l[0])
#plt.axes('equal')
#plt.show
                    
                
        
        
    

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

    
    
    
