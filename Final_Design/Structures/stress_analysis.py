import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from math import pi, sqrt
import copy
from typing import Union
from Final_Design.aerodynamics.utils import readAeroLoads

#################################################################
### SOME LINES NEED TO BE REVIEWED IF QC SWEEP IS IMPLEMENTED ###
############## NEEDS REVIEW FOR SYMMETRIC AIRFOIL ###############
##### ENSURE DATA FOR AIRFOIL IS GIVEN FROM TE -> TOP -> BOT ####
######### SPAR Iyy IGNORED DUE TO THIN WALL ASSUMPTION ##########
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



class MatProps():
    def __init__(self, sigma_y, E, poisson, rho, sigma_comp=None, name=None, alpha=0.8, n=0.6):
        # Density (rho) in g/cc, NOT kg/m3!
        self.sigma_y = sigma_y
        if sigma_comp == None:
            self.sigma_comp = sigma_y
        else:
            self.sigma_comp = sigma_comp
        self.E = E
        self.poisson = poisson
        self.rho = rho*1000
        self.alpha = alpha
        self.n = n
        if name != None:
            self.name = name



class SparProp:
    def __init__(self, t=0.01, material=MatProps(1,1,1,1)):
        self.t        = t
        self.material = material
        self.h        = 0
        self.Ixx      = self.calc_Ixx(t)
        self.Iyy      = self.calc_Iyy(t)
    
    def calc_Ixx(self, t):
        return (t*self.h**3)/12
    
    def calc_Iyy(self, t):
        return (self.h*t**3)/12
    
    def calc_A(self, t):
        return t*self.h
    
class Spar:
    def __init__(self, xpos, zpos, h, spar_prop=SparProp):
        self.properties = copy.deepcopy(spar_prop)
        self.xpos = xpos
        self.zpos = zpos
        self.h    = h



class Stringer:
    def __init__(self, xpos, zpos, stringer_instance):
        self.properties = copy.deepcopy(stringer_instance)
        self.xpos = xpos
        self.zpos = zpos
        
    





class StiffenedWing(WingPlanform):
    def __init__(self, n, z_stiff, j_stiff, n_string_u, n_string_l,
                 spar_le_loc, spar_te_loc, spar):
        super().__init__()
        self.stringers_u_lst = [list(zeros) for zeros in np.zeros((n, n_string_u))]
        # self.stringers_u_lst = [[0] * n_string_u] * n                   #stringer type on upper sheet (class)
        self.stringers_l_lst = [[0] * n_string_l] * n                   #stringer type on lower sheet (class)
        self.spar_lst        = [[0] * 2] * n
        
        self.n_string_u      = n_string_u                               #number of stringers on upper sheet
        self.n_string_l      = n_string_l                               #number of stringers on lower sheet
        self.spar_le_loc     = spar_le_loc                              #x/c location of leading edge spar
        self.spar_te_loc     = spar_te_loc                              #x/c location of trailing edge spar
        self.cross_sections  = self.create_cross_sections(n)[0]         #list of cross sections x,z coordinates (output: [[[x1,x2],[z1,z2]],...])
        self.ylst            = self.create_cross_sections(n)[1]         #list of all y coordinates of the cross sections
        self.clst            = self.create_cross_sections(n)[2]         #list of all chord lengths of cross sections
        
        self.string_z_u      = self.position_stringers(n_string_u, n_string_l)[0]
        self.string_x_u      = self.position_stringers(n_string_u, n_string_l)[1]
        self.string_z_l      = self.position_stringers(n_string_u, n_string_l)[2]
        self.string_x_l      = self.position_stringers(n_string_u, n_string_l)[3]
        self.x_box_u         = self.position_stringers(n_string_u, n_string_l)[4]
        self.z_box_u         = self.position_stringers(n_string_u, n_string_l)[5]
        self.x_box_l         = self.position_stringers(n_string_u, n_string_l)[6]
        self.z_box_l         = self.position_stringers(n_string_u, n_string_l)[7]
        
        self.spar_x_loc_le   = self.calc_spar_centroid_loc()[0]
        self.spar_z_loc_le   = self.calc_spar_centroid_loc()[1]
        self.spar_x_loc_te   = self.calc_spar_centroid_loc()[2]
        self.spar_z_loc_te   = self.calc_spar_centroid_loc()[3]
        self.spar_h_le       = self.calc_spar_centroid_loc()[4]
        self.spar_h_te       = self.calc_spar_centroid_loc()[5]
        
        
        

        for i in range(len(self.cross_sections)):
            stringers_u_cs   = self.stringers_u_lst[i]
            stringers_l_cs   = self.stringers_l_lst[i]
            
            stringers_u_x_cs = self.string_x_u[i]
            stringers_u_z_cs = self.string_z_u[i]
            stringers_l_x_cs = self.string_x_l[i]
            stringers_l_z_cs = self.string_z_l[i]
            
            spar_x_loc_le_cs = self.spar_x_loc_le[i]
            spar_z_loc_le_cs = self.spar_z_loc_le[i]
            spar_x_loc_te_cs = self.spar_x_loc_te[i]
            spar_z_loc_te_cs = self.spar_z_loc_te[i]
            spar_h_le_cs     = self.spar_h_le[i]
            spar_h_te_cs     = self.spar_h_te[i]
            
            
            spars_cs         = self.spar_lst[i]
            spars_cs[0]      = Spar(spar_x_loc_le_cs, spar_z_loc_le_cs, spar_h_le_cs)
            spars_cs[1]      = Spar(spar_x_loc_te_cs, spar_z_loc_te_cs, spar_h_te_cs)

            for j in range(len(stringers_u_cs)):
                stringers_u_cs[j] = Stringer(stringers_u_x_cs[j], stringers_u_z_cs[j], j_stiff)
            for k in range(len(stringers_l_cs)):
                stringers_l_cs[k] = Stringer(stringers_l_x_cs[k], stringers_l_z_cs[k], z_stiff)
        
        


    
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
    
    
        
    def position_stringers(self, n_string_u, n_string_l):   ############################################
                                                            ### CREATE WB CS AND POSITIONS STRINGERS ###
                                                            ### EQUIDISTANCE ALONG UP AND LOWER SEC- ###
                                                            ### TIONS, RETURNS LIST OF
        l_ucs  = []
        l_lcs  = []
        l_u    = []
        l_l    = []
        l_lst  = []
        l_cum  = []
        
        string_z_u = []
        string_x_u = []
        
        string_z_l = []
        string_x_l = []
            
        
        
        for cross_section in self.cross_sections:
            l_u_cs    = []
            l_l_cs    = []
            l_lst_cs  = []
            l_cum_cs  = [0]
            
            
            for i in range(len(self.cross_sections[0][0])-1):
                x0 = cross_section[0][i]
                z0 = cross_section[1][i]
                x1 = cross_section[0][i+1]
                z1 = cross_section[1][i+1]
                l = np.sqrt((z1-z0)**2+(x1-x0)**2)
                l_lst_cs.append(l)
                l_cum_cs.append(sum(l_lst_cs[:i+1]))
                
                if z0 >=0 and z1 >=0:
                    l_u_cs.append(l)
                    
                elif z0 <=0 and z0 <=0 or z0>=0 and z1<=0:
                    l_l_cs.append(l)
            
            
            l_ucs.append(l_u_cs)
            l_lcs.append(l_l_cs)
            l_lst.append(l_lst_cs)
            
            
            l_cum.append(l_cum_cs)
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
            
        
            l_cum_lst = l_cum[j]
            
            string_z_u_cs = []
            string_x_u_cs = []
            
            string_z_l_cs = []
            string_x_l_cs = []
            
            
            for string_posu in string_pos_u_cs:
                for m in range(len(l_cum_lst)-1):
                    if l_cum_lst[m] <= string_posu  <= l_cum_lst[m+1]:
                        sint  = (cross_section[1][m+1]-cross_section[1][m])/(l_cum_lst[m+1]-l_cum_lst[m])
                        cost = (cross_section[0][m+1]-cross_section[0][m])/(l_cum_lst[m+1]-l_cum_lst[m])
                        
                        string_zi = (string_posu-l_cum_lst[m])*sint+cross_section[1][m]
                        string_xi = (string_posu-l_cum_lst[m])*cost+cross_section[0][m]
                        
                        string_z_u_cs.append(string_zi)
                        string_x_u_cs.append(string_xi)
                        
                        break
            for string_posl in string_pos_l_cs:
                for m in range(len(l_cum_lst)-1):
                    if l_cum_lst[m] <= string_posl  <= l_cum_lst[m+1]:
                        sint  = (cross_section[1][m+1]-cross_section[1][m])/(l_cum_lst[m+1]-l_cum_lst[m])
                        cost = (cross_section[0][m+1]-cross_section[0][m])/(l_cum_lst[m+1]-l_cum_lst[m])
                        
                        string_zi = (string_posl-l_cum_lst[m])*sint+cross_section[1][m]
                        string_xi = (string_posl-l_cum_lst[m])*cost+cross_section[0][m]
                        
                        string_z_l_cs.append(string_zi)
                        string_x_l_cs.append(string_xi)
                        
                        break
            
        
            string_z_u.append(string_z_u_cs)
            string_x_u.append(string_x_u_cs)
            string_z_l.append(string_z_l_cs)
            string_x_l.append(string_x_l_cs)
            
        return string_z_u, string_x_u, string_z_l, string_x_l, x_box_u, z_box_u, x_box_l, z_box_l
    
    
    def calc_spar_centroid_loc(self):
        
        spar_x_loc_le = []
        spar_z_loc_le = []
        
        spar_x_loc_te = []
        spar_z_loc_te = []
        
        h_spar_le     = []
        h_spar_te     = []
        
        for i in range(len(self.cross_sections)):
            x_box_ui   = self.x_box_u[i]
            z_box_ui   = self.z_box_u[i]
            
            x_box_li   = self.x_box_l[i]
            z_box_li   = self.z_box_l[i]
            
            spar_x_loci = min(min(x_box_ui), min(x_box_li))
            
            if spar_x_loci == min(x_box_ui):
                ind_le = x_box_ui.index(spar_x_loci)
                spar_x_loc_le.append(spar_x_loci)
                
            else:
                ind_le = x_box_li.index(spar_x_loci)
                spar_x_loc_le.append(spar_x_loci)
                
            spar_x_loci_te = max(max(x_box_ui), max(x_box_li))
            
            if spar_x_loci_te == max(x_box_ui):
                ind_te = x_box_ui.index(spar_x_loci_te)
                spar_x_loc_te.append(spar_x_loci_te)
                
            else:
                ind_te = x_box_li.index(spar_x_loci_te)
                spar_x_loc_te.append(spar_x_loci_te)
            
            spar_z_loc_lei = z_box_ui[ind_le]+z_box_li[x_box_li.index(min(x_box_li))]
            h_spar_lei     = abs(z_box_ui[ind_le])+abs(z_box_li[x_box_li.index(min(x_box_li))])
            
            spar_z_loc_tei = z_box_ui[ind_te]+z_box_li[x_box_li.index(max(x_box_li))]
            h_spar_tei     = abs(z_box_ui[ind_te])+abs(z_box_li[x_box_li.index(max(x_box_li))])
            
            
            
            spar_z_loc_le.append(spar_z_loc_lei)
            h_spar_le.append(h_spar_lei)
            
            spar_z_loc_te.append(spar_z_loc_tei)
            h_spar_te.append(h_spar_tei)
            
            
            
        return spar_x_loc_le, spar_z_loc_le, spar_x_loc_te, spar_z_loc_te, h_spar_le, h_spar_te
    
            
            
            
            
        
        
        
        
        
        

class J_Stringer:
    name = "J_Stringer"
    def __init__(self, Le, material, stiff_ratio=0.5, ts=0.001, t1=0.001, t2=0.001, t3=0.001, t4=0.001,
                 b1=0.005, b2=0.006, b3=0.005, b4=0.005):
        # Geometry defined following https://puu.sh/FSdDm/df0747d2c5.png, gray area included in 2_alt
        # Reference coordinate system is x positive to the left, y positive upwards, origin bottom right.
        self.Le = Le

        self.material = material

        self.ts = ts

        self.t1 = t1
        self.t2 = t2
        self.t3 = t3
        self.t4 = t4
        self.b1 = b1
        self.b2 = b2
        self.b3 = b3
        self.b4 = b4
        self.area1 = self.t1*self.b1
        self.area2 = self.t2*self.b2
        self.area2_alt = self.area2 + self.t2*(self.t1 + self.t4)
        self.area3 = self.t3*self.b3
        self.area4 = self.t4*self.b4
        self.total_area = self.area1 + self.area2_alt + self.area3 + self.area4

        self.xbar, self.ybar = self.calc_centroid()
        self.Ixx, self.Iyy = self.momentofinertia()

        self.sigma_cc = self.crippling_stress(self.material)

        self.sigma_cr = self.calc_critical_stress(self.Le)


    def calc_centroid(self):
        # Returns xbar, ybar
        Qx = self.area1 * (self.b2+self.t1*0.5+self.t3) + self.area2_alt * (self.t1+self.b2+self.t4)*0.5 + \
             self.area4 * self.t4*0.5 + self.area3 * self.t3*0.5
        Qy = self.area1 * (self.b4 + self.t2 + self.b1*0.5) + self.area2_alt * (self.b4 + self.t2*0.5) + \
             self.area4 * self.b4*0.5 + self.area3 * (self.b4 + self.t2 + self.b3*0.5)
        return Qy/self.total_area, Qx/self.total_area

    def crippling_partial(self, t, b, C, material: MatProps):
        return material.alpha * (C/material.sigma_comp * pi*pi*material.E*t*t/(b*b*12*(1-material.poisson)))**(1-material.n) * material.sigma_comp

    def crippling_stress(self, material: MatProps):
        return (self.crippling_partial(self.t1, self.b1, 0.425, material) * self.area1 +
                self.crippling_partial(self.t2, self.b2, 0.400, material) * self.area2 +
                self.crippling_partial(self.t4, self.b4, 0.425, material) * self.area4 +
                self.crippling_partial(self.t3, self.b3, 0.425, material) * self.area3) / \
               (self.total_area-self.area2_alt+self.area2)

    def momentofinertia(self):
        Ixx = (self.b1*self.t1**3+self.t2*(self.b2+self.t1+self.t4)**3+self.b4*self.t4**3+self.b3*self.t3**3)/12 + \
              self.area1*(self.b2+self.t4+self.t1*0.5-self.ybar)**2 + \
              self.area2_alt * (self.area2_alt/(2*self.t2)-self.ybar)**2 + \
              self.area4*(self.t4*0.5-self.ybar)**2 + self.area3*(self.t3*0.5-self.ybar)**2
        Iyy = (self.t1*self.b1**3+(self.b2+self.t1+self.t4)*self.t2**3+self.t4+self.b4**3+self.t3*self.b3**3) + \
              self.area1*(self.b4+self.t2+self.b1*0.5-self.xbar)**2 + self.area2*(self.b4+self.t2*0.5-self.xbar)**2 + \
              self.area4*(self.b4*0.5-self.xbar)**2 + self.area3*(self.b4+self.t2+self.b3*0.5-self.xbar)**2
        return Ixx, Iyy

    def calc_critical_stress(self, Le):
        self.slendercrit = sqrt(2*pi*pi*self.material.E/self.sigma_cc)
        rhoxx = sqrt(self.Ixx/self.total_area)
        rhoyy = sqrt(self.Iyy/self.total_area)
        self.slender = Le/min(rhoxx, rhoyy)
        if self.slender < self.slendercrit:
            sigma_cr = self.sigma_cc * (1 - self.sigma_cc * self.slender * self.slender / (4*pi*pi*self.material.E))
        else:
            sigma_cr = pi*pi*self.material.E / (self.slender*self.slender)

        if sigma_cr < self.material.sigma_comp:
            return sigma_cr
        else:
            return self.material.sigma_comp


class Z_Stringer(J_Stringer):
    name = "Z_Stringer"
    def __init__(self, Le, material, stiff_ratio=0.5, ts=0.001, t1=0.001, t2=0.001, t4=0.001,
                 b1=0.005, b2=0.006, b4=0.005, t3=None, b3=None):
        # Geometry defined following https://puu.sh/FSdDm/df0747d2c5.png without the t3 and b3 elements
        super().__init__(Le=Le, material=material, stiff_ratio=stiff_ratio,
                         ts=ts, t1=t1, t2=t2, t4=t4, b1=b1, b2=b2, b4=b4)

        self.total_area = self.area1 + self.area2_alt + self.area4

        self.sigma_cc=self.crippling_stress(self.material)*self.total_area
        self.sigma_cr = self.calc_critical_stress(self.Le)

    def calc_centroid(self):
        # Returns xbar, ybar
        Qx = self.area1 * (self.b2+self.t1*0.5+self.t3) + self.area2_alt * (self.t1+self.b2+self.t4)*0.5 + \
             self.area4 * self.t4*0.5
        Qy = self.area1 * (self.b4 + self.t2 + self.b1*0.5) + self.area2_alt * (self.b4 + self.t2*0.5) + \
             self.area4 * self.b4*0.5
        return Qy/self.total_area, Qx/self.total_area

    def crippling_stress(self, material: MatProps):
        return (self.crippling_partial(self.t1, self.b1, 0.425, material) * self.area1 +
                self.crippling_partial(self.t2, self.b2, 0.400, material) * self.area2 +
                self.crippling_partial(self.t4, self.b4, 0.425, material) * self.area4) / \
               (self.total_area-self.area2_alt+self.area2)

    def momentofinertia(self):
        Ixx = (self.b1*self.t1**3+self.t2*(self.b2+self.t1+self.t4)**3+self.b4*self.t4**3)/12 + \
              self.area1*(self.b2+self.t4+self.t1*0.5-self.ybar)**2 + \
              self.area2_alt * (self.area2_alt/(2*self.t2)-self.ybar)**2 + \
              self.area4*(self.t4*0.5-self.ybar)**2
        Iyy = (self.t1*self.b1**3+(self.b2+self.t1+self.t4)*self.t2**3+self.t4+self.b4**3)/12 + \
              self.area1*(self.b4+self.t2+self.b1*0.5-self.xbar)**2 + self.area2*(self.b4+self.t2*0.5-self.xbar)**2 + \
              self.area4*(self.b4*0.5-self.xbar)**2
        return Ixx, Iyy
    






aluminium=MatProps(sigma_y=450000000, E=72400000000, poisson=0.33, rho=2.87, name="AA2024", alpha=0.8, n=0.6)
z_stiff = Z_Stringer(1, aluminium)
j_stiff = J_Stringer(1, aluminium)
wing = StiffenedWing(5, z_stiff, j_stiff, 10, 7, 0.25, 0.75, SparProp)



    
    
    
