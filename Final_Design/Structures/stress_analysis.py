import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from math import pi, sqrt
import copy
from Final_Design.aerodynamics.utils import readAeroLoads
from scipy import interpolate



#################################################################
### SOME LINES NEED TO BE REVIEWED IF QC SWEEP IS IMPLEMENTED ###
#################################################################
############## NEEDS REVIEW FOR SYMMETRIC AIRFOIL ###############
#################################################################
##### ENSURE DATA FOR AIRFOIL IS GIVEN FROM TE -> TOP -> BOT ####
#################################################################
######### SPAR Iyy IGNORED DUE TO THIN WALL ASSUMPTION ##########
#################################################################
### A/C WEIGHT AND AIR DENSITY MANUALLY ENTERED IN V_MAXG AND ###
######## V_MING FUNCTIONS (ALSO MAX AND MIN LOAD FACTORS) #######
#################################################################
### BATTERY VOLUME HARD CODED IN CALC_BATT_ENDPOINT FUNCTION ####
#################################################################
######### AERODYNAMIC STUFF HARD CODED IN VMIN AND VMAX #########
#################################################################





class MatProps():
    def __init__(self, sigma_y, E, poisson, rho, G, sigma_comp=None, name=None, alpha=0.8, n=0.6):
        # Density (rho) in g/cc, NOT kg/m3!
        self.sigma_y = sigma_y
        if sigma_comp == None:
            self.sigma_comp = sigma_y
        else:
            self.sigma_comp = sigma_comp
        self.E = E
        self.poisson = poisson
        self.rho = rho*1000
        self.G   = G
        self.alpha = alpha
        self.n = n
        if name != None:
            self.name = name



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
    



class WingPlanform:
    def __init__(self):
        self.b          = 6.2643*2               #m
        self.S          = 20                     #m2
        self.qc_sweep   = 0                      #rad
        self.taper      = 1.4
        self.airfoil    = 'naca4415'
        self.cr         = 1.98                   #m
        self.ct         = 0.79
        self.y_MAC      = (2/3)*self.cr*((1+self.taper+self.taper**2)/(1+self.taper))



aluminum=MatProps(sigma_y=450000000, E=72400000000, poisson=0.33, rho=2.87, G=28E9, name="AA2024", alpha=0.8, n=0.6)


class SparProp:
    def __init__(self, t=0.01, material=aluminum):
        self.t        = t
        self.material = material
    
    
class Spar:
    def __init__(self, xpos, zpos, h, t=0.01, material=aluminum):
        self.t        = t
        self.material = material
        self.xpos     = xpos
        self.zpos     = zpos
        self.h        = h
        self.Ixx      = self.calc_Ixx(self.t)
        self.Izz      = self.calc_Izz(self.t)
        self.A        = self.calc_A(self.t)
    
    
    def calc_Ixx(self, t):
        return (t*self.h**3)/12
    
    def calc_Izz(self, t):
        return (self.h*t**3)/12
    
    def calc_A(self, t):
        return t*self.h


class Stringer:
    def __init__(self, xpos, zpos, stringer_instance):
        self.properties = copy.deepcopy(stringer_instance)
        self.xpos = xpos
        self.zpos = zpos
        
class Skin:
    def __init__(self, t=0.003, material=aluminum):
        self.t = t
        self.material = material
        
    





class StiffenedWing(WingPlanform):
    def __init__(self, n, z_stiff, j_stiff, n_string_u, n_string_l,
                 spar_le_loc, spar_te_loc, spar, skin, batt = False):
        super().__init__()
        self.stringers_u_lst = [list(zeros) for zeros in np.zeros((n, n_string_u))]
        self.stringers_l_lst = [list(zeros) for zeros in np.zeros((n, n_string_l))]
        self.spar_lst        = [list(zeros) for zeros in np.zeros((n,2))]
        self.skin            = skin
        
        self.batt            = batt
        
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
        
        self.x_bar           = self.calc_centroid_wb()[0]
        self.z_bar           = self.calc_centroid_wb()[1]
        
        self.Ixx_lst         = self.calc_MoI_cs()[0]
        self.Izz_lst         = self.calc_MoI_cs()[1]
        self.Ixz_lst         = self.calc_MoI_cs()[2]
        
        self.Am_lst          = self.calc_wb_A()
        self.V_wb            = self.calc_half_wing_volume()
        
        self.x_sc_cg_lst     = self.calc_sc_loc()[0]
        self.x_sc_lst        = self.calc_sc_loc()[1]
        
        self.v_maxg          = self.calc_v_maxg()[0]
        self.ylst_maxL       = self.calc_v_maxg()[1]
        self.L_dist_maxg     = self.calc_v_maxg()[2]
        self.Lmax_tot        = self.calc_v_maxg()[3]
        self.cp_maxL_lst     = self.calc_v_maxg()[4]
        self.D_dist_maxg     = self.calc_v_maxg()[5]
        self.Dmax_tot        = self.calc_v_maxg()[6]
        
        self.v_ming          = self.calc_v_ming()[0]
        self.ylst_minL       = self.calc_v_ming()[1]
        self.L_dist_ming     = self.calc_v_ming()[2]
        self.Lmin_tot        = self.calc_v_ming()[3]
        self.cp_minL_lst     = self.calc_v_ming()[4]
        self.D_dist_ming     = self.calc_v_ming()[5]
        self.Dmin_tot        = self.calc_v_ming()[6]
        
        #INTERPOLATE LIFT TO GET LIFT AND CP AT ALL CROSS SECTIONS
        
        self.Lmax_cs_lst    = list(interpolate.interp1d(self.ylst_maxL, self.L_dist_maxg, kind='slinear', fill_value='extrapolate')(self.ylst))
        self.Dmax_cs_lst    = list(interpolate.interp1d(self.ylst_maxL, self.D_dist_maxg, kind='slinear', fill_value='extrapolate')(self.ylst))
        self.cp_maxL_int    = list(interpolate.interp1d(self.ylst_maxL, self.cp_maxL_lst, kind='slinear', fill_value='extrapolate')(self.ylst))
        self.cp_maxL        = interpolate.interp1d(self.ylst_maxL, self.cp_maxL_lst, kind='slinear', fill_value='extrapolate')([self.y_MAC])[0]
        
        self.xsc_at_yMAC    = interpolate.interp1d(self.ylst, self.x_sc_lst, kind = 'slinear', fill_value='extrapolate')([self.y_MAC])[0]
        
        self.Lmin_cs_lst    = list(interpolate.interp1d(self.ylst_minL, self.L_dist_ming, kind='slinear', fill_value='extrapolate')(self.ylst))
        self.Dmin_cs_lst    = list(interpolate.interp1d(self.ylst_minL, self.L_dist_ming, kind='slinear', fill_value='extrapolate')(self.ylst))
        self.cp_minL_int    = list(interpolate.interp1d(self.ylst_minL, self.cp_minL_lst, kind='slinear', fill_value='extrapolate')(self.ylst))
        self.cp_minL        = interpolate.interp1d(self.ylst_minL, self.cp_minL_lst, kind='slinear', fill_value='extrapolate')([self.y_MAC])[0]
        
        self.xcp_maxg_lst   = self.calc_xloc_cp()[0]
        self.xcp_ming_lst   = self.calc_xloc_cp()[1]
        
        if batt == False:
        
            self.T_lst_maxg     = self.calc_torque_all_cs()[0]
            self.T_lst_ming     = self.calc_torque_all_cs()[1]
            self.T_root_maxg    = self.calc_torque_all_cs()[2]
            self.T_root_ming    = self.calc_torque_all_cs()[3]
            
            self.Vz_maxg_lst    = self.calc_Vz_dist()[0]
            self.Vz_ming_lst    = self.calc_Vz_dist()[1]
            self.Vz_maxg_root   = self.calc_Vz_dist()[2]
            self.Vz_ming_root   = self.calc_Vz_dist()[3]
            
            self.Vx_maxg_lst    = self.calc_Vx_dist()[0]
            self.Vx_ming_lst    = self.calc_Vx_dist()[1]
            self.Vx_maxg_root   = self.calc_Vx_dist()[2]
            self.Vx_ming_root   = self.calc_Vx_dist()[3]
            
            self.Mx_maxg_lst    = self.calc_Mx_dist()[0]
            self.Mx_ming_lst    = self.calc_Mx_dist()[1]
            self.Mx_maxg_root   = self.calc_Mx_dist()[2]
            self.Mx_ming_root   = self.calc_Mx_dist()[3]
            
            self.Mz_maxg_lst    = self.calc_Mz_dist()[0]
            self.Mz_ming_lst    = self.calc_Mz_dist()[1]
            self.Mz_maxg_root   = self.calc_Mz_dist()[2]
            self.Mz_ming_root   = self.calc_Mz_dist()[3]
        
        if self.batt == True:
            self.y_batt         = self.calc_batt_endpoint_y()[0]
            self.V_batt         = self.calc_batt_endpoint_y()[1]
            
            self.batt_load_lst  = self.create_batt_loads_cs()
            
            self.T_lst_maxg     = self.calc_torque_all_cs()[0]
            self.T_lst_ming     = self.calc_torque_all_cs()[1]
            self.T_root_maxg    = self.calc_torque_all_cs()[2]
            self.T_root_ming    = self.calc_torque_all_cs()[3]
            
            self.Vz_maxg_lst    = self.calc_Vz_dist()[0]
            self.Vz_ming_lst    = self.calc_Vz_dist()[1]
            self.Vz_maxg_root   = self.calc_Vz_dist()[2]
            self.Vz_ming_root   = self.calc_Vz_dist()[3]
            
            self.Vx_maxg_lst    = self.calc_Vx_dist()[0]
            self.Vx_ming_lst    = self.calc_Vx_dist()[1]
            self.Vx_maxg_root   = self.calc_Vx_dist()[2]
            self.Vx_ming_root   = self.calc_Vx_dist()[3]
            
            self.Mx_maxg_lst    = self.calc_Mx_dist()[0]
            self.Mx_ming_lst    = self.calc_Mx_dist()[1]
            self.Mx_maxg_root   = self.calc_Mx_dist()[2]
            self.Mx_ming_root   = self.calc_Mx_dist()[3]
            
            self.Mz_maxg_lst    = self.calc_Mz_dist()[0]
            self.Mz_ming_lst    = self.calc_Mz_dist()[1]
            self.Mz_maxg_root   = self.calc_Mz_dist()[2]
            self.Mz_ming_root   = self.calc_Mz_dist()[3]
            
        
        
    
        
        
        
        
        
        
        


    
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
    
    def calc_centroid_wb(self):                             #CALCULATES CENTROID OF WB FOR ALL CROSS SECTIONS
        cross_sections = self.cross_sections
        
        x_bar = []
        z_bar = []
        
        for i in range(len(cross_sections)):
            
            xA = []
            zA = []
            A  = []
            
            
            stringers_x_u_cs = self.string_x_u[i]
            stringers_z_u_cs = self.string_z_u[i]
            stringers_x_l_cs = self.string_x_l[i]
            stringers_z_l_cs = self.string_z_l[i]
            
            
            xA.append(self.spar_lst[i][0].A*self.spar_x_loc_le[i])
            xA.append(self.spar_lst[i][1].A*self.spar_x_loc_te[i])
            zA.append(self.spar_lst[i][0].A*self.spar_z_loc_le[i])
            zA.append(self.spar_lst[i][1].A*self.spar_z_loc_te[i])
            A.append(self.spar_lst[i][0].A)
            A.append(self.spar_lst[i][1].A)
            
            
            for j in range(len(stringers_x_u_cs)):
                x_co   = stringers_x_u_cs[j]
                z_co   = stringers_z_u_cs[j]
                A_st_u = self.stringers_u_lst[i][j].properties.total_area
                xA.append(x_co*A_st_u)
                zA.append(z_co*A_st_u)
                A.append(A_st_u)
            
            for k in range(len(stringers_x_l_cs)):
                x_co   = stringers_x_l_cs[k]
                z_co   = stringers_z_l_cs[k]
                A_st_l = self.stringers_l_lst[i][k].properties.total_area
                xA.append(x_co*A_st_l)
                zA.append(z_co*A_st_l)
                A.append(A_st_l)
                
            
            
            x_bar_cs = sum(xA)/sum(A)
            z_bar_cs = sum(zA)/sum(A)
            
            x_bar.append(x_bar_cs)
            z_bar.append(z_bar_cs)
        
        return x_bar, z_bar

    def calc_MoI_cs(self):
        
        Ixx = []
        Izz = []
        Ixz = []
        
        for i in range(len(self.cross_sections)):
            x2A = []
            z2A = []
            xzA = []
            
            Ixx_spar_le = self.spar_lst[i][0].Ixx
            Izz_spar_le = self.spar_lst[i][0].Izz
            stein_le_xx = self.spar_lst[i][0].A * ((self.spar_lst[i][0].zpos-self.z_bar[i])**2)
            stein_le_zz = self.spar_lst[i][0].A * ((self.spar_lst[i][0].xpos-self.x_bar[i])**2)
            Ixx_cont_le = Ixx_spar_le+stein_le_xx
            Izz_cont_le = Izz_spar_le+stein_le_zz
            Ixz_cont_le = (self.spar_lst[i][0].zpos-self.z_bar[i])*(self.spar_lst[i][0].xpos-self.x_bar[i])*self.spar_lst[i][0].A
            
            Ixx_spar_te = self.spar_lst[i][1].Ixx
            Izz_spar_te = self.spar_lst[i][1].Izz
            stein_te_xx = self.spar_lst[i][1].A * ((self.spar_lst[i][1].zpos-self.z_bar[i])**2)
            stein_te_zz = self.spar_lst[i][1].A * ((self.spar_lst[i][1].xpos-self.x_bar[i])**2)
            Ixx_cont_te = Ixx_spar_te+stein_te_xx
            Izz_cont_te = Izz_spar_te+stein_te_zz
            Ixz_cont_te = (self.spar_lst[i][1].zpos-self.z_bar[i])*(self.spar_lst[i][1].xpos-self.x_bar[i])*self.spar_lst[i][1].A
            
            
            
                        
            for j in range(len(self.stringers_u_lst[i])):
                x = self.stringers_u_lst[i][j].xpos-self.x_bar[i]
                z = self.stringers_u_lst[i][j].zpos-self.z_bar[i]
                A = self.stringers_u_lst[i][j].properties.total_area
                x2A.append(x**2*A)
                z2A.append(z**2*A)
                xzA.append(x*z*A)
            
            for k in range(len(self.stringers_l_lst[i])):
                x = self.stringers_l_lst[i][k].xpos-self.x_bar[i]
                z = self.stringers_l_lst[i][k].zpos-self.z_bar[i]
                A = self.stringers_l_lst[i][k].properties.total_area
                x2A.append(x**2*A)
                z2A.append(z**2*A)
                xzA.append(x*z*A)
            
            Ixx.append(sum(x2A)+Ixx_cont_le+Ixx_cont_te)
            Izz.append(sum(z2A)+Izz_cont_le+Izz_cont_te)
            Ixz.append(sum(xzA)+Ixz_cont_le+Ixz_cont_te)
        
        return Ixx, Izz, Ixz
    
    def calc_wb_A(self):
        
        Am = []
        
        for i in range(len(self.cross_sections)):
            Am_cs = []
            x_box_ui = self.x_box_u[i]
            z_box_ui = self.z_box_u[i]
            x_box_li = self.x_box_l[i]
            z_box_li = self.x_box_l[i]
            
            for j in range(len(x_box_ui)-1):
                A_ub = abs(z_box_ui[j+1]*(x_box_ui[j+1]-x_box_ui[j]))
                A_lb = abs(z_box_ui[j]*(x_box_ui[j+1]-x_box_ui[j]))
                A_av = (A_ub+A_lb)/2
                Am_cs.append(A_av)
            
            for k in range(len(x_box_li)-1):
                A_ub = abs(z_box_li[k+1]*(x_box_li[k+1]-x_box_li[k]))
                A_lb = abs(z_box_li[k]*(x_box_li[k+1]-x_box_li[k]))
                A_av = (A_ub+A_lb)/2
                Am_cs.append(A_av)
            
            Am.append(sum(Am_cs))
        
        return Am
    
    def calc_sc_loc(self):
        
        d           = []
        l_bt        = []
        qb_cs       = []
        moment_all  = []
        qs_0        = []
        dqbl_lst    = []
        
        x_sc_cg_lst = []
        x_sc_lst    = []
        
        
        for i in range(len(self.cross_sections)):
            
            d_cs = []
            
            l_bt_cs   = []
            
            stringers_u_cs  = self.stringers_u_lst[i]
            stringers_l_cs  = self.stringers_l_lst[i]
            
            spars_cs        = self.spar_lst[i]
            
            x_cg = self.x_bar[i]
            z_cg = self.z_bar[i]
            
            
            for j in range(len(stringers_u_cs)-1):
                grad = (stringers_u_cs[j].zpos-stringers_u_cs[j+1].zpos)/(stringers_u_cs[j].xpos-stringers_u_cs[j+1].xpos)  #THIS IS CORRECT WAY BC OF DATA ORGANISATION
                c    = stringers_u_cs[j].zpos-grad*stringers_u_cs[j].xpos
                A    = -(stringers_u_cs[j].zpos-stringers_u_cs[j+1].zpos)
                B    = (stringers_u_cs[j].xpos-stringers_u_cs[j+1].xpos)
                C    = -(stringers_u_cs[j].xpos-stringers_u_cs[j+1].xpos)*c
                di   = abs(A*x_cg+B*z_cg+C)/np.sqrt(A**2+B**2)
                d_cs.append(di)
                
                l_bt_i = np.sqrt((stringers_u_cs[j].xpos-stringers_u_cs[j+1].xpos)**2+(stringers_u_cs[j].zpos-stringers_u_cs[j+1].zpos)**2)
                l_bt_cs.append(l_bt_i)
                
            
            d_cs.append(x_cg-spars_cs[0].xpos)
            l_bt_cs.append(spars_cs[0].h)
            
            for k in range(len(stringers_l_cs)-1):
                grad = (stringers_l_cs[k+1].zpos-stringers_l_cs[k].zpos)/(stringers_l_cs[k+1].xpos-stringers_l_cs[k].xpos)
                c    = stringers_l_cs[k].zpos-grad*stringers_l_cs[k].xpos
                A    = -(stringers_l_cs[k+1].zpos-stringers_l_cs[k].zpos)
                B    = (stringers_l_cs[k+1].xpos-stringers_l_cs[k].xpos)
                C    = -(stringers_l_cs[k+1].xpos-stringers_l_cs[k].xpos)
                di   = abs(A*x_cg+B*z_cg+C)/np.sqrt(A**2+B**2)
                d_cs.append(di)
                
                l_bt_i = np.sqrt((stringers_l_cs[k].xpos-stringers_l_cs[k+1].xpos)**2+(stringers_l_cs[k].zpos-stringers_l_cs[k+1].zpos)**2)
                l_bt_cs.append(l_bt_i)
            
            d_cs.append(spars_cs[1].xpos-x_cg)
            l_bt_cs.append(spars_cs[1].h)
            
            
            
            d.append(d_cs)
            l_bt.append(l_bt_cs)
            
            
            
            qb  = [0]
            Vz  = 1
            
            Ixx = self.Ixx_lst[i]
            Izz = self.Izz_lst[i]
            Ixz = self.Ixz_lst[i]
            
            for l in range(len(stringers_u_cs)-2):
                qbi = -(Vz*Izz)/(Ixx*Izz-Ixz**2)*stringers_u_cs[l+1].zpos*stringers_u_cs[l+1].properties.total_area + (Vz*Ixz)/(Ixx*Izz-Ixz**2)*stringers_u_cs[l+1].xpos*stringers_u_cs[l+1].properties.total_area + qb[l]
                qb.append(qbi)
            
            
            qb_sp_le = -(Vz*Izz)/(Ixx*Izz-Ixz**2)*stringers_u_cs[-1].zpos*stringers_u_cs[-1].properties.total_area + (Vz*Ixz)/(Ixx*Izz-Ixz**2)*stringers_u_cs[-1].xpos*stringers_u_cs[-1].properties.total_area + qb[-1]
            qb.append(qb_sp_le)
            
            
            for m in range(len(stringers_l_cs)-1):
                qbi = -(Vz*Izz)/(Ixx*Izz-Ixz**2)*stringers_l_cs[m+1].zpos*stringers_l_cs[m+1].properties.total_area + (Vz*Ixz)/(Ixx*Izz-Ixz**2)*stringers_l_cs[m+1].xpos*stringers_l_cs[m+1].properties.total_area + qb[l+m+2]
                qb.append(qbi)
            
            qb_sp_te = -(Vz*Izz)/(Ixx*Izz-Ixz**2)*stringers_l_cs[-1].zpos*stringers_l_cs[-1].properties.total_area + (Vz*Ixz)/(Ixx*Izz-Ixz**2)*stringers_l_cs[-1].xpos*stringers_l_cs[-1].properties.total_area + qb[-1]
            qb.append(qb_sp_te)
            
            qb_cs.append(qb)
            
            F = []
            moment_cs = []
            
            for n in range(len(qb)):
                q  = qb[n]
                l  = l_bt_cs[n]
                Fn = q*l
                F.append(Fn)
                dn = d_cs[n]
                moment = Fn*dn
                moment_cs.append(moment)
                
            moment_all.append(moment_cs)
                
            qbGt_cs = []
            lGt_cs  = []
            
            for o in range(len(qb[:self.n_string_u-1])):
                qbo = qb[o]
                G   = self.skin.material.G
                t   = self.skin.t
                l_s = l_bt_cs[o]
                
                qbGt_cs.append((qbo/(G*t))*l_s)
                lGt_cs.append((1/(G*t))*l_s)
            
            qb_le = qb[self.n_string_u-1]
            G_le  = spars_cs[0].material.G
            t_le  = spars_cs[0].t
            l_le  = l_bt_cs[self.n_string_u-1]
            
            qbGt_cs.append((qb_le/(G_le*t_le))*l_le)
            lGt_cs.append((1/(G_le*t_le))*l_le)
            
            for p in range(len(qb[self.n_string_u:-1])):
                qbp = qb[p]
                G   = self.skin.material.G
                t   = self.skin.t
                l_s = l_bt_cs[p]
                
                qbGt_cs.append((qbp/(G*t))*l_s)
                lGt_cs.append((1/(G*t))*l_s)
                
            qb_te = qb[-1]
            G_te  = spars_cs[1].material.G
            t_te  = spars_cs[1].t
            l_te  = l_bt_cs[-1]
            
            qbGt_cs.append((qb_te/(G_te*t_te))*l_te)
            lGt_cs.append((1/(G_te*t_te))*l_te)
            
            qs_0_cs = -(sum(qbGt_cs)/sum(lGt_cs))
            qs_0.append(qs_0_cs)
            
            
            
            Amqs_0 = 2*self.Am_lst[i]*qs_0_cs
            
            dqbl_lst_cs = []
            
            for r in range(len(qb)):
                dr = d_cs[r]
                qbr = qb[r]
                lr = l_bt_cs[r]
                dqbl = dr*qbr*lr
                dqbl_lst_cs.append(dqbl)
                
            dqbl_lst.append(dqbl_lst_cs)
            
            x_sc_cg = Amqs_0 + sum(dqbl_lst_cs)
            x_sc_cg_lst.append(x_sc_cg)
            

            x_sc = self.x_bar[i]+x_sc_cg
            
            
            x_sc_lst.append(x_sc)
        
            
        return x_sc_cg_lst, x_sc_lst
    
    def calc_y_lift_loc(self):
        
        ylst_aero, cl_lst, cd_lst, cp_lst = readAeroLoads(5)
        cl_lst[-1] = 0
        cd_lst[-1] = 0
        
        yA = []
        A  = []
        
        for i in range(len(ylst_aero)-1):
            y = (ylst_aero[i]+ylst_aero[i+1])/2
            A_ub = cl_lst[i+1]*(ylst_aero[i+1]-ylst_aero[i])
            A_lb = cl_lst[i]*(ylst_aero[i+1]-ylst_aero[i])
            A_av = (A_ub+A_lb)/2
            
            yA.append(y*A_av)
            A.append(A_av)
        
        return sum(yA)/sum(A)
    
    def calc_v_maxg(self):
        
        ylst_aero, cl_lst, cd_lst, cp_lst = readAeroLoads(5)
        cl_lst[-1] = 0
        cd_lst[-1] = 0
        A_lst = []
        
        for i in range(len(ylst_aero)-1):
            y = ylst_aero[i+1]-ylst_aero[i]
            A_ub = y*cl_lst[i+1]
            A_lb = y*cl_lst[i]
            A_av = (A_ub+A_lb)/2
            A_lst.append(A_av)
            
        v_maxg = np.sqrt((4.45*750*9.80665)/(1.225*self.S*sum(A_lst)))
        
        L_lst = [j * 0.5*1.225*v_maxg**2*self.S for j in cl_lst]
        D_lst = [l * 0.5*1.225*v_maxg**2*self.S for l in cd_lst]
        A_L = []
        A_D = []
        
        for k in range(len(ylst_aero)-1):
            y = ylst_aero[k+1]-ylst_aero[k]
            A_ub = y*L_lst[k+1]
            A_lb = y*L_lst[k]
            A_av = (A_ub+A_lb)/2
            
            A_ub_D = y*D_lst[k+1]
            A_lb_D = y*D_lst[k]
            A_av_D = (A_ub_D+A_lb_D)/2
            
            A_L.append(A_av)
            A_D.append(A_av_D)
            
        return v_maxg, ylst_aero, L_lst ,sum(A_L), cp_lst, D_lst, sum(A_D)
    
    def calc_v_ming(self):
        
        ylst_aero, cl_lst, cd_lst, cp_lst = readAeroLoads(-10)
        cl_lst[-1] = 0
        cd_lst[-1] = 0
        A_lst = []
        
        for i in range(len(ylst_aero)-1):
            y = ylst_aero[i+1]-ylst_aero[i]
            A_ub = y*cl_lst[i+1]
            A_lb = y*cl_lst[i]
            A_av = (A_ub+A_lb)/2
            A_lst.append(A_av)
            
        v_ming = np.sqrt((-2.45*750*9.80665)/(1.225*self.S*sum(A_lst)))
        
        L_lst = [j * 0.5*1.225*v_ming**2*self.S for j in cl_lst]
        D_lst = [l * 0.5*1.225*v_ming**2*self.S for l in cd_lst]
        A_L = []
        A_D = []
        
        for k in range(len(ylst_aero)-1):
            y = ylst_aero[k+1]-ylst_aero[k]
            A_ub = y*L_lst[k+1]
            A_lb = y*L_lst[k]
            A_av = (A_ub+A_lb)/2
            
            A_ub_D = y*D_lst[k+1]
            A_lb_D = y*D_lst[k]
            A_av_D = (A_ub_D+A_lb_D)/2
            
            A_L.append(A_av)
            A_D.append(A_av_D)
            
        return v_ming, ylst_aero, L_lst ,sum(A_L), cp_lst, D_lst, sum(A_D)
    
    def interpolate_splines(self, ylst_cs, ylst_aero, L_lst):
        
        grad_lst = []
        
        for i in range(len(ylst_aero)-1):
            grad = (L_lst[i+1]-L_lst[i])/(ylst_aero[i+1]-ylst_aero[i])
            grad_lst.append(grad)
        
        L_int_lst = []
        
        for j in range(len(ylst_cs)):
            y_cs = ylst_cs[j]
            
            for k in range(len(ylst_aero)-1):
                if ylst_aero[k] <= y_cs <= ylst_aero[k+1]:
                    grad = grad_lst[k]
                    
                    dL = grad*(y_cs-ylst_aero[k])
                    L_int = L_lst[k]+dL
                    L_int_lst.append(L_int)
                    
        
        return L_int_lst
    
    def calc_half_wing_volume(self):
        
        Vi = []
        
        for i in range(len(self.Am_lst)-1):
            dy = self.ylst[i+1]-self.ylst[i]
            V_ub = self.Am_lst[i]*dy
            V_lb = self.Am_lst[i+1]*dy
            V_av = (V_ub+V_lb)*dy
            Vi.append(V_av)
        
        return sum(Vi)
    
    def calc_batt_endpoint_y(self):
        
        V_batt = 131/2/1000  #L->m3
        Vi     = []
        y_batt = []
        i=0
        while sum(Vi) < V_batt:
            
            y  = self.ylst[i+1]
            dy = self.ylst[i+1]-self.ylst[i]
            V_ub = self.Am_lst[i]*dy
            V_lb = self.Am_lst[i+1]*dy
            V_av = (V_ub+V_lb)*dy
            Vi.append(V_av)
            y_batt.append(y)
            i += 1
            
        return y_batt[-1], sum(Vi)
            
            
        
    def calc_xloc_cp(self):
        
        x_cp_maxg = []
        x_cp_ming = []
        
        for i in range(len(self.ylst)):
            x_le = min(self.cross_sections[i][0])
            x_cp_maxg_i = x_le+self.cp_maxL_int[i]*self.clst[i]
            x_cp_ming_i = x_le+self.cp_minL_int[i]*self.clst[i]
            
            x_cp_maxg.append(x_cp_maxg_i)
            x_cp_ming.append(x_cp_ming_i)
            
        return x_cp_maxg, x_cp_ming
            
    
    def create_batt_loads_cs(self):
        
        m_batt = 200                    #kg
        W_batt = m_batt*9.80665
        w = W_batt/2/self.y_batt        #N/m
        
        batt_load_lst = []
        
        for i in range(len(self.cross_sections)):
            if self.ylst[i] < self.y_batt:
                batt_load_lst.append(w)
            else:
                batt_load_lst.append(0)
        
        return batt_load_lst
        
        
    
    def calc_torque_all_cs(self):
        
        T_lst_maxg = []
        T_lst_ming = []
        
        ylst_interm = []
        L_maxg_interm_lst = []
        L_ming_interm_lst = []
        W_batt_interm = []
        
        xcp_maxg_interm = []
        xcp_ming_interm = []
        
        
        if self.batt == False:
            
            for i in range(1, len(self.cross_sections)+1):
                ylst_interm.append(self.ylst[-i])
                
                Li_maxg = self.Lmax_cs_lst[-i]
                Li_ming = self.Lmin_cs_lst[-i]
                
                L_maxg_interm_lst.append(Li_maxg)
                L_ming_interm_lst.append(Li_ming)
                
                xcp_maxg_i = self.xcp_maxg_lst[-i]
                xcp_ming_i = self.xcp_ming_lst[-i]
                
                xcp_maxg_interm.append(xcp_maxg_i)
                xcp_ming_interm.append(xcp_ming_i)
                
                if len(T_lst_maxg) != 0:
                    L_maxg_toti = self.integrate(ylst_interm, L_maxg_interm_lst)
                    L_ming_toti = self.integrate(ylst_interm, L_ming_interm_lst)
                    
                    xcp_maxg = self.integrate(ylst_interm, xcp_maxg_interm)/(ylst_interm[0]-ylst_interm[i-1])
                    xcp_ming = self.integrate(ylst_interm, xcp_maxg_interm)/(ylst_interm[0]-ylst_interm[i-1])
                    
                    xsc_i   = self.x_sc_lst[-i]
                    
                    T_maxg_toti = L_maxg_toti*(xcp_maxg-xsc_i)
                    T_ming_toti = L_ming_toti*(xcp_ming-xsc_i)
                    
                    T_lst_maxg.append(T_maxg_toti)
                    T_lst_ming.append(T_ming_toti)
                
                if len(T_lst_maxg) == 0:
                    T_lst_maxg.append(0)
                    T_lst_ming.append(0)
                
                
            T_root_maxg = -T_lst_maxg[-1]
            T_root_ming = -T_lst_ming[-1]
            
            T_lst_maxg[-1] = T_lst_maxg[-1]+T_root_maxg
            T_lst_ming[-1] = T_lst_ming[-1]+T_root_ming
            
            
            T_lst_maxg.reverse()
            T_lst_ming.reverse()
        
        if self.batt == True:
            for i in range(1, len(self.cross_sections)+1):
                ylst_interm.append(self.ylst[-i])
                
                Li_maxg = self.Lmax_cs_lst[-i]
                Li_ming = self.Lmin_cs_lst[-i]
                W_batti = self.batt_load_lst[-i]
                
                L_maxg_interm_lst.append(Li_maxg)
                L_ming_interm_lst.append(Li_ming)
                W_batt_interm.append(W_batti)
                
                xcp_maxg_i = self.xcp_maxg_lst[-i]
                xcp_ming_i = self.xcp_ming_lst[-i]
                
                xcp_maxg_interm.append(xcp_maxg_i)
                xcp_ming_interm.append(xcp_ming_i)
                
                if len(T_lst_maxg) != 0:
                    L_maxg_toti = self.integrate(ylst_interm, L_maxg_interm_lst)
                    L_ming_toti = self.integrate(ylst_interm, L_ming_interm_lst)
                    W_batt_toti = self.integrate(ylst_interm, W_batt_interm)
                    
                    xcp_maxg = self.integrate(ylst_interm, xcp_maxg_interm)/(ylst_interm[0]-ylst_interm[i-1])
                    xcp_ming = self.integrate(ylst_interm, xcp_maxg_interm)/(ylst_interm[0]-ylst_interm[i-1])
                    
                    xsc_i   = self.x_sc_lst[-i]
                    
                    T_maxg_toti = L_maxg_toti*(xcp_maxg-xsc_i)+W_batt_toti*(self.x_bar[-1]-xsc_i)
                    T_ming_toti = L_ming_toti*(xcp_ming-xsc_i)+W_batt_toti*(self.x_bar[-1]-xsc_i)
                    
                    T_lst_maxg.append(T_maxg_toti)
                    T_lst_ming.append(T_ming_toti)
                
                if len(T_lst_maxg) == 0:
                    T_lst_maxg.append(0)
                    T_lst_ming.append(0)
                
                
            T_root_maxg = -T_lst_maxg[-1]
            T_root_ming = -T_lst_ming[-1]
            
            T_lst_maxg[-1] = T_lst_maxg[-1]+T_root_maxg
            T_lst_ming[-1] = T_lst_ming[-1]+T_root_ming
            
            
            T_lst_maxg.reverse()
            T_lst_ming.reverse()
        return T_lst_maxg, T_lst_ming, T_root_maxg, T_root_ming
    
    def calc_Vz_dist(self):
        
        Vz_maxg = []
        Vz_ming = []
        
        ylst_interm  = []
        L_lst_interm_maxg = []
        L_lst_interm_ming = []
        W_batt_interm     = []
        
        if self.batt == False:
            for i in range(1, len(self.cross_sections)+1):
                
                Li_maxg = self.Lmax_cs_lst[-i]
                Li_ming = self.Lmin_cs_lst[-i]
                
                ylst_interm.append(self.ylst[-i])
                
                L_lst_interm_maxg.append(Li_maxg)
                L_lst_interm_ming.append(Li_ming)
                
                
                if len(Vz_maxg) != 0:
                    L_toti_maxg = self.integrate(ylst_interm, L_lst_interm_maxg)
                    L_toti_ming = self.integrate(ylst_interm, L_lst_interm_ming)
                    
                    Vz_loc_maxg = -L_toti_maxg
                    Vz_loc_ming = -L_toti_ming
                    
                    Vz_maxg.append(Vz_loc_maxg)
                    Vz_ming.append(Vz_loc_ming)
                    
                if len(Vz_maxg) == 0:
                    Vz_maxg.append(-Li_maxg)
                    Vz_ming.append(-Li_ming)
                    
            Vz_root_maxg = -Vz_maxg[-1]
            Vz_root_ming = -Vz_ming[-1]
            
            Vz_maxg[-1] = Vz_maxg[-1]+Vz_root_maxg
            Vz_ming[-1] = Vz_ming[-1]+Vz_root_ming
            
            Vz_maxg.reverse()
            Vz_ming.reverse()
            
            
        if self.batt == True:
            for i in range(1, len(self.cross_sections)+1):
                
                Li_maxg = self.Lmax_cs_lst[-i]
                Li_ming = self.Lmin_cs_lst[-i]
                W_batti = self.batt_load_lst[-i]
                
                ylst_interm.append(self.ylst[-i])
                
                L_lst_interm_maxg.append(Li_maxg)
                L_lst_interm_ming.append(Li_ming)
                W_batt_interm.append(W_batti)
                
                
                if len(Vz_maxg) != 0:
                    L_toti_maxg = self.integrate(ylst_interm, L_lst_interm_maxg)
                    L_toti_ming = self.integrate(ylst_interm, L_lst_interm_ming)
                    W_batt_toti = self.integrate(ylst_interm, W_batt_interm)
                    
                    Vz_loc_maxg = -L_toti_maxg+W_batt_toti
                    Vz_loc_ming = -L_toti_ming+W_batt_toti
                    
                    Vz_maxg.append(Vz_loc_maxg)
                    Vz_ming.append(Vz_loc_ming)
                    
                if len(Vz_maxg) == 0:
                    Vz_maxg.append(-Li_maxg+W_batti)
                    Vz_ming.append(-Li_ming+W_batti)
                    
            Vz_root_maxg = -Vz_maxg[-1]
            Vz_root_ming = -Vz_ming[-1]
            
            Vz_maxg[-1] = Vz_maxg[-1]+Vz_root_maxg
            Vz_ming[-1] = Vz_ming[-1]+Vz_root_ming
            
            Vz_maxg.reverse()
            Vz_ming.reverse()
            
            
        return Vz_maxg, Vz_ming, Vz_root_maxg, Vz_root_ming
    
    def calc_Vx_dist(self):
                
        Vx_maxg = []
        Vx_ming = []
        
        ylst_interm  = []
        D_lst_interm_maxg = []
        D_lst_interm_ming = []
        
        
        for i in range(1, len(self.cross_sections)+1):
            
            Di_maxg = self.Dmax_cs_lst[-i]
            Di_ming = self.Dmin_cs_lst[-i]
            
            ylst_interm.append(self.ylst[-i])
            
            D_lst_interm_maxg.append(Di_maxg)
            D_lst_interm_ming.append(Di_ming)
            
            
            if len(Vx_maxg) != 0:
                D_toti_maxg = self.integrate(ylst_interm, D_lst_interm_maxg)
                D_toti_ming = self.integrate(ylst_interm, D_lst_interm_ming)
                
                Vx_loc_maxg = D_toti_maxg
                Vx_loc_ming = D_toti_ming
                
                Vx_maxg.append(Vx_loc_maxg)
                Vx_ming.append(Vx_loc_ming)
                
            if len(Vx_maxg) == 0:
                Vx_maxg.append(-Di_maxg)
                Vx_ming.append(-Di_ming)
                
        Vx_root_maxg = -Vx_maxg[-1]
        Vx_root_ming = -Vx_ming[-1]
        
        Vx_maxg[-1] = Vx_maxg[-1]+Vx_root_maxg
        Vx_ming[-1] = Vx_ming[-1]+Vx_root_ming
        
        Vx_maxg.reverse()
        Vx_ming.reverse()
        
        return Vx_maxg, Vx_ming, Vx_root_maxg, Vx_root_ming
    
    
    def calc_Mx_dist(self):
        Mx_maxg = []
        Mx_ming = []
        
        ylst_interm       = []
        L_lst_maxg_interm = []
        L_lst_ming_interm = []
        W_batt_interm     = []
        
        if self.batt == False:
            for i in range(1, len(self.Vz_maxg_lst)+1):
                ylst_interm.append(self.ylst[-i])
                L_lst_maxg_interm.append(self.Lmax_cs_lst[-i])
                L_lst_ming_interm.append(self.Lmin_cs_lst[-i])
                
                if len(Mx_maxg) != 0:
                    d_maxg = self.calc_force_loc(ylst_interm, L_lst_maxg_interm)
                    d_ming = self.calc_force_loc(ylst_interm, L_lst_ming_interm)
                    
                    L_maxg_toti = self.integrate(ylst_interm, L_lst_maxg_interm)
                    L_ming_toti = self.integrate(ylst_interm, L_lst_ming_interm)
                    
                    Mx_maxg_i = L_maxg_toti*d_maxg
                    Mx_ming_i = L_ming_toti*d_ming
                    
                    Mx_maxg.append(Mx_maxg_i)
                    Mx_ming.append(Mx_ming_i)
                    
                if len(Mx_maxg) == 0:
                    Mx_maxg.append(0)
                    Mx_ming.append(0)
                
            Mx_root_maxg = -Mx_maxg[-1]
            Mx_root_ming = -Mx_ming[-1]
            
            Mx_maxg[-1]  = Mx_maxg[-1]+Mx_root_maxg
            Mx_ming[-1]  = Mx_ming[-1]+Mx_root_ming
            
            Mx_maxg.reverse()
            Mx_ming.reverse()
        
        if self.batt == True:
            for i in range(1, len(self.Vz_maxg_lst)+1):
                ylst_interm.append(self.ylst[-i])
                L_lst_maxg_interm.append(self.Lmax_cs_lst[-i])
                L_lst_ming_interm.append(self.Lmin_cs_lst[-i])
                W_batt_interm.append(self.batt_load_lst[-i])
                
                if len(Mx_maxg) != 0:
                    d_maxg = self.calc_force_loc(ylst_interm, L_lst_maxg_interm)
                    d_ming = self.calc_force_loc(ylst_interm, L_lst_ming_interm)                    
                    d_batt = self.calc_force_loc(ylst_interm, W_batt_interm)
                    
                    L_maxg_toti = self.integrate(ylst_interm, L_lst_maxg_interm)
                    L_ming_toti = self.integrate(ylst_interm, L_lst_ming_interm)
                    W_batt_toti = self.integrate(ylst_interm, W_batt_interm)
                    
                    Mx_maxg_i = L_maxg_toti*d_maxg-W_batt_toti*d_batt
                    Mx_ming_i = L_ming_toti*d_ming-W_batt_toti*d_batt
                    
                    Mx_maxg.append(Mx_maxg_i)
                    Mx_ming.append(Mx_ming_i)
                
                if len(Mx_maxg) == 0:
                    Mx_maxg.append(0)
                    Mx_ming.append(0)
            
            Mx_root_maxg = -Mx_maxg[-1]
            Mx_root_ming = -Mx_ming[-1]
            
            Mx_maxg[-1]  = Mx_maxg[-1]+Mx_root_maxg
            Mx_ming[-1]  = Mx_ming[-1]+Mx_root_ming
            
            Mx_maxg.reverse()
            Mx_ming.reverse()
        
        return Mx_maxg, Mx_ming, Mx_root_maxg, Mx_root_ming
    
    def calc_Mz_dist(self):
        Mx_maxg = []
        Mx_ming = []
        
        ylst_interm       = []
        D_lst_maxg_interm = []
        D_lst_ming_interm = []
        
        for i in range(1, len(self.Vz_maxg_lst)+1):
            ylst_interm.append(self.ylst[-i])
            D_lst_maxg_interm.append(self.Dmax_cs_lst[-i])
            D_lst_ming_interm.append(self.Dmin_cs_lst[-i])
            
            if len(Mx_maxg) != 0:
                d_maxg = self.calc_force_loc(ylst_interm, D_lst_maxg_interm)
                d_ming = self.calc_force_loc(ylst_interm, D_lst_ming_interm)
                
                D_maxg_toti = self.integrate(ylst_interm, D_lst_maxg_interm)
                D_ming_toti = self.integrate(ylst_interm, D_lst_ming_interm)
                
                Mx_maxg_i = D_maxg_toti*d_maxg
                Mx_ming_i = D_ming_toti*d_ming
                
                Mx_maxg.append(Mx_maxg_i)
                Mx_ming.append(Mx_ming_i)
                
            if len(Mx_maxg) == 0:
                Mx_maxg.append(0)
                Mx_ming.append(0)
            
        Mx_root_maxg = -Mx_maxg[-1]
        Mx_root_ming = -Mx_ming[-1]
        
        Mx_maxg[-1]  = Mx_maxg[-1]+Mx_root_maxg
        Mx_ming[-1]  = Mx_ming[-1]+Mx_root_ming
        
        Mx_maxg.reverse()
        Mx_ming.reverse()
        
        return Mx_maxg, Mx_ming, Mx_root_maxg, Mx_root_ming
    
    def integrate(self, ylst, L_lst):
        A_lst = []
        for i in range(len(ylst)-1):
            dy = abs(ylst[i+1]-ylst[i])
            
            A_ub = L_lst[i+1]*dy
            A_lb = L_lst[i]*dy
            A_av = (A_ub+A_lb)/2
            
            A_lst.append(A_av)
            
        return sum(A_lst)
    
    def calc_force_loc(self, ylst_aero, cl_lst):
        
        yA = []
        A  = []
        
        for i in range(len(ylst_aero)-1):
            y = (ylst_aero[i]+ylst_aero[i+1])/2
            A_ub = cl_lst[i+1]*(ylst_aero[i+1]-ylst_aero[i])
            A_lb = cl_lst[i]*(ylst_aero[i+1]-ylst_aero[i])
            A_av = (A_ub+A_lb)/2
            
            yA.append(y*A_av)
            A.append(A_av)
            
        if sum(A) != 0:
            return sum(yA)/sum(A)
        
        if sum(A) == 0:
            return 0
    
        
    
    
    
        
                    
            
            
            
    
    
                
                
                    
               
        
          
        
        
        
    
    
        
        


aluminium=MatProps(sigma_y=450000000, E=72400000000, poisson=0.33, rho=2.87, G=28E9, name="AA2024", alpha=0.8, n=0.6)
z_stiff = Z_Stringer(1, aluminium)
j_stiff = J_Stringer(1, aluminium)
skin    = Skin()
wing = StiffenedWing(100, z_stiff, j_stiff, 50, 50, 0.25, 0.75, Spar,skin)