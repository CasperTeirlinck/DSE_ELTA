try:
    from plateclass import *
    from stringerclass import *
    from loadinfo import *
except:
    from Structures.plateclass import *
    from Structures.stringerclass import *
    from Structures.loadinfo import *
import copy
from scipy.interpolate import interp1d
from numpy import pi
import numpy as np
from math import sqrt, cos, sin
from matplotlib import pyplot as plt
import matplotlib.patches as patches

def get_C(bovert):
    if bovert < 40:
        return 4.0
    elif bovert > 110:
        return 6.98
    else:
        return 4.0 + 2.98 * (bovert - 40) / 70

class EmptyCrossSectionProps:
    def __init__(self):
        self.E = 70000000000
        self.internal_area = 0.0
        self.mass = 0.0
        self.cost = 0.0
        self.l = 0.0
        self.total_area = 0.0
        plates = []

class EmptyCrossSection:
    def __init__(self):
        self.properties = EmptyCrossSectionProps()
        self.ypos = 0.0
        self.xpos = 0.0
        self.zpos = 0.0
        self.rotation = 0.0
        self.Ixx_global = 0.0000000001
        self.Izz_global = 0.0000000001
        self.Ixz_global = 0.0000000001
        self.xbar_global = 0.0
        self.zbar_global = 0.0


class CircCrossSection:
    def __init__(self, length, r, ts, material):
        self.l = length
        self.material = material

        self.radius = r
        self.thickness = ts
        self.circumference = 2*pi*r

        self.total_area = self.circumference*ts

        self.xbar, self.ybar = r, r
        self.Ixx, self.Iyy, self.Ixy = pi*ts*r*r*r, pi*ts*r*r*r, 0

        self.E = self.material.E

        self.mass = self.total_area * self.l * self.material.rho

        self.internal_area = pi*r*r
        self.tresca_yield = self.material.sigma_y*self.material.sigma_y / 3.

        self.stiffeners = []
        self.we2s = []
        self.stringerpitches = []

        self.sigma_cr = self.calc_critical_stress()
        self.calc_sigma_ult_tensile()


    def calc_we(self):
        for stringer in self.stiffeners:
            C = get_C(self.circumference/self.l)
            self.we2s.append(self.thickness * sqrt(C*pi*pi*self.material.E/(12*(1-self.material.poisson**2)*stringer.properties.sigma_cc)))

    def calc_stringerpitch(self):
        for index, stringer in enumerate(self.stiffeners):
            if index == 0:
                if len(self.stiffeners) != 1:
                    stringerpitch = self.stiffeners[index+1].xpos - stringer.xpos
                else:
                    stringerpitch = self.circumference*0.5
            elif index == len(self.stiffeners)-1:
                stringerpitch = stringer.xpos-self.stiffeners[index-1].xpos
            else:
                stringerpitch = (self.stiffeners[index+1].xpos - self.stiffeners[index-1].xpos)/2
            self.stringerpitches.append(stringerpitch)

    def calc_critical_stress(self):
        rhoxx = sqrt(self.Ixx / self.total_area)
        rhoyy = sqrt(self.Iyy / self.total_area)
        slender = self.l / min(rhoxx, rhoyy)
        sigma_cr_circ = pi * pi * self.material.E / (slender * slender)
        stiffenerforces = 0
        totalstiffenerarea = 0
        if sum(self.we2s) > self.circumference:
            f_plate = 0
            f_stiff = [self.circumference / len(self.stiffeners) for _ in range(len(self.stiffeners))]
        else:
            f_plate = self.circumference - sum(self.we2s)
            f_stiff = self.we2s
        for we, stiff in zip(f_stiff, self.stiffeners):
            stiffenerforces += stiff.properties.sigma_cr * (stiff.properties.total_area + we * self.thickness)
            totalstiffenerarea += stiff.properties.total_area
        return (sigma_cr_circ * self.thickness * f_plate + stiffenerforces) / (totalstiffenerarea + self.thickness * self.circumference)

    def calc_sigma_ult_tensile(self):
        size = len(self.stiffeners)+1
        if size == 1:
            self.sigma_ult_tensile = self.material.sigma_y
            self.stresses = np.array([1./self.total_area])
            self.averagestress = 1./self.total_area
            self.ultimatestress = self.material.sigma_y
            self.tresca_yield = self.ultimatestress*self.ultimatestress / 3
        else:
            forcematrix = np.zeros((size, size))
            forcematrix[0,:] = 1
            forcematrix[1:,0] = self.l/(self.total_area*self.material.E)
            for idx, stringer in enumerate(self.stiffeners):
                forcematrix[idx+1, idx+1] = -stringer.properties.Le/(stringer.properties.total_area*stringer.properties.material.E)
            bvector = np.zeros(size)
            bvector[0] = 1.
            try:
                fractions = np.linalg.solve(forcematrix, bvector)
            except np.linalg.LinAlgError:
                print(forcematrix, bvector)
                print("Error in ult stress calculation in a circular section, matrix is singular!")

            self.stresses = np.copy(fractions)
            self.stresses[0] = self.stresses[0]/(self.total_area)
            for idx, stringer in enumerate(self.stiffeners):
                self.stresses[idx+1] = self.stresses[idx+1]/(stringer.properties.total_area)
            self.averagestress = np.sum(fractions)/self.total_area
            ult_stresses = []
            for stringer, fraction in zip(self.stiffeners, fractions):
                ult_stresses.append(stringer.properties.material.sigma_y/fraction)
            ult_stresses = np.array([ult_stresses])
            self.ultimatestress = np.min(ult_stresses)
            self.tresca_yield = self.ultimatestress*self.ultimatestress / 3

    def addStringer(self, stringer: Stringer):
        self.stiffeners.append(stringer)

    def recalculate_component_effects(self):
        self.calc_total_area()
        self.xbar, self.ybar = self.calc_centroid()
        self.Ixx, self.Iyy, self.Ixy = self.momentofinertia()
        self.E = self.youngsmodulus()
        self.calc_mass()
        self.calc_cost()
        self.calc_we()
        self.calc_stringerpitch()
        self.sigma_cr = self.calc_critical_stress()
        self.calc_sigma_ult_tensile()


    def calc_total_area(self):
        self.total_area = 2*pi*self.radius*self.thickness
        for stiff in self.stiffeners:
            self.total_area += stiff.properties.total_area

    def calc_centroid(self):
        self.Qx, self.Qy = 2*pi*self.radius*self.radius*self.thickness, 2*pi*self.radius*self.radius*self.thickness
        for stiff in self.stiffeners:
            self.Qx += stiff.ybar_global * stiff.properties.total_area
            self.Qy += stiff.xbar_global * stiff.properties.total_area
        return self.Qy/self.total_area, self.Qx/self.total_area

    def momentofinertia(self):
        Ixx = pi*self.thickness*self.radius*self.radius*self.radius + self.radius*2*pi*self.thickness*(self.radius-self.ybar)**2
        Iyy = pi*self.thickness*self.radius*self.radius*self.radius + self.radius*2*pi*self.thickness*(self.radius-self.xbar)**2
        Ixy = self.radius*2*pi*self.thickness*(self.radius-self.ybar)*(self.radius-self.xbar)
        for stiff in self.stiffeners:
            Ixx += stiff.Ixx_global + stiff.properties.total_area*(stiff.ybar_global-self.ybar)**2
            Iyy += stiff.Iyy_global + stiff.properties.total_area*(stiff.xbar_global-self.xbar)**2
            Ixy += stiff.Ixy_global + stiff.properties.total_area*(stiff.ybar_global-self.ybar)*(stiff.xbar_global-self.xbar)
        return Ixx, Iyy, Ixy

    def youngsmodulus(self):
        E = self.material.E*self.radius*2*pi*self.thickness/self.total_area
        for stiff in self.stiffeners:
            E += stiff.properties.material.E*stiff.properties.total_area/self.total_area
        return E

    def calc_mass(self):
        self.mass = self.radius*2*pi*self.thickness*self.l*self.material.rho
        for stiff in self.stiffeners:
            self.mass += stiff.properties.mass
        return self.mass

    def calc_cost(self):
        self.cost = self.radius*2*pi*self.thickness*self.l*self.material.rho*self.material.cost
        for stiff in self.stiffeners:
            self.cost += stiff.properties.cost
        return self.cost

    def plot_section(self, show=True):
        fig, ax = plt.subplots(1)
        ax.set_xlim(-0.3, 0.4)
        ax.set_ylim(-0.2, 0.6)
        stringerplot = plt.Circle((self.xbar, self.ybar), radius=self.radius, color='r', fill=False)
        ax.add_artist(stringerplot)
        for stiff in self.stiffeners:
            x, y = stiff.xpos, stiff.ypos
            stringerplot = plt.Circle((x, y), 0.005, color='b', fill=True)
            ax.add_artist(stringerplot)
        if show:
            plt.show()


class CrossSection:
    def __init__(self, length):
        self._l = length

        self.plates = []
        self.stiffeners = []

        self.total_area = 0
        self.internal_area = 0
        self.Qx, self.Qy = None, None
        self.xbar, self.ybar = None, None
        self.Ixx, self.Iyy, self.Ixy = None, None, None

        self.E = None

        self.mass = 0


    @property
    def l(self):
        return self._l

    @l.setter
    def l(self, val):
        self._l = val
        self.recalculate_component_effects()

    def resetData(self, length):
        self.__init__(length)

    def addStringers(self, *stringers: Stringer):
        for stringer in stringers:
            self.stiffeners.append(stringer)
        # self.recalculate_component_effects()

    def addStringer(self, stringer: Stringer):
        self.stiffeners.append(stringer)
        # self.recalculate_component_effects()

    def addPanels(self, *panels: Panel):
        for panel in panels:
            panel.properties.a = self.l
            self.plates.append(panel)
        # self.recalculate_component_effects()

    def calc_total_area(self):
        self.total_area = 0
        for stiff in self.stiffeners:
            self.total_area += stiff.properties.total_area
        for plate in self.plates:
            self.total_area += plate.properties.total_area
        self.internal_area = 0
        xmin, xmax, ymin, ymax = None, None, None, None
        for plate in self.plates:
            if xmin is None or plate.xpos < xmin:
                xmin = plate.xpos
            if xmax is None or plate.xpos > xmax:
                xmax = plate.xpos
            if ymin is None or plate.ypos < ymin:
                ymin = plate.ypos
            if ymax is None or plate.ypos > ymax:
                ymax = plate.ypos
            if plate.rotation != 0 and abs(plate.rotation-0.5*pi) > 0.001 and abs(plate.rotation-pi) > 0.001 and abs(plate.rotation-1.5*pi) > 0.001:
                self.internal_area -= plate.properties.b * plate.properties.b * 0.5
        self.internal_area += (xmax-xmin)*(ymax-ymin)

    def calc_mass(self):
        self.mass = 0
        for stiff in self.stiffeners:
            self.mass += stiff.properties.mass
        for plate in self.plates:
            self.mass += plate.properties.mass
        return self.mass

    def calc_cost(self):
        self.cost = 0
        for stiff in self.stiffeners:
            self.cost += stiff.properties.cost
        for plate in self.plates:
            self.cost += plate.properties.cost
        return self.cost

    def recalculate_component_effects(self):
        self.calc_total_area()
        self.xbar, self.ybar = self.calc_centroid()
        self.Ixx, self.Iyy, self.Ixy = self.momentofinertia()
        self.E = self.youngsmodulus()
        self.calc_mass()
        self.calc_cost()
        self.calc_Qs()

    def calc_centroid(self):
        self.Qx, self.Qy = 0, 0
        for stiff in self.stiffeners:
            self.Qx += stiff.ybar_global * stiff.properties.total_area
            self.Qy += stiff.xbar_global * stiff.properties.total_area
        for plate in self.plates:
            self.Qx += plate.ybar_global * plate.properties.total_area
            self.Qy += plate.xbar_global * plate.properties.total_area
        return self.Qy/self.total_area, self.Qx/self.total_area

    def momentofinertia(self):
        Ixx, Iyy, Ixy = 0., 0., 0.
        for stiff in self.stiffeners:
            Ixx += stiff.Ixx_global + stiff.properties.total_area*(stiff.ybar_global-self.ybar)**2
            Iyy += stiff.Iyy_global + stiff.properties.total_area*(stiff.xbar_global-self.xbar)**2
            Ixy += stiff.Ixy_global + stiff.properties.total_area*(stiff.ybar_global-self.ybar)*(stiff.xbar_global-self.xbar)
        for plate in self.plates:
            Ixx += plate.Ixx_global + plate.properties.total_area*(plate.ybar_global-self.ybar)**2
            Iyy += plate.Iyy_global + plate.properties.total_area*(plate.xbar_global-self.xbar)**2
            Ixy += plate.Ixy_global + plate.properties.total_area*(plate.ybar_global-self.ybar)*(plate.xbar_global-self.xbar)
        return Ixx, Iyy, Ixy

    def youngsmodulus(self):
        E = 0
        for stiff in self.stiffeners:
            E += stiff.properties.material.E*stiff.properties.total_area/self.total_area
        for plate in self.plates:
            E += plate.properties.material.E*plate.properties.total_area/self.total_area
        return E

    def calc_Qs(self):
        Qx = 0
        Qy = 0
        self.Qxs = []
        self.Qys = []
        self.Q_xs = []
        self.Q_ys = []
        self.Q_xs = []
        self.Q_ys = []
        for plate in self.plates:
            Qx = [plate.properties.area1*(plate.pos_i1_global[1]-self.ybar),
                  plate.properties.area2*(plate.pos_i2_global[1]-self.ybar),
                  plate.properties.area3*(plate.pos_i3_global[1]-self.ybar)]
            Qy = [plate.properties.area1*(plate.pos_i1_global[0]-self.xbar),
                  plate.properties.area2*(plate.pos_i2_global[0]-self.xbar),
                  plate.properties.area3*(plate.pos_i3_global[0]-self.xbar)]
            self.Qxs.extend(Qx)
            self.Qys.extend(Qy)
            self.Q_xs.extend([plate.pos_i1_global[0], plate.pos_i2_global[0], plate.pos_i3_global[0]])
            self.Q_ys.extend([plate.pos_i1_global[1], plate.pos_i2_global[1], plate.pos_i3_global[1]])
        self.Qxs = np.asarray(self.Qxs)
        self.Qys = np.asarray(self.Qys)
        self.Q_xs = np.asarray(self.Q_xs)
        self.Q_ys = np.asarray(self.Q_ys)

    def plot_section(self, show=True):
        fig, ax = plt.subplots(1)
        ax.set_xlim(-0.3, 0.4)
        ax.set_ylim(-0.2, 0.6)
        for stiff in self.stiffeners:
            stringerplot = plt.Circle((-stiff.xpos, stiff.ypos), 0.002, color='g', fill=True)
            ax.add_artist(stringerplot)
        for idx, plate in enumerate(self.plates):
            br = (-plate.xpos, plate.ypos)
            rect = patches.Rectangle(br, -plate.properties.b, plate.properties.ts, 360-plate.rotation*180/pi, linewidth=1, edgecolor='r',facecolor='r')
            ax.add_patch(rect)
            # brpoint = plt.Circle(br, 0.005, color='g', fill=True)
            # ax.add_artist(brpoint)
            for stiff in plate.properties.stringers:
                x, y = -stiff.xpos*cos(plate.rotation)+ stiff.ypos*sin(plate.rotation), \
                       stiff.xpos*sin(plate.rotation)+ stiff.ypos*cos(plate.rotation)
                x -= plate.xpos
                y += plate.ypos
                stringerplot = plt.Circle((x, y), 0.005, color='b', fill=True)
                ax.add_artist(stringerplot)
        if show:
            plt.show()



class FuselageSection:
    def __init__(self, ypos, xpos, zpos, section_instance, rotation=0):
        self.properties = copy.deepcopy(section_instance)
        self.xpos = xpos
        self.ypos = ypos
        self.zpos = zpos
        self.rotation = rotation

        self.xbar_global, self.zbar_global = -self.properties.xbar*cos(rotation) + self.properties.ybar*sin(rotation), \
                                             self.properties.xbar*sin(rotation) + self.properties.ybar*cos(rotation)
        self.zbar_global += zpos
        self.xbar_global += xpos

        self.Ixx_global, self.Izz_global, self.Ixz_global = translate_mmoi(
            self.properties.Ixx, self.properties.Iyy, -self.properties.Ixy, rotation)



class Fuselage:
    def __init__(self, framelocs):
        self.sections = []

        self.framelocs = framelocs

        self.ys = []
        self.El = []
        self.Ixx = []
        self.Izz = []
        self.Ixz = []
        self.xbars = []
        self.zbars = []
        self.xbar, self.ybar, self.zbar = None, None, None

        self.mass = None

        self.cost = None

    def addSection(self, section: FuselageSection):
        self.sections.append(section)

    def calc_total_mass(self):
        self.mass = 0
        for section in self.sections:
            self.mass += section.properties.mass
            if isinstance(section.properties, CircCrossSection):
                self.mass += section.properties.internal_area * 0.75 * 0.002 * section.properties.material.rho
            elif isinstance(section.properties, CrossSection):
                self.mass += section.properties.internal_area * 0.75 * 0.002 * section.properties.plates[0].properties.material.rho
        return self.mass

    def calc_total_cost(self):
        self.cost = 0
        for section in self.sections:
            self.cost += section.properties.cost
            if isinstance(section.properties, CircCrossSection):
                self.cost += section.properties.internal_area * 0.75*0.002*section.properties.material.rho*section.properties.material.cost
            elif isinstance(section.properties, CrossSection):
                self.cost += section.properties.internal_area * 0.75 * 0.002 * section.properties.plates[0].properties.material.rho * section.properties.plates[0].properties.material.cost
        return self.cost

    def calc_cg(self):
        self.Qyy = 0
        self.Qz = 0
        self.Qx = 0
        for section, y in zip(self.sections, self.framelocs):
            self.Qyy += section.properties.mass * (y - section.properties.l * 0.5)
            if isinstance(section.properties, CircCrossSection):
                self.Qyy += section.properties.internal_area * 0.75 * 0.002 * section.properties.material.rho * y
            elif isinstance(section.properties, CrossSection):
                self.Qyy += section.properties.internal_area * 0.75 * 0.002 * section.properties.plates[0].properties.material.rho * y
            self.Qz += section.properties.mass * section.xbar_global
            self.Qx += section.properties.mass * section.zbar_global
        self.ybar = self.Qyy/self.mass
        self.xbar, self.zbar = self.Qz/self.mass, self.Qx/self.mass

    def recalculate_lists(self):
        self.ys = []
        self.E = []
        self.Ixx = []
        self.Izz = []
        self.Ixz = []
        self.xbars = []
        self.zbars = []
        self.areas = []
        self.sections_amount = 0
        for section in self.sections:
            if not isinstance(section, EmptyCrossSection):
                self.sections_amount += 1
                self.ys.append(section.ypos)
                self.E.append(section.properties.E)
                self.Ixx.append(section.Ixx_global)
                self.Izz.append(section.Izz_global)
                self.Ixz.append(section.Ixz_global)
                self.xbars.append(section.xbar_global)
                self.zbars.append(section.zbar_global)
                self.areas.append(section.properties.internal_area)
        self.ys = np.asarray(self.ys)
        self.E = np.asarray(self.E)
        self.Ixx, self.Izz, self.Ixz = np.asarray(self.Ixx), np.asarray(self.Izz), np.asarray(self.Ixz)
        self.xbars, self.zbars = np.asarray(self.xbars), np.asarray(self.zbars)
        self.areas = np.asarray(self.areas)
        self.calc_total_mass()
        self.calc_total_cost()
        self.calc_cg()
