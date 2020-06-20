from plateclass import *
from stringerclass import *
from loadinfo import *
import copy
from scipy.interpolate import interp1d
from numpy import pi
import numpy as np
from math import sqrt, cos, sin
from matplotlib import pyplot as plt
import matplotlib.patches as patches


class CircCrossSection:
    def __init__(self, length, r, ts, material):
        self.l = length
        self.material = material

        self.radius = r
        self.thickness = ts

        self.total_area = 2*pi*r*ts

        self.xbar, self.ybar = r, r
        self.Ixx, self.Iyy, self.Ixy = pi*ts*r*r*r, pi*ts*r*r*r, 0

        self.E = self.material.E

        self.mass = self.total_area * self.l * self.material.rho

        self.internal_area = pi*r*r
        self.tresca_yield = self.material.sigma_y*self.material.sigma_y / 3.

        rhoxx = sqrt(self.Ixx/self.total_area)
        rhoyy = sqrt(self.Iyy/self.total_area)
        slender = self.l/min(rhoxx, rhoyy)

        self.sigma_cr = pi*pi*self.material.E / (slender*slender)

        self.stiffeners = []

    def addStringer(self, stringer: Stringer):
        self.stiffeners.append(stringer)

    def recalculate_component_effects(self):
        self.calc_total_area()
        self.xbar, self.ybar = self.calc_centroid()
        self.Ixx, self.Iyy, self.Ixy = self.momentofinertia()
        self.E = self.youngsmodulus()
        self.calc_mass()

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

    def recalculate_component_effects(self):
        self.calc_total_area()
        self.xbar, self.ybar = self.calc_centroid()
        self.Ixx, self.Iyy, self.Ixy = self.momentofinertia()
        self.E = self.youngsmodulus()
        self.calc_mass()
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
    def __init__(self):
        self.sections = []

        self.ys = []
        self.El = []
        self.Ixx = []
        self.Izz = []
        self.Ixz = []
        self.xbars = []
        self.zbars = []

        self.mass = None

    def addSection(self, section: FuselageSection):
        self.sections.append(section)

    def calc_total_mass(self):
        self.mass = 0
        for section in self.sections:
            self.mass += section.properties.mass
        return self.mass

    def recalculate_lists(self):
        self.ys = []
        self.E = []
        self.Ixx = []
        self.Izz = []
        self.Ixz = []
        self.xbars = []
        self.zbars = []
        self.areas = []
        for section in self.sections:
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
