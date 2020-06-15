from plateclass import *
from stringerclass import *
from loadinfo import *
import copy

class CrossSection:
    def __init__(self, length):
        self.l = length

        self.plates = []
        self.stiffeners = []

        self.total_area = 0

        self.xbar, self.ybar = None, None
        self.Ixx, self.Iyy, self.Ixy = None, None, None

        self.Eu, self.El = None, None

        self.mass = self.calc_mass()

    def addStringers(self, *stringers: Stringer):
        for stringer in stringers:
            stringer.properties.Le = self.l
            self.stiffeners.append(stringer)
        self.recalculate_component_effects()

    def addPanels(self, *panels: Sheet):
        for panel in panels:
            panel.properties.a = self.l
            self.plates.append(panel)
        self.recalculate_component_effects()

    def calc_total_area(self):
        self.total_area = 0
        for stiff in self.stiffeners:
            self.total_area += stiff.properties.total_area
        for plate in self.plates:
            self.total_area += plate.properties.total_area

    def calc_mass(self):
        mass = 0
        for stiff in self.stiffeners:
            mass += stiff.properties.mass
        for plate in self.plates:
            mass += plate.properties.mass
        return mass

    def recalculate_component_effects(self):
        self.calc_total_area()
        self.xbar, self.ybar = self.calc_centroid()
        self.Ixx, self.Iyy, self.Ixy = self.momentofinertia()
        self.Eu, self.El = self.youngsmodulus()
        self.mass = self.calc_mass()

    def calc_centroid(self):
        Qx, Qy = 0, 0
        for stiff in self.stiffeners:
            Qx += stiff.ybar_global * stiff.properties.total_area
            Qy += stiff.xbar_global * stiff.properties.total_area
        for plate in self.plates:
            Qx += plate.ybar_global * plate.properties.total_area
            Qy += plate.xbar_global * plate.properties.total_area
        return Qy/self.total_area, Qx/self.total_area

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
        Eu = 0
        El_num = 1
        for stiff in self.stiffeners:
            Eu += stiff.properties.material.E*stiff.properties.total_area/self.total_area
            El_num *= stiff.properties.material.E
        for plate in self.plates:
            Eu += plate.properties.material.E*plate.properties.total_area/self.total_area
            El_num *= plate.properties.material.E
        return Eu, El_num/Eu

class FuselageSection:
    def __init__(self, ypos, section_instance, rotation=0):
        self.properties = copy.deepcopy(section_instance)
        self.ypos = ypos
        self.rotation = rotation

        self.zbar = self.properties.ybar
        self.

        self.Ixx_global, self.Izz_global, self.Ixz_global = translate_mmoi(
            self.properties.Ixx, self.properties.Iyy, -self.properties.Ixy, rotation)



class Fuselage:
    def __init__(self):
        self.sections = []

    def addSection(self, section: FuselageSection):
        self.sections.append(section)
        self.recalculate_EIs()

    def recalculate_EIs(self):
        El = []
        Eu = []
        Ixx = []
        Iyy = []
        Ixy = []
        for section in self.sections:
            El.append(section.El)