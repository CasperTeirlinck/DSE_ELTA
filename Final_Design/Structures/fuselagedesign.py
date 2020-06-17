from plateclass import *
from stringerclass import *
from loadinfo import *
import copy
from scipy.interpolate import interp1d
from numpy import pi
import numpy as np
from math import sqrt, cos, sin
from matplotlib import pyplot as plt

class CircCrossSection:
    def __init__(self, length, r, ts, material):
        self.l = length
        self.material = material

        self.total_area = 2*pi*r*ts

        self.xbar, self.ybar = r, r
        self.Ixx, self.Iyy, self.Ixy = pi*ts*2*r*r*r, pi*ts*2*r*r*r, 0

        self.E = self.material.E

        self.mass = self.total_area * self.material.rho


class CrossSection:
    def __init__(self, length):
        self._l = length

        self.plates = []
        self.stiffeners = []

        self.total_area = 0

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
            # print(plate.rotation)
            # print("Base", plate.Ixx_global)
            # print("Steiner", plate.properties.total_area*(plate.ybar_global-self.ybar)**2)
            # print()
            Ixx += plate.Ixx_global + plate.properties.total_area*(plate.ybar_global-self.ybar)**2
            Iyy += plate.Iyy_global + plate.properties.total_area*(plate.xbar_global-self.xbar)**2
            # print("Base", plate.Ixy_global)
            # print("Steiner", plate.properties.total_area*(plate.ybar_global-self.ybar)*(plate.xbar_global-self.xbar))
            # print("Rotation", plate.rotation)
            # print()
            Ixy += plate.Ixy_global + plate.properties.total_area*(plate.ybar_global-self.ybar)*(plate.xbar_global-self.xbar)
        return Ixx, Iyy, Ixy

    def youngsmodulus(self):
        E = 0
        for stiff in self.stiffeners:
            E += stiff.properties.material.E*stiff.properties.total_area/self.total_area
        for plate in self.plates:
            E += plate.properties.material.E*plate.properties.total_area/self.total_area
        return E


class FuselageSection:
    def __init__(self, ypos, xpos, zpos, section_instance, rotation=0):
        self.properties = copy.deepcopy(section_instance)
        self.ypos = ypos
        self.rotation = rotation

        self.xbar_global, self.zbar_global = -self.properties.xbar*cos(rotation) - self.properties.ybar*sin(rotation), \
                                             -self.properties.xbar*sin(rotation) + self.properties.ybar*cos(rotation)
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

    def addSection(self, section: FuselageSection):
        self.sections.append(section)

    def recalculate_lists(self):
        self.ys = []
        self.E = []
        self.Ixx = []
        self.Izz = []
        self.Ixz = []
        self.xbars = []
        self.zbars = []
        for section in self.sections:
            self.ys.append(section.ypos)
            self.E.append(section.properties.E)
            self.Ixx.append(section.Ixx_global)
            self.Izz.append(section.Izz_global)
            self.Ixz.append(section.Ixz_global)
            self.xbars.append(section.xbar_global)
            self.zbars.append(section.zbar_global)
        self.ys = np.asarray(self.ys)
        self.E = np.asarray(self.E)
        self.Ixx, self.Izz, self.Ixz = np.asarray(self.Ixx), np.asarray(self.Izz), np.asarray(self.Ixz)
        self.xbars, self.zbars = np.asarray(self.xbars), np.asarray(self.zbars)

def create_compound_segment(material, length, ts, height, width, r_top, r_bot):
    cs = CrossSection(length)
    wbot = width-2*r_bot
    wtop = width-2*r_top
    hleft = height-r_bot-r_top
    hright = hleft
    xs = [r_bot, r_bot+wbot, width, width, r_top+wtop, r_top, 0., 0.]
    ys = [0., 0., r_bot, r_bot+hleft, height, height, r_bot+hright, r_bot]
    widths = [wbot, sqrt(2)*r_bot, hleft, sqrt(2)*r_top, wtop, sqrt(2)*r_top, hright, sqrt(2)*r_bot]
    for idx, (x, y, w_sheet) in enumerate(zip(xs, ys, widths)):
        angle = idx*pi*0.25
        sheet = Sheet(length, w_sheet, ts, material, mode=2)
        pan = Panel(x, y, sheet, rotation=angle)
        cs.addPanels(pan)
    return cs

def create_round_segment(material, length, ts, height, width=None, r_top=None, r_bot=None):
    cs = CircCrossSection(length, height*0.5, ts, material)
    return cs


def generate_fuselage(data, sparlocs, zshift, base_material, circ_material=None, stiff_material=None, base_ts=0.002):
    if circ_material is None:
        circ_material = base_material
    if stiff_material is None:
        stiff_material = base_material
    assert(len(sparlocs)==len(zshift))
    fus = Fuselage()
    xref = 0.0
    ys = data[:,0]
    interpw = interp1d(ys, data[:,1])
    interph = interp1d(ys, data[:,2])
    interprt = interp1d(ys, data[:,3])
    interprb = interp1d(ys, data[:,4])

    for idx, (y, zref) in enumerate(zip(sparlocs, zshift)):
        if idx == len(sparlocs)-1:
            length = 2.2-y
        else:
            length = y-sparlocs[idx+1]
        if y < 6.099:
            width = interpw(y)
            cs = create_compound_segment(base_material, length, base_ts, interph(y), width, interprt(y), interprb(y))
        else:
            width = interph(y)
            cs = create_round_segment(circ_material, length, base_ts, interph(y), width, interprt(y), interprb(y))

        fus.addSection(FuselageSection(y, xref+width/2, zref, cs))

    fus.recalculate_lists()
    return fus, sparlocs

def make_default_fuselage():
    data = np.genfromtxt("fuselage_data.csv", delimiter=",", skip_header=1)
    data *= 0.001

    aluminium = MatProps(sigma_y=450000000, E=72400000000, poisson=0.33, rho=2.87, name="AA2024", alpha=0.8,
                         n=0.6)
    carbon = MatProps(sigma_y=600000000, E=70000000000, poisson=0.1, rho=1.60, sigma_comp=570000000,
                      name="carbonfibre")
    sparlocs = np.array([2.2, 2.3, 2.5, 2.755, 3.9, 5.9, 6., 6.05, 6.1, 6.5, 7., 8.5, 9., 9.4])[::-1]
    zfunc = lambda y: 0
    zshift = [zfunc(y) for y in sparlocs]
    fus, sparlocs = generate_fuselage(data, sparlocs, zshift, base_material=aluminium, circ_material=carbon)
    return fus, sparlocs

def analyze_fuselage(sparlocs, fus, bcs, *loads):
    EIx = fus.E*fus.Ixx
    EIz = fus.E*fus.Izz
    # plt.plot(sparlocs, EIx)
    # plt.show()
    system = LoadCase(y_coordinates=np.asarray(sparlocs), EIx=EIx, EIz=EIz)
    ymin = min(sparlocs)
    flag = False
    for load in loads:
        if type(load) == Force or type(load) == Moment:
            system.add_loads(load)
            if load.position[1] < ymin:
                ymin = load.position[1]
                flag = True
        elif type(load) == dict:
            system.add_loads(*load.values())
            for single_load in load.values():
                # print(single_load)
                if single_load.position[1] < ymin:
                    ymin = single_load.position[1]
                    flag = True
        else:
            print("Unrecognized load:", load)
    for bc in bcs:
        # print(bc)
        system.add_boundary_condition(**bc)
        if bc["y"] < ymin:
            ymin = bc["y"]
            flag = True
    if flag == True:
        sparlocs = np.concatenate((sparlocs, np.array([ymin])))
        EIx = np.concatenate((EIx, np.array([EIx[-1]])))
        EIz = np.concatenate((EIz, np.array([EIz[-1]])))
        system.change_geometry(sparlocs, EIx, EIz)

    system.solve_bcs()

    sparlocs = np.concatenate
    sx, sz = system.get_shearforce(sparlocs)
    mx, mz = system.get_moment(sparlocs)
    ny = system.get_normalforce(sparlocs)
    totalforce, totalmoment, totalforce_res, totalmoment_res = \
        np.array([0.0,0.0,0.0]), np.array([0.0,0.0,0.0]), np.array([0.0,0.0,0.0]), np.array([0.0,0.0,0.0])
    for force in system.forces:
        totalforce += force.vector
        totalmoment += np.cross(force.position, force.vector)
    for force in system.reactionforces:
        totalforce_res += force.vector
        totalmoment_res += np.cross(force.position, force.vector)
    for moment in system.moments:
        totalmoment += moment.vector
    for moment in system.reactionmoments:
        totalmoment_res += moment.vector
    plt.subplot(131)
    plt.plot(sparlocs, ny)
    plt.subplot(132)
    plt.plot(sparlocs, sx)
    plt.plot(sparlocs, sz)
    plt.subplot(133)
    plt.plot(sparlocs, mx)
    plt.plot(sparlocs, mz)
    plt.show()

if __name__ == "__main__":
    fus, sparlocs = make_default_fuselage()
    y_origin=1.835
    bc1 = {'y' : 0.0+y_origin, 'defl_x' : 0.0}
    bc2 = {'y' : 0.0+y_origin, 'defl_z' : 0.0}
    bc3 = {'y' : 0.0+y_origin, 'angle_x' : 0.0}
    bc4 = {'y' : 0.0+y_origin, 'angle_z' : 0.0}
    bc5 = {'y' : 0.0+y_origin, 'defl_y' : 0.0}
    bc6 = {'y' : 0.0+y_origin, 'angle_y' : 0.0}
    bcs = [bc1, bc2, bc3, bc4, bc5, bc6]
    loads = {}
    # loads['f1'] = Force(xpos=0.0, ypos=0.5, zpos=0.0, xmag=1., ymag=0.0, zmag=0.0)
    # loads['f2'] = Force(xpos=0.0, ypos=0.75, zpos=0.0, xmag=-50.0, ymag=0.0, zmag=10.0)
    # loads['f3'] = Force(xpos=0.0, ypos=1.0, zpos=0.0, xmag=2.0, ymag=0.0, zmag=0.0)
    # loads['m1'] = Moment(xpos=0.0, ypos=0.8, zpos=0.0, xmag=0.0, ymag=0.0, zmag=50.)
    # loads['m1'] = Moment(xpos=0.0, ypos=1.0, zpos=0.0, xmag=0.0, ymag=0.0, zmag=-25.)
    # loads["fhtail"] = Force(xpos=0.0, ypos=8.9, zpos=0.0, xmag=0, zmag=-700, ymag=70)
    loads["fvtail"] = Force(xpos=0.0, ypos=8.9, zpos=0.0, xmag=300, ymag=30, zmag=0)
    analyze_fuselage(sparlocs, fus, bcs, loads)
    # print(fus.Ixx)
    # print(fus.Izz)
    # print(fus.Ixz)
    # print(fus.xbars)
    # print(fus.zbars)
    # print(fus.E)
    # print(fus.Ixx*fus.E)
