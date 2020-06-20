import numpy as np
from math import pi, sqrt, cos, sin
import itertools
from matplotlib import pyplot as plt
from stringerclass import *
import copy


def get_C(bovert):
    if bovert < 40:
        return 4.0
    elif bovert > 110:
        return 6.98
    else:
        return 4.0 + 2.98 * (bovert - 40) / 70


def make_stringers(positions, stringerinstance, rotation=0):
    stringers = []
    for pos in list(positions):
        stringers.append(Stringer(*pos, stringerinstance, rotation=rotation))
    stringers = sorted(stringers, key=lambda x: x.xpos) # Sort stringers by x coordinate
    return stringers


def make_panel(panel, *stringers):
    for stringer in stringers:
        panel.addStringer(stringer)
    panel.recalculate_stringer_effects()
    return None


class Sheet:
    def __init__(self, a, b, ts, material: MatProps, mode=1):
        # Rotation in radians, clockwise positive. Keep 0 for correct orientation for bottom plate.
        self.material = material
        self.mode = mode

        self.ts = ts # Thickness of the panel
        self._a = a # Length of the panel
        self.b = b # Width of the panel
        # self.rotation=rotation

        self.stringers = []
        self.we2s = []
        self.stringerpitches = []

        self.total_area = None # Total area of the stiffened plate cross-section
        self.calc_total_area()

        self.sigma_cr = self.calc_critical_stress(mode=self.mode)
        self.calc_sigma_ult_tensile()

        self.stresses = None # Stresses on plate and subsequently the stringers for a unit tensile load
        self.averagestress = None # Average stress on stiffened panel for a unit tensile load
        self.ultimatestress = None # Ultimate tensile stress that the panel can take before yielding
        self.calc_sigma_ult_tensile()

        self.Qx = None
        self.Qy = None

        self.xbar, self.ybar = self.calc_centroid()

        self.Ixx, self.Iyy, self.Ixy = self.momentofinertia()

        self.mass = self.get_mass()

    @property
    def a(self):
        return self._a

    @a.setter
    def a(self, val):
        self._a = val
        for stringer in self.stringers:
            stringer.properties.Le = val
        self.recalculate_stringer_effects()

    def construct_panel(self, *stringers: Stringer):
        for stringer in stringers:
            self.stringers.append(stringer)
        self.recalculate_stringer_effects()

    def recalculate_stringer_effects(self):
        self.calc_we()
        self.calc_stringerpitch()
        self.sigma_cr = self.calc_critical_stress(mode=self.mode)
        self.calc_total_area()
        self.calc_sigma_ult_tensile()
        self.xbar, self.ybar = self.calc_centroid()
        self.Ixx, self.Iyy, self.Ixy = self.momentofinertia()
        self.mass = self.get_mass()
        self.idealize()

    def calc_we(self):
        for stringer in self.stringers:
            C = get_C(self.b/self._a)
            self.we2s.append(self.ts * sqrt(C*pi*pi*self.material.E/(12*(1-self.material.poisson**2)*stringer.properties.sigma_cc)))

    def calc_stringerpitch(self):
        for index, stringer in enumerate(self.stringers):
            if index == 0:
                if len(self.stringers) != 1:
                    stringerpitch = self.stringers[index+1].xpos - stringer.xpos
                else:
                    stringerpitch = self.b*0.5
            elif index == len(self.stringers)-1:
                stringerpitch = stringer.xpos-self.stringers[index-1].xpos
            else:
                stringerpitch = (self.stringers[index+1].xpos - self.stringers[index-1].xpos)/2
            self.stringerpitches.append(stringerpitch)

    def calc_critical_stress(self, mode):
        # mode determines the type of support of the panel. mode 1 is CCSS, 2 is SSSS,
        # 3 is CCCC, 4 is SSCC, 5 is SSCS, 6 is CCFF, 7 is SSFC, 8 is SSFF, 9 is SSFS
        # If required mode is anything else, the program asks for manual input. Look at
        # https://puu.sh/FSyYi/82a36e2870.png to find the correct value for C.
        aoverb = self._a/self.b
        # if aoverb < 3:
        if aoverb < 0:
            C = input("Manual input required because a/b < 3, look at https://puu.sh/FSyYi/82a36e2870.png for data. "
                      "a/b = {}".format(aoverb))
        else:
            if mode == 1 or mode == 2:
                C = 4.0
            elif mode == 3 or mode == 4:
                C = 6.98
            elif mode == 5:
                C = 6.98
            elif mode == 6:
                C = 3.72
            elif mode == 7:
                if round(self.material.poisson, 1) != 0.3:
                    C = 1.33
                else:
                    C = 1.28
            elif mode == 8:
                C = 0.92
            elif mode == 9:
                if round(self.material.poisson, 1) != 0.3:
                    C = 0.456
                else:
                    C = 0.425
            else:
                C = input("Manual input required because mode is unknown, look at https://puu.sh/FSyYi/82a36e2870.png "
                          "for data. a/b = {}".format(aoverb))
        self.sigma_cr_skin = C * pi*pi*self.material.E*(self.ts/self.b)**2/(12*(1-self.material.poisson**2))
        stiffenerforces = 0
        totalstiffenerarea = 0
        if sum(self.we2s) > self.b:
            f_plate = 0
            f_stiff = [self.b/len(self.stringers) for _ in range(len(self.stringers))]
        else:
            f_plate = self.b-sum(self.we2s)
            f_stiff = self.we2s
        for we, stiff in zip(f_stiff, self.stringers):
            stiffenerforces += stiff.properties.sigma_cr*(stiff.properties.total_area+we*self.ts)
            totalstiffenerarea += stiff.properties.total_area
        return (self.sigma_cr_skin * self.ts * f_plate + stiffenerforces) / (totalstiffenerarea+ self.ts * self.b)

    def calc_total_area(self):
        self.total_area = self.ts * self.b
        for stringer in self.stringers:
            self.total_area += stringer.properties.total_area

    def calc_sigma_ult_tensile(self):
        size = len(self.stringers)+1
        if size == 1:
            self.sigma_ult_tensile = self.material.sigma_y
            self.stresses = np.array([1./self.total_area])
            self.averagestress = 1./self.total_area
            self.ultimatestress = self.material.sigma_y
            self.tresca_yield = self.ultimatestress*self.ultimatestress / 3
        else:
            forcematrix = np.zeros((size, size))
            forcematrix[0,:] = 1
            forcematrix[1:,0] = self._a/(self.total_area*self.material.E)
            for idx, stringer in enumerate(self.stringers):
                forcematrix[idx+1, idx+1] = -stringer.properties.Le/(stringer.properties.total_area*stringer.properties.material.E)
            bvector = np.zeros(size)
            bvector[0] = 1.
            try:
                fractions = np.linalg.solve(forcematrix, bvector)
            except np.linalg.LinAlgError:
                print(forcematrix, bvector)
                print("Error in ult stress calculation in a panel, matrix is singular!")

            self.stresses = np.copy(fractions)
            self.stresses[0] = self.stresses[0]/(self.total_area)
            for idx, stringer in enumerate(self.stringers):
                self.stresses[idx+1] = self.stresses[idx+1]/(stringer.properties.total_area)
            self.averagestress = np.sum(fractions)/self.total_area
            ult_stresses = []
            for stringer, fraction in zip(self.stringers, fractions):
                ult_stresses.append(stringer.properties.material.sigma_y/fraction)
            ult_stresses = np.array([ult_stresses])
            self.ultimatestress = np.min(ult_stresses)
            self.tresca_yield = self.ultimatestress*self.ultimatestress / 3

    def get_mass(self) -> float:
        mass = self.ts * self._a * self.b * self.material.rho
        for stringer in self.stringers:
            mass += stringer.properties.mass
        return mass

    def calc_centroid(self):
        if self.total_area is None:
            self.calc_total_area()
        self.Qx = 0.5*self.ts * self.ts * self.b
        self.Qy = 0.5*self.b * self.ts * self.b
        for stringer in self.stringers:
            self.Qx += stringer.ybar_global*stringer.properties.total_area
            self.Qy += stringer.xbar_global*stringer.properties.total_area
        return self.Qy/self.total_area, self.Qx/self.total_area

    def momentofinertia(self):
        self.Ixxown = self.b*self.ts*self.ts*self.ts/12 + self.b*self.ts*(self.ts*0.5-self.ybar)**2
        self.Iyyown = self.ts*self.b*self.b*self.b/12 + self.b*self.ts*(self.b*0.5-self.xbar)**2
        self.Ixyown = self.b*self.ts*(self.ts*0.5-self.ybar)*(self.b*0.5-self.xbar)
        # Ixx, Iyy, Ixy = translate_mmoi(Ixx, Iyy, Ixy, self.rotation)
        Ixx = self.Ixxown
        Iyy = self.Iyyown
        Ixy = self.Ixyown
        for stringer in self.stringers:
            Ixx += stringer.Ixx_global + stringer.properties.total_area*(stringer.ybar_global-self.ybar)**2
            Iyy += stringer.Iyy_global + stringer.properties.total_area*(stringer.xbar_global-self.xbar)**2
            Ixy += stringer.Ixy_global + stringer.properties.total_area*(stringer.xbar_global-self.xbar)*(stringer.ybar_global-self.ybar)
        return Ixx, Iyy, Ixy

    def print_stiffening_effect(self):
        for we, stiff in zip(self.we2s, self.stringers):
            print(we, stiff)
            print("x={}, y={}, we={}, sig_cr_stiff={}".format(stiff.xpos, stiff.ypos, we, stiff.properties.sigma_cr))

    def idealize(self):
        self.area1 = self.area2 = self.area3 = self.total_area/3
        self.pos_i1 = (self.b, 0.0)
        self.pos_i2 = (self.b/2, 0.0)
        self.pos_i3 = (0.0, 0.0)
        for stiff in self.stringers:
            if stiff.xpos > self.xbar*1.5:
                self.area1 += stiff.properties.total_area
            elif stiff.xpos < self.xbar*0.5:
                self.area3 += stiff.properties.total_area
            else:
                self.area2 += stiff.properties.total_area


class Panel:
    def __init__(self, xpos, ypos, panel_instance, rotation=0):
        self.properties = copy.deepcopy(panel_instance)
        self.xpos = xpos
        self.ypos = ypos
        self.rotation = rotation
        # Centroid and mmoi with respect to reference coordinate system (using xpos, ypos, rotation)
        # x assumed positive to the left, y positive upwards
        self.xbar_global, self.ybar_global = self.properties.xbar*cos(rotation) - self.properties.ybar*sin(rotation), \
                                             self.properties.xbar*sin(rotation) + self.properties.ybar*cos(rotation)
        self.xbar_global += xpos
        self.ybar_global += ypos

        self.pos_i1_global = [self.properties.pos_i1[0]*cos(rotation) - self.properties.pos_i1[1]*sin(rotation), \
                             self.properties.pos_i1[0]*sin(rotation) + self.properties.pos_i1[1]*cos(rotation)]
        self.pos_i2_global = [self.properties.pos_i2[0]*cos(rotation) - self.properties.pos_i2[1]*sin(rotation), \
                             self.properties.pos_i2[0]*sin(rotation) + self.properties.pos_i2[1]*cos(rotation)]
        self.pos_i3_global = [self.properties.pos_i3[0]*cos(rotation) - self.properties.pos_i3[1]*sin(rotation), \
                             self.properties.pos_i3[0]*sin(rotation) + self.properties.pos_i3[1]*cos(rotation)]

        self.pos_i1_global[0], self.pos_i1_global[1] = self.pos_i1_global[0] + xpos, self.pos_i1_global[1] + ypos
        self.pos_i2_global[0], self.pos_i2_global[1] = self.pos_i2_global[0] + xpos, self.pos_i2_global[1] + ypos
        self.pos_i3_global[0], self.pos_i3_global[1] = self.pos_i3_global[0] + xpos, self.pos_i3_global[1] + ypos
        self.ipos = [self.pos_i1_global, self.pos_i2_global, self.pos_i3_global]

        self.Ixx_global, self.Iyy_global, self.Ixy_global = translate_mmoi(
            self.properties.Ixx, self.properties.Iyy, self.properties.Ixy, rotation)



def find_best_material():
    Le, a, b, ts = 0.3, 0.3, 0.6, 0.0009
    materials = {}
    materials['alu2024'] = MatProps(sigma_y=450000000, E=72400000000, poisson=0.33, rho=2.87, name="AA2024", alpha=0.8,
                                    n=0.6)  # http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA2024T81
    materials['alu5052'] = MatProps(sigma_y=255000000, E=72300000000, poisson=0.33, rho=2.68, name="AA5052", alpha=0.8,
                                    n=0.6)  # http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA5052H38
    materials['alu6061'] = MatProps(sigma_y=455000000, E=69000000000, poisson=0.33, rho=2.70, name="AA6061", alpha=0.8,
                                    n=0.6)  # http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA6061T913
    materials['alu6063'] = MatProps(sigma_y=295000000, E=69000000000, poisson=0.33, rho=2.70, name="AA6063", alpha=0.8,
                                    n=0.6)  # http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA6063T835
    materials['alu7050'] = MatProps(sigma_y=490000000, E=71700000000, poisson=0.33, rho=2.83, name="AA7050", alpha=0.8,
                                    n=0.6)  # http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA7050T765
    materials['alu7075'] = MatProps(sigma_y=503000000, E=71700000000, poisson=0.33, rho=2.81, name="AA7075", alpha=0.8,
                                    n=0.6)  # http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA7075T6
    materials['carbonfibre'] = MatProps(sigma_y=600000000, E=70000000000, poisson=0.1, rho=1.60, sigma_comp=570000000,
                                        name="carbonfibre")  # http://www.performance-composites.com/carbonfibre/mechanicalproperties_2.asp
    materials['glassfibre'] = MatProps(sigma_y=440000000, E=25000000000, poisson=0.2, rho=1.90, sigma_comp=425000000,
                                       name="glassfibre")  # http://www.performance-composites.com/carbonfibre/mechanicalproperties_2.asp

    functions = [J_Stringer, Z_Stringer]
    # functions = []

    l = [materials.values(), functions, materials.values()]

    best_ratio, best_names = 0, []

    for stringermat, stringertype, panelmat in itertools.product(*l):
        # print(stringermat, stringertype, panelmat)

        templatestringer = stringertype(Le=0.4, material=stringermat)
        stringers = make_stringers([(6, 0), (-2, 0), (3, 0)], templatestringer)

        panel = Sheet(a=a, b=b, ts=ts, material=panelmat)
        panel.construct_panel(*stringers)

        mass, sigma_cr = panel.mass, panel.sigma_cr
        ratio = sigma_cr / mass

        if ratio > best_ratio:
            best_ratio = ratio
            best_names = [stringermat.name, stringertype.name, panelmat.name]

        print("stiffmat={}, stringertype={}, panelmat={}".format(stringermat.name, stringertype.name, panelmat.name))
        print("mass={} kg, sig_cr={} MPa, sig_cr/mass={} MPa/kg\n".format(mass, round(sigma_cr / 1000000, 2),
                                                                          round(sigma_cr / mass, 1)))

    print("Best sig_cr/mass={} MPa/kg; {}".format(best_ratio, best_names))

if __name__ == "__main__":
    Le, a, b, ts = 0.5, 0.5, 0.2, 0.001
    rotation = 0
    angle = rotation*pi/180
    aluminium = MatProps(sigma_y=503000000, E=71700000000, poisson=0.33, rho=2.81, name="AA7075", alpha=0.8, n=0.6)
    # locs = [(0.0, 0.0), (0.1, 0.0), (0.2, 0.0), (0.3, 0.0)]
    x, y = 0.0, 0.0
    w_sheet = 0.2
    Ns = 0
    locs = [(x + w_sheet * cos(angle) * n / (Ns+1), y + w_sheet * sin(angle) * n / (Ns+1)) for n in range(1, (Ns+1))]
    # print(locs)
    stringers = make_stringers(locs, Z_Stringer(Le=Le, material=aluminium))
    panel = Sheet(a=a, b=b, ts=ts, material=aluminium)
    panel.construct_panel(*stringers)
    stiff_panel = Panel(0.0, 0.0, panel, rotation=rotation)
    # print(stiff_panel.properties.we2s)
    # print(stiff_panel.properties.we2s[0]*Ns)
    print("\n\nFinal critical stress=",stiff_panel.properties.sigma_cr)
