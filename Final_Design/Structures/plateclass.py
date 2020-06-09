import numpy as np
from math import pi, sqrt
import itertools
from matplotlib import pyplot as plt
from stringerclass import *


def get_C(bovert):
    if bovert < 40:
        return 4.0
    elif bovert > 110:
        return 6.98
    else:
        return 4.0 + 2.98 * (bovert - 40) / 70


def make_stringers(positions, stringerinstance):
    stringers = []
    for pos in list(positions):
        stringers.append(Stringer(*pos, stringerinstance))
    stringers = sorted(stringers, key=lambda x: x.xpos) # Sort stringers by x coordinate
    return stringers



def make_panel(panel, *args):
    for arg in args:
        panel.addStringer(arg)
    panel.recalculate_stringer_effects()
    return None


class Sheet:
    def __init__(self, a, b, ts, material: MatProps, mode=1):
        self.material = material
        self.mode = mode

        self.ts = ts # Thickness of the panel
        self.a = a # Length of the panel
        self.b = b # Width of the panel

        self.stringers = []
        self.we2s = []
        self.stringerpitches = []

        self.total_area = None # Total area of the stiffened plate cross-section
        self.calc_total_area()

        self.sigma_cr = self.calc_critical_stress(mode=self.mode)
        self.sigma_ult_tensile = self.calc_sigma_ult_tensile()

        self.stresses = None # Stresses on plate and subsequently the stringers for a unit tensile load
        self.averagestress = None # Average stress on stiffened panel for a unit tensile load
        self.ultimatestress = None # Ultimate tensile stress that the panel can take before yielding
        self.calc_sigma_ult_tensile()

    def addStringer(self, stringer: Stringer):
        self.stringers.append(stringer)

    def recalculate_stringer_effects(self):
        self.calc_we()
        self.calc_stringerpitch()
        self.calc_critical_stress(mode=self.mode)
        self.calc_total_area()
        self.calc_sigma_ult_tensile()

    def calc_we(self):
        for stringer in self.stringers:
            C = get_C(self.b/self.a)
            self.we2s.append(self.ts * sqrt(C*pi*pi*self.material.E/(12*(1-self.material.poisson**2)*stringer.properties.sigma_cc)))

    def calc_stringerpitch(self):
        for index, stringer in enumerate(self.stringers):
            if index == 0:
                stringerpitch = stringers[index+1].xpos - stringer.xpos
            elif index == len(self.stringers)-1:
                stringerpitch = stringer.xpos-stringers[index-1].xpos
            else:
                stringerpitch = (stringers[index+1].xpos - stringers[index-1].xpos)/2
            self.stringerpitches.append(stringerpitch)

    def calc_critical_stress(self, mode):
        # mode determines the type of support of the panel. mode 1 is CCSS, 2 is SSSS,
        # 3 is CCCC, 4 is SSCC, 5 is SSCS, 6 is CCFF, 7 is SSFC, 8 is SSFF, 9 is SSFS
        # If required mode is anything else, the program asks for manual input. Look at
        # https://puu.sh/FSyYi/82a36e2870.png to find the correct value for C.
        aoverb = self.a/self.b
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
        for we, stiff in zip(self.we2s, self.stringers):
            stiffenerforces += stiff.properties.sigma_cr*(stiff.properties.total_area+we*self.ts)
            totalstiffenerarea += stiff.properties.total_area
        return (self.sigma_cr_skin * self.ts * (self.b-sum(self.we2s)) + stiffenerforces) / (totalstiffenerarea+ self.ts + self.b)

    def calc_total_area(self):
        self.total_area = self.ts * self.b
        for stringer in stringers:
            self.total_area += stringer.properties.total_area

    def calc_sigma_ult_tensile(self):
        size = len(stringers)+1
        forcematrix = np.zeros((size, size))
        forcematrix[0,:] = 1
        forcematrix[1:,0] = self.a/(self.ts*self.b*self.material.E)
        for idx, stringer in enumerate(stringers):
            forcematrix[idx+1, idx+1] = -stringer.properties.Le/(stringer.properties.total_area*stringer.properties.material.E)
        # print(forcematrix)
        bvector = np.zeros(size)
        bvector[0] = 1.
        # print(bvector)
        fractions = np.linalg.solve(forcematrix, bvector)
        # print(fractions)
        self.stresses = np.copy(fractions)
        self.stresses[0] = self.stresses[0]/(self.ts*self.b)
        for idx, stringer in enumerate(stringers):
            self.stresses[idx+1] = self.stresses[idx+1]/(stringer.properties.total_area)
        self.averagestress = np.sum(fractions)/self.total_area
        ult_stresses = []
        for stringer, fraction in zip(self.stringers, fractions):
            # print("asdf", fraction)
            ult_stresses.append(stringer.properties.material.sigma_y/fraction)
        # print("This", ult_stresses)
        # self.ultimatestress = min([stringer.properties.material.sigma_y/fraction for stringer, fraction in zip(self.stringers, fractions)])
        # print(self.ultimatestress)

    def get_mass(self) -> float:
        mass = self.ts * self.a * self.b * self.material.rho
        for stringer in self.stringers:
            mass += stringer.properties.total_area*stringer.properties.Le * stringer.properties.material.rho
        return mass

    def print_stiffening_effect(self):
        for we, stiff in zip(self.we2s, self.stringers):
            print(we, stiff)
            print("x={}, y={}, we={}, sig_cr_stiff={}".format(stiff.xpos, stiff.ypos, we, stiff.properties.sigma_cr))


if __name__ == "__main__":

    Le, a, b, ts = 0.3, 0.3, 0.6, 0.0009
    materials = {}
    materials['alu2024'] = MatProps(sigma_y=450000000, E=72400000000, poisson=0.33, rho=2.87, name="AA2024", alpha=0.8, n=0.6) # http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA2024T81
    materials['alu5052'] = MatProps(sigma_y=255000000, E=72300000000, poisson=0.33, rho=2.68, name="AA5052", alpha=0.8, n=0.6) # http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA5052H38
    materials['alu6061'] = MatProps(sigma_y=455000000, E=69000000000, poisson=0.33, rho=2.70, name="AA6061", alpha=0.8, n=0.6) # http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA6061T913
    materials['alu6063'] = MatProps(sigma_y=295000000, E=69000000000, poisson=0.33, rho=2.70, name="AA6063", alpha=0.8, n=0.6) # http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA6063T835
    materials['alu7050'] = MatProps(sigma_y=490000000, E=71700000000, poisson=0.33, rho=2.83, name="AA7050", alpha=0.8, n=0.6) # http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA7050T765
    materials['alu7075'] = MatProps(sigma_y=503000000, E=71700000000, poisson=0.33, rho=2.81, name="AA7075", alpha=0.8, n=0.6) # http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA7075T6
    materials['carbonfibre'] = MatProps(sigma_y=600000000, E=70000000000, poisson=0.1, rho=1.60, sigma_comp=570000000, name="carbonfibre") # http://www.performance-composites.com/carbonfibre/mechanicalproperties_2.asp
    materials['glassfibre'] = MatProps(sigma_y=440000000, E=25000000000, poisson=0.2, rho=1.90, sigma_comp=425000000, name="glassfibre") # http://www.performance-composites.com/carbonfibre/mechanicalproperties_2.asp

    functions = [J_Stringer, Z_Stringer]
    # functions = []

    l = [materials.values(), functions, materials.values()]

    best_ratio, best_names = 0, []

    for stringermat, stringertype, panelmat in itertools.product(*l):
        # print(stringermat, stringertype, panelmat)
        templatestringer = stringertype(Le=0.4, material=stringermat)
        stringers = make_stringers([(6,0),(-2,0),(3,0)], templatestringer)
        panel = Sheet(a=a, b=b, ts = ts, material=panelmat)
        make_panel(panel, *stringers)

        mass, sigma_cr = panel.get_mass(), panel.sigma_cr
        ratio = sigma_cr/mass

        if ratio > best_ratio:
            best_ratio = ratio
            best_names = [stringermat.name, stringertype.name, panelmat.name]

        print("stiffmat={}, stringertype={}, panelmat={}".format(stringermat.name, stringertype.name, panelmat.name))
        print("mass={} kg, sig_cr={} MPa, sig_cr/mass={} MPa/kg\n".format(mass, round(sigma_cr/1000000, 2), round(sigma_cr/mass, 1)))

    print("Best sig_cr/mass={} MPa/kg; {}".format(best_ratio, best_names))

    # panel.print_stiffening_effect()
    # print(panel.stringers[0].properties.Ixx)
