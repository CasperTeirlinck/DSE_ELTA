from math import pi, sqrt, sin, cos
import copy

class MatProps():
    def __init__(self, sigma_y, E, poisson, rho, cost, sigma_comp=None, name=None, alpha=0.8, n=0.6):
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
        self.cost = cost


class J_Stringer():
    name = "J_Stringer"
    def __init__(self, Le, material: MatProps, t1=0.001, t2=0.001, t3=0.001, t4=0.001,
                 b1=0.005, b2=0.006, b3=0.005, b4=0.005):
        # Geometry defined following https://puu.sh/FSdDm/df0747d2c5.png, gray area included in 2_alt
        # Reference coordinate system is x positive to the left, y positive upwards, origin bottom right.
        self._Le = Le
        self.material = material

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
        self.total_area = self.calc_total_area()

        self.Qx = None
        self.Qy = None

        self.xbar, self.ybar = self.calc_centroid()
        self.Ixx, self.Iyy, self.Ixy = self.momentofinertia()

        self.sigma_cc = self.crippling_stress()

        self.sigma_cr = self.calc_critical_stress(self._Le)

        self.mass = self.calc_mass()

        self.cost = self.calc_cost()

    @property
    def Le(self):
        return self._Le

    @Le.setter
    def Le(self, val):
        self._Le = val
        self.sigma_cr = self.calc_critical_stress(self._Le)
        self.mass = self.calc_mass()
        self.cost = self.calc_cost()

    def calc_mass(self):
        return self._Le*self.total_area*self.material.rho

    def calc_cost(self):
        return self.mass*self.material.cost

    def calc_total_area(self):
        return self.area1 + self.area2_alt + self.area3 + self.area4

    def calc_centroid(self):
        # Returns xbar, ybar
        self.Qx = self.area1 * (self.b2+self.t1*0.5+self.t3) + self.area2_alt * (self.t1+self.b2+self.t4)*0.5 + \
             self.area4 * self.t4*0.5 + self.area3 * self.t3*0.5
        self.Qy = self.area1 * (self.b4 + self.t2 + self.b1*0.5) + self.area2_alt * (self.b4 + self.t2*0.5) + \
             self.area4 * self.b4*0.5 + self.area3 * (self.b4 + self.t2 + self.b3*0.5)
        return self.Qy/self.total_area, self.Qx/self.total_area

    def crippling_partial(self, t, b, C, material: MatProps):
        return material.alpha * (C/material.sigma_comp * pi*pi*material.E*t*t/(b*b*12*(1-material.poisson)))**(1-material.n) * material.sigma_comp

    def crippling_stress(self):
        return (self.crippling_partial(self.t1, self.b1, 0.425, self.material) * self.area1 +
                self.crippling_partial(self.t2, self.b2, 0.400, self.material) * self.area2 +
                self.crippling_partial(self.t4, self.b4, 0.425, self.material) * self.area4 +
                self.crippling_partial(self.t3, self.b3, 0.425, self.material) * self.area3) / \
               (self.total_area-self.area2_alt+self.area2)

    def momentofinertia(self):
        Ixx = (self.b1*self.t1**3+self.t2*(self.b2+self.t1+self.t4)**3+self.b4*self.t4**3+self.b3*self.t3**3)/12 + \
              self.area1*(self.b2+self.t4+self.t1*0.5-self.ybar)**2 + \
              self.area2_alt * (self.area2_alt/(2*self.t2)-self.ybar)**2 + \
              self.area4*(self.t4*0.5-self.ybar)**2 + self.area3*(self.t3*0.5-self.ybar)**2
        Iyy = (self.t1*self.b1**3+(self.b2+self.t1+self.t4)*self.t2**3+self.t4*self.b4**3+self.t3*self.b3**3) + \
              self.area1*(self.b4+self.t2+self.b1*0.5-self.xbar)**2 + self.area2*(self.b4+self.t2*0.5-self.xbar)**2 + \
              self.area4*(self.b4*0.5-self.xbar)**2 + self.area3*(self.b4+self.t2+self.b3*0.5-self.xbar)**2
        Ixy = self.area1*(self.b2+self.t4+self.t1*0.5-self.ybar)*(self.b4+self.t2+self.b1*0.5-self.xbar) + \
              self.area2_alt * (self.area2_alt/(2*self.t2)-self.ybar)*(self.b4+self.t2*0.5-self.xbar) + \
              self.area4 * (self.t4*0.5-self.ybar)*(self.b4*0.5-self.xbar) + \
              self.area3 * (self.t3*0.5-self.ybar)*(self.b4+self.t2+self.b3*0.5-self.xbar)
        return Ixx, Iyy, Ixy

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
        super().__init__(Le=Le, material=material, t1=t1, t2=t2, t4=t4, b1=b1, b2=b2, b4=b4)

    def calc_total_area(self):
        return self.area1 + self.area2_alt + self.area4

    def calc_centroid(self):
        # Returns xbar, ybar
        Qx = self.area1 * (self.b2+self.t1*0.5+self.t3) + self.area2_alt * (self.t1+self.b2+self.t4)*0.5 + \
             self.area4 * self.t4*0.5
        Qy = self.area1 * (self.b4 + self.t2 + self.b1*0.5) + self.area2_alt * (self.b4 + self.t2*0.5) + \
             self.area4 * self.b4*0.5
        return Qy/self.total_area, Qx/self.total_area

    def crippling_stress(self):
        return (self.crippling_partial(self.t1, self.b1, 0.425, self.material) * self.area1 +
                self.crippling_partial(self.t2, self.b2, 0.400, self.material) * self.area2 +
                self.crippling_partial(self.t4, self.b4, 0.425, self.material) * self.area4) / \
               (self.total_area-self.area2_alt+self.area2)

    def momentofinertia(self):
        Ixx = (self.b1*self.t1**3+self.t2*(self.b2+self.t1+self.t4)**3+self.b4*self.t4**3)/12 + \
              self.area1*(self.b2+self.t4+self.t1*0.5-self.ybar)**2 + \
              self.area2_alt * (self.area2_alt/(2*self.t2)-self.ybar)**2 + \
              self.area4*(self.t4*0.5-self.ybar)**2
        Iyy = (self.t1*self.b1**3+(self.b2+self.t1+self.t4)*self.t2**3+self.t4+self.b4**3)/12 + \
              self.area1*(self.b4+self.t2+self.b1*0.5-self.xbar)**2 + self.area2*(self.b4+self.t2*0.5-self.xbar)**2 + \
              self.area4*(self.b4*0.5-self.xbar)**2
        Ixy = self.area1*(self.b2+self.t4+self.t1*0.5-self.ybar)*(self.b4+self.t2+self.b1*0.5-self.xbar) + \
              self.area2_alt * (self.area2_alt/(2*self.t2)-self.ybar)*(self.b4+self.t2*0.5-self.xbar) + \
              self.area4 * (self.t4*0.5-self.ybar)*(self.b4*0.5-self.xbar)
        return Ixx, Iyy, Ixy


def translate_mmoi(ixx, iyy, ixy, rotation):
    # Rotation in radians
    iuu = 0.5*(ixx + iyy) + 0.5*(ixx - iyy)*cos(2*rotation) - ixy*sin(2*rotation)
    ivv = 0.5*(ixx + iyy) - 0.5*(ixx - iyy)*cos(2*rotation) + ixy*sin(2*rotation)
    iuv = 0.5*(ixx - iyy)*sin(2*rotation) + ixy*cos(2*rotation)
    return iuu, ivv, iuv


class Stringer:
    def __init__(self, xpos, ypos, stringer_instance, rotation=0):
        # Rotation in radians, clockwise positive. Keep 0 for correct orientation for bottom plate.
        # xpos and ypos is where the bottom right corner of the stringer is placed.
        self.properties = copy.deepcopy(stringer_instance)
        self.xpos = xpos
        self.ypos = ypos
        self.rotation = rotation
        # Centroid and mmoi with respect to reference coordinate system (using xpos, ypos, rotation)
        # origin of a panel is bottom right of plate, x positive to the left, y positive upwards
        self.xbar_global, self.ybar_global = self.properties.xbar*cos(rotation) - self.properties.ybar*sin(rotation), \
                                             self.properties.xbar*sin(rotation) + self.properties.ybar*cos(rotation)
        self.xbar_global += xpos
        self.ybar_global += ypos

        self.Ixx_global, self.Iyy_global, self.Ixy_global = translate_mmoi(
            self.properties.Ixx, self.properties.Iyy, self.properties.Ixy, rotation)

if __name__ == "__main__":
    print(translate_mmoi(1,2,0,0*pi/180))
    aluminium = MatProps(sigma_y=450000000, E=72400000000, poisson=0.33, rho=2.87, name="AA2024", alpha=0.8, n=0.6, cost=1.6)

    kwargs = {'t1':0.001, 't2':0.001, 't3':0.001, 't4':0.001, 'b1':0.01, 'b2':0.02, 'b3':0.01, 'b4':0.01}


    def constructstringerarg(modifier, material):
        return {'material': material, 't1': 0.0005 * modifier, 't2': 0.0005 * modifier, 't3': 0.0005 * modifier,
                't4': 0.0005 * modifier, 'b1': 0.005 * modifier, 'b2': 0.015 * modifier, 'b3': 0.003 * modifier,
                'b4': 0.003 * modifier}
    kwargs = constructstringerarg(1.0, aluminium)

    j = J_Stringer(Le=0.4, **kwargs)
    z = Z_Stringer(Le=0.4, **kwargs)
    print(z.total_area)
    # stringer = Stringer(0.0, 1.0, j)

    print(j.total_area)
    # print(stringer.properties.total_area)
    # print(stringer.properties.material.sigma_y)
    #
    # factor = sqrt(1.8/2.3)
    # for key in kwargs:
    #     kwargs[key] *= factor
    # j = J_Stringer(Le=0.4, material=aluminium, **kwargs)
    # s = Stringer(0.1, 4.5, j, 0.5*pi)


