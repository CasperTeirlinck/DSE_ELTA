

from math import pi, sqrt
import copy


class MatProps():
    def __init__(self, sigma_y, E, poisson, alpha=0.8, n=0.6):
        self.sigma_y = sigma_y
        self.E = E
        self.poisson = poisson
        self.alpha = alpha
        self.n = n


class J_Stringer():
    def __init__(self, Le, material: MatProps, stiff_ratio=0.5, ts=0.001, t1=0.001, t2=0.001, t3=0.001, t4=0.001,
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
        return material.alpha * (C/material.sigma_y * pi*pi*material.E*t*t/(b*b*12*(1-material.poisson)))**(1-material.n) * material.sigma_y

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

        if sigma_cr < self.material.sigma_y:
            return sigma_cr
        else:
            return self.material.sigma_y


class Z_Stringer(J_Stringer):
    # def __init__(self, Le, yieldstress, E, poisson, alpha=0.8, n=0.6, stiff_ratio=0.5, ts=0.001,
    #              t1=0.001, t2=0.001, t3=0.001, t4=0.001, b1=0.005, b2=0.006, b3=0.005, b4=0.005):
    def __init__(self, Le, material, stiff_ratio=0.5, ts=0.001, t1=0.001, t2=0.001, t4=0.001,
                 b1=0.005, b2=0.006, b4=0.005):
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


class Stringer:
    def __init__(self, xpos, ypos, stringer_instance):
        self.properties = copy.deepcopy(stringer_instance)
        self.xpos = xpos
        self.ypos = ypos


if __name__ == "__main__":
    aluminium = MatProps(sigma_y=450000000, E=72000000000, poisson=0.3, alpha=0.8, n=0.6)
    j = J_Stringer(Le=0.4, material=aluminium)
    print("\nZ_stringer\n")
    z = Z_Stringer(Le=0.4, material=aluminium)


