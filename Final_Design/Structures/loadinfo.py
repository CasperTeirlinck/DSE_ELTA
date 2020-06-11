import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
from functools import partial
from typing import Union
from matplotlib import pyplot as plt

try:
    from tqdm import tqdm
    # Install this module for a low-impact progress bar!
    loadbar = lambda x: tqdm(x)
except:
    loadbar = lambda x: x


def step_function_twoarg(scalar, target):
    if target > scalar:
        return target-scalar
    else:
        return 0.0


def step_function_moment_twoarg(scalar, target):
    if target > scalar:
        return 1.0
    else:
        return 0.0


def step_function(scalar):
    return partial(step_function_twoarg, scalar)

def step_function_moment(scalar):
    return partial(step_function_moment_twoarg, scalar)


class Force:
    def __init__(self, xpos, ypos, zpos, xmag=0., ymag=0., zmag=0.):
        self.position = np.array([xpos, ypos, zpos])
        self.vector = np.array([xmag, ymag, zmag])
        self.magnitude = np.abs(self.vector)


class Moment(Force):
    def __init__(self, xpos, ypos, zpos, xmag=0., ymag=0., zmag=0.):
        super().__init__(xpos=xpos, ypos=ypos, zpos=zpos, xmag=xmag, ymag=ymag, zmag=zmag)
        self.magnitude = None


class LoadCase:
    def __init__(self, y_coordinates=None, EIx=None, EIz=None):

        self.forces = []
        self.moments = []

        self.reactionforces = []
        self.reactionmoments = []

        self.totalforces = [self.forces, self.reactionforces]
        self.totalmoments = [self.moments, self.reactionmoments]

        self.C1 = None
        self.C2 = None

        self.boundary_conditions = []

        self.y_samples = y_coordinates
        self.EIx = EIx
        self.EIz = EIz

    def change_geometry(self, y_coordinates, EIx, EIz):
        self.y_samples = y_coordinates
        self.EIx = EIx
        self.EIz = EIz

    def add_boundary_condition(self, y=0., defl_x=None, defl_z=None, angle_x=None, angle_z=None):
        constraints = (defl_x != None) + (defl_z != None) + (angle_x != None) + (angle_z != None)
        if constraints > 1:
            print("Only one constraint allowed per boundary condition, add them separately!")
            return
        if constraints == 0:
            print("Not a valid boundary condition!")
        bc = {"y": y, "defl_x": defl_x, "defl_z": defl_z, "angle_x": angle_x, "angle_z": angle_z}
        self.boundary_conditions.append(bc)
        return

    def remove_all_boundary_conditions(self):
        self.boundary_conditions = []

    def add_loads(self, *loads: Union[Force, Moment]):
        for load in loads:
            if type(load) == Force:
                self.forces.append(load)
            else:
                self.moments.append(load)

    def calc_reactions(self):
        reactionforcex, reactionforcez = 0.0, 0.0
        for force in self.forces:
            reactionforcex -= force.vector[0]
            reactionforcez -= force.vector[2]
        reactionmomentx, reactionmomentz = 0.0, 0.0
        for moment in self.moments:
            reactionmomentx -= moment.vector[0]
            reactionmomentz -= moment.vector[2]
        for bc in self.boundary_conditions:
            if bc["angle_x"] != None:
                self.reactionmoments.append(Moment(xpos=0.0, ypos=bc["y"], zpos=0.0, zmag=reactionmomentz))
            if bc["defl_x"] != None:
                self.reactionforces.append(Force(xpos=0.0, ypos=bc["y"], zpos=0.0, xmag=reactionforcex))
            if bc["angle_z"] != None:
                self.reactionmoments.append(Moment(xpos=0.0, ypos=bc["y"], zpos=0.0, xmag=reactionmomentx))
            if bc["defl_z"] != None:
                self.reactionforces.append(Force(xpos=0.0, ypos=bc["y"], zpos=0.0, zmag=reactionforcez))

    def do_bc_checks(self):
        if len(self.boundary_conditions) != 4:
            print("Too many or too few boundary conditions. Required amount of b.c. is 4, current amount is {}".format(
                len(self.boundary_conditions)))
            return
        xs, zs = 0, 0
        for bc in self.boundary_conditions:
            if bc["defl_x"] != None or bc["angle_x"] != None:
                xs += 1
            elif bc["defl_z"] != None or bc["angle_z"] != None:
                zs += 1
        if xs != 2 or zs != 2:
            print("Over- or underconstrained system! Required boundary conditions in every direction is 2, currently x={}, z={}".format(xs, zs))
            return

        if self.EIx is None or self.EIz is None:
            print("No geometry specified, use the 'change_geometry()' method or add during initialization!")
            return


    def solve_bcs(self, printing=False):
        self.calc_reactions()
        self.do_bc_checks()

        matrix = np.zeros((4,4))
        bvect = np.zeros(4)
        for idx, bc in enumerate(self.boundary_conditions):
            angles = self.angle_bare(y=bc["y"])
            deflections = self.deflection_bare(y=bc["y"])
            if bc["angle_x"] != None:
                matrix[0,:] = np.array([1, 0, 0, 0])
                bvect[0] = bc["angle_x"] - angles[0]
            if bc["defl_x"] != None:
                matrix[1,:] = np.array([bc["y"], 1, 0, 0])
                bvect[1] = bc["defl_x"] - deflections[0]
            if bc["angle_z"] != None:
                matrix[2,:] = np.array([0, 0, 1, 0])
                bvect[2] = bc["angle_z"] - angles[1]
            if bc["defl_z"] != None:
                matrix[3,:] = np.array([0, 0, bc["y"], 1])
                bvect[3] = bc["defl_z"] - deflections[1]
        c1c2 = np.linalg.solve(matrix, bvect)
        self.C1 = np.array([c1c2[0],c1c2[2]])
        self.C2 = np.array([c1c2[1],c1c2[3]])
        return

    def calc_coefficient(self, EI, y, y_application, numerator, degree_of_integration=1):
        if np.all(EI==EI[0]):
            interpolant = lambda val: numerator(val)
            multiplier = 1/(EI[0]**degree_of_integration)
        else:
            interpolant = lambda val: numerator(val)/(interp1d(self.y_samples, EI)(val))
            multiplier = 1.0
        if degree_of_integration == 2:
            return multiplier * quad(lambda t: quad(func=interpolant, a=y_application, b=t)[0], a=y_application, b=y)[0]
        else:
            return multiplier * quad(func=interpolant, a=y_application, b=y)[0]

    def angle_bare(self, y):
        # point = self.origin_y - y
        angle = np.array([0., 0.])
        tempforces = [item for sublist in self.totalforces for item in sublist]
        tempmoments = [item for sublist in self.totalmoments for item in sublist]
        for force in tempforces:
            coefficient = step_function(force.position[1])
            angle[0] += force.vector[0] * self.calc_coefficient(self.EIz, y, force.position[1], coefficient)
            angle[1] += force.vector[2] * self.calc_coefficient(self.EIx, y, force.position[1], coefficient)
        for moment in tempmoments:
            coefficient = step_function_moment(moment.position[1])
            angle[0] += moment.vector[2] * self.calc_coefficient(self.EIz, y, moment.position[1], coefficient)
            angle[1] += moment.vector[0] * self.calc_coefficient(self.EIx, y, moment.position[1], coefficient)
        del tempforces, tempmoments, coefficient
        return angle

    def calc_angle(self, y_in):
        try:
            angles = np.zeros((len(y_in),2))
            for idx, y in loadbar(enumerate(y_in)):
                angles[idx,:] = self.angle_bare(y=y) + self.C1
        except TypeError:
            angles = self.angle_bare(y=y_in) + self.C1
        return angles

    def deflection_bare(self, y):
        # point = self.origin_y - y
        deflection = np.array([0., 0.])
        tempforces = [item for sublist in self.totalforces for item in sublist]
        tempmoments = [item for sublist in self.totalmoments for item in sublist]
        for force in tempforces:
            coefficient = step_function(force.position[1])
            deflection[0] += force.vector[0] * self.calc_coefficient(self.EIz, y, force.position[1], coefficient, 2)
            deflection[1] += force.vector[2] * self.calc_coefficient(self.EIx, y, force.position[1], coefficient, 2)
        for moment in tempmoments:
            coefficient = step_function_moment(moment.position[1])
            deflection[0] += moment.vector[2] * self.calc_coefficient(self.EIz, y, moment.position[1], coefficient, 2)
            deflection[1] += moment.vector[0] * self.calc_coefficient(self.EIx, y, moment.position[1], coefficient, 2)
        del tempforces, tempmoments, coefficient
        return deflection

    def calc_deflection(self, y_in):
        try:
            deflections = np.zeros((len(y_in),2))
            for idx, y in loadbar(enumerate(y_in)):
                deflections[idx,:] = self.deflection_bare(y=y) + self.C2 + self.C1*y
        except TypeError:
            deflections = self.deflection_bare(y=y_in) + self.C2 + self.C1*y_in
        return deflections


def run_tests():#
    EIz = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])*2
    ys = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    system = LoadCase(y_coordinates=ys, EIx=EIz, EIz=EIz)

    loads = {}
    loads['f1'] = Force(xpos=0.0, ypos=1.0, zpos=0.0, xmag=1, ymag=0.0, zmag=0.0)
    # loads['f2'] = Force(xpos=0.0, ypos=0.5, zpos=0.0, xmag=-1000.0, ymag=0.0, zmag=0.0)
    loads['m1'] = Moment(xpos=0.0, ypos=1.0, zpos=0.0, xmag=0.0, ymag=0.0, zmag=-1.0)
    system.add_loads(*loads.values())

    system.add_boundary_condition(y=0.0, defl_x=0)
    system.add_boundary_condition(y=0.0, defl_z=0)
    system.add_boundary_condition(y=0.0, angle_x=0)
    system.add_boundary_condition(y=0.0, angle_z=0)
    system.solve_bcs()

    totalforce, totalmoment = np.array([0.0,0.0,0.0]), np.array([0.0,0.0,0.0])
    tempforces = [item for sublist in system.totalforces for item in sublist]
    tempmoments = [item for sublist in system.totalmoments for item in sublist]
    for force in tempforces:
        totalforce += force.vector
    for moment in tempmoments:
        totalmoment += moment.vector
    assert(np.all(totalforce==0.0))
    assert(np.all(totalmoment==0.0))

    def exact_deflection(yout):
        return yout*yout*(3.-yout)/(6.*2.) + yout*yout/(2.*2.)

    def exact_angle(length):
        return length*length/(2.*2.) - length / 2.

    yout = np.linspace(0,1.0,51)
    defl_exact = exact_deflection(yout)
    defl_out = system.calc_deflection(yout)
    angle_out = system.calc_angle(yout)
    angle_exact = exact_angle(np.max(yout))
    difference = defl_exact - defl_out[:, 0]
    print(defl_exact[-1])
    print(defl_out[-1, 0])
    print(difference[-1])
    plt.plot(yout, difference)
    plt.plot(yout, defl_exact)
    plt.plot(yout, defl_out[:, 0])
    plt.legend(['difference', 'exact', 'calculated'])
    plt.show()

    print("All tests passed.")



if __name__ == "__main__":
    run_tests()
    direction=0  # 0 for x, 1 for z

    # EIx = np.array([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0])
    EIz = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])*2
    ys = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

    # f1 = Force(xpos=0.0, ypos=0.1, zpos=0.0, xmag=2.0, ymag=0., zmag=1.0)
    # m1 = Moment(xpos=0.0, ypos=0.3, zpos=0.0, xmag=-1.0, ymag=0., zmag=10.0)

    loads = {}
    # f1 = Force(xpos=0.0, ypos=1.0, zpos=0.0, xmag=1000, ymag=0.0, zmag=0.0)
    # m1 = Moment(xpos=0.0, ypos=1.0, zpos=0.0, xmag=0.0, ymag=0.0, zmag=1.0)
    loads['f1'] = Force(xpos=0.0, ypos=1.0, zpos=0.0, xmag=1, ymag=0.0, zmag=0.0)
    loads['f2'] = Force(xpos=0.0, ypos=0.5, zpos=0.0, xmag=-1000.0, ymag=0.0, zmag=0.0)
    loads['m1'] = Moment(xpos=0.0, ypos=1.0, zpos=0.0, xmag=0.0, ymag=0.0, zmag=-1.0)

    system = LoadCase(y_coordinates=ys, EIx=EIz, EIz=EIz)

    system.add_loads(*loads.values())

    system.add_boundary_condition(y=0.05, defl_x=0)
    system.add_boundary_condition(y=0.05, defl_z=0)
    system.add_boundary_condition(y=0.05, angle_x=0)
    system.add_boundary_condition(y=0.05, angle_z=0)
    system.solve_bcs(printing=True)

    print(system.C1[direction])
    print(system.C2[direction])

    yout = np.linspace(0,1.0,51)
    # yout = 0.05
    defl_out = system.calc_deflection(yout)
    angle_out = system.calc_angle(yout)

    # print(angle_out)

    # plt.plot(yout, defl_out[:,direction])
    # plt.plot(yout, angle_out[:,direction])
    # plt.legend(['Deflection', 'Slope'], loc='best')
    # plt.show()







