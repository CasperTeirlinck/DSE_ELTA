import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
from functools import partial
from typing import Union
from matplotlib import pyplot as plt
from time import *
import cProfile

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
    return np.vectorize(partial(step_function_twoarg, scalar))

def step_function_moment(scalar):
    return np.vectorize(partial(step_function_moment_twoarg, scalar))


def interpolator(EIfunc, scalar, moment, target):
    if target < scalar:
        return 0.0
    else:
        if moment:
            return 1.0/EIfunc(target)
        else:
            return (target-scalar)/EIfunc(target)

def interpolant(EIfunc, scalar, moment):
    return partial(interpolator, EIfunc, scalar, moment)


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
    def __init__(self, y_coordinates=None, EIx=None, EIz=None, areas=None, shearcentre_z=None, shearcentre_x=None):

        self.forces = []
        self.moments = []

        self.reactionforces = []
        self.reactionmoments = []

        self.C1 = None
        self.C2 = None

        self.boundary_conditions = []

        self.y_samples = y_coordinates
        self.EIx = EIx
        self.EIz = EIz
        self.EIxfunc = interp1d(self.y_samples, EIx)
        self.EIzfunc = interp1d(self.y_samples, EIz)

        self.areas = areas
        # self.shearcentre_x = shearcentre_x
        self.shearcentre_z = shearcentre_z
        self.areasfunc = interp1d(self.y_samples, areas)
        # self.shearcentre_x_func = interp1d(self.y_samples, shearcentre_x)
        self.shearcentre_z_func = interp1d(self.y_samples, shearcentre_z)

    def change_geometry(self, y_coordinates, EIx, EIz):
        self.y_samples = y_coordinates
        self.EIx = EIx
        self.EIz = EIz
        self.EIxfunc = interp1d(self.y_samples, EIx)
        self.EIzfunc = interp1d(self.y_samples, EIz)

    def change_areas(self, areas, shearcentre_z, y_locs=None):
        if y_locs is None:
            y_locs = self.y_samples
        self.areas = areas
        # self.shearcentre_x = shearcentre_x
        self.shearcentre_z = shearcentre_z
        self.areasfunc = interp1d(y_locs, areas)
        # self.shearcentre_x_func = interp1d(y_locs, shearcentre_x)
        self.shearcentre_z_func = interp1d(y_locs, shearcentre_z)

    def add_boundary_condition(self, x=0., y=0., z=0., defl_x=None, defl_y=None, defl_z=None, angle_x=None, angle_y=None, angle_z=None):
        constraints = (defl_x != None) + (defl_y != None) + (defl_z != None) + (angle_x != None) + (angle_y != None) + (angle_z != None)
        if constraints > 1:
            print("Only one constraint allowed per boundary condition, add them separately!")
            return
        if constraints == 0:
            print("Not a valid boundary condition!")
        bc = {"x": x, "y": y, "z": z, "defl_x": defl_x, "defl_y": defl_y, "defl_z": defl_z, "angle_x": angle_x, "angle_y": angle_y, "angle_z": angle_z}
        self.boundary_conditions.append(bc)
        return

    def remove_all_boundary_conditions(self):
        self.boundary_conditions = []

    def remove_all_loads(self):
        self.forces, self.moments = [], []
        self.reactionforces, self.reactionmoments = [], []

    def add_loads(self, *loads: Union[Force, Moment]):
        for load in loads:
            if type(load) == Force:
                self.forces.append(load)
            else:
                self.moments.append(load)

    def calc_resultant_force(self):
        reactionforcex, reactionforcey, reactionforcez = 0.0, 0.0, 0.0
        for force in self.forces:
            reactionforcex += force.vector[0]
            reactionforcey += force.vector[1]
            reactionforcez += force.vector[2]
        return reactionforcex, reactionforcey, reactionforcez

    def calc_resultant_moment(self, origin=np.array([0., 0., 0.])):
        reactionmomentx, reactionmomenty, reactionmomentz = 0.0, 0.0, 0.0
        for moment in self.moments:
            reactionmomentx += moment.vector[0]
            reactionmomenty += moment.vector[1]
            reactionmomentz += moment.vector[2]
        vect_r_m = np.array([reactionmomentx, reactionmomenty, reactionmomentz])
        for force in self.forces:
            vect_r_m += np.cross(force.position-origin, force.vector)
        reactionmomentx, reactionmomenty, reactionmomentz = vect_r_m[0], vect_r_m[1], vect_r_m[2]
        return reactionmomentx, reactionmomenty, reactionmomentz

    def check_equilibrium(self):
        totalforce, totalmoment, totalforce_res, totalmoment_res = \
            np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.0])
        for force in self.forces:
            totalforce += force.vector
            totalmoment += np.cross(force.position, force.vector)
        for force in self.reactionforces:
            totalforce_res += force.vector
            totalmoment_res += np.cross(force.position, force.vector)
        for moment in self.moments:
            totalmoment += moment.vector
        for moment in self.reactionmoments:
            totalmoment_res += moment.vector
        assert(np.all(totalforce+totalforce_res<0.0000000001))
        assert(np.all(totalmoment+totalmoment_res<0.0000000001))

    def _single_yout_torque_shear(self, vect_mag, position, shearcentre):
        return vect_mag*(position-shearcentre)

    def _vectorized_numpy_torque_shear(self, vect_mag, position):
        return np.vectorize(partial(self._single_yout_torque_shear, vect_mag, position))

    def _always_returns_number(self, number, vector): # Don't ask me why, but this works, just leave it
        return number

    def _vectorized_number(self, number):
        return np.vectorize(partial(self._always_returns_number, number))

    def get_torque_shearflow(self, yout):
        firsty = [yout.item(0)]
        yout = [y+0.000000001 for y in yout[1:]]
        firsty.extend(yout)
        yout = np.asarray(firsty)
        shearflow = []
        shearflowfunc = []
        s_centres = self.shearcentre_z_func(yout)
        areas = self.areasfunc(yout)
        for force in self.forces:
            shearflow.append(self._vectorized_numpy_torque_shear(force.vector[0], force.position[2]))
            shearflowfunc.append(step_function_moment(force.position[1]))
        for force in self.reactionforces:
            shearflow.append(self._vectorized_numpy_torque_shear(force.vector[0], force.position[2]))
            shearflowfunc.append(step_function_moment(force.position[1]))
        for moment in self.moments:
            shearflow.append(self._vectorized_number(moment.vector[1]))
            shearflowfunc.append(step_function_moment(moment.position[1]))
        for moment in self.reactionmoments:
            shearflow.append(self._vectorized_number(moment.vector[1]))
            shearflowfunc.append(step_function_moment(moment.position[1]))
        for i in range(len(shearflow)):
            shearflow[i] = shearflow[i](s_centres)*shearflowfunc[i](yout)
        shearflow = np.sum(shearflow, axis=0)
        shearflow = shearflow * 0.5 / areas
        return shearflow

    def get_normalforce(self, yout):
        yout = np.array([y+0.000000001 for y in yout])
        normaly = []
        normalyfunc = []
        for force in self.forces:
            normaly.append(force.vector[1])
            normalyfunc.append(step_function_moment(force.position[1]))
        for force in self.reactionforces:
            normaly.append(force.vector[1])
            normalyfunc.append(step_function_moment(force.position[1]))
        for i in range(len(normaly)):
            normaly[i] = normaly[i]*normalyfunc[i](yout)
        return np.sum(normaly, axis=0)

    def get_shearforce(self, yout):
        yout = np.array([y+0.000001 for y in yout])
        shearx, shearz = [], []
        shearxfunc, shearzfunc = [], []
        for force in self.forces:
            shearx.append(force.vector[0])
            shearxfunc.append(step_function_moment(force.position[1]))
            shearz.append(force.vector[2])
            shearzfunc.append(step_function_moment(force.position[1]))
        for force in self.reactionforces:
            shearx.append(force.vector[0])
            shearxfunc.append(step_function_moment(force.position[1]))
            shearz.append(force.vector[2])
            shearzfunc.append(step_function_moment(force.position[1]))
        for i in range(len(shearx)):
            shearx[i] = shearx[i]*shearxfunc[i](yout)
        for i in range(len(shearz)):
            shearz[i] = shearz[i]*shearzfunc[i](yout)
        sx = np.sum(shearx, axis=0)
        sz = np.sum(shearz, axis=0)
        return sx, sz

    def get_moment(self, yout):
        yout = np.array([y+0.000001 for y in yout])
        momentx, momentz = [], []
        momentxfunc, momentzfunc = [], []
        for force in self.forces:
            momentx.append(-1*force.vector[2])
            momentxfunc.append(step_function(force.position[1]))
            momentz.append(force.vector[0])
            momentzfunc.append(step_function(force.position[1]))
        for force in self.reactionforces:
            momentx.append(-1*force.vector[2])
            momentxfunc.append(step_function(force.position[1]))
            momentz.append(force.vector[0])
            momentzfunc.append(step_function(force.position[1]))
        for moment in self.moments:
            momentx.append(moment.vector[0])
            momentxfunc.append(step_function_moment(moment.position[1]))
            momentz.append(moment.vector[2])
            momentzfunc.append(step_function_moment(moment.position[1]))
        for moment in self.reactionmoments:
            momentx.append(moment.vector[0])
            momentxfunc.append(step_function_moment(moment.position[1]))
            momentz.append(moment.vector[2])
            momentzfunc.append(step_function_moment(moment.position[1]))
        for i in range(len(momentx)):
            momentx[i] = momentx[i]*momentxfunc[i](yout)
        for i in range(len(momentz)):
            momentz[i] = momentz[i]*momentzfunc[i](yout)
        mx = np.sum(momentx, axis=0)
        mz = np.sum(momentz, axis=0)
        return mx, mz


    def do_bc_checks(self):
        if len(self.boundary_conditions) != 6:
            print("Too many or too few boundary conditions. Required amount of b.c. is 6, current amount is {}".format(
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
        self.do_bc_checks()

        self.matrix = np.zeros((10,10))
        self.bvect = np.zeros(10)
        for bc in self.boundary_conditions:
            angles = self.angle_bare(y=bc["y"], forces=self.forces, moments=self.moments)
            deflections = self.deflection_bare(y=bc["y"], forces=self.forces, moments=self.moments)
            if bc["angle_x"] != None:
                angle_x_bc = bc
                self.matrix[0,:4] = np.array([1, 0, 0, 0])
                self.bvect[0] = bc["angle_x"] - angles[0]
            if bc["defl_x"] != None:
                defl_x_bc = bc
                self.matrix[1,:4] = np.array([bc["y"], 1, 0, 0])
                self.bvect[1] = bc["defl_x"] - deflections[0]
            if bc["angle_z"] != None:
                angle_z_bc = bc
                self.matrix[2,:4] = np.array([0, 0, 1, 0])
                self.bvect[2] = bc["angle_z"] - angles[1]
            if bc["defl_z"] != None:
                defl_z_bc = bc
                self.matrix[3,:4] = np.array([0, 0, bc["y"], 1])
                self.bvect[3] = bc["defl_z"] - deflections[1]
            if bc["defl_y"] != None:
                defl_y_bc = bc
            if bc["angle_y"] != None:
                angle_y_bc = bc

        self.matrix[0, 5] = self.calc_coefficient(interpolant(self.EIzfunc, defl_x_bc["y"], False), angle_x_bc["y"], defl_x_bc["y"])
        self.matrix[0, 9] = -(defl_y_bc["x"] - angle_x_bc["x"])*self.calc_coefficient(interpolant(self.EIzfunc, defl_y_bc["y"], True), angle_x_bc["y"], defl_y_bc["y"])
        self.matrix[1, 4] = self.calc_coefficient(interpolant(self.EIzfunc, angle_x_bc["y"], True), defl_x_bc["y"], angle_x_bc["y"])
        self.matrix[1, 9] = -(defl_y_bc["x"] - defl_x_bc["x"]) *self.calc_coefficient(interpolant(self.EIzfunc, defl_y_bc["y"], True),  defl_x_bc["y"], defl_y_bc["y"], 2)
        self.matrix[2, 7] = self.calc_coefficient(interpolant(self.EIxfunc, defl_z_bc["y"], False), angle_z_bc["y"], defl_z_bc["y"])
        self.matrix[2, 9] = -(defl_y_bc["z"] - angle_z_bc["z"])*self.calc_coefficient(interpolant(self.EIxfunc, defl_y_bc["y"], True), angle_z_bc["y"], defl_y_bc["y"])
        self.matrix[3, 6] = self.calc_coefficient(interpolant(self.EIxfunc, angle_z_bc["y"], True), defl_z_bc["y"], angle_z_bc["y"])
        self.matrix[3, 9] = -(defl_y_bc["z"] - defl_z_bc["z"]) *self.calc_coefficient(interpolant(self.EIxfunc, defl_y_bc["y"], True),  defl_z_bc["y"], defl_y_bc["y"], 2)

        res_f_x, res_f_y, res_f_z = self.calc_resultant_force()
        res_m_x, res_m_y, res_m_z = self.calc_resultant_moment(np.array([defl_z_bc["x"], defl_z_bc["y"], defl_z_bc["z"]]))
        # res_m_x, res_m_z = self.calc_resultant_moment(angle_z_bc["y"], angle_x_bc["y"])
        self.bvect[4:] = np.array([-res_f_x, -res_m_x, -res_f_z, -res_m_z, -res_f_y, -res_m_y])
        self.matrix[4, 5], self.matrix[6,7], self.matrix[8,9] = 1.0, 1.0, 1.0 # Force equilibrium lines
        self.matrix[5, 4:] = np.array([0., 0., 1., defl_z_bc["y"]-angle_z_bc["y"], 0., -(defl_y_bc["z"]-angle_z_bc["z"])])
        self.matrix[7, 4:] = np.array([1., -(defl_x_bc["y"]-angle_x_bc["y"]), 0., 0., 0., defl_y_bc["x"]-angle_x_bc["x"]])
        self.matrix[9, 4:] = np.array([0., defl_x_bc["z"]-angle_y_bc["z"], 0., -(defl_z_bc["x"]-angle_y_bc["x"]), 1., 0.])

        self.solution = np.linalg.solve(self.matrix, self.bvect)
        self.C1 = np.array([self.solution[0],self.solution[2]])
        self.C2 = np.array([self.solution[1],self.solution[3]])
        self.reactionmoments.append(Moment(angle_x_bc["x"], angle_x_bc["y"], angle_x_bc["z"], zmag=self.solution[4]))
        self.reactionforces.append(Force(defl_x_bc["x"], defl_x_bc["y"], defl_x_bc["z"], xmag=self.solution[5]))
        self.reactionmoments.append(Moment(angle_z_bc["x"], angle_z_bc["y"], angle_z_bc["z"], xmag=self.solution[6]))
        self.reactionforces.append(Force(defl_z_bc["x"], defl_z_bc["y"], defl_z_bc["z"], zmag=self.solution[7]))
        self.reactionmoments.append(Moment(angle_z_bc["x"], angle_z_bc["y"], angle_z_bc["z"], ymag=self.solution[8]))
        self.reactionforces.append(Force(defl_z_bc["x"], defl_z_bc["y"], defl_z_bc["z"], ymag=self.solution[9]))
        self.check_equilibrium()
        return

    def calc_coefficient(self, interpolant, y, y_application, degree_of_integration=1):
        if degree_of_integration == 2:
            return quad(lambda t: quad(func=interpolant, a=y_application, b=t, epsabs=1.49e-05, epsrel=1.49e-05)[0], a=y_application, b=y, epsabs=1.49e-06, epsrel=1.49e-06)[0]
        else:
            return quad(func=interpolant, a=y_application, b=y, epsabs=1.49e-05, epsrel=1.49e-05)[0]

    def angle_bare(self, y, forces, moments):
        angle = np.array([0., 0.])
        for force in forces:
            angle[0] += force.vector[0] * self.calc_coefficient(interpolant(self.EIzfunc, force.position[1], False), y, force.position[1])
            angle[1] += force.vector[2] * self.calc_coefficient(interpolant(self.EIxfunc, force.position[1], False), y, force.position[1])
        for moment in moments:
            angle[0] += moment.vector[2] * self.calc_coefficient(interpolant(self.EIzfunc, moment.position[1], True), y, moment.position[1])
            angle[1] -= moment.vector[0] * self.calc_coefficient(interpolant(self.EIxfunc, moment.position[1], True), y, moment.position[1])
        return angle
        # return np.multiply(angle,np.array([1,-1]))

    def calc_angle(self, y_in):
        forces = list(np.append(self.forces, self.reactionforces))
        moments = list(np.append(self.moments, self.reactionmoments))
        try:
            angles = np.zeros((len(y_in),2))
            for idx, y in loadbar(enumerate(y_in)):
                angles[idx,:] = self.angle_bare(y=y, forces=forces, moments=moments) + self.C1
        except TypeError:
            angles = self.angle_bare(y=y_in, forces=forces, moments=moments) + self.C1
        return angles

    def deflection_bare(self, y, forces, moments):
        deflection = np.array([0., 0.])
        for force in forces:
            deflection[0] += force.vector[0] * self.calc_coefficient(interpolant(self.EIzfunc, force.position[1], False), y, force.position[1], 2)
            deflection[1] += force.vector[2] * self.calc_coefficient(interpolant(self.EIxfunc, force.position[1], False), y, force.position[1], 2)
        for moment in moments:
            deflection[0] += moment.vector[2] * self.calc_coefficient(interpolant(self.EIzfunc, moment.position[1], True), y, moment.position[1], 2)
            deflection[1] -= moment.vector[0] * self.calc_coefficient(interpolant(self.EIxfunc, moment.position[1], True), y, moment.position[1], 2)
        return deflection

    def calc_deflection(self, y_in):
        forces = list(np.append(self.forces, self.reactionforces))
        moments = list(np.append(self.moments, self.reactionmoments))
        try:
            deflections = np.zeros((len(y_in),2))
            for idx, y in loadbar(enumerate(y_in)):
                deflections[idx,:] = self.deflection_bare(y=y, forces=forces, moments=moments) + self.C2 + self.C1*y
        except TypeError:
            deflections = self.deflection_bare(y=y_in, forces=forces, moments=moments) + self.C2 + self.C1*y_in
        return deflections
        # return np.multiply(deflections,np.array([1,-1]))


def run_tests():
    global loadbar
    loadbar = lambda x: x
    EIfactor = 7.34
    EIz = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])*EIfactor
    ys = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    system = LoadCase(y_coordinates=ys, EIx=EIz, EIz=EIz)

    loads = {}
    testload = 3
    if testload == 1 or testload == 3:
        loads['end_f'] = Force(xpos=0.0, ypos=1.0, zpos=0.0, xmag=1.0, ymag=0.0, zmag=1.0)
    if testload == 2 or testload == 3:
        loads['end_m'] = Moment(xpos=0.0, ypos=1.0, zpos=0.0, xmag=3.0, ymag=0.0, zmag=3.0)
    system.add_loads(*loads.values())

    offset = 0.0
    system.add_boundary_condition(y=0.0+offset, defl_x=0)
    system.add_boundary_condition(y=0.0+offset, defl_y=0)
    system.add_boundary_condition(y=0.0+offset, defl_z=0)
    system.add_boundary_condition(y=0.0+offset, angle_x=0)
    system.add_boundary_condition(y=0.0+offset, angle_y=0)
    system.add_boundary_condition(y=0.0+offset, angle_z=0)
    system.solve_bcs()  # Includes test for equilibrium

    def exact_deflection(yout, length=1.0):
        factor = direction*2 - 1
        if testload == 1:
            return loads['end_f'].vector[0]*yout*yout*(3.-yout)/(6.*EIfactor)*length*length*length
        elif testload == 2:
            return factor * loads['end_m'].vector[2]*yout*yout/(2.*EIfactor)*length*length
        elif testload == 3:
            return loads['end_f'].vector[0]*yout*yout*(3.-yout)/(6.*EIfactor)*length*length*length + \
                   factor * loads['end_m'].vector[2] * yout*yout/(2.*EIfactor)*length*length

    def exact_angle(length=1.0):
        factor = direction*2 - 1
        if testload == 1:
            return loads['end_f'].vector[0]*length * length / (2. * EIfactor)
        elif testload == 2:
            return factor * loads['end_m'].vector[2]*length / EIfactor
        elif testload == 3:
            return loads['end_f'].vector[0]*length*length/(2.*EIfactor) + \
                   factor * loads['end_m'].vector[2]* length / EIfactor

    for direction in [0,1]:
        yout_amount = 6
        yout = np.linspace(0.0001,1.0,yout_amount)
        defl_exact = exact_deflection(yout)
        defl_out = system.calc_deflection(yout)
        angle_out = system.calc_angle(yout)
        angle_exact = exact_angle(1.0)
        assert(np.all(defl_out[:, direction]-defl_exact)<0.00000000001)
        assert((angle_out[-1, direction]-angle_exact)<0.00000000001)
    assert(np.all(system.calc_deflection([-0.05])<0.00000000001))
    system.remove_all_boundary_conditions()
    system.remove_all_loads()
    loads = {}
    loads['f1'] = Force(xpos=0.0, ypos=0.5, zpos=0.0, xmag=1., ymag=0.0, zmag=0.0)
    loads['f2'] = Force(xpos=0.0, ypos=0.75, zpos=0.0, xmag=-50.0, ymag=0.0, zmag=10.0)
    loads['f3'] = Force(xpos=0.0, ypos=1.0, zpos=0.0, xmag=2.0, ymag=0.0, zmag=0.0)
    loads['m1'] = Moment(xpos=0.0, ypos=0.8, zpos=0.0, xmag=0.0, ymag=0.0, zmag=50.)
    loads['m1'] = Moment(xpos=0.0, ypos=1.0, zpos=0.0, xmag=0.0, ymag=0.0, zmag=-25.)
    system.add_loads(*loads.values())
    offset = 0.0
    system.add_boundary_condition(y=0.0+offset, defl_x=0)
    system.add_boundary_condition(y=0.0+offset, defl_y=0)
    system.add_boundary_condition(y=0.0+offset, defl_z=0)
    system.add_boundary_condition(y=0.0+offset, angle_x=0)
    system.add_boundary_condition(y=0.0+offset, angle_y=0)
    system.add_boundary_condition(y=0.0+offset, angle_z=0)
    system.solve_bcs()
    yout = np.array([0.525, 0.785])
    sx, sz = system.get_shearforce(yout)
    mx, mz = system.get_moment(yout)
    assert(np.all((sx-np.array([48, -2]))<0.00000000001))
    assert(np.all((mz-np.array([14.7, 25.43]))<0.001))
    assert(np.all((sz-np.array([-10, 0]))<0.00000000001))
    assert(np.all((mx-np.array([-2.25, 0]))<0.0001))
    print("Tests passed succesfully")


if __name__ == "__main__":
    run_tests()
    # system = LoadCase([0.0, 0.5, 1.0], [2, 2, 2], [2, 2, 2])







