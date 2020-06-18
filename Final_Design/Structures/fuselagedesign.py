from plateclass import *
from stringerclass import *
from loadinfo import *
from Fuselageclass import *
import copy
from scipy.interpolate import interp1d
from numpy import pi
import numpy as np
from math import sqrt, cos, sin
from matplotlib import pyplot as plt


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
        # print(sheet.xbar, sheet.ybar, "|||||||||", x, y, angle)
        pan = Panel(x, y, sheet, rotation=angle)
        # print(pan.xbar_global, pan.ybar_global)
        cs.addPanels(pan)
        # print("=======================")
    # print()
    cs.recalculate_component_effects()
    return cs

def create_round_segment(material, length, ts, height, width=None, r_top=None, r_bot=None):
    cs = CircCrossSection(length, height*0.5, ts, material)
    return cs


def generate_fuselage(data, sparlocs, zshift, ts, base_material, circ_material=None, stiff_material=None):
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
            cs = create_compound_segment(base_material, length, ts(y), interph(y), width, interprt(y), interprb(y))
        else:
            width = interph(y)
            cs = create_round_segment(circ_material, length, ts(y), interph(y), width, interprt(y), interprb(y))

        fus.addSection(FuselageSection(y, xref+width*0.5, zref, cs))

    fus.recalculate_lists()
    return fus, sparlocs

def make_default_fuselage(skint):
    data = np.genfromtxt("fuselage_data.csv", delimiter=",", skip_header=1)
    data *= 0.001

    aluminium = MatProps(sigma_y=450000000, E=72400000000, poisson=0.33, rho=2.87, name="AA2024", alpha=0.8,
                         n=0.6)
    carbon = MatProps(sigma_y=600000000, E=70000000000, poisson=0.1, rho=1.60, sigma_comp=570000000,
                      name="carbonfibre")
    sparlocs = np.concatenate((np.arange(22, 60, 1)*0.1, np.array([5.95, 6.0, 6.05]), np.arange(61, 95, 1)*0.1))[::-1]
    zfunc = lambda y: 0.0022*y*y*y-0.0458*y*y+0.3535*y-1.10606383
    zshift = [zfunc(y) for y in sparlocs]
    fus, sparlocs = generate_fuselage(data, sparlocs, zshift, skint, base_material=aluminium, circ_material=carbon)
    return fus, sparlocs

def analyse_fuselage(sparlocs, fus, bcs, *loads):
    EIx = fus.E*fus.Ixx
    EIz = fus.E*fus.Izz
    areas = fus.areas
    shearcentres = fus.zbars
    system = LoadCase(y_coordinates=np.asarray(sparlocs), EIx=EIx, EIz=EIz, areas=areas, shearcentre_z = shearcentres)
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
    old_sparlocs = sparlocs
    if flag == True:
        sparlocs = np.concatenate((sparlocs, np.array([ymin])))
        EIx = np.concatenate((EIx, np.array([EIx[-1]])))
        EIz = np.concatenate((EIz, np.array([EIz[-1]])))
        system.change_geometry(sparlocs, EIx, EIz)
    system.solve_bcs()

    yout = np.concatenate((sparlocs, np.array([sparlocs[-1]-0.000003, 0.])))
    areas = np.concatenate((areas, np.array([areas[-1], areas[-1], areas[-1]])))
    shearcentres = np.concatenate((shearcentres, np.array([shearcentres[-1], shearcentres[-1], shearcentres[-1]])))
    print("in-function:", min(sparlocs), max(sparlocs))
    print(len(areas), len(shearcentres), len(yout))
    system.change_areas(areas, shearcentres, y_locs=yout)
    sx, sz = system.get_shearforce(yout)
    mx, mz = system.get_moment(yout)
    ny = system.get_normalforce(yout)
    print("yout:", min(yout), max(yout))
    shearflow = system.get_torque_shearflow(yout)
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
    print(totalforce, totalforce_res)
    print(totalmoment, totalmoment_res)
    return yout, old_sparlocs, fus, sx, sz, mx, mz, ny, shearflow

def n_stress_moments(x, z, mx, mz, ixx, izz, ixz):
    return ((mx*izz-mz*ixz)*z + (mz*ixx-mx*ixz)*x)/(ixx*izz-ixz*ixz)

def n_stress_int_n(area, ny):
    return ny/area

def analyse_forces(yout, sparlocs, fus, sx, sz, mx, mz, ny, shearflow):
    maxstresses = []
    minstresses = []
    critloc_mins = []
    critloc_maxs = []
    normalstresses = []
    for sparloc, section in zip(sparlocs, fus.sections):
        critlen = sparloc - section.properties.l
        nyi = ny[np.abs(yout - critlen) < 0.001]
        normalstresses.append((nyi/section.properties.total_area).item())
        sxi = sx[np.abs(yout - critlen) < 0.001]
        szi = sz[np.abs(yout - critlen) < 0.001]
        mxi = mx[np.abs(yout - critlen) < 0.001]
        mzi = mz[np.abs(yout - critlen) < 0.001]
        maxstress = np.array([0])
        minstress = np.array([0])
        critloc_min = (0,0)
        critloc_max = (0,0)
        if isinstance(section.properties, CircCrossSection):
            for i in range(8):
                x, z = section.xbar_global, section.zbar_global
                angle = i*pi*0.25
                x += cos(angle)*section.properties.radius
                z += sin(angle)*section.properties.radius
                stress = n_stress_moments(x, z, mxi, mzi, section.Ixx_global, section.Izz_global, section.Ixz_global)
                if stress < minstress:
                    minstress = stress
                    critloc_min = (x, z)
                if stress > maxstress:
                    maxstress = stress
                    critloc_max = (x, z)
        else:
            for plate in section.properties.plates:
                x, z = -plate.xbar_global + section.xpos, plate.ybar_global + section.zpos
                for k in range(2):
                    i = k*2-1
                    length = i * plate.properties.b*0.5
                    (xi, zi) = (x + length*cos(plate.rotation), z - length*sin(plate.rotation))
                    stress = n_stress_moments(xi, zi, mxi, mzi, section.Ixx_global, section.Izz_global, section.Ixz_global)
                    if stress < minstress:
                        minstress = stress
                        critloc_min = (xi, zi)
                    if stress > maxstress:
                        maxstress = stress
                        critloc_max = (xi, zi)
                    if stress < -plate.properties.sigma_cr:
                        print("Buckling occurs at y={}, x={}, z={}".format(sparloc, xi, zi))

            for stiff in section.properties.stiffeners:
                x, z = -stiff.xbar_global + section.xpos, stiff.ybar_global + section.zpos
                stress = n_stress_moments(x, z, mxi, mzi, section.Ixx_global, section.Izz_global, section.Ixz_global)
                if stress < minstress:
                    minstress = stress
                    critloc_min = (x, z)
                if stress > maxstress:
                    maxstress = stress
                    critloc_max = (x, z)
        maxstresses.append(maxstress.item())
        critloc_maxs.append(critloc_max)
        minstresses.append(minstress.item())
        critloc_mins.append(critloc_min)
    maxstresses = np.array(maxstresses)
    minstresses = np.asarray(minstresses)
    normalstresses = np.asarray(normalstresses)
    critloc_maxs = np.asarray(critloc_maxs)
    critloc_mins = np.asarray(critloc_mins)
    return maxstresses, critloc_maxs, minstresses, critloc_mins, normalstresses


def skinthickness(y):
    if y > 6.099:
        return 0.001
    else:
        return 0.002


if __name__ == "__main__":
    t = skinthickness
    fus, sparlocs = make_default_fuselage(t)
    y_origin=1.835
    bc1 = {'y' : 0.0+y_origin, 'defl_x' : 0.0}
    bc2 = {'y' : 0.0+y_origin, 'defl_z' : 0.0}
    bc3 = {'y' : 0.0+y_origin, 'angle_x' : 0.0}
    bc4 = {'y' : 0.0+y_origin, 'angle_z' : 0.0}
    bc5 = {'y' : 0.0+y_origin, 'defl_y' : 0.0}
    bc6 = {'y' : 0.0+y_origin, 'angle_y' : 0.0}
    bcs = [bc1, bc2, bc3, bc4, bc5, bc6]
    loads = {}
    loads["fhtail"] = Force(xpos=0.0, ypos=8.9, zpos=0.0, xmag=0, zmag=-8144, ymag=815)
    loads["fvtail"] = Force(xpos=0.0, ypos=8.9, zpos=0.0, xmag=408, ymag=41, zmag=0)

    yout, old_sparlocs, fus, sx, sz, mx, mz, ny, shearflow = analyse_fuselage(sparlocs, fus, bcs, loads)

    sig_max, loc_max, sig_min, loc_min, sig_norm = analyse_forces(yout, old_sparlocs, fus, sx, sz, mx, mz, ny, shearflow)
    print(fus.mass)
    print(fus.zbars)
    plt.plot(old_sparlocs, sig_min)
    plt.plot(old_sparlocs, sig_max)
    plt.show()

