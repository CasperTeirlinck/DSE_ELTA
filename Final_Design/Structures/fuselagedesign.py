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


def calc_tresca(normal, shear):
    return np.sqrt(normal*normal+3*shear*shear)

def n_stress_moments(x, z, mx, mz, ixx, izz, ixz):
    return ((mx*izz-mz*ixz)*z + (mz*ixx-mx*ixz)*x)/(ixx*izz-ixz*ixz)

def n_stress_int_n(area, ny):
    return ny/area

def shear_stress(sx, sz, Qx, Qz, ixx, izz, ixz):
    return -(sz*izz-sx*ixz)/(ixx*izz-ixz*ixz)*Qx-(sx*ixx-sz*ixz)/(ixx*izz-ixz*ixz)*Qz

def shearflow_to_stress(q, ts):
    return q*ts

def construct_locs(x, w_sheet, Ns, y=0.0, angle=0.0, ts=0.0):
    return [(w_sheet * n / (Ns+1), ts) for n in range(1, Ns+1)]

def construct_locs_circ(r, Ns, ts=0.0):
    locs = [(r-(r-ts)*sin(n*2*pi/Ns), r+(r-ts)*cos(n*2*pi/Ns)) for n in range(Ns)]
    angles = [pi+n*2*pi/Ns for n in range(Ns)]
    return locs, angles


def create_compound_segment(material, length, ts, height, width, r_top, r_bot, stringerlist=[], stringerargs=None, n_longs=0, longargs=None):
    cs = CrossSection(length)
    wbot = width-2*r_bot
    wtop = width-2*r_top
    hleft = height-r_bot-r_top
    hright = hleft
    xs = [r_bot, r_bot+wbot, width, width, r_top+wtop, r_top, 0., 0.]
    ys = [0., 0., r_bot, r_bot+hleft, height, height, r_bot+hright, r_bot]
    if n_longs != 0:
        longs_locs = [(r_bot+ts, ts), (r_bot+wbot-ts, ts), (width-ts, r_bot-ts),(width-ts, r_bot+hleft-ts),
                      (r_top+wtop-ts,height-ts), (r_top+ts, height-ts), (ts, r_bot+hright-ts), (ts, r_bot+ts)]
        angles = [-pi*0.25, pi*0.25, pi*0.25, pi*0.75, pi*0.75, pi*1.25, pi*1.25, -pi*0.25]
    widths = [wbot, sqrt(2)*r_bot, hleft, sqrt(2)*r_top, wtop, sqrt(2)*r_top, hright, sqrt(2)*r_bot]
    j = J_Stringer(Le=length, **stringerargs)
    z = Z_Stringer(Le=length, **stringerargs)
    for idx, (x, y, w_sheet, n_stiff) in enumerate(zip(xs, ys, widths, stringerlist)):
        angle = idx*pi*0.25
        sheet = Sheet(length, w_sheet, ts, material, mode=2)
        locs = construct_locs(x, w_sheet, int(round(n_stiff)), y=y, angle=angle, ts=sheet.ts)
        if idx == 0 or idx == 1 or idx == 7:
            template_stringer = z
        else:
            template_stringer = j

        stringers = make_stringers(locs, template_stringer)
        sheet.construct_panel(*stringers)
        # print(angle*180/pi)
        pan = Panel(x, y, sheet, rotation=angle)
        cs.addPanels(pan)
    if n_longs != 0:
        template_longeron = J_Stringer(Le=length, **longargs)
        for loc, angle in zip(longs_locs, angles):
            longeron = make_stringers([loc], template_longeron, angle)
            cs.addStringer(*longeron)
    cs.recalculate_component_effects()
    return cs

def create_round_segment(material, length, ts, height, width=None, r_top=None, r_bot=None, n_stiff_circ=0, circstringerargs=None, n_longs=0, longargs=None):
    cs = CircCrossSection(length, height*0.5, ts, material)
    if n_stiff_circ != 0:
        locs, angles = construct_locs_circ(height*0.5, n_stiff_circ, ts)
        for loc, angle in zip(locs, angles):
            stringer = make_stringers([loc], Z_Stringer(Le=length, **circstringerargs), rotation=angle)
            cs.addStringer(*stringer)
    cs.recalculate_component_effects()
    return cs


def generate_fuselage(data, sparlocs, zshift, ts, base_material, circ_material=None, stringerlist=[], stringerargs=None, n_stiff_circ=0, circstringerargs=None, n_longs=0, longeronsarg=None):
    if circ_material is None:
        circ_material = base_material
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
        if length == 0.0:
            length = 0.00001
        if y < 6.099:
            width = interpw(y)
            cs = create_compound_segment(base_material, length, ts(y), interph(y), width, interprt(y), interprb(y), stringerlist, stringerargs, n_longs=n_longs, longargs=longeronsarg)
        else:
            width = interph(y)
            cs = create_round_segment(circ_material, length, ts(y), interph(y), width, interprt(y), interprb(y), n_stiff_circ=n_stiff_circ, circstringerargs=circstringerargs, n_longs=n_longs, longargs=longeronsarg)

        fus.addSection(FuselageSection(y, xref+width*0.5, zref, cs))

    fus.recalculate_lists()
    return fus, sparlocs

def make_default_fuselage(skint, stringerlist=[], stringerargs=None, materials=None, sparlocs=None, n_stiff_circ=0, circstringerargs=None, n_longs=0, longargs=None):
    data = np.genfromtxt("fuselage_data.csv", delimiter=",", skip_header=1)
    data *= 0.001
    if materials is None:
        aluminium = MatProps(sigma_y=450000000, E=72400000000, poisson=0.33, rho=2.87, name="AA2024", alpha=0.8,
                             n=0.6)
        carbon = MatProps(sigma_y=600000000, E=70000000000, poisson=0.1, rho=1.60, sigma_comp=570000000,
                          name="carbonfibre")
        materials = [aluminium, carbon]
    if sparlocs is None:
        sparlocs = np.concatenate((np.arange(22, 60, 1)*0.1, np.array([5.95, 6.0, 6.05]), np.arange(61, 95, 1)*0.1))[::-1]
    zfunc = lambda y: 0.0022*y*y*y-0.0458*y*y+0.3535*y-1.10606383
    zshift = [zfunc(y) for y in sparlocs]
    fus, sparlocs = generate_fuselage(data, sparlocs, zshift, skint, base_material=materials[0], circ_material=materials[1], stringerlist=stringerlist, stringerargs=stringerargs, n_stiff_circ=n_stiff_circ, circstringerargs=circstringerargs, n_longs=n_longs, longeronsarg=longargs)
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
                if single_load.position[1] < ymin:
                    ymin = single_load.position[1]
                    flag = True
        else:
            print("Unrecognized load:", load)
    for bc in bcs:
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
    yout = np.concatenate((sparlocs, np.array([sparlocs[-1]-0.000003, 0.])))
    areas = np.concatenate((areas, np.array([areas[-1], areas[-1], areas[-1]])))
    shearcentres = np.concatenate((shearcentres, np.array([shearcentres[-1], shearcentres[-1], shearcentres[-1]])))
    system.change_areas(areas, shearcentres, y_locs=yout)

    system.solve_bcs()

    sx, sz = system.get_shearforce(yout)
    mx, mz = system.get_moment(yout)
    ny = system.get_normalforce(yout)
    shearflow = system.get_torque_shearflow(yout)

    return yout, old_sparlocs, fus, sx, sz, mx, mz, ny, shearflow


def analyse_forces(yout, sparlocs, fus, sx, sz, mx, mz, ny, q, printing=True):
    flag = np.zeros(8)
    maxstresses = []
    minstresses = []
    shearstressesmax = []
    critloc_mins = []
    critloc_maxs = []
    critloc_shears = []
    critloc_trescas = []
    normalstresses = []
    trescastresses = []
    for sparloc, section in zip(sparlocs, fus.sections):
        critlen = sparloc - section.properties.l
        nyi = ny[np.abs(yout - critlen) < 0.001]
        normalstresses.append((nyi/section.properties.total_area).item())
        sxi = sx[np.abs(yout - critlen) < 0.001]
        szi = sz[np.abs(yout - critlen) < 0.001]
        mxi = mx[np.abs(yout - critlen) < 0.001]
        mzi = mz[np.abs(yout - critlen) < 0.001]
        qi = q[np.abs(yout-critlen) < 0.001]
        maxstress = np.array([0])
        minstress = np.array([0])
        shearstressmax = np.array([0])
        trescamax = np.array([0])
        critloc_min = (0,0)
        critloc_max = (0,0)
        critloc_shear = (0,0)
        critloc_tresca = (0,0)
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
                shearstress = shearflow_to_stress(qi, section.properties.thickness)
                Q = 2*section.properties.thickness*section.properties.radius*section.properties.radius
                shearstress += abs(shear_stress(sxi, szi, Q, Q, section.Ixx_global, section.Izz_global, section.Ixz_global))
                if shearstress > shearstressmax:
                    shearstressmax = shearstress
                    critloc_shear = (x, z)
                tresca = calc_tresca(max(stress, -stress), shearstress)
                if tresca > section.properties.tresca_yield:
                    flag[0] = 1
                    if printing:
                        print("Tresca yield exceeded at y={}, x={}, z={}; stress={}, max_stress={}".format(sparloc, x, z, tresca, section.properties.tresca_yield))
                if tresca > trescamax:
                    trescamax = tresca
                    critloc_tresca = (x, z)
                if abs(minstress)>abs(section.properties.sigma_cr):
                    flag[1] = 1
                    if printing:
                        print("Buckling occurs at y={}, x={}, z={}; stress={}, sig_cr={}".format(sparloc, x, z, stress,
                                                                                             -section.properties.sigma_cr))

            if maxstress > section.properties.material.sigma_y:
                flag[2] = 1
                if printing:
                    print("Ultimate stress exceeded at y={}, x={}, z={}; stress={}, ult_stress={}"
                          .format(sparloc, x, z, maxstress, section.properties.material.sigma_y))
        else:
            for plate in section.properties.plates:
                for k in range(3):
                    (x, z) = plate.ipos[k]
                    xi = -x + section.xpos
                    zi = z + section.zpos
                    stress = n_stress_moments(xi, zi, mxi, mzi, section.Ixx_global, section.Izz_global, section.Ixz_global)
                    if stress < minstress:
                        minstress = stress
                        critloc_min = (xi, zi)
                    if stress > maxstress:
                        maxstress = stress
                        critloc_max = (xi, zi)
                    if stress < -plate.properties.sigma_cr:
                        flag[3] = 1
                        if printing:
                            print("Buckling occurs at y={}, x={}, z={}; stress={}, sig_cr={}".format(sparloc, xi, zi, stress, -plate.properties.sigma_cr))

                    Qx = np.sum(section.properties.Qxs[np.where(section.properties.Q_ys >= z-0.001)])
                    Qy = np.sum(section.properties.Qys[np.where(section.properties.Q_xs >= x-0.001)])
                    shearstress = shearflow_to_stress(qi, plate.properties.ts)
                    shearstress += abs(shear_stress(sxi, szi, Qx, Qy, section.Ixx_global, section.Izz_global, section.Ixz_global))
                    if shearstress > shearstressmax:
                        shearstressmax = shearstress
                        critloc_shear = (xi, zi)
                    tresca = calc_tresca(max(stress, -stress), shearstress)
                    if tresca > plate.properties.tresca_yield:
                        flag[4] = 1
                        if printing:
                            print("Tresca yield exceeded at y={}, x={}, z={}; stress={}, max_stress={}"
                                  .format(sparloc, xi, zi, tresca, section.properties.tresca_yield))
                    if tresca > trescamax:
                        trescamax = tresca
                        critloc_tresca = (xi, zi)
                    if maxstress > plate.properties.ultimatestress:
                        flag[5] = 1
                        if printing:
                            print("Ultimate stress exceeded at y={}, x={}, z={}; stress={}, ult_stress={}"
                                  .format(sparloc, xi, zi, stress, plate.ultimatestress))
            for stiff in section.properties.stiffeners:
                x, z = -stiff.xbar_global + section.xpos, stiff.ybar_global + section.zpos
                stress = n_stress_moments(x, z, mxi, mzi, section.Ixx_global, section.Izz_global, section.Ixz_global)
                if stress < minstress:
                    minstress = stress
                    critloc_min = (x, z)
                if stress > maxstress:
                    maxstress = stress
                    critloc_max = (x, z)
                if stress < 0 and stress < -stiff.properties.sigma_cr:
                    flag[6] = 1
                    if printing:
                        print("Buckling of stiffener at y={}, x={}, z={}; stress={}, sigma_cr={}".format(sparloc, x, z, stress, -stiff.properties.sigma_cr))
                if maxstress > stiff.properties.material.sigma_y:
                    flag[7] = 1
                    if printing:
                        print("Ultimate stress exceeded at y={}, x={}, z={}; stress={}, ult_stress={}"
                              .format(sparloc, x, z, maxstress, stiff.properties.materials.sigma_y))
        maxstresses.append(maxstress.item())
        critloc_maxs.append(critloc_max)
        minstresses.append(minstress.item())
        critloc_mins.append(critloc_min)
        shearstressesmax.append(shearstressmax.item())
        critloc_shears.append(critloc_shear)
        trescastresses.append(tresca)
        critloc_trescas.append(critloc_tresca)
    maxstresses = np.asarray(maxstresses)
    minstresses = np.asarray(minstresses)
    normalstresses = np.asarray(normalstresses)
    shearstressesmax = np.asarray(shearstressesmax)
    trescastresses = np.asarray(trescastresses)
    critloc_maxs = np.asarray(critloc_maxs)
    critloc_mins = np.asarray(critloc_mins)
    critloc_shears = np.asarray(critloc_shears)
    critloc_trescas = np.asarray(critloc_trescas)
    return flag, maxstresses, critloc_maxs, minstresses, critloc_mins, normalstresses, shearstressesmax, critloc_shears, trescastresses, critloc_trescas


def skinthickness(modifier, y):
    if y > 6.099:
        return 0.001*modifier
    else:
        return 0.001*modifier

def constructstringerarg(modifier, material):
    return {'material':material, 't1':0.0005*modifier, 't2':0.0005*modifier, 't3':0.0005*modifier, 't4':0.0005*modifier, 'b1':0.005*modifier, 'b2':0.015*modifier, 'b3':0.003*modifier, 'b4':0.003*modifier}

def constructlongeronarg(modifier, material):
    return {'material':material, 't1':0.0005*modifier, 't2':0.0005*modifier, 't3':0.0005*modifier, 't4':0.0005*modifier, 'b1':0.005*modifier, 'b2':0.01*modifier, 'b3':0.005*modifier, 'b4':0.005*modifier}



# def design_fuselage(v):


if __name__ == "__main__":
    from materials import materials
    from functools import partial
    t = np.vectorize(partial(skinthickness, 1))
    Wbat = 256*9.81
    mats = materials()
    stringerargs = constructstringerarg(1.0, mats['alu2024'])
    circstringerargs = constructstringerarg(1.0, mats['alu2024'])
    longargs = constructlongeronarg(4.0, mats['carbonfibre'])
    n_stringers = np.ones(8)*5
    n_stiff_circ = 6
    n_longs = 8

    fus, sparlocs = make_default_fuselage(t, stringerlist=n_stringers, stringerargs=stringerargs, materials=[mats['alu7075'], mats['carbonfibre']], n_stiff_circ=n_stiff_circ, circstringerargs=circstringerargs, n_longs=n_longs, longargs=longargs)

    # print(fus.sections[-2].properties.plates[0].properties.tresca_yield)

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
    loads["battery1"] = Force(xpos=0.0, ypos=2.3, zpos=0.0, xmag=0, ymag=0, zmag=-Wbat*0.5*1.5)
    loads["battery2"] = Force(xpos=0.0, ypos=2.7, zpos=0.0, xmag=0, ymag=0, zmag=-Wbat*0.5*1.5)

    yout, old_sparlocs, fus, sx, sz, mx, mz, ny, shearflow = analyse_fuselage(sparlocs, fus, bcs, loads)

    flag, sig_max, loc_max, sig_min, loc_min, sig_norm, shear, loc_shear, tresca, loc_tresca = analyse_forces(yout, old_sparlocs, fus, sx, sz, mx, mz, ny, shearflow, printing=True)
    # print(fus.mass)
    # print(np.where(np.abs(np.asarray(old_sparlocs)-6.0)<0.15))

    secs = np.asarray(fus.sections)
    secs = secs[np.where(np.abs(np.asarray(old_sparlocs)-6.0)<0.25)]
    # for sec in secs[1:3]:
    for sec in secs[1:-4]:
        sec.properties.plot_section(show=True)
    plt.plot(old_sparlocs, fus.Ixx)
    # plt.plot(yout, mx)
    # plt.plot(yout, sz)
    # plt.plot(yout, mz)
    # plt.plot(yout, sx)
    # plt.plot(old_sparlocs, sig_max)
    plt.show()

