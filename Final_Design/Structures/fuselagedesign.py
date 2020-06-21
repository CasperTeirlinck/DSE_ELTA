try:
    from plateclass import *
    from stringerclass import *
    from loadinfo import *
    from Fuselageclass import *
except ModuleNotFoundError:
    from Structures.plateclass import *
    from Structures.stringerclass import *
    from Structures.loadinfo import *
    from Structures.Fuselageclass import *

import copy
from scipy.interpolate import interp1d
from numpy import pi
import numpy as np
from math import sqrt, cos, sin
from matplotlib import pyplot as plt
from functools import partial
from itertools import product
from Final_Design.new_variables import NewVariables


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

def construct_locs(w_sheet, Ns, ts=0.0):
    return [(w_sheet * n / (Ns+1), ts) for n in range(1, Ns+1)]

def construct_locs_circ(r, Ns, ts=0.0):
    locs = [(r-(r-ts)*sin(n*2*pi/Ns), r+(r-ts)*cos(n*2*pi/Ns)) for n in range(Ns)]
    angles = [pi+n*2*pi/Ns for n in range(Ns)]
    return locs, angles


def change_parameter(minmax, parameter, increase, amount=0.1):
    if isinstance(parameter, int):
        if increase:
            if parameter != minmax[1]:
                return parameter + 1, True
            else:
                return minmax[1], False
        else:
            if parameter != minmax[0]:
                return parameter - 1, True
            else:
                return minmax[0], False
    elif isinstance(parameter, float):
        if increase:
            if minmax[1]-parameter > parameter*amount:
                return parameter * (1.0 + amount), True
            else:
                return minmax[1], False
        else:
            if parameter-minmax[0] > parameter*amount:
                return parameter * (1.0 - amount), True
            else:
                return minmax[0], False


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
        locs = construct_locs(w_sheet, int(round(n_stiff)), ts=sheet.ts)
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
    if n_longs != 0:
        r = height*0.5
        locs = [(r-(r-ts)*sin(n*2*pi/4+pi*0.25), r+(r-ts)*cos(n*2*pi/4+pi*0.25)) for n in range(4)]
        angles = [pi+n*2*pi/4+pi*0.25 for n in range(4)]
        template_longeron = J_Stringer(Le=length, **longargs)
        for loc, angle in zip(locs, angles):
            longeron = make_stringers([loc], template_longeron, rotation=angle)
            cs.addStringer(*longeron)
    cs.recalculate_component_effects()
    return cs


def generate_fuselage(data, framelocs, zshift, ts, base_material, circ_material=None, stringerlist=[], stringerargs=None, n_stiff_circ=0, circstringerargs=None, n_longs=0, longeronsarg=None):
    if circ_material is None:
        circ_material = base_material
    assert(len(framelocs)==len(zshift))
    fus = Fuselage(framelocs)
    xref = 0.0
    ys = data[:,0]
    interpw = interp1d(ys, data[:,1])
    interph = interp1d(ys, data[:,2])
    interprt = interp1d(ys, data[:,3])
    interprb = interp1d(ys, data[:,4])
    for idx, (y, zref) in enumerate(zip(framelocs, zshift)):
        if idx == len(framelocs)-1:
            length = 2.2-y
        else:
            length = y-framelocs[idx+1]
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
    return fus

def make_fuselage(skint, stringerlist=[], stringerargs=None, materials=None, framelocs=None, n_stiff_circ=0, circstringerargs=None, n_longs=0, longargs=None):
    data = np.genfromtxt("fuselage_data.csv", delimiter=",", skip_header=1)
    data *= 0.001
    if materials is None:
        aluminium = MatProps(sigma_y=450000000, E=72400000000, poisson=0.33, rho=2.87, name="AA2024", alpha=0.8,
                             n=0.6)
        carbon = MatProps(sigma_y=600000000, E=70000000000, poisson=0.1, rho=1.60, sigma_comp=570000000,
                          name="carbonfibre")
        materials = [aluminium, carbon]
    if framelocs is None:
        framelocs = np.concatenate((np.arange(22, 60, 1)*0.1, np.array([5.95, 6.0, 6.05]), np.arange(61, 95, 1)*0.1))[::-1]
    zfunc = lambda y: 0.0022*y*y*y-0.0458*y*y+0.3535*y-1.10606383
    zshift = [zfunc(y) for y in framelocs]
    fus = generate_fuselage(data, framelocs, zshift, skint, base_material=materials[0], circ_material=materials[1], stringerlist=stringerlist, stringerargs=stringerargs, n_stiff_circ=n_stiff_circ, circstringerargs=circstringerargs, n_longs=n_longs, longeronsarg=longargs)
    return fus, framelocs

def analyse_fuselage(framelocs, fus, bcs, *loads, yout=None):
    EIx = fus.E*fus.Ixx
    EIz = fus.E*fus.Izz
    areas = fus.areas
    shearcentres = fus.zbars
    system = LoadCase(y_coordinates=np.asarray(framelocs), EIx=EIx, EIz=EIz, areas=areas, shearcentre_z = shearcentres)
    ymin = min(framelocs)
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
    old_framelocs = framelocs
    if flag == True:
        framelocs = np.concatenate((framelocs, np.array([ymin])))
        EIx = np.concatenate((EIx, np.array([EIx[-1]])))
        EIz = np.concatenate((EIz, np.array([EIz[-1]])))
        system.change_geometry(framelocs, EIx, EIz)
    if yout is None:
        yout = np.concatenate((framelocs, np.array([framelocs[-1]-0.000003, 0.])))
    areas = np.concatenate((areas, np.array([areas[-1], areas[-1], areas[-1]])))
    shearcentres = np.concatenate((shearcentres, np.array([shearcentres[-1], shearcentres[-1], shearcentres[-1]])))
    system.change_areas(areas, shearcentres, y_locs=yout)

    system.solve_bcs()

    sx, sz = system.get_shearforce(yout)
    mx, mz = system.get_moment(yout)
    ny = system.get_normalforce(yout)
    shearflow = system.get_torque_shearflow(yout)

    return yout, old_framelocs, fus, sx, sz, mx, mz, ny, shearflow


def analyse_forces(yout, framelocs, fus, sx, sz, mx, mz, ny, q, printing=True):
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
    for frameloc, section in zip(framelocs, fus.sections):
        critlen = frameloc - section.properties.l
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
                        print("Tresca yield exceeded at y={}, x={}, z={}; Circular section; stress={}, max_stress={}".format(frameloc, x, z, tresca, section.properties.tresca_yield))
                if tresca > trescamax:
                    trescamax = tresca
                    critloc_tresca = (x, z)
                if abs(minstress)>abs(section.properties.sigma_cr):
                    flag[1] = 1
                    if printing:
                        print("Buckling occurs at y={}, x={}, z={}; Circular section; stress={}, sig_cr={}".format(frameloc, x, z, stress,
                                                                                             -section.properties.sigma_cr))

            if maxstress > section.properties.material.sigma_y:
                flag[2] = 1
                if printing:
                    print("Ultimate stress exceeded at y={}, x={}, z={}; Circular section; stress={}, ult_stress={}"
                          .format(frameloc, x, z, maxstress, section.properties.material.sigma_y))
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
                            print("Buckling occurs at y={}, x={}, z={}; Combined section; stress={}, sig_cr={}".format(frameloc, xi, zi, stress, -plate.properties.sigma_cr))

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
                            print("Tresca yield exceeded at y={}, x={}, z={}; Combined section; stress={}, max_stress={}"
                                  .format(frameloc, xi, zi, tresca, section.properties.tresca_yield))
                    if tresca > trescamax:
                        trescamax = tresca
                        critloc_tresca = (xi, zi)
                    if maxstress > plate.properties.ultimatestress:
                        flag[5] = 1
                        if printing:
                            print("Ultimate stress exceeded at y={}, x={}, z={}; Combined section; stress={}, ult_stress={}"
                                  .format(frameloc, xi, zi, stress, plate.ultimatestress))
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
                        print("Buckling of stiffener at y={}, x={}, z={}; Longerons; stress={}, sigma_cr={}".format(frameloc, x, z, stress, -stiff.properties.sigma_cr))
                if maxstress > stiff.properties.material.sigma_y:
                    flag[7] = 1
                    if printing:
                        print("Ultimate stress exceeded at y={}, x={}, z={}; Longerons; stress={}, ult_stress={}"
                              .format(frameloc, x, z, maxstress, stiff.properties.materials.sigma_y))
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

def save_fuselage(v):
    t = np.vectorize(partial(skinthickness, v.skin_t)) if v.skin_t_func is None else np.vectorize(
        partial(v.skin_t_func, v.skin_t))
    stringerargs = constructstringerarg(v.stringermod, v.stringermat)
    circstringerargs = constructstringerarg(v.circstringermod, v.circstringermat)
    longargs = constructlongeronarg(v.longeronmod, v.longeronmat)
    if v.framelocs is None:
        framelocs = np.concatenate(
            (np.arange(22, 60, 1) * 0.1, np.array([5.95, 6.0, 6.05]), np.arange(61, 95, 1) * 0.1))[::-1]
    else:
        framelocs = v.framelocs
    fus, framelocs = make_fuselage(t, stringerlist=v.n_stringers, stringerargs=stringerargs,
                                  materials=[v.material_regular, v.material_circular], framelocs=framelocs,
                                  n_stiff_circ=v.n_stiff_circ, circstringerargs=circstringerargs, n_longs=v.n_longs,
                                  longargs=longargs)
    v.framelocs = framelocs
    v.fuselage = fus
    return v

def size_parameters(v):
    fus = v.fuselage
    framelocs = v.framelocs
    y_origin = 1.835 if v.xcg_min is None else v.xcg_min
    bc1 = {'y': 0.0 + y_origin, 'defl_x': 0.0}
    bc2 = {'y': 0.0 + y_origin, 'defl_z': 0.0}
    bc3 = {'y': 0.0 + y_origin, 'angle_x': 0.0}
    bc4 = {'y': 0.0 + y_origin, 'angle_z': 0.0}
    bc5 = {'y': 0.0 + y_origin, 'defl_y': 0.0}
    bc6 = {'y': 0.0 + y_origin, 'angle_y': 0.0}
    bcs = [bc1, bc2, bc3, bc4, bc5, bc6]
    
    Wbat = v.W_batt
    loads = {}
    # VTO, VL, CL_h, XMAC_v, MAC_v, VvV, S_v, CD
    if v.CLh_TO is not None:
        if v.CLh_L*v.VL*v.VL < v.CLh_TO*v.VTO*v.VTO:
            clh = v.CLh_L
            vel = v.VL
        else:
            clh = v.CLh_TO
            vel = v.VTO
    else:
        clh = v.CLh_L
    h_lift = 0.5 * v.rhoTO * vel * vel * v.VhV * v.VhV * v.Sh * clh * 1.5 * v.n_ult
    v_lift = 0.5 * v.rhoTO * vel * vel * v.VhV * v.VhV * v.Sv * v.CnB*20*pi/180*1.5*v.n_ult
    loads["fhtail"] = Force(xpos=0.0, ypos=v.XMAC_h + 0.25 * v.MAC_h, zpos=0.0, xmag=0, zmag=h_lift, ymag=h_lift * 0.2)
    loads["fvtail"] = Force(xpos=0.0, ypos=v.XMAC_h + 0.25 * v.MAC_h, zpos=0.0, xmag=v_lift, ymag=v_lift * 0.2,
                            zmag=0)
    loads["battery1"] = Force(xpos=0.0, ypos=v.cockpitbulkhead + v.batteryoffset, zpos=0.0, xmag=0, ymag=0,
                              zmag=-Wbat * 0.5 * 1.5)
    loads["battery2"] = Force(xpos=0.0, ypos=v.cockpitbulkhead + v.batteryoffset + v.batterywidth, zpos=0.0, xmag=0,
                              ymag=0, zmag=-Wbat * 0.5 * 1.5)
    v.loads = loads
    yout, old_framelocs, fus, sx, sz, mx, mz, ny, shearflow = analyse_fuselage(framelocs, fus, bcs, loads, yout=v.yout)
    flag, sig_max, loc_max, sig_min, loc_min, sig_norm, shear, loc_shear, tresca, loc_tresca = \
        analyse_forces(yout, old_framelocs, fus, sx, sz, mx, mz, ny, shearflow, printing=True)
    return flag

def design_fuselage(v: NewVariables, maxiterations=20):
    if v._max_fuselage_iterations is None:
        maxits = maxiterations
    else:
        maxits = v._max_fuselage_iterations

    if maxits == 1:
        if v.fuselage is None:
            return save_fuselage(v)
        else:
            return v
    else:
        v = save_fuselage(v)
        flags = size_parameters(v)
        parameters = {"framesamount": (3, 12), "skin_t": (0.5, 10.0), "stringermod": (1.0, 5.0),
                      "circstringermod": (1.0, 5.0), "longeronmod": (1.0, 6.0), "n_stiff": (0, 12),
                      "n_stiff_circ": (0, 20), "batteryoffset": (0.05, 0.4), "batterywidth": (0.15, 0.5)}

        if sum(flags) > 0: # Check for not failing under load
            counter = 0
            while sum(flags) > 0 and counter < maxits:
                counter += 1
                skinchange = False
                stiffchange = False
                longchange = False
                if flags[0] == 1:
                    v.skin_t, skinchange = change_parameter(parameters["skin_t"], v.skin_t, True)
                if flags[1] == 1:
                    if not skinchange:
                        v.skin_t, skinchange = change_parameter(parameters["skin_t"], v.skin_t, True)
                if flags[2] == 1:
                    v.n_stiff_circ, _ = change_parameter(parameters["n_stiff_circ"], v.n_stiff_circ, True)
                if flags[3] == 1:
                    if not skinchange:
                        v.skin_t, skinchange = change_parameter(parameters["skin_t"], v.skin_t, True)
                        print("Changing skin, done={}, value={}".format(skinchange, v.skin_t))
                    if not skinchange: # Double, in case minimum or maximum is already reached.
                        v.n_stiff, stiffchange = change_parameter(parameters["n_stiff"], v.n_stiff, True)
                        print("Changing stiffeners, done={}, amount={}".format(skinchange, v.n_stiff))
                if flags[4] == 1:
                    v.n_stiff, stiffchange = change_parameter(parameters["n_stiff"], v.n_stiff, True)
                if flags[5] == 1:
                    if not stiffchange:
                        v.n_stiff, stiffchange = change_parameter(parameters["n_stiff"], v.n_stiff, True)
                if flags[6] == 1:
                    v.stringermod, stiffchange = change_parameter(parameters["stringermod"], v.stringermod, True)
                    v.longeronmod, longchange = change_parameter(parameters["longeronmod"], v.longeronmod, True)
                if flags[7] == 1:
                    if not stiffchange:
                        v.stringermod, _ = change_parameter(parameters["stringermod"], v.stringermod, True)
                    if not longchange:
                        v.longeronmod, _ = change_parameter(parameters["longeronmod"], v.longeronmod, True)
                v = save_fuselage(v)
                flags = size_parameters(v)
            if sum(flags) > 0:
                print(
                    "Structures convergence failed. From now on, fuselage will not be redesigned, only weight will be recalculated following the input. Different input parameters should be chosen, or more iterations should be allowed.")
                v._max_fuselage_iterations = 1
                return v
        # Minimize mass within constraints
        vnew = copy.deepcopy(v)
        increaselist = [True, False]
        increase_results = []
        basemass = v.Wfus_aft
        for parameter, minmax in parameters.items():
            baseval = getattr(vnew, parameter)
            newmasses = []
            no_failure = []
            for increase in increaselist:
                newval, changed = change_parameter(minmax, baseval, increase)
                if changed:
                    setattr(vnew, parameter, newval)
                    vnew = save_fuselage(vnew)
                    flags = size_parameters(v)
                    if sum(flags) == 0:
                        no_failure.append(True)
                    else:
                        no_failure.append(False)
                    newmasses.append(vnew.Wfus_aft)
                else:
                    no_failure.append(None)
                setattr(vnew, parameter, baseval)
            if no_failure[0] == True and newmasses[0] < basemass:
                if no_failure[1] == True and newmasses[1] < basemass:
                    if newmasses[0] > newmasses[1]:
                        increase_results.append(True)
                        setattr(vnew, parameter, change_parameter(minmax, baseval, True)[0])
                    else:
                        increase_results.append(False)
                        setattr(vnew, parameter, change_parameter(minmax, baseval, False)[0])
                else:
                    increase_results.append(True)
                    setattr(vnew, parameter, change_parameter(minmax, baseval, True)[0])
            elif no_failure[1] == True and newmasses[1] < basemass:
                increase_results.append(False)
                setattr(vnew, parameter, change_parameter(minmax, baseval, False)[0])
            else:
                increase_results.append(None)
        return save_fuselage(vnew)






if __name__ == "__main__":
    from materials import materials
    t = np.vectorize(partial(skinthickness, 1))
    Wbat = 256*9.81
    mats = materials()
    stringerargs = constructstringerarg(1.0, mats['alu2024'])
    circstringerargs = constructstringerarg(1.0, mats['alu2024'])
    longargs = constructlongeronarg(4.0, mats['carbonfibre'])
    n_stringers = np.ones(8)*5
    n_stiff_circ = 6
    n_longs = 8

    fus, framelocs = make_fuselage(t, stringerlist=n_stringers, stringerargs=stringerargs, materials=[mats['alu7075'], mats['carbonfibre']], n_stiff_circ=n_stiff_circ, circstringerargs=circstringerargs, n_longs=n_longs, longargs=longargs)

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

    yout, old_framelocs, fus, sx, sz, mx, mz, ny, shearflow = analyse_fuselage(framelocs, fus, bcs, loads)

    flag, sig_max, loc_max, sig_min, loc_min, sig_norm, shear, loc_shear, tresca, loc_tresca = analyse_forces(yout, old_framelocs, fus, sx, sz, mx, mz, ny, shearflow, printing=True)
    # print(fus.mass)
    # print(np.where(np.abs(np.asarray(old_framelocs)-6.0)<0.15))

    secs = np.asarray(fus.sections)
    secs = secs[np.where(np.abs(np.asarray(old_framelocs)-6.0)<0.25)]
    # for sec in secs[1:3]:
    for sec in secs[1:-4]:
        sec.properties.plot_section(show=True)
    plt.plot(old_framelocs, fus.Ixx)
    # plt.plot(yout, mx)
    # plt.plot(yout, sz)
    # plt.plot(yout, mz)
    # plt.plot(yout, sx)
    # plt.plot(old_framelocs, sig_max)
    plt.show()

