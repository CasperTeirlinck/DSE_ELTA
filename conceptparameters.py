import numpy as np


def conceptparameters(number):
    x_cg_pass = 0.97
    x_cg_batt = 1.96
    x_cg_fuselage = 1.965
    if number == 1:
        n_engines = 2
        wing_mounted = False
        t_tail = True
        ducted = True
        lowwing = True
        strutted_wing = False
        l_fuselage = 4.41
        d_fuselage = 1.27
        spinner_length = 0.365
        bulkhead_loc = 0.7
        cabin_length = 1.3
        tailheight = 1.15
    elif number == 2:
        n_engines = 1
        wing_mounted = False
        t_tail = False
        ducted = False
        lowwing = True
        strutted_wing = False
        l_fuselage = 4.3
        d_fuselage = 1.2
        spinner_length = 0.365
        bulkhead_loc = 1.095
        cabin_length = 1.3
        tailheight = 0.91
    elif number == 3:
        n_engines = 1
        wing_mounted = False
        t_tail = False
        ducted = False
        lowwing = True
        strutted_wing = False
        l_fuselage = 4.3
        d_fuselage = 1.2
        spinner_length = 0.365
        bulkhead_loc = 1.095
        cabin_length = 1.3
        tailheight =0.91
    elif number == 4:
        n_engines = 2
        wing_mounted = True
        t_tail = False
        ducted = False
        lowwing = False
        strutted_wing = True
        l_fuselage = 5.945
        d_fuselage = 1.2
        spinner_length = 0.365
        bulkhead_loc = 0.875
        cabin_length = 1.73
        tailheight =1.2
    elif number == 5:
        n_engines = 12
        wing_mounted = False
        t_tail = False
        ducted = True
        lowwing = False
        strutted_wing = True
        l_fuselage = 5.945
        d_fuselage = 1.2
        spinner_length = 0.365
        bulkhead_loc = 0.875
        cabin_length = 1.73
        tailheight = 1.2
    return (number, n_engines, wing_mounted, t_tail, x_cg_pass, x_cg_batt, x_cg_fuselage, ducted, lowwing, strutted_wing,
            l_fuselage, d_fuselage, spinner_length, bulkhead_loc, cabin_length, tailheight)