import numpy as np

""" QUANTITY DISCOUNT FACTOR """ # (learning curve)

F_exp = 0.9 # [-] experience effectiveness ref: gadesign
N = 100 # [-] nr of units produced

QDF = F_exp**(1.4427 * np.log(N))

""" MAN_HOURS """ # DAPCA-IV Method ref GADesign

# ENGINEERING
W_af = # [kg] airframe weight
N = 100 # [-] nr of planned units produced over 5y span
V_h = # [ktas] max level airspeed
F_cert = (0.67 + 1)/2 # [-] 0.67 for LSA, 1 for CS23 ==> take avg
f_comp = 0 # [-] composite fraction
F_comp = 1 + f_comp
F_cf = # [-] 1.03 for compl. flaps, 1 for simple flaps
F_press = 1 # [-] 1.03 for pressurized cabin, else 1

H_eng = 0.0396 * W_af**0.791 * V_h**1.526 * N**0.183 * F_cert * F_cf * F_comp * F_press

# TOOLING
Q_m # [#a/c per month] estimated prod rate (= N/60 for 60 months/5years)
F_taper = 1 # [-] 0.95 for no taper, 1 for tapered
F_cf = # [-] 1.02 for compl. flaps, 1 for simple flaps
F_press = 1 # [-] 1.01 for pressurized, else 1

H_tool = 1.0032 * W_af**0.764 * V_h**0.899 * N**0.178 * Q_m**0.066 * F_taper * F_cf * F_comp * F_press

# MANUFACTURING
F_cert = (0.75 + 1)/2 # [-] 0.75 for LSA, 1 fpr CS23 => avg
F_cf = # [-] 1.01 for compl. flaps, 1 for simple flaps
F_comp = 1 + 0.25*f_comp

H_mfg = 9.6613 * W_af**0.74 * V_h**0.543 * N**0.524 * F_cert * F_cf * F_comp

""" COSTS """

# ENGINEERING
R_eng = 100 # [curr/h] 
CPI_2012 = 1 # [-] consumer price index wrt. 2012

C_eng = 2.0969 * H_eng * R_eng * CPI_2012

# DEVELOPMENT SUPPORT (overhead)
N_p = 1 # [-] nr of prototypes
F_cert = (0.5 + 1)/2 # [-] 0.5 for LSA, 1 for CS23 +> agv
F_cf = # [-] 1.01 for compl. flaps, 1 for simple flaps
F_comp = 1 + 0.5*f_comp
F_press = 1 # [-] 1.03 for pressurized, else 1

C_dev = 0.06458 * W_af**0.873 * V_h**1.89 * N_p**0.346 * CPI_2012 * F_cert * F_cf * F_comp * F_press

# FLIGHT TESTS
F_cert = (10 + 1)/2 # [-] 10 for LSA, 1 for CS23

C_ft = 0.009646 * W_af**1.16 * V_h**1.3718 * N_p**1.281 * CPI_2012 * F_cert

# TOOLING
R_tool = 60 # [curr/h]

C_tool = 2.0969 * H_tool * R_tool * CPI_2012

# MANUFACTURING
R_mfg = 50 # [curr/h]

C_mfg = 2.0969 * H_mfg * R_mfg * CPI_2012

# QC
F_cert = (0.5 + 1)/2 # [-] 0.5 for LSA, 1 for CS23
F_comp = 1 + 0.5*f_comp

C_qc = 0.13 * C_mfg * F_cert * F_comp

# MATERIALS
F_cert = (0.75 + 1)/2 # [-] 0.75 fro LSA, 1 for CS23
F_cf = # [-] 1.02 for compl. flaps, 1 for simple flaps
F_press = 1 # [-] 1.01 for press, else 1

C_mat = 24.896 * W_af**0.689 * V_h**0.624 * N**0.792 * CPI_2012 * F_cert * F_cf * F_press

# CERTIFICATION
C_cert = C_eng + C_dev + C_ft + C_tool

""" PROPULSION """
C_prop = 0

""" AVIONICS """ # ref: TE package cost estimation
C_te = 10395 # concept 1
# C_te = 29128 # concept 1
# C_te = 20482 # concept 1
# C_te = 8778 # concept 1
# C_te = 4543 # concept 1

""" TOTAL COST """
C_lg = -7500 # fixed lg discout ref: gadesign
C_tot = C_cert + C_mfg + C_qc + C_mat + C_prop + C_te + C_lg