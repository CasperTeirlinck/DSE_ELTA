from Midterm_Design.variables import *

variables = CurrentVariables()

#""" QUANTITY DISCOUNT FACTOR """ # (learning curve)
#
#F_exp = 0.9 # [-] experience effectiveness ref: gadesign
#N = 100 # [-] nr of units produced
#
#QDF = F_exp**(1.4427 * np.log(N))
#
#


#Input parameter that influence the final trade-off:
""" PROPULSION """

#C_prop = 13500 + 5000 #eng + prop for concept 1
#C_prop = 7000 + 3500 #eng + prop for concept 2
C_prop =  9000 + 5000 #eng + prop for concept 3
#C_prop = 10500 + 6000 #eng + prop for concept 4
#C_prop = 13500 + 10000 #eng + prop for concept 5
bat = 500*51 #500 euro per kWh

""" AVIONICS """ # ref: TE package cost estimation
#C_te = 10395 # concept 1
#C_te = 29128 # concept 2
C_te = 20482 # concept 3
#C_te = 8778 # concept 4
#C_te = 4543 # concept 5


""" OTHER """
W_af = 250*2.20462262 # [lbs] airframe weight SO NOT (O)EW!!!!!!!!!!!!!!!!!!!!!!!!!
V_h = 120 # [ktas] max level airspeed
N = 650 # [-] nr of planned units produced over 5y span














""" MAN_HOURS """ # DAPCA-IV Method ref: gadesign
# ENGINEERING
#W_af = 450*2.20462262 # [lbs] airframe weight
#V_h = 110 # [ktas] max level airspeed
#N = 1100 # [-] nr of planned units produced over 5y span
F_cert = (0.67 + 1)/2 # [-] 0.67 for LSA, 1 for CS23 ==> take avg
f_comp = 0 # [-] composite fraction
F_comp = 1 + f_comp
F_cf = 1 # [-] 1.03 for compl. flaps, 1 for simple flaps
F_press = 1 # [-] 1.03 for pressurized cabin, else 1

H_eng = 0.0396 * W_af**0.791 * V_h**1.526 * N**0.183 * F_cert * F_cf * F_comp * F_press
#print(f'HOURS: eng: {H_eng}')

# TOOLING
Q_m = N/60 # [#a/c per month] estimated prod rate (= N/60 for 60 months/5years)
F_taper = 1 # [-] 0.95 for no taper, 1 for tapered
F_cf = 1 # [-] 1.02 for compl. flaps, 1 for simple flaps
F_press = 1 # [-] 1.01 for pressurized, else 1

H_tool = 1.0032 * W_af**0.764 * V_h**0.899 * N**0.178 * Q_m**0.066 * F_taper * F_cf * F_comp * F_press
#print(f'HOURS: tooling: {H_tool}')

# MANUFACTURING
F_cert = (0.75 + 1)/2 # [-] 0.75 for LSA, 1 fpr CS23 => avg
F_cf = 1 # [-] 1.01 for compl. flaps, 1 for simple flaps
F_comp = 1 + 0.25*f_comp

H_mfg = 9.6613 * W_af**0.74 * V_h**0.543 * N**0.524 * F_cert * F_cf * F_comp
#print(f'HOURS: mnfg: {H_mfg}')

""" COSTS """

# ENGINEERING
R_eng = 90 # [curr/h]
CPI_2012 = 1.1376 # [-] consumer price index 2020 wrt. 2012

C_eng = 2.0969 * H_eng * R_eng * CPI_2012
#print(f'eng: {C_eng/N}')

# DEVELOPMENT SUPPORT (overhead)
N_p = 2 # [-] nr of prototypes
F_cert =(0.5 + 1)/2 # [-] 0.5 for LSA, 1 for CS23 => agv
F_cf = 1 # [-] 1.01 for compl. flaps, 1 for simple flaps
F_comp = 1 + 0.5*f_comp
F_press = 1 # [-] 1.03 for pressurized, else 1

C_dev = 0.06458 * W_af**0.873 * V_h**1.89 * N_p**0.346 * CPI_2012 * F_cert * F_cf * F_comp * F_press
#print(f'dev: {C_dev}')

# FLIGHT TESTS
F_cert = (10 + 1)/2 # [-] 10 for LSA, 1 for CS23

C_ft = 0.009646 * W_af**1.16 * V_h**1.3718 * N_p**1.281 * CPI_2012 * F_cert
#print(f'flighttest: {C_ft}')

# TOOLING
R_tool = 60 # [curr/h]

C_tool = 2.0969 * H_tool * R_tool * CPI_2012
#print(f'tooling: {C_tool/N}')

# MANUFACTURING
R_mfg = 50 # [curr/h]

C_mfg = 2.0969 * H_mfg * R_mfg * CPI_2012
print(f'manufacturing per a/c: {C_mfg/N}')

# QC
F_cert =  (0.5 + 1)/2 # [-] 0.5 for LSA, 1 for CS23
F_comp = 1 + 0.5*f_comp

C_qc = 0.13 * C_mfg * F_cert * F_comp
print(f'QC per a/c: {C_qc/N}')

# MATERIALS
F_cert = (0.75 + 1)/2 # [-] 0.75 fro LSA, 1 for CS23
F_cf = 1 # [-] 1.02 for compl. flaps, 1 for simple flaps
F_press = 1 # [-] 1.01 for press, else 1

C_mat = 24.896 * W_af**0.689 * V_h**0.624 * N**0.792 * CPI_2012 * F_cert * F_cf * F_press
print(f'Materials/equip per a/c: {C_mat/N}')

# CERTIFICATION
C_cert = C_eng + C_dev + C_ft + C_tool
print(f'TOT Certifiction: {C_cert}')


""" TOTAL COST """
C_lg = -7500 # fixed lg discout ref: gadesign
C_tot = C_cert + C_mfg + C_qc + C_mat + N*C_prop + N*C_te + N*C_lg +N*bat

print(f'TOT unit cost to break-even in 5 years: {C_tot/N}')



'''
fixed costs: 15 000 000 certification
variable costs:
    engine:         14 000
    avionics:       20 000
    battery:        25 500
    manufacturing:  67 000
    QC:              6 500
    materials:      10 500
    +++++++++++++++++++++++++++++++++++
                   143 500 EUR per A/C
'''
