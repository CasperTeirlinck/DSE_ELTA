import numpy as np

n_engines = 1
if n_engines == 1:
    a = -0.25
    b = 1.14
else:
    a = -0.90
    b = 1.32
C1 = -0.2
C2 = 0.08
C3 = 0.05
C4 = -0.05
C5 = 0.27

hp_to_W = 745.699872
kts_to_fps = 1.68780986
# kts_to_fps = 0.514444444 # kts to m/s
kg_to_pound = 2.20462262
N_to_lbf = 1/4.44822162
g0 = 9.80665


def oew_fraction(WTO, A, powerloading, wingloading, Vmax_kts):
    print(1/(1/powerloading*hp_to_W/kg_to_pound))
    return a + b*(WTO*N_to_lbf)**C1*A**C2*(1/powerloading*hp_to_W/N_to_lbf)**C3*(N_to_lbf*wingloading/0.09290304)**C4*(Vmax_kts*kts_to_fps)**C5

def ductweight(variables):
    # AR_duct = 5.0  D_fan/l_fan
    V_duct = np.pi*(variables.D_fan/variables.AR_duct)*(((variables.D_fan/2)+0.005)**2 - (variables.D_fan/2)**2)
    return (V_duct*variables.rho_ductmaterial)*1000


if __name__ == "__main__":
    WTO = 750*g0
    A = 12
    powerloading = 0.1218
    wingloading = 465
    Vmax_kts = 120

    #print(oew_fraction(WTO, A, powerloading, wingloading, Vmax_kts))

    print("mass duct:", ductweight(2))