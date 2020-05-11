import numpy as np
from math import sqrt, pi

range_m = 250000
endurance_s = 9000
g0 = 9.80665 # m/s2
rho0 = 1.225 # kg/m3

def range_battery(CD0, e, A, WTO, Especific, eff_total):
    # Returns the battery mass. If WTO in kg, then in kg, if WTO in N, then in N.
    # TODO: Add taxi, climb and descent battery mass
    CL = sqrt(pi * A * e * CD0)
    CD = 2 * CD0
    LoverD = CL/CD
    battery_weight_fraction = range_m / LoverD * g0 / Especific / eff_total
    return WTO * battery_weight_fraction

def endurance_battery(CD0, e, A, WTO, Especific, eff_total, wingloading):
    # Returns the battery weight in N. WTO required to be in N.
    # TODO: Add taxi, climb and descent battery mass
    S = WTO / wingloading
    CD = 4 * CD0
    CL = sqrt(3 * pi * CD0 * e * A)
    V_loiter = sqrt(WTO / (0.5 * rho0 * S * CL))
    total_energy_loiter = endurance_s * 0.5 * rho0 * V_loiter**3 * S * CD
    return total_energy_loiter / (Especific * eff_total) * g0



if __name__ == "__main__":
    CD0 = 0.011 # -
    e = 0.83 # -
    A = 8 # -
    WTO = 542.4761973 * g0 # N
    Especific = 900000 # J/kg Li-ion from Maarten
    eff_motor = 0.95 # -, PM motor from Maarten
    eff_battery = 0.8 # -, lower estimate for Li-ion
    eff_total = eff_motor * eff_battery
    wingloading = 592.289 # N/m2

    print(range_battery(CD0, e, A, WTO, Especific, eff_total)/g0)
    print(endurance_battery(CD0, e, A, WTO, Especific, eff_total, wingloading)/g0)