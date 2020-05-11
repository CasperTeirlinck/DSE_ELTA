import numpy as np
from math import sqrt, pi

range_m = 250000
endurance_s = 9000
g0 = 9.80665 # m/s2
rho0 = 1.225 # kg/m3
hcr = 914.4 # m
Vy = 2 # m/s
taxitime = 330 # s
taxipower = 0.05 # -
climbtime = hcr/Vy + 30 # s, added time for take-off
climbpower = 1 # -
# Descent assumed to be included in endurance and range


def range_battery(CD0, e, A, WTO, Especific, eff_total, powerloading):
    # Returns the battery weight in N. WTO required to be in N.
    CL = sqrt(pi * A * e * CD0)
    CD = 2 * CD0
    total_energy_range = range_m / (CL/CD)
    total_energy_range += taxipower / powerloading * taxitime * 2 # Taxi before TO and after L
    total_energy_range += climbpower / powerloading * climbtime
    return total_energy_range / (Especific * eff_total) * WTO * g0


def endurance_battery(CD0, e, A, WTO, Especific, eff_total, powerloading, wingloading):
    # Returns the battery weight in N. WTO required to be in N.
    CD = 4 * CD0
    CL = sqrt(3 * pi * CD0 * e * A)
    V_loiter = CL/CD / sqrt(0.5 * rho0 * CL**3 / CD**2 / wingloading)
    total_energy_loiter = endurance_s * 0.5 * rho0 * V_loiter**3 * CL * CD / wingloading
    total_energy_loiter += taxipower / powerloading * taxitime * 2 # Taxi before TO and after L
    total_energy_loiter += climbpower / powerloading * climbtime
    return total_energy_loiter / (Especific * eff_total) * WTO * g0



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
    powerloading = 0.1 # N/W

    print(range_battery(CD0, e, A, WTO, Especific, eff_total, powerloading)/g0)
    print(endurance_battery(CD0, e, A, WTO, Especific, eff_total, powerloading, wingloading)/g0)