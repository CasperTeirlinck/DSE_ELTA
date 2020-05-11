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


def range_power_per_WTO(CD0, e, A, range_m):
    CL, CD = sqrt(pi * A * e * CD0), 2 * CD0
    return range_m / (CL/CD)


def endurance_power_per_WTO(CD0, e, A, endurance_s, wingloading):
    CL, CD = sqrt(3 * pi * CD0 * e * A), 4 * CD0
    V_loiter = CL/CD / sqrt(0.5 * rho0 * CL**3 / CD**2 / wingloading)
    return endurance_s * 0.5 * rho0 * V_loiter**3 * CL * CD / wingloading


def taxi_power_per_WTO(powerloading):
    return taxipower * taxitime / powerloading


def climb_power_per_WTO(powerloading):
    return climbpower * climbtime / powerloading


def flight_profile_energy_per_WTO(CD0, e, A, Especific, eff_total, powerloading, wingloading, taxi=2, climb=1,
                                  endurance_s=9000, range_m=250000):
    total_energy = 0
    if taxi != 0:
        for _ in range(taxi):
            total_energy += taxi_power_per_WTO(powerloading)
    if climb != 0:
        for _ in range(climb):
            total_energy += climb_power_per_WTO(powerloading)
    if endurance_s != 0:
        total_energy += endurance_power_per_WTO(CD0, e, A, endurance_s, wingloading)
    if range_m != 0:
        total_energy += range_power_per_WTO(CD0, e, A, range_m)
    return total_energy * g0 / (Especific * eff_total)

if __name__ == "__main__":
    CD0 = 0.011 # -
    e = 0.83 # -
    A = 8 # -
    WTO = 730 * g0 # N
    Especific = 900000 # J/kg Li-ion from Maarten
    eff_motor = 0.95 # -, PM motor from Maarten
    eff_battery = 0.8 # -, lower estimate for Li-ion
    eff_total = eff_motor * eff_battery
    wingloading = 592.289 # N/m2
    powerloading = 0.1 # N/W

    print(flight_profile_energy_per_WTO(CD0, e, A, Especific, eff_total, powerloading, wingloading, range_m=0)*WTO/g0)