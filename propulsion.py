import numpy as np
from math import sqrt, pi
from variables import *

range_m = 250000
endurance_s = 9000
# endurance_s *= 2/2.5
g0 = 9.80665 # m/s2
rho0 = 1.225 # kg/m3
hcr = 914.4 # m
Vy = 2 # m/s
taxitime = 330 # s
taxipower = 0.05 # -
climbtime = hcr/Vy + 30 # s, added time for take-off
climbpower = 1 # -
# Descent assumed to be included in endurance and range


def range_power_per_WTO(variables, range_m):
    CL, CD = sqrt(pi * variables.A * variables.e * variables.CD0), 2 * variables.CD0
    return range_m / (CL/CD)


def endurance_power_per_WTO(variables, endurance_s):
    CL, CD = sqrt(3 * pi * variables.CD0 * variables.e * variables.A), 4 * variables.CD0
    V_loiter = CL/CD / sqrt(0.5 * rho0 * CL**3 / CD**2 / variables.WS)
    return endurance_s * 0.5 * rho0 * V_loiter**3 * CD / variables.WS


def taxi_power_per_WTO(WP):
    return taxipower * taxitime / WP


def climb_power_per_WTO(WP):
    return climbpower * climbtime / WP


def flight_profile_energy_per_WTO(variables: CurrentVariables, taxi=2, climb=1,
                                  endurance_s=9000., range_m=250000.):
    total_energy = 0
    if taxi != 0:
        for _ in range(taxi):
            total_energy += taxi_power_per_WTO(variables.WP)
    if climb != 0:
        for _ in range(climb):
            total_energy += climb_power_per_WTO(variables.WP)
    if endurance_s != 0:
        total_energy += endurance_power_per_WTO(variables, endurance_s)
    if range_m != 0:
        total_energy += range_power_per_WTO(variables, range_m)
    return total_energy * g0 / (variables.Especif_bat * variables.eff_tot_prop)


def motor_mass_volume(variables):
    """
    NOTE IT SETS THE MOTOR WEIGHT IN NEWTONS IN THE VARIABLE CLASS!!!
    :return: The mass of the electric motor in kg
    """
    Pmax = (variables.WTO / variables.WP) / 1000    # kW
    motor_mass = Pmax/variables.motor_spec_mass     # kg
    motor_volume = Pmax/variables.motor_spec_volume # L

    # Setting the motor weight in variables class
    variables.Wmotor = motor_mass * 9.80665

    return variables





if __name__ == "__main__":
    # Cfes = [0.0055, 0.0045]
    # Cfe = Cfes[n_engines-1] # -
    # SwetS = 3.7 # -
    # variables.CD0 = Cfe * SwetS # -
    # # CD0 = 0.015
    # e = 0.83 # -
    # A = 12 # -0.2677
    # WTO_target = 750 * g0 # N
    # Especific = 900000 # J/kg Li-ion from Maarten
    # eff_motor = 0.95 # -, PM motor from Maarten
    # eff_battery = 0.8 # -, lower estimate for Li-ion
    # eff_total = eff_motor * eff_battery
    # wingloading = 465 # N/m2
    # powerloading = 0.1218 # N/W
    # m_PL = 200 # kg
    # rho_battery = 1800000 # J/L

    n_engines = 1
    variables = CurrentVariables(n_engines=n_engines)

    # bat_range = flight_profile_energy_per_WTO(variables, range_m=variables.range_m, endurance_s=0, climb=1)
    bat_endurance = flight_profile_energy_per_WTO(variables, range_m=0, endurance_s=variables.endurance_s*0.5, climb=1)
    Wbat = bat_endurance*variables.WTO
    # Pmax = (variables.WTO / variables.WP) / 1000  # kW
    # E_bat = (max(bat_endurance, bat_range) * variables.WTO * (variables.Especif_bat/9.81))/1000
    #
    # print("start\n===================================================")
    # print("Num engines            =  ", n_engines)
    # print("Pmax [kW]              =  ", Pmax)
    # print("E battery [kJ]         =  ", E_bat)
    # print("Wbat/WTO range         =  ", round(bat_range, 4))
    # print("Wbat/WTO endurance     =  ", round(bat_endurance, 4))
    # print("Battery volume [L]     =  ", (E_bat*1000)/variables.rho_bat)
    # print("Battery capacity [kWh] =  ", (max(bat_endurance, bat_range) * (variables.WTO/9.81) * 250)/1000)
    # print("Motor weight [kg]      =  ", Pmax/2.5)
    # print("Motor volume [L]       =  ", Pmax/7)
    # print("Motor*bat efficiency   =  ", variables.eff_tot_prop)
    # print("Testing the motor size func")
    # print("old:", variables.Wmotor)
    # motor_mass_volume(variables)
    # print("new:", variables.Wmotor, "[N] =", variables.Wmotor/9.80665, "[kg]")
    # print("===================================================\nend")






