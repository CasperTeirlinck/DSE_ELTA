"""
This file calculates the total amount of energy required by the aircraft --> for battery sizing
Author: Matthijs van Ede
"""
from math import sqrt, pi

def range_power_per_WTO(variables, range_m):
    CL, CD = sqrt(pi * variables.A * variables.e * variables.CD0), 2 * variables.CD0
    return range_m / (CL/CD)


def endurance_power_per_WTO(variables, endurance_s):
    rho0 = 1.225
    CL, CD = sqrt(3 * pi * variables.CD0 * variables.e * variables.A), 4 * variables.CD0
    V_loiter = CL/CD / sqrt(0.5 * rho0 * CL**3 / CD**2 / variables.WS)
    return endurance_s * 0.5 * rho0 * V_loiter**3 * CD / variables.WS


def taxi_power_per_WTO(WP, taxipower = 0.05, taxitime = 330):
    """
    :param WP: Design point W/P value
    :param taxipower: fraction of maximum power, i.e. power setting
    :param taxitime: time to taxi in seconds
    :return: energy required for single taxi phase in J
    """
    return taxipower * taxitime / WP


def climb_power_per_WTO(WP, climbpower = 1.0, climbtime = (914.4/2 + 30)):
    """
    :param WP: Design point W/P value
    :param climbpower: fraction of maximum power, i.e. power setting
    :param climbtime: time to climb in seconds; 3000ft/2m/s + 30 seconds reserve
    :return: energy required for single climb phase in J
    """
    return climbpower * climbtime / WP


def bms_engine(variables, taxi=2, climb=1, endurance_s=9000., range_m=250000.):
    """
    :param variables: variables class, containing all the variables suchs as A,e,CDo etc.
    :param taxi: amount of taxi phases in flight profile, nominal = 2
    :param climb: amount of full power climb phases in flight profile, nominal = 1
    :param endurance_s: endurance required in seconds, set to 0 if calculating for range
    :param range_m: range required in meters, set to 0 if calculating for endurance
    :return: batter mass fraction to ONLY power the engine
    """

    # Simple check that checks whether we ar overdesigning
    if (endurance_s != 0) and (range_m != 0):
        print("Whatch out you are calucalting the BMS for both range and endurance, check pp_energy_calculation.py")

    total_energy = 0

    # Taxi phase
    if taxi != 0:
        for _ in range(taxi):
            total_energy += taxi_power_per_WTO(variables.WP)

    # Climb phase
    if climb != 0:
        for _ in range(climb):
            total_energy += climb_power_per_WTO(variables.WP)

    # Cruise phase
    if endurance_s != 0:
        total_energy += endurance_power_per_WTO(variables, endurance_s)
    if range_m != 0:
        total_energy += range_power_per_WTO(variables, range_m)

    # Returning battery mass fraction
    return total_energy * g0 / (variables.Especif_bat * variables.eff_tot_prop)


def bms_avionics():
    pass


def bms_total():
    pass

