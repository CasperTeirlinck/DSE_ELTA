"""
This file calculates the total amount of energy required by the aircraft --> for battery sizing
Author: Matthijs van Ede
"""
from math import sqrt, pi

def E_cruise_range(variables, range_m):
    """
    Returning the required energy for cruising based on range requirement
    :param variables:
    :param range_m: range in meters
    :return: energy required for cruising based on range [J]
    """
    CL, CD = sqrt(pi * variables.A * variables.e * variables.CD0), 2 * variables.CD0
    E = (range_m * variables.WTO_range)/(CL / CD)
    return E


def E_cruise_endurance(variables, endurance_s):
    """
    :param variables:
    :param endurance_s: endurance in seconds
    :return: the energy for cruising based on range [J]
    """
    CL, CD = sqrt(3 * pi * variables.CD0 * variables.e * variables.A), 4 * variables.CD0
    V_loiter = CL/CD / sqrt(0.5 * variables.rho0 * CL**3 / CD**2 / variables.WS)
    E = endurance_s*(1/2)*variables.rho0*(V_loiter**3)*variables.S*CD
    return E


def E_taxi(P_max, taxipower = 0.05, taxitime = 330):
    """
    :param P_max: maximum power Rik knows this number
    :param taxipower: fraction of maximum power, i.e. power setting
    :param taxitime: time to taxi in seconds
    :return: energy required for single taxi phase in [J]
    """
    E = P_max*taxipower*taxitime
    return E


def E_climb(P_max, climbpower = 1.0, climbtime = (914.4/2 + 30)):
    """
    :param P_max: maximum power Rik knows this number
    :param climbpower: fraction of maximum power, i.e. power setting
    :param climbtime: time to climb in seconds; 3000ft/2m/s + 30 seconds reserve
    :return: energy required for single climb phase in [J]
    """
    E = P_max*climbpower*climbtime
    return E


def energy_engine(variables, taxi=2, climb=1, endurance_s=9000., range_m=250000.):
    """
    This function brings all functions described above together
    :param variables: variables class, containing all the variables suchs as A,e,CDo etc.
    :param taxi: amount of taxi phases in flight profile, nominal = 2
    :param climb: amount of full power climb phases in flight profile, nominal = 1
    :param endurance_s: endurance required in seconds, set to 0 if calculating for range
    :param range_m: range required in meters, set to 0 if calculating for endurance
    :return: total energy required for the engine only based on the flight profile in [J]
    """

    # Simple check that checks whether we ar overdesigning
    if (endurance_s != 0) and (range_m != 0):
        print("Whatch out you are calucalting the BMS for both range and endurance, check pp_energy_calculation.py")

    total_energy = 0

    # Taxi phase
    if taxi != 0:
        for _ in range(taxi):
            total_energy += E_taxi(variables.P_max)

    # Climb phase
    if climb != 0:
        for _ in range(climb):
            total_energy += E_climb(variables.P_max)

    # Cruise phase
    if endurance_s != 0:
        total_energy += E_cruise_endurance(variables, endurance_s)
    if range_m != 0:
        total_energy += E_cruise_range(variables, range_m)

    # Returning total energy required
    return total_energy/variables.eff_tot_prop


## Calculating the energy required by the avionics
def energy_avionics(variables):
    """
    This function calculates the required energy for the avionics
    :return: Energy required for the avionics, note that this is independent of the take-off weight
    """
    # Power [W] requirements by avionics
    # See spreadsheet "List of EBS items"
    pfd_mfd = 250
    autopilot = 14
    intercom = 8.4
    portible_instrument_panel = 12
    engine_computer = 6.1
    lights = 28

    total_power = pfd_mfd + autopilot + intercom + portible_instrument_panel + engine_computer + lights
    total_time = 2.5*3600 + 2*330 + (914.4/2 + 30) # 2.5h endurance + 2*taxi phase + climb phase
    avionics_energy = total_power * total_time
    return avionics_energy/variables.eff_batt


if __name__ ==  "__main__":
    ## Testing
    # THESE ARE TEST VARIABLES!
    class Test_variables_pp:
        def __init__(self):
            # Pipistrel alpha validation
            self.CD0 = 0.0280
            self.A = 11.6
            self.e = 0.83
            self.rho0 = 1.225
            self.WP = 0.131
            self.WS = 434
            self.Especif_bat = 900000
            self.rho_bat = 500 * 3600
            self.eff_tot_prop = 0.95* 0.88
            self.eff_batt = 0.95
            self.g0 = 9.80665
            self.WTO_endurance = 550*9.81
            self.WTO_range = self.WTO_endurance
            self.range_m = 0.0
            self.endurance_s = 1.5*3600
            self.P_max = 60*1000
            self.S = 9.51

    test_v = Test_variables_pp()
    print(energy_engine(test_v, endurance_s = test_v.endurance_s, range_m= 0.0 )/3600/1000, "[kWh]")
    print("Pipistrel website says 21 kWh :D :D :D")



