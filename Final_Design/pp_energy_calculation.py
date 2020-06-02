"""
This file calculates the total amount of energy required by the aircraft --> for battery sizing
Author: Matthijs van Ede
"""
from math import sqrt, pi

## All of the bms calculations for the engine regarding the flight profile
def rangepower_per_WTO(variables, range_m):
    CL, CD = sqrt(pi * variables.A * variables.e * variables.CD0), 2 * variables.CD0
    return range_m / (CL/CD)


def endurancepower_per_WTO(variables, endurance_s):
    CL, CD = sqrt(3 * pi * variables.CD0 * variables.e * variables.A), 4 * variables.CD0
    V_loiter = CL/CD / sqrt(0.5 * variables.rho0 * CL**3 / CD**2 / variables.WS)
    return endurance_s * 0.5 * variables.rho0 * V_loiter**3 * CD / variables.WS


def taxipower_per_WTO(WP, taxipower = 0.05, taxitime = 330):
    """
    :param WP: Design point W/P value
    :param taxipower: fraction of maximum power, i.e. power setting
    :param taxitime: time to taxi in seconds
    :return: energy required for single taxi phase in J
    """
    return taxipower * taxitime / WP


def climbpower_per_WTO(WP, climbpower = 1.0, climbtime = (914.4/2 + 30)):
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
            total_energy += taxipower_per_WTO(variables.WP)

    # Climb phase
    if climb != 0:
        for _ in range(climb):
            total_energy += climbpower_per_WTO(variables.WP)

    # Cruise phase
    if endurance_s != 0:
        total_energy += endurancepower_per_WTO(variables, endurance_s)
    if range_m != 0:
        total_energy += rangepower_per_WTO(variables, range_m)

    # Returning battery mass fraction
    return total_energy * variables.g0 / (variables.Especif_bat * variables.eff_tot_prop)


## Calculating the energy required by the avionics
def energy_avionics():
    """
    :return: Energy required for the avionics, note that this is independent of the take-off weight
    """
    # Power [W] requirements by avionics
    # TODO: Check these power numbers, and add the other avionics ones available
    pfd_mfd = 250 # Garmin g1000 https://www.safeflightintl.com/downloads/g1000specsheet.pdf
    autopilot = 0 # I DUNNO YET
    intercom = 16 # I took transmit power for now https://buy.garmin.com/nl-NL/NL/p/102764#specs
    portible_instrument_panel = 12 # I took an Ipad Pro for now

    total_power = pfd_mfd + autopilot + intercom + portible_instrument_panel
    total_time = 2.5*3600 + 2*330 + (914.4/2 + 30) # 2.5h endurance + 2*taxi phase + climb phase
    avionics_energy = total_power * total_time
    return avionics_energy


## Calculating the actual weight and volume of the batteries
def battery_mass_total(variables):
    """
    :param variables: variables class with all the variables
    :return: the total battery weight, the total battery volume
    """
    # TODO: Check if Especif_bat and rho_bat are in Joules!
    # TODO: Check that a distinction is made in WTO_endurance and in WTO_range
    # TODO: Check that every weight used is in Newtons!

    # take-off weight for endurance and range flight (100 or 200 kg of payload)
    WTO_endurance = variables.WTO_endurance
    WTO_range = variables.WTO_range

    # calculation of the battery mass fractions
    bms_endurance = bms_engine(variables, taxi = 2, climb = 1, endurance_s = variables.endurance_s, range_m = 0.0)
    bms_range = bms_engine(variables, taxi = 2, climb = 1, endurance_s = 0.0, range_m = variables.range_m)

    # caculation of the battery weights in newtons
    W_bat_endurance = bms_endurance * WTO_endurance
    W_bat_range = bms_range * WTO_range

    # total battery weight in newtons
    W_bat = max(W_bat_endurance, W_bat_range) + (energy_avionics()/variables.Especif_bat)

    # total battery volume in liters
    E_bat = max(W_bat_endurance, W_bat_range)/variables.g0 * variables.Especif_bat # battery energy in J
    V_bat = E_bat/variables.rho_bat

    return W_bat, V_bat


## Testing
# THESE ARE TEST VARIABLES!
class Test_variables_pp:
    def __init__(self):
        self.CD0            = 0.0280
        self.A              = 12
        self.e              = 0.83
        self.rho0           = 1.225
        self.WP             = 0.131
        self.WS             = 434
        self.Especif_bat    = 900000
        self.rho_bat        = 500*3600
        self.eff_tot_prop   = 0.72
        self.g0             = 9.80665
        self.WTO_endurance  = 9978.38
        self.WTO_range      = self.WTO_endurance - (100*self.g0)
        self.range_m        = 250000
        self.endurance_s    = 2.5*3600

class Test_alpha_electro:
    def __init__(self):
        self.CD0            = 0.0280
        self.A              = 10.5**2/9.51
        self.e              = 0.83
        self.rho0           = 1.225
        self.WP             = (550*9.80665)/(60*1000)
        self.WS             = (550*9.80665)/(9.51)
        self.Especif_bat    = 900000
        self.rho_bat        = 500*3600
        self.eff_tot_prop   = 0.72
        self.g0             = 9.80665
        self.WTO_endurance  = 550*9.80665
        self.WTO_range      = 0
        self.range_m        = 0
        self.endurance_s    = 3600 + (20*60)


class Test_eflyer2:
    def __init__(self):
        self.CD0            = 0.0280
        self.A              = 12**2/12.0
        self.e              = 0.83
        self.rho0           = 1.225
        self.WP             = (862*9.80665)/(90*1000)
        self.WS             = (862*9.80665)/12.0
        self.Especif_bat    = 900000
        self.rho_bat        = 500*3600
        self.eff_tot_prop   = 0.72
        self.g0             = 9.80665
        self.WTO_endurance  = 862*9.80665
        self.WTO_range      = 0
        self.range_m        = 0
        self.endurance_s    = 3.5*3600


if __name__ ==  "__main__":
    test_v = Test_variables_pp()
    W, V = battery_mass_total(test_v)

    print("Our aircraft:\n------------------------------------------")
    print("Battery weight [N]   =", W)
    print("Battery weight [kg]  =", W/test_v.g0)
    print("Battery volume [L]   =", V)
    print("Battery volume [m^3] =", V*0.001)
    print("------------------------------------------")

    test_v_2 = Test_alpha_electro()
    W2, V2 = battery_mass_total(test_v_2)

    print("\nAlpha electro:\n------------------------------------------")
    print("DATA: batt mass [kg] = 114")
    print("Battery weight [N]   =", W2)
    print("Battery weight [kg]  =", W2/test_v_2.g0)
    print("Battery volume [L]   =", V2)
    print("Battery volume [m^3] =", V2*0.001)
    print("------------------------------------------")




