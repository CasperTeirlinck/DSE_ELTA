"""
Main battery sizing file
outputting battery mass, volume and composition
Author: Matthijs
"""
import pp_energy_calculation as ppe
import batteries as bat


# This class describes all the variables required for the battery sizing, the values are OLD values from the midterm report.
class Variables_battery:
    def __init__(self):
        # For the energy required
        self.CD0 = 0.0250
        self.A = 10.1
        self.e = 0.85
        self.rho0 = 1.225
        self.WP = 0.122
        self.WS = 508
        self.Especif_bat = 900000
        self.rho_bat = 500 * 3600
        self.eff_tot_prop = 0.95*0.88
        self.eff_batt = 0.95
        self.g0 = 9.80665
        self.WTO_endurance = 809.95*9.81
        self.WTO_range = self.WTO_endurance  # THIS IS OFCOURSE WRONG BUT I DON'T HAVE A BETTER VALUE FOR NOW
        self.range_m = 250000
        self.endurance_s = 2.5 * 3600
        self.P_max = self.WTO_endurance / self.WP
        self.S = self.WTO_endurance / self.WS

        # For the battery composition
        self.V_req_bat1 = 400    # Voltage required for engine subsytem [V]
        self.V_req_bat2 = 28     # Voltage required for avionics subsytem [V]
        self.I_req_bat1 = 189.75 # Current required for engine subsytem [A]
        self.I_req_bat2 = 6.34   # Current required for avionics subsytem [A]
        self.DOD = 90            # Depth of discharge in %


# Function to print results
def bat_print_func(N_bat1, N_bat2,  E_engine, E_avionics, m_bat1 = 0.0, V_bat1 = 0.0, m_bat2 = 0.0, V_bat2 = 0.0, E_prod1 = 0.0, E_prod2 = 0.0):
    print("Total energy:", (E_avionics + E_engine)/3600/1000, "[kWh] "
          "\nMidterm calc: Total energy/(E_specific * DOD (=90%)):", (E_avionics + E_engine)/(900000*0.90), "[kg]\n"
          "\n------Battery 1------\n"
          "Number of cells in series:", N_bat1[0], "\nNumber of cells in parallel:", N_bat1[1], "\nTotal number of cells:", N_bat1[2],
          "Energy required :", E_engine/3600/1000, "[kWh]\n"
          "Energy produced :", E_prod1, "[kWh]\n"
          "Battery 1 mass cells  :", m_bat1, "[kg]\n"
          "Battery 1 volume:", V_bat1, "[L]\n"
          "\n------Battery 2------\n"
          "Number of cells in series:", N_bat2[0], "\nNumber of cells in parallel:", N_bat2[1], "\nTotal number of cells:", N_bat2[2],
          "\nEnergy required :", E_avionics/3600/1000, "[kWh]\n"
          "Energy produced :", E_prod2, "[kWh]\n"
          "Battery 2 mass cells  :", m_bat2, "[kg]\n"
          "Battery 2 volume:", V_bat2, "[L]\n"
          "\nTotal mass:", m_bat1+m_bat2, "* 1.15 =", (m_bat1+m_bat2)*1.15, "[kg]\n")
    return None



def main_bat(variables):
    E_engine_range = ppe.energy_engine(variables, endurance_s = 0.0, range_m = 250000.0) # Energy for engine in [J] for range req.
    E_engine_endur = ppe.energy_engine(variables, endurance_s = 9000.0, range_m = 0.0) # Energy for engine in [J] for endurance req.
    E_engine = max(E_engine_range, E_engine_endur) # Energy for engine in [J]
    E_avionics = ppe.energy_avionics(variables) # Energy for avionics in [J]
    # TODO: add energy for cooling system

    # Just for information which requirement requires the most energy
    if E_engine == E_engine_endur: print("Enudance uses most energy")
    else: print("Range uses most energy")

    # Total energy required from the batteries
    E_total = E_engine + E_avionics

    # Sizing the battery
    # Bat 1 = the battery for the engine subsystem
    # Bat 2 = the battery for the avioncis --> see Electrical block diagram
    # Different cell types
    cell1 = "Sam_21700_50E"
    cell2 = "Sam_21700_40T"
    cell3 = "Pan_18650B"

    # Data belonging to each cell type
    cell_data = bat.cell_type(name = cell1)

    # Calculating battery composition (# of cells in parallel and series)
    bat_1 = bat.num_of_cells(cell_data, E_req = E_engine/3600, V_req = variables.V_req_bat1, I_req = variables.I_req_bat1, DoD = variables.DOD)
    bat_2 = bat.num_of_cells(cell_data, E_req=E_avionics/3600, V_req = variables.V_req_bat2, I_req = variables.I_req_bat2, DoD = variables.DOD)

    # Calculating mass and volume based on cells only
    m_bat_1 = bat_1[2] * cell_data.mass
    m_bat_2 = bat_2[2] * cell_data.mass
    v_bat_1 = bat_1[2] * cell_data.volume * 1000 # from m^3 to L
    v_bat_2 = bat_2[2] * cell_data.volume * 1000 # from m^3 to L

    # Calculating the actual energy produced due to rounding in the cell composition
    E_bat1_produced = ((variables.DOD/100)*bat_1[1] * cell_data.C_Ah * bat_1[0]*cell_data.V_nom)/1000
    E_bat2_produced = ((variables.DOD/100)*bat_2[1] * cell_data.C_Ah * bat_2[0]*cell_data.V_nom)/1000

    # Print function to print the outcomes
    bat_print_func(bat_1, bat_2, E_engine, E_avionics, m_bat1 = m_bat_1, V_bat1 = v_bat_1, m_bat2 = m_bat_2, V_bat2 = v_bat_2, E_prod1 = E_bat1_produced, E_prod2 = E_bat2_produced)
    return None


# Testing block
if __name__ == "__main__":
    test_v = Variables_battery()
    main_bat(test_v)