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
        self.CD0 = 0.0280
        self.A = 12
        self.e = 0.83
        self.rho0 = 1.225
        self.WP = 0.131
        self.WS = 434
        self.Especif_bat = 900000
        self.rho_bat = 500 * 3600
        self.eff_tot_prop = 0.72
        self.g0 = 9.80665
        self.WTO_endurance = 9978.38
        self.WTO_range = self.WTO_endurance # THIS IS OFCOURSE WRONG BUT I DON'T HAVE A BETTER VALUE FOR NOW
        self.range_m = 250000
        self.endurance_s = 2.5 * 3600
        self.P_max = self.WTO_endurance / self.WP
        self.S = self.WTO_endurance / self.WS

        # For the battery composition
        self.V_req_bat1 = 400    # Voltage required for engine subsytem [V]
        self.V_req_bat2 = 28     # Voltage required for avionics subsytem [V]
        self.I_req_bat1 = 189.75 # Current required for engine subsytem [A]
        self.I_req_bat2 = 6.34   # Current required for avionics subsytem [A]
        self.DOD = 80            # Depth of discharge in %


# Function to print results
def bat_print_func(N_bat1, N_bat2,  E_engine, E_avionics, m_bat1 = 0.0, V_bat1 = 0.0, m_bat2 = 0.0, V_bat2 = 0.0):
    print("Total energy:", (E_avionics + E_engine)/3600/1000, "[kWh]\n")
    print("------Battery 1------")
    print("Number of cells in series:", N_bat1[0], "\nNumber of cells in parallel:", N_bat1[1], "\nTotal number of cells:", N_bat1[2])
    print("Energy required :", E_engine/3600/1000, "[kWh]")
    print("Battery 1 mass  :", m_bat1, "[kg]")
    print("Battery 1 volume:", V_bat1, "[L]\n")
    print("------Battery 2------")
    print("Number of cells in series:", N_bat2[0], "\nNumber of cells in parallel:", N_bat2[1], "\nTotal number of cells:", N_bat2[2])
    print("Energy required :", E_avionics/3600/1000, "[kWh]")
    print("Battery 2 mass  :", m_bat2, "[kg]")
    print("Battery 2 volume:", V_bat2, "[L]\n")
    return None



def main_bat(variables):
    E_engine_range = ppe.energy_engine(variables, endurance_s = 0.0, range_m = 250000.0) # Energy for engine in [J] for range req.
    E_engine_endur = ppe.energy_engine(variables, endurance_s = 9000.0, range_m = 0.0) # Energy for engine in [J] for endurance req.
    E_engine = max(E_engine_range, E_engine_endur) # Energy for engine in [J]
    E_avionics = ppe.energy_avionics() # Energy for avionics in [J]
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

    # Calculating battery composition
    bat_1 = bat.num_of_cells(cell_data, E_req = E_engine/3600, V_req = variables.V_req_bat1, I_req = variables.I_req_bat1, DoD = variables.DOD)
    bat_2 = bat.num_of_cells(cell_data, E_req=E_avionics/3600, V_req = variables.V_req_bat2, I_req = variables.I_req_bat2, DoD = variables.DOD)

    bat_print_func(bat_1, bat_2, E_engine, E_avionics)
    return None

# Testing block
if __name__ == "__main__":
    test_v = Variables_battery()
    main_bat(test_v)