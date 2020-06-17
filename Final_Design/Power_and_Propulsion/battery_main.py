"""
Main battery sizing file
outputting battery mass, volume and composition
Author: Matthijs
"""
import pp_energy_calculation as ppe
import batteries as batc


# This class describes all the variables required for the battery sizing, the values are OLD values from the midterm report.
class Variables_battery:
    def __init__(self):
        # For the energy required
        self.CD0 = 0.0250
        self.A = 10.1
        self.e = 0.85
        self.rho0 = 1.225
        self.eff_tot_prop = 0.95*0.88
        self.eff_batt = 0.95
        self.g0 = 9.80665
        self.WTO_endurance = 809.95*9.81
        self.WTO_range = self.WTO_endurance  # THIS IS OFCOURSE WRONG BUT I DON'T HAVE A BETTER VALUE FOR NOW
        self.range_m = 250000
        self.endurance_s = 2.5 * 3600
        self.P_max = 65*1000 #[W]
        self.S = 15.6

        # For the battery composition
        self.V_req_bat = 400    # Voltage required for engine subsytem [V]
        self.I_req_bat = 189.75 # Current required for engine subsytem [A]
        self.DOD = 90            # Depth of discharge in %


def main_bat(variables):
    E_engine_range = ppe.energy_engine(variables, endurance_s = 0.0, range_m = variables.range_m) # Energy for engine in [J] for range req.
    E_engine_endur = ppe.energy_engine(variables, endurance_s = variables.endurance_s, range_m = 0.0) # Energy for engine in [J] for endurance req.
    E_engine = max(E_engine_range, E_engine_endur) # Energy for engine in [J]
    E_avionics = ppe.energy_avionics(variables) # Energy for avionics in [J]

    # Just for information which requirement requires the most energy
    if E_engine == E_engine_endur: print("Enudance uses most energy")
    else: print("Range uses most energy")

    # Total energy required from the batteries
    E_total = E_engine + E_avionics

    # Sizing the battery
    # Different cell types
    cell1 = "Sam_21700_50E"
    cell2 = "Sam_21700_40T"
    cell3 = "Pan_18650B"

    # Data belonging to each cell type
    cell_data = batc.cell_type(name = cell1)

    # Calculating battery composition (# of cells in parallel and series)
    bat = batc.num_of_cells(cell_data, E_req = E_total/3600, V_req = variables.V_req_bat, I_req = variables.I_req_bat, DoD = variables.DOD)

    # Calculating mass and volume based on cells only
    m_bat_cells = bat[2] * cell_data.mass
    v_bat_cells = bat[2] * cell_data.volume * 1000 # from m^3 to L

    # Actual mass
    m_bat = m_bat_cells * 1.15

    # Calculating the actual energy produced due to rounding in the cell composition
    E_bat_produced = ((variables.DOD/100)*bat[1] * cell_data.C_Ah * bat[0]*cell_data.V_nom)/1000

    print("Mass bat [kg]:", m_bat)
    print("Cells volume [L]:", v_bat_cells)
    print("E required kWh:", E_total/(3600*1000))
    print("E produced kWh:", E_bat_produced)
    return None


# Testing block
if __name__ == "__main__":
    test_v = Variables_battery()
    main_bat(test_v)