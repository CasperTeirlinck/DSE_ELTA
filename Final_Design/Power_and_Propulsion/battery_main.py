"""
Main battery sizing file
outputting battery mass, volume and composition
Author: Matthijs
"""
import Power_and_Propulsion.pp_energy_calculation as ppe
import Power_and_Propulsion.batteries as batc
# import pp_energy_calculation as ppe
# import batteries as batc
import numpy as np

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
        self.V_req_batt = 400    # Voltage required for engine subsytem [V]
        self.I_req_batt = 189.75 # Current required for engine subsytem [A]
        self.DoD = 90            # Depth of discharge in %

        self.batt_cell_diameter   = 0.0211                                                          # Cell diameter [m]
        self.batt_cell_length     = 0.0707                                                          # Cell length [m]
        self.batt_cell_volume     = np.pi * (self.batt_cell_diameter / 2) ** 2 * self.batt_cell_length # [m^3]
        self.batt_cell_mass       = 0.0687                                                          # Cell mass [kg]
        self.batt_cell_V_nom      = 3.6                                                             # Nominal voltage [V]
        self.batt_cell_V_max      = 4.2                                                             # Maximum voltage [V]
        self.batt_cell_V_cut      = 2.5                                                             # Cut-off voltage [V]
        self.batt_cell_I_max      = 9.8                                                             # Maximum discharge current [A]
        self.batt_cell_C_Ah       = 5.0                                                             # Capacity [Ah]
        self.batt_cell_C_Wh       = self.batt_cell_C_Ah * self.batt_cell_V_nom                      # Capacity [Wh]
        self.batt_cell_E_spec     = self.batt_cell_C_Wh / self.batt_cell_mass                       # Capacity per kg of weight [Wh/kg]
        self.batt_cell_E_vol_spec = self.batt_cell_C_Wh / self.batt_cell_volume                     # Capacity per unit volume [Wh/m^3]
        self.batt_cell_P          = self.batt_cell_I_max * self.batt_cell_V_nom




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
    # cell_data = batc.cell_type(name = cell1) #OLD code

    # Calculating battery composition (# of cells in parallel and series)
    bat = batc.num_of_cells(variables, E_req = E_total/3600)

    # Calculating mass and volume based on cells only
    m_bat_cells = bat[2] * variables.batt_cell_mass
    v_bat_cells = bat[2] * variables.batt_cell_volume * 1000 # from m^3 to L

    # Actual mass
    m_bat = m_bat_cells * 1.15

    # Calculating the actual energy produced due to rounding in the cell composition
    E_bat_produced = ((variables.DoD/100)*bat[1] * variables.batt_cell_C_Ah * bat[0]*variables.batt_cell_V_nom)/1000

    # print("Mass bat [kg]:", m_bat)
    # print("Cells volume [L]:", v_bat_cells)
    # print("E required kWh:", E_total/(3600*1000))
    # print("E produced kWh:", E_bat_produced)

    variables.W_batt = m_bat*variables.g0
    variables.v_batt = v_bat_cells
    variables.batt_E_prod = E_bat_produced*1000
    variables.batt_Ns = bat[0]
    variables.batt_Np = bat[1]
    variables.batt_N = bat[2]
    return variables


# Testing block
if __name__ == "__main__":
    test_v = Variables_battery()
    main_bat(test_v)