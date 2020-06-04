"""
Author: Matthijs
"""
from math import pi, ceil

class cell_type:
    # This class basically is a keeper of all the data per cell just use the names to switch between the cells
    def __init__(self, name = ""):
        if name == "":
            print("ERROR: You forgot to specify the cell type! --> see batteries.py\nAsk Matthijs")

        # https://www.imrbatteries.com/content/samsung_50E.pdf
        if name == "Sam_21700_50E":
            self.diameter   = 0.0211 # Cell diameter [m]
            self.length     = 0.0707 # Cell length [m]
            self.volume     = pi * (self.diameter/2)**2 * self.length
            self.mass       = 0.0687 # Cell mass [kg]
            self.V_nom      = 3.6 # Nominal voltage [V]
            self.V_max      = 4.2 # Maximum voltage [V]
            self.V_cut      = 2.5 # Cut-off voltage [V]
            self.I_max      = 9.8 # Maximum discharge current [A]
            self.C_Ah       = 5.0 # Capacity [Ah]
            self.C_Wh       = self.C_Ah * self.V_nom # Capacity [Wh]
            self.E_spec     = self.C_Wh / self.mass # Capacity per kg of weight [Wh/kg]
            self.E_vol_spec = self.C_Wh / self.volume # Capacity per unit volume [Wh/m^3]
            self.P          = self.I_max * self.V_nom # Maximum power [W]

        # https://www.imrbatteries.com/content/samsung_40T.pdf
        if name == "Sam_21700_40T":
            self.diameter   = 0.0211 # Cell diameter [m]
            self.length     = 0.0704 # Cell length [m]
            self.volume     = pi * (self.diameter/2)**2 * self.length
            self.mass       = 0.0669 # Cell mass [kg]
            self.V_nom      = 3.6 # Nominal voltage [V]
            self.V_max      = 4.2 # Maximum voltage [V]
            self.V_cut      = 2.5 # Cut-off voltage [V]
            self.I_max      = 35 # Maximum discharge current [A]
            self.C_Ah       = 4.0 # Capacity [Ah]
            self.C_Wh       = self.C_Ah * self.V_nom # Capacity [Wh]
            self.E_spec     = self.C_Wh / self.mass # Capacity per kg of weight [Wh/kg]
            self.E_vol_spec = self.C_Wh / self.volume # Capacity per unit volume [Wh/m^3]
            self.P          = self.I_max * self.V_nom # Maximum power [W]

        # https://www.imrbatteries.com/content/panasonic_ncr18650b-2.pdf
        if name == "Pan_18650B":
            self.diameter   = 0.0185 # Cell diameter [m]
            self.length     = 0.0653 # Cell length [m]
            self.volume     = pi * (self.diameter/2)**2 * self.length
            self.mass       = 0.0475 # Cell mass [kg]
            self.V_nom      = 3.6 # Nominal voltage [V]
            self.V_max      = 4.2 # Maximum voltage [V]
            self.V_cut      = 2.5 # Cut-off voltage [V]
            self.I_max      = 4.87 # Maximum discharge current [A]
            self.C_Ah       = 3.350 # Capacity [Ah]
            self.C_Wh       = self.C_Ah * self.V_nom # Capacity [Wh]
            self.E_spec     = self.C_Wh / self.mass # Capacity per kg of weight [Wh/kg]
            self.E_vol_spec = self.C_Wh / self.volume # Capacity per unit volume [Wh/m^3]
            self.P          = self.I_max * self.V_nom # Maximum power [W]


def num_of_cells(cell_data, E_req, V_req, I_req, DoD, eta_batt_load = 1.0):
    """
    This functions calculates the number of cells required in ceries as well as in parallel
    :param cell_data: data about cell --> see cell_type class
    :param E_req: Required energy from battery [Wh]
    :param V_req: Required voltage from battery [V]
    :param I_req: Required current from battery [A]
    :param DoD: Dept of discharge for the battery [%] (e.g. 80%)
    :param eta_batt_load: the total efficiency of the connection between battery and load
    :return: list; number of cells in series, parallel and total
    """
    # total number of cells required for the required energy
    Num_cells_energy = E_req/cell_data.C_Wh

    # Required capacity
    C_req = E_req/(V_req * eta_batt_load * (DoD/100))

    # Number of cells is series for voltage
    Ns = V_req/cell_data.V_nom
    # Rounding to the smallest integer larger than Ns
    Ns = ceil(Ns)

    # Number of cells in parallel for capacity
    Np = C_req/cell_data.C_Ah
    Np = ceil(Np)

    # Check if enough current is generated
    if Np*cell_data.I_max < I_req:
        print("ERROR: battery cannont provide enough current, I_req should be lowered, ask Matthijs")
        return None

    # Check if enough cells are aivalable for the energy requirement
    if Num_cells_energy > Np*Ns:
        print("ERROR: number of cells required for energy is larger that Np*Ns, ask Matthijs")

    return [Ns, Np, Ns*Np]


# Main block
if __name__ == "__main__":
    cell1 = "Sam_21700_50E"
    cell2 = "Sam_21700_40T"
    cell3 = "Pan_18650B"

    cell = cell_type(name = cell1)

    # Test run
    lol = num_of_cells(cell, E_req = 52*1000, V_req = 400, I_req = 189.75, DoD = 90)
    print("Ns :", lol[0])
    print("Np :", lol[1])
    print("N  :", lol[2])









