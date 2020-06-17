import numpy as np




class NewVariables:
    def __init__(self):
        self.init_general()
        self.init_aerodynamics()
        self.init_weight()
        self.init_propulsion()
        self.init_cs()

    def init_general(self):
        self.WTO = None
        self.Woew_classII = None
        self.Wbat = None
        self.WPL = None

    def init_aerodynamics(self):
        self.a = 100

    def init_weight(self):
        self.b = 200

    def init_propulsion(self):
        # Sizing
        self.W_batt       = None           # Battery weight in Newtons
        self.v_batt       = None           # Battery volume in liters
        self.W_motor      = 30 * 9.80665   # Motor weight in Newtons
        self.W_shaft      = 4.48 * 9.80665 # Engine shaft weight in Newtons
        self.W_prop       = 12 * 9.80665   # Propeller weight in Newtons
        self.P_max        = 65 * 1000      # Maximum power produced by the engine in W

        # Battery characteristics
        self.batt_E_prod  = None           # Energy produced by the battery in Wh
        self.batt_Ns      = None           # Number of battery cells in series
        self.batt_Np      = None           # Number of battery cells in parallel
        self.batt_N       = None           # Number of battery cells
        self.eff_batt     = 0.95           # Battery efficiency
        self.eff_prop     = 0.85           # Propeller efficiency
        self.eff_tot_prop = 0.95 * 0.85    # Total propulsive efficiency
        self.V_req_batt   = 400            # Required voltage in Volts
        self.I_req_batt   = 189.75         # Required current in Amps
        self.DoD          = 90             # Depth of discharge of the battery

        # Requirements
        self.range_m      = 250*1000       # Range requirement in meters
        self.endurance_s  = 2.5*3600       # Endurance requirement in seconds

        # Battery cell data of the Samsung 21700 50E
        self.batt_cell_diameter   = 0.0211                                                          # Cell diameter [m]
        self.batt_cell_length     = 0.0707                                                          # Cell length [m]
        self.batt_cell_volume     = pi * (self.batt_cell_diameter / 2) ** 2 * self.batt_cell_length # [m^3]
        self.batt_cell_mass       = 0.0687                                                          # Cell mass [kg]
        self.batt_cell_V_nom      = 3.6                                                             # Nominal voltage [V]
        self.batt_cell_V_max      = 4.2                                                             # Maximum voltage [V]
        self.batt_cell_V_cut      = 2.5                                                             # Cut-off voltage [V]
        self.batt_cell_I_max      = 9.8                                                             # Maximum discharge current [A]
        self.batt_cell_C_Ah       = 5.0                                                             # Capacity [Ah]
        self.batt_cell_C_Wh       = self.batt_cell_C_Ah * self.batt_cell_V_nom                      # Capacity [Wh]
        self.batt_cell_E_spec     = self.batt_cell_C_Wh / self.batt_cell_mass                       # Capacity per kg of weight [Wh/kg]
        self.batt_cell_E_vol_spec = self.batt_cell_C_Wh / self.batt_cell_volume                     # Capacity per unit volume [Wh/m^3]
        self.batt_cell_P          = self.batt_cell_I_max * self.batt_cell_V_nom                     # Maximum power [W]

    def init_cs(self):
        pass

    ###################################

    def update_WTO(self):
        self.WTO = self.Woew_classII + self.Wbat + self.WPL


if __name__ == "__main__":
    v = NewVariables()
    print(v.a)
    print(v.b)
