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
        pass

    def init_cs(self):
        pass

    ###################################

    def update_WTO(self):
        self.WTO = self.Woew_classII + self.Wbat + self.WPL


if __name__ == "__main__":
    v = NewVariables()
    print(v.a)
    print(v.b)
