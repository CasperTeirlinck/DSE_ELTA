class MatProps():
    def __init__(self, sigma_y, E, poisson, rho, cost, sigma_comp=None, name=None, alpha=0.8, n=0.6):
        # Density (rho) in g/cc, NOT kg/m3!
        self.sigma_y = sigma_y
        if sigma_comp == None:
            self.sigma_comp = sigma_y
        else:
            self.sigma_comp = sigma_comp
        self.E = E
        self.poisson = poisson
        self.rho = rho*1000
        self.alpha = alpha
        self.n = n
        self.cost = cost
        if name != None:
            self.name = name

def materials():
    materials = {}
    materials['alu2024'] = MatProps(sigma_y=450000000, E=72400000000, poisson=0.33, rho=2.87, name="AA2024", alpha=0.8,
                                    n=0.6, cost=1.6)  # http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA2024T81
    materials['alu5052'] = MatProps(sigma_y=255000000, E=72300000000, poisson=0.33, rho=2.68, name="AA5052", alpha=0.8,
                                    n=0.6, cost=1.5)  # http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA5052H38
    materials['alu6061'] = MatProps(sigma_y=455000000, E=69000000000, poisson=0.33, rho=2.70, name="AA6061", alpha=0.8,
                                    n=0.6, cost=1.5)  # http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA6061T913
    materials['alu6063'] = MatProps(sigma_y=295000000, E=69000000000, poisson=0.33, rho=2.70, name="AA6063", alpha=0.8,
                                    n=0.6, cost=1.5)  # http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA6063T835
    materials['alu7050'] = MatProps(sigma_y=490000000, E=71700000000, poisson=0.33, rho=2.83, name="AA7050", alpha=0.8,
                                    n=0.6, cost=1.7)  # http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA7050T765
    materials['alu7075'] = MatProps(sigma_y=503000000, E=71700000000, poisson=0.33, rho=2.81, name="AA7075", alpha=0.8,
                                    n=0.6, cost=1.7)  # http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA7075T6
    materials['carbonfibre'] = MatProps(sigma_y=600000000, E=70000000000, poisson=0.1, rho=1.60, sigma_comp=570000000,
                                        name="carbonfibre", cost=40)  # http://www.performance-composites.com/carbonfibre/mechanicalproperties_2.asp
    materials['glassfibre'] = MatProps(sigma_y=440000000, E=25000000000, poisson=0.2, rho=1.90, sigma_comp=425000000,
                                       name="glassfibre", cost=20)  # http://www.performance-composites.com/carbonfibre/mechanicalproperties_2.asp
    return materials