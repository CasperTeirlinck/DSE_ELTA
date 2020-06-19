from Midterm_Design.variables import *
import math as m


class Test():
    def __init__(self):
        self.fus_height = 1.5
        self.tail_height = 0.75
        self.l_tail = 3
        self.xcg_aft = -1
        self.fuselagecg_z = self.fus_height*0.6
        self.x_lemac = -1
        self.z_lemac = 0.2
        self.ct = 0.6
        self.b = 12.2
        self.MAC = 1.1
        self.dihedral = 1
        self.WTO = 750*9.81
        self.propclear = 0.23
        self.tailtiplength = 3.325
        self.lowwing = True


def size_gear(variables):
    tailstrike = m.radians(15.75)  # rad
    frw_cg_angle = m.radians(5.25) # rad

    if variables.lowwing:
        z_lemac = variables.propclear
    else:
        z_lemac = variables.propclear+variables.fus_height

    h_cg = (variables.fus_height + variables.propclear)-variables.fuselagecg_z
    h_tail = (variables.fus_height + variables.propclear)-variables.tail_height

    l0 = variables.tailtiplength + variables.xcg_aft
    dist0 = l0 + m.tan(frw_cg_angle) * (h_cg - h_tail)
    z_landinggear = dist0/(1/m.tan(tailstrike) + m.tan(frw_cg_angle))+h_tail
    x_landinggear = -(z_landinggear - (h_cg - h_tail))*m.tan(frw_cg_angle)+variables.xcg_aft + (h_cg - h_tail) * m.tan(frw_cg_angle)

    x_aft_tip_wing = variables.x_lemac - variables.MAC*0.25 - variables.ct*0.75
    z_aft_tip_wing = z_lemac + 0.15*variables.ct-variables.b*0.5*m.tan(m.radians(variables.dihedral))
    angle = m.atan((z_aft_tip_wing-z_landinggear)/(-x_aft_tip_wing+x_landinggear))
    # print(angle)
    if angle > 0 and angle < m.radians(15.75):
        print("Wing tip strike! Tell Max to change code to account for that further. He did not expect this to occur so did not waste time implementing it.")

    fnoverw = 0.12
    d_ng = (variables.xcg_aft - x_landinggear) / fnoverw
    x_nosegear = x_landinggear+d_ng
    # print(variables.zcg/((x_nosegear - variables.xcg_aft)*m.tan(m.radians(57.75))))
    y_landinggear = d_ng * m.tan(m.asin(variables.fuselagecg_z/((x_nosegear - variables.xcg_aft)*m.tan(m.radians(57.75)))))

    variables.x_maingear, variables.y_maingear, variables.z_maingear, variables.x_nosegear = \
        -x_landinggear, y_landinggear, z_landinggear, -x_nosegear

    noseload = fnoverw*variables.WTO # in N
    mainload = 0.5*(variables.WTO - noseload)

    # Tire dimentions in meters
    diameternose, diametermain = 0.01*5.1*noseload**0.349, 0.01*5.1*mainload**0.349
    widthnose, widthmain = 0.01*2.3*noseload**0.312, 0.01*2.3*mainload**0.312

    variables.maingeardiameter = diametermain
    variables.maingearwidth = widthmain
    variables.nosegeardiameter = diameternose
    variables.nosegearwidth = widthnose

    return variables


if __name__ == "__main__":
    varis = Test()
    v = CurrentVariables()
    v.fuselagecg_z = 1.3497224763553965
    v.xcg_aft = 2.103121271382179
    # size_gear(varis)
    size_gear(v)