# Fuselage directional stability
def k_beta(lcg,lf,hfmax):
    return 0.3*lcg/lf + 0.75*hfmax/lf - 0.105

def fuselage_stability(lcg,lf,hfmax,Sfs,Sw,bw,hf1,hf2,bf1,bf2):
    kB = k_beta(lcg,lf,hfmax)
    return -kB*Sfs*lf/(Sw*bw) * (hf1/hf2)**(1/2) * (bf2/bf1)**(1/3)

# Propeller directional stability
def propeller_stability(lp1,lp2,Dp1,Dp2,Sw,bw,Bp):
    sum = lp1*Dp1**2/(Sw*bw) + lp2*Dp2**2/(Sw*bw)
    return -0.053*Bp * sum

# Directional stability
def directional_stability(lcg,lf,hfmax,Sfs,Sw,bw,hf1,hf2,bf1,bf2,lp1,lp2,Dp1,Dp2,Bp,CnBi):
    CnBf = fuselage_stability(lcg,lf,hfmax,Sfs,Sw,bw,hf1,hf2,bf1,bf2)
    CnBp = propeller_stability(lp1,lp2,Dp1,Dp2,Sw,bw,Bp)
    return CnBf + CnBi + CnBp

# Vertical tail sizing
def verticaltail_sizing(variables):
    lcg = variables.xcg_fuselage        # [m]       Fuselage centre of gravity location
    lf = variables.fuselagelength       # [m]       Fuselage length
    hfmax = variables.fuselageheight    # [m]       Maximum fuselage height
    Sfs = variables.Sfs                 # [m2]      Fuselage lateral surface area
    hf1 = variables.hf1                 # [m]       Fuselage nose height
    hf2 = variables.hf2                 # [m]       Fuselage tail height
    bf1 = variables.bf1                 # [m]       Fuselage nose width
    bf2 = variables.bf2                 # [m]       Fuselage tail width
    Bp = variables.Bp                   # [-]       Number of blades per porpeller
    xcg = variables.xcg_max
    xp1 = variables.xp1
    xp2 = variables.xp2
    lp1 = xcg - xp1
    lp2 = xcg - xp2                     # [m]       Distance 1st propeller plane - aircraft centre of gravity
    Dp1 = variables.Dp1                 # [m2]      1st propeller disk diameter
    Dp2 = variables.Dp2                 # [m2]      1st propeller disk diameter
    Sw = variables.S                    # [m2]      Wing surface area
    bw = variables.b                    # [m]       Wing span
    CnBi = variables.CnBi               # [-]       Wing configuration stability component
    xcg = variables.xcg_max             # [m]       Centre of gravity location
    xtail = variables.xtail             # [m]       Tail location
    lv = xtail - xcg                    # [m]       Vertical tail arm

    CnB = directional_stability(lcg,lf,hfmax,Sfs,Sw,bw,hf1,hf2,bf1,bf2,lp1,lp2,Dp1,Dp2,Bp,CnBi)
    print('---- Vertical tail sizing -----')
    print('Directional stability coefficient =',round(CnB,4))

    SvlvSb = float(input('Vertical tail size = '))

    Sv = SvlvSb*(Sw*bw)/lv

    variables.CnB = CnB
    variables.Sv = Sv
    return variables
