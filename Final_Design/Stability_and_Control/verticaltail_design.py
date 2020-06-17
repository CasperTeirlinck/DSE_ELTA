from math import sqrt

lcg = 2.3           # [m]       Fuselage centre of gravity location
lf = 9.420          # [m]       Fuselage length
hfmax = 1.213       # [m]       Maximum fuselage height
Sfs = 5.258         # [m2]      Fuselage lateral surface area
hf1 = 1.146         # [m]       Fuselage nose height
hf2 = 0.306         # [m]       Fuselage tail height
bf1 = 0.960         # [m]       Fuselage nose width
bf2 = 0.243         # [m]       Fuselage tail width
Bp = 3              # [-]       Number of blades per porpeller
lp1 = 2             # [m]       Distance 1st propeller plane - aircraft centre of gravity
lp2 = 2             # [m]       Distance 1st propeller plane - aircraft centre of gravity
Dp1 = 2             # [m2]      1st propeller disk diameter
Dp2 = 2             # [m2]      1st propeller disk diameter
Sw = 15.1           # [m2]      Wing surface area
bw = sqrt(Sw*10.1)  # [m]       Wing span
CnBi = 0.024        # [-]       Wing configuration stability component
lv = 9              # [m]

# Fuselage directional stability
def k_beta(lcg,lf,hfmax):
    return 0.3*lcg/lf + 0.75*hfmax/lf - 0.105

def fuselage_stability(lcg,lf,hfmax,Sfs,Sw,bw,hf1,hf2,bf1,bf2):
    kB = k_beta(lcg,lf,hfmax)
    return -kB*Sfs*lf/(Sw*bw) * (hf1/hf2)**0.5 * (bf2/bf1)**0.5

# Propeller directional stability
def propeller_stability(lp1,lp2,Dp1,Dp2,Sw,bw,Bp):
    sum = lp1*Dp1**2/(Sw*bw) + lp2*Dp2**2/(Sw*bw)
    return -0.053*Bp * sum

# Directional stability
def directional_stability(lcg,lf,hfmax,Sfs,Sw,bw,hf1,hf2,bf1,bf2,lp1,lp2,Dp1,Dp2,Bp,CnBi):
    CnBf = fuselage_stability(lcg,lf,hfmax,Sfs,Sw,bw,hf1,hf2,bf1,bf2)
    CnBp = propeller_stability(lp1,lp2,Dp1,Dp2,Sw,bw,Bp)
    return CnBf + CnBi + CnBp

def verticaltail_sizing():
    CnB = directional_stability(lcg,lf,hfmax,Sfs,Sw,bw,hf1,hf2,bf1,bf2,lp1,lp2,Dp1,Dp2,Bp,CnBi)

    print('Directional stability coefficient =',round(CnB,2))

    SvlvSb = float(input('Vertical tail size = '))

    return SvlvSb*(Sw*bw)/lv
