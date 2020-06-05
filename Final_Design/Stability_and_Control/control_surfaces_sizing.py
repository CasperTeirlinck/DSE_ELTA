'''
This script performs the sizing of the control surfaces.
Author: Bob
'''

# Inputs
b1 = 6.36       # [m]       Aileron start
b2 = 7.81       # [m]       Aileron end
da = 20         # [deg]     Aileron deflection angle

Clp = -(cla + cd0)*Cr*b/(24*Sw) * (1 + 3*taperw)

Clda = clda*CR/(Sw*bw) * ((b2**2-b1**2)+4*(taperw-1)/(3*bw)*(b2**3-b1**3))
