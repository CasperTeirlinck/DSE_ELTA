import matplotlib.pyplot as plt
import numpy as np
from variables import *

variables = CurrentVariables()

# SPEEDS: all eas !!
# assumtion: VEAS == VIAS
V_s = variables.Vs_eas
V_c = variables.V_eas
V_d = variables.Vmax
# MANOUVER
n_max = 3.8 # ref: CS-VLA 337
n_min = -1.5 # ref: CS-VLA 337
# GUST
dCLda = 2*np.pi #[1/rad]
rho = variables.rhocruise
rho0 = variables.rho0
S = variables.S  # m2
W = variables.WTO  # [N]
M = W/9.81 # [kg]
MAC = variables.MAC #[m]
U_c = 15.24  # [m/s] gust velocity ref: CS-VLA 333
U_d = 7.62  # [m/s] gust velocity ref: CS-VLA 333
K_c = ( 0.5*rho0 * dCLda * ( 0.88*(2*M/S)/(rho*dCLda*MAC) )/( 5.3 + (2*M/S)/(rho*dCLda*MAC) ) * U_c ) / (W/S) # ref: CS-VLA
K_d = ( 0.5*rho0 * dCLda * ( 0.88*(2*M/S)/(rho*dCLda*MAC) )/( 5.3 + (2*M/S)/(rho*dCLda*MAC) ) * U_d ) / (W/S) # ref: CS-VLA
# K_c = dCLda * 0.5*rho * (S / W) * U_c # ref: tudeft
# K_d = dCLda * 0.5*rho * (S / W) * U_d # ref: tudeft

""" MANOUVER ENVELOPE """
# STALL
plt.plot(
    np.arange(0, V_s * np.sqrt(n_max), 0.1),
    (np.arange(0, V_s * np.sqrt(n_max), 0.1) / V_s) ** 2,
    linestyle="--",
    color="black",
    linewidth=1,
    label="Manouver Limits",
)
plt.plot(
    np.arange(0, V_s * np.sqrt(-n_min), 0.1),
    -((np.arange(0, V_s * np.sqrt(-n_min), 0.1) / V_s) ** 2),
    linestyle="--",
    color="black",
    linewidth=1,
)

plt.plot([V_s, V_s], [0, 1], linestyle="-", color="black", linewidth=1.5)
plt.annotate("VS", (V_s, 0), xytext=(5, 5), textcoords="offset points", ha="left")
plt.plot(
    [V_s * np.sqrt(-n_min), V_s * np.sqrt(-n_min)],
    [n_min, 0],
    linestyle="-",
    color="black",
    linewidth=1.5,
)
# TOP
plt.plot(
    [V_s * np.sqrt(n_max), V_d],
    [n_max, n_max],
    linestyle="--",
    color="black",
    linewidth=1,
)
# BOTTOM
plt.plot(
    [V_s * np.sqrt(-n_min), V_c],
    [n_min, n_min],
    linestyle="--",
    color="black",
    linewidth=1,
)
plt.plot([V_c, V_d], [n_min, 0], linestyle="--", color="black", linewidth=1)
# DIVE
plt.plot([V_d, V_d], [0, n_max], linestyle="--", color="black", linewidth=1)

""" GUST ENVELOPE """
# CRUISE
plt.plot(
    np.arange(0, V_c, 0.1),
    1 + K_c * np.arange(0, V_c, 0.1),
    linestyle="-.",
    color="black",
    linewidth=1,
    label="Gust Limits",
)
plt.plot(
    np.arange(0, V_c, 0.1),
    1 - K_c * np.arange(0, V_c, 0.1),
    linestyle="-.",
    color="black",
    linewidth=1,
)
# DIVE
plt.plot(
    np.arange(0, V_d, 0.1),
    1 + K_d * np.arange(0, V_d, 0.1),
    linestyle="-.",
    color="black",
    linewidth=1,
)
plt.plot(
    np.arange(0, V_d, 0.1),
    1 - K_d * np.arange(0, V_d, 0.1),
    linestyle="-.",
    color="black",
    linewidth=1,
)
# CONNECTION
a = (K_c*V_c - K_d*V_d)/(V_c - V_d)
b_top = 1 + (K_c - a)*V_c
b_bottom = 1 - (K_c - a)*V_c
plt.plot(
    np.arange(V_c, V_d, 0.1),
    a * np.arange(V_c, V_d, 0.1) + b_top,
    linestyle="-.",
    color="black",
    linewidth=1,
)
plt.plot(
    np.arange(V_c, V_d, 0.1),
    -a * np.arange(V_c, V_d, 0.1) + b_bottom,
    linestyle="-.",
    color="black",
    linewidth=1,
)

plt.plot(
    [V_c, V_c],
    [1 - K_c * V_c, 1 + K_c * V_c],
    linestyle="-",
    color="black",
    linewidth=1.5,
)
plt.annotate("VC", (V_c, 0), xytext=(5, 5), textcoords="offset points", ha="left")

plt.plot(
    [V_d, V_d],
    [1 - K_d * V_d, 1 + K_d * V_d],
    linestyle="-",
    color="black",
    linewidth=1.5,
)
plt.annotate("VD", (V_d, 0), xytext=(5, 5), textcoords="offset points", ha="left")

""" COMBINED ENVELOPE: from left CW """
if True:
    # stall
    plt.plot(
        np.arange(V_s, V_s * np.sqrt(n_max), 0.1),
        (np.arange(V_s, V_s * np.sqrt(n_max), 0.1) / V_s) ** 2,
        linestyle="-",
        color="cornflowerblue",
        linewidth=2,
        label='Combined Envelope'
    )
    # nmax
    # plt.plot(
    #     [V_s * np.sqrt(n_max), (n_max - 1)/K_c],
    #     [n_max, n_max],
    #     linestyle="-",
    #     color="cornflowerblue",
    #     linewidth=2,
    # )
    plt.plot(
        [V_s * np.sqrt(n_max), V_c],
        [n_max, n_max],
        linestyle="-",
        color="cornflowerblue",
        linewidth=2,
    )
    plt.plot(
        [V_c, V_c],
        [n_max, a * V_c + b_top],
        linestyle="-",
        color="cornflowerblue",
        linewidth=2,
    )
    # gust
    # plt.plot(
    #     np.arange((n_max - 1)/K_c, V_c, 0.1),
    #     1 + K_c * np.arange((n_max - 1)/K_c, V_c, 0.1),
    #     linestyle="-",
    #     color="cornflowerblue",
    #     linewidth=2,    
    # )
    plt.plot(
        np.arange(V_c, (n_max - b_top)/a, 0.1),
        a * np.arange(V_c, (n_max - b_top)/a, 0.1) + b_top,
        linestyle="-",
        color="cornflowerblue",
        linewidth=2,
    )
    # nmax
    plt.plot(
        [(n_max - b_top)/a, V_d],
        [n_max, n_max],
        linestyle="-",
        color="cornflowerblue",
        linewidth=2,
    )
    # dive
    plt.plot(
        [V_d, V_d],
        [1 - K_d * V_d, n_max],
        linestyle="-",
        color="cornflowerblue",
        linewidth=2,
    )
    # gust
    _V = ( (n_min*V_d)/(V_c-V_d) + b_bottom )/( n_min/(V_c-V_d) + a )
    # plt.plot(
    #     np.arange(_V, V_d, 0.1),
    #     -a * np.arange(_V, V_d, 0.1) + b_bottom,
    #     linestyle="-",
    #     color="cornflowerblue",
    #     linewidth=2,
    # )
    # manouver
    # plt.plot(
    #     np.arange(V_c, _V, 0.1),
    #     n_min/(V_c-V_d) * np.arange(V_c, _V, 0.1) - (n_min*V_d)/(V_c-V_d),
    #     linestyle="-",
    #     color="cornflowerblue",
    #     linewidth=2,
    # )
    plt.plot(
        [V_s * np.sqrt(-n_min), V_c],
        [n_min, n_min],
        linestyle="-",
        color="cornflowerblue",
        linewidth=2,
    )

""" FINAL PLOT """
if __name__ == "__main__":
    plt.plot([0, 300], [0, 0], color="black")
    plt.ylim(-6, 8)
    plt.xlim(0, 80)
    plt.ylabel("Load Factor n [-]")
    plt.xlabel("Speed V EAS [m/s]")
    plt.legend(loc='upper left')
    plt.tight_layout()
    plt.show()
