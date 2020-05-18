import matplotlib.pyplot as plt
import numpy as np

# SPEEDS: all eas !!
V_s = 80
V_c = 200
V_d = 250
# MANOUVER
n_max = 2.5
n_min = -1.2
# GUST
dCLda = 0.25
rho = 1.225
S = 14  # m2
W = 1200  # [N]
U_c = 5  # gust velocity
U_d = 3  # gust velocity
K_c = dCLda * (rho / 2) * (S / W) * U_c
K_d = dCLda * (rho / 2) * (S / W) * U_d

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
plt.plot(
    [V_s * np.sqrt(n_max), (n_max - 1)/K_c],
    [n_max, n_max],
    linestyle="-",
    color="cornflowerblue",
    linewidth=2,
)
# gust
plt.plot(
    np.arange((n_max - 1)/K_c, V_c, 0.1),
    1 + K_c * np.arange((n_max - 1)/K_c, V_c, 0.1),
    linestyle="-",
    color="cornflowerblue",
    linewidth=2,    
)
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
plt.plot(
    np.arange(_V, V_d, 0.1),
    -a * np.arange(_V, V_d, 0.1) + b_bottom,
    linestyle="-",
    color="cornflowerblue",
    linewidth=2,
)
# manouver
plt.plot(
    np.arange(V_c, _V, 0.1),
    n_min/(V_c-V_d) * np.arange(V_c, _V, 0.1) - (n_min*V_d)/(V_c-V_d),
    linestyle="-",
    color="cornflowerblue",
    linewidth=2,
)
plt.plot(
    [V_s * np.sqrt(-n_min), V_c],
    [n_min, n_min],
    linestyle="-",
    color="cornflowerblue",
    linewidth=2,
)

""" FINAL PLOT """
plt.plot([0, 300], [0, 0], color="black")
plt.ylim(-2, 4)
plt.xlim(0, 300)
plt.ylabel("Load Factor n [-]")
plt.xlabel("Speed V EAS [kts]")
plt.legend()
plt.tight_layout()
plt.show()
