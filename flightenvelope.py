import matplotlib.pyplot as plt
import numpy as np

# all eas !!
V_s = 80
V_c = 200
V_d = 250

""" MANOUVER ENVELOPE """
n_max = 2.5
n_min = -1
# STALL
plt.plot(np.arange(0, V_s*np.sqrt(n_max), 0.1), (np.arange(0, V_s*np.sqrt(n_max), 0.1)/V_s)**2, linestyle='--', color='black', linewidth=1)
plt.plot(np.arange(0, V_s*np.sqrt(-n_min), 0.1), -(np.arange(0, V_s*np.sqrt(-n_min), 0.1)/V_s)**2, linestyle='--', color='black', linewidth=1)

plt.plot([V_s, V_s], [-1, 1], linestyle='-', color='black', linewidth=1.5)
plt.annotate('VS', (V_s, 0), xytext=(5,5), textcoords='offset points', ha='left')
# TOP
plt.plot([V_s*np.sqrt(n_max), V_d], [n_max, n_max], linestyle='--', color='black', linewidth=1)
# BOTTOM
plt.plot([V_s*np.sqrt(-n_min), V_c], [n_min, n_min], linestyle='--', color='black', linewidth=1)
plt.plot([V_c, V_d], [n_min, 0], linestyle='--', color='black', linewidth=1)
# DIVE
plt.plot([V_d, V_d], [0, n_max], linestyle='--', color='black', linewidth=1)

""" GUST ENVELOPE """
dCLda = 0.3
rho = 1.225
S = 14 # m2
W = 1200 # [N]
U_c = 5 # gust velocity
U_d = 3 # gust velocity
K_c = dCLda*(rho/2)*(S/W)*U_c
K_d = dCLda*(rho/2)*(S/W)*U_d
# CRUISE
plt.plot(np.arange(0, V_c, 0.1), 1+K_c*np.arange(0, V_c, 0.1), linestyle='-.', color='black', linewidth=1)
plt.plot(np.arange(0, V_c, 0.1), 1-K_c*np.arange(0, V_c, 0.1), linestyle='-.', color='black', linewidth=1)
# DIVE
plt.plot(np.arange(0, V_d, 0.1), 1+K_d*np.arange(0, V_d, 0.1), linestyle='-.', color='black', linewidth=1)
plt.plot(np.arange(0, V_d, 0.1), 1-K_d*np.arange(0, V_d, 0.1), linestyle='-.', color='black', linewidth=1)

plt.plot(np.linspace(V_c, V_d, 2), np.linspace(1+K_c*V_c, 1+K_d*V_d, 2), linestyle='-.', color='black', linewidth=1)
plt.plot(np.linspace(V_c, V_d, 2), np.linspace(1-K_c*V_c, 1-K_d*V_d, 2), linestyle='-.', color='black', linewidth=1)

plt.plot([V_c, V_c], [1-K_c*V_c, 1+K_c*V_c], linestyle='-', color='black', linewidth=1.5)
plt.annotate('VC', (V_c, 0), xytext=(5,5), textcoords='offset points', ha='left')

plt.plot([V_d, V_d], [1-K_d*V_d, 1+K_d*V_d], linestyle='-', color='black', linewidth=1.5)
plt.annotate('VD', (V_d, 0), xytext=(5,5), textcoords='offset points', ha='left')

""" FINAL PLOT """
plt.plot([0, 300], [0, 0], color='black')
plt.ylim(-2, 4)
plt.xlim(0, 300)
plt.ylabel('Load Factor n [-]')
plt.xlabel('Speed V [kts]')
plt.tight_layout()
plt.show()