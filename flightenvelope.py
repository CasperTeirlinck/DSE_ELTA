import matplotlib.pyplot as plt
import numpy as np

VS = 80
VA = 165
VC = 250
VD = 300

n_VA = { 'mnvr': [3.8, -1.5],   'gst': [None, None] }
n_VC = { 'mnvr': [3.8, -1.5],   'gst': [5.5, -3.5] }
n_VD = { 'mnvr': [3.8, 0],      'gst': [3.5, -1] }

""" MANOUVER ENVELOPE """
plt.plot(np.arange(0, VA, 0.1), n_VA['mnvr'][0]/VA**2 * np.arange(0, VA, 0.1)**2, linestyle='--', color='black', linewidth=1)
plt.plot([VA, VD], [n_VA['mnvr'][0], n_VD['mnvr'][0]], linestyle='--', color='black', linewidth=1)
plt.plot([VC, VD], [n_VC['mnvr'][1], n_VD['mnvr'][1]], linestyle='--', color='black', linewidth=1)

""" GUST ENVELOPE """
plt.plot([0, VD], [1, n_VD['gst'][0]], linestyle='-.', color='black', linewidth=1)
plt.plot([0, VC], [1, n_VC['gst'][0]], linestyle='-.', color='black', linewidth=1)
plt.plot([VC, VD], [n_VC['gst'][0], n_VD['gst'][0]], linestyle='-.', color='black', linewidth=1)

plt.plot([0, VD], [1, n_VD['gst'][1]], linestyle='-.', color='black', linewidth=1)
plt.plot([0, VC], [1, n_VC['gst'][1]], linestyle='-.', color='black', linewidth=1)
plt.plot([VC, VD], [n_VC['gst'][1], n_VD['gst'][1]], linestyle='-.', color='black', linewidth=1)

plt.plot([VD, VD], [min(n_VD['mnvr'][1], n_VD['gst'][1]), max(n_VD['mnvr'][0], n_VD['gst'][0])], linestyle='--', color='black', linewidth=1)


plt.plot([0, 350], [0, 0], color='black')
plt.ylim(-4, 6)
plt.xlim(0, 350)
plt.ylabel('Load Factor n [-]')
plt.xlabel('Speed V [kts]')
plt.tight_layout()
plt.show()