# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 08:53:43 2019

@author: JingQIN
"""

import matplotlib.pyplot as plt
import numpy as np

r = np.load("Quaternion_DCM_rotation_precision.npz")

fig, ax = plt.subplots()
ax.semilogx(r['X1'], r['MEANDIFF_Q'], 'k--', label = 'Quaternion', color = 'red')
ax.semilogx(r['X1'], r['MEANDIFF_DCM'], 'k--', label = 'DCM', color = 'green')
ax.set_xlabel("1/Num_Full_Rotations, [1/degrees]")
ax.set_ylabel("Differences")
ax.set_title("Mean Differences of DCM and Quaternion, semilog_plot")
legend = ax.legend(loc = 2)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_Precision_semilog_x.png', dpi = 100)
del r