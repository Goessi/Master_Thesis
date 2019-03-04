# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 09:56:28 2018

test on quaternion

@author: JingQIN
"""

import matplotlib.pyplot as plt
from Quaternion import Quaternion
import functions as f
import numpy as np

print("-------Quaternion-----------------------")
(X1, X2, MEANDIFF, MINDIFF, MAXDIFF, TIME) = f.Quaternion_rotation_precision(300, 1000, np.pi, np.pi, np.pi)
np.savez("Quaternion_rotation_precision.npz", X1 = X1, X2 = X2, MEANDIFF = MEANDIFF, MINDIFF = MINDIFF, MAXDIFF = MAXDIFF, TIME = TIME)

fig, ax = plt.subplots()
ax.plot(X2, TIME, 'o--', label = 'time', color = 'red')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Time in seconds")
ax.set_title("Time in seconds(Quaternion)_plot")
legend = ax.legend(loc = 2, bbox_to_anchor = (1.05, 1.0),borderaxespad = 0.)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('Time_Quaternion1.png', dpi = 100)

print("-------DCM-----------------------")
(X1, X2, MEANDIFF, MINDIFF, MAXDIFF, TIME) = f.DCM_rotation_precision(300, 1000, np.pi, np.pi, np.pi)
np.savez("DCM_rotation_precision.npz", X1 = X1, X2 = X2, MEANDIFF = MEANDIFF, MINDIFF = MINDIFF, MAXDIFF = MAXDIFF, TIME = TIME)

fig, ax = plt.subplots()
ax.plot(X2, TIME, 'o--', label = 'time', color = 'green')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Time in seconds")
ax.set_title("Time in seconds(DCM)_plot")
legend = ax.legend(loc = 2, bbox_to_anchor = (1.05, 1.0),borderaxespad = 0.)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('Time_DCM1.png', dpi = 100)



