# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 17:32:27 2019

@author: JingQIN
"""

import matplotlib.pyplot as plt
#from Quaternion import Quaternion
#import functions as f
import numpy as np
#import csv
plt.tick_params(labelsize=28)
plt.rcParams.update({'font.size': 28})
r = np.load("Quaternion_rotation_precision.npz")
k = np.load("DCM_rotation_precision.npz")
print(np.mean(r['TIME']))
print(np.mean(k['TIME']))
fig, ax = plt.subplots()
ax.plot(r['X2'], r['TIME'], 'o--', label = 'Quaternion', color = 'red')
ax.plot(k['X2'], k['TIME'], 'o--', label = 'DCM', color = 'green')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Time in seconds")
ax.set_title("Time in seconds Quaternion and DCM")
legend = ax.legend(loc = 1)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('Time.png', dpi = 100)

r = np.load("Quaternion_rotation_precision.npz")
k = np.load("DCM_rotation_precision.npz")
fig, ax = plt.subplots()
ax.plot(r['X2'], r['TIME']/k['TIME'], '--', label = 'Quaternion/DCM', color = 'blue')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Time Ratios, Quaternion/DCM")
ax.set_title("Time Ratios")
legend = ax.legend(loc = 4)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('Time Ratios.png', dpi = 100)
del r
del k

r = np.load("Quaternion_DCM_rotation_precision.npz")

fig, ax = plt.subplots()
plt.tick_params(labelsize=28)
plt.rcParams.update({'font.size': 28})
ax.semilogx(r['X1'], r['MEANDIFF_DCM'], 'k--', label = 'DCM', color = 'green')
ax.semilogx(r['X1'], r['MEANDIFF_Q'], 'k--', label = 'Quaternion', color = 'red',alpha=0.7)
ax.set_xlabel("1/Num_Full_Rotations, [1/degrees]")
ax.set_ylabel("Differences")
ax.set_title("Mean Differences of DCM and Quaternion")
legend = ax.legend(loc = 1)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_Precision_semilog_x.png', dpi = 100)
del r