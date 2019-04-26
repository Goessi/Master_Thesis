# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 08:09:44 2019

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
print(np.mean(r['TOTALTIME']))
print(np.mean(k['TOTALTIME']))
print(np.mean(r['TIME'])/np.mean(k['TIME']))

r = np.load("Quaternion_rotation_precision.npz")
k = np.load("DCM_rotation_precision.npz")
fig, ax = plt.subplots()
ax.plot(r['X2'], r['TIME'], 'o--', label = 'Quaternion', color = 'red')
ax.plot(k['X2'], k['TIME'], 'o--', label = 'DCM', color = 'green')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Time in seconds")
ax.set_title("Time in seconds Quaternion and DCM")
legend = ax.legend(loc = 1)
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('Timev2.png', dpi = 100)

fig, ax = plt.subplots()
ax.plot(r['X2'], r['TIME']/k['TIME'], '--', label = 'Quaternion/DCM', color = 'blue')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Time Ratios, Quaternion/DCM")
ax.set_title("Time Ratios")
legend = ax.legend(loc = 4)
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('Time Ratiosv2.png', dpi = 100)
del r
del k