# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 15:21:00 2019

@author: JingQIN
"""
import matplotlib.pyplot as plt
#from Quaternion import Quaternion
#import functions as f
import numpy as np
#import csv
r = np.load("Quaternion_rotation_precision.npz")
k = np.load("DCM_rotation_precision.npz")
j = np.load("Quaternion_DCM_rotation_precision.npz")
fig, ax = plt.subplots()
ax.plot(r['X2'], r['TIME'], 'o--', label = 'Quaternion', color = 'red')
ax.plot(k['X2'], k['TIME'], 'o--', label = 'DCM', color = 'green')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Time in seconds")
ax.set_title("Time in seconds Quaternion and DCM")
legend = ax.legend(loc = 2, bbox_to_anchor = (1.05, 1.0),borderaxespad = 0.)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('Time.png', dpi = 100)

fig, ax = plt.subplots()
ax.plot(r['X2'], r['TIME']/k['TIME'], '--', label = 'Quaternion', color = 'blue')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Time Ratios, Quaternion/DCM")
ax.set_title("Time Ratios")
legend = ax.legend(loc = 2, bbox_to_anchor = (1.05, 1.0),borderaxespad = 0.)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('Time Ratios.png', dpi = 100)

fig, ax = plt.subplots()
ax.loglog(j['X1'], j['MEANDIFF_Q'], 'k--', label = 'Quaternion', color = 'red')
ax.loglog(j['X1'], j['MEANDIFF_DCM'], 'k--', label = 'DCM', color = 'green')
ax.set_xlabel("1/Num_Full_Rotations, [1/degrees]")
ax.set_ylabel("Differences")
ax.set_title("Mean Differences of DCM and Quaternion, loglog_plot")
legend = ax.legend(loc = 2)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_Precision_loglog.png', dpi = 100)

fig, ax = plt.subplots()
ax.semilogy(j['X1'], j['MEANDIFF_Q'], 'k--', label = 'Quaternion', color = 'red')
ax.semilogy(j['X1'], j['MEANDIFF_DCM'], 'k--', label = 'DCM', color = 'green')
ax.set_xlabel("1/Num_Full_Rotations, [1/degrees]")
ax.set_ylabel("Differences")
ax.set_title("Mean Differences of DCM and Quaternion, semilog_plot")
legend = ax.legend(loc = 2)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_Precision_semilog.png', dpi = 100)

fig, ax = plt.subplots()
ax.semilogx(j['X1'], j['MEANDIFF_Q'], 'k--', label = 'Quaternion', color = 'red')
ax.semilogx(j['X1'], j['MEANDIFF_DCM'], 'k--', label = 'DCM', color = 'green')
ax.set_xlabel("1/Num_Full_Rotations, [1/degrees]")
ax.set_ylabel("Differences")
ax.set_title("Mean Differences of DCM and Quaternion, semilog_plotx")
legend = ax.legend(loc = 2)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_Precision_semilogx.png', dpi = 100)
del r
del k
del j