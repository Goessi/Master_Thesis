# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 16:56:15 2019

@author: JingQIN
"""
import matplotlib.pyplot as plt
import numpy as np

r = np.load("DCM_Quaternion_rotator_check0.npz")

fig, ax = plt.subplots()
ax.plot(r['X2'], r['DIAGONAL_CHECK_DCM'], 'k--', label = 'DCM', color = 'green')
ax.plot(r['X2'], r['DIAGONAL_CHECK_Q'], 'k--', label = 'Quaternion', color = 'red')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Diagonal_Differences")
ax.set_title("DCM and Quaternion in diagonal check")
legend = ax.legend(loc = 2)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_diagonal_check_biubiubiu.png', dpi = 100)

fig, ax = plt.subplots()
ax.plot(r['X2'], r['OFF_DIAGONAL_CHECK_DCM'], 'k--', label = 'DCM', color = 'green')
ax.plot(r['X2'], r['OFF_DIAGONAL_CHECK_Q'], 'k--', label = 'Quaternion', color = 'red')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Off_Diagonal_Differences")
ax.set_title("DCM and Quaternion in off diagonal check")
legend = ax.legend(loc = 2)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_off_diagonal_check_biubiubiu.png', dpi = 100)
del r

r = np.load("DCM_Quaternion_rotator_check1.npz")

fig, ax = plt.subplots()
ax.plot(r['X2'], r['ORTHONORMALITY_ROW_DCM'], 'k--', label = 'DCM', color = 'green')
ax.plot(r['X2'], r['ORTHONORMALITY_ROW_Q'], 'k--', label = 'Quaternion', color = 'red')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Orthonormality_check_row")
ax.set_title("Orthonormality_check_row")
legend = ax.legend(loc = 2)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_orthonormality_row_check_biubiubiu.png', dpi = 100)

fig, ax = plt.subplots()
ax.plot(r['X2'], r['ORTHONORMALITY_COL_DCM'], 'k--', label = 'DCM', color = 'green')
ax.plot(r['X2'], r['ORTHONORMALITY_COL_Q'], 'k--', label = 'Quaternion', color = 'red')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Orthonormality_check_col")
ax.set_title("Orthonormality_check_col")
legend = ax.legend(loc = 2)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_orthonormality_col_check_biubiubiu.png', dpi = 100)
del r

r = np.load("DCM_Quaternion_rotator_check2.npz")
fig, ax = plt.subplots()
ax.plot(r['X2'], r['DIFF_CHECK_DCM'], 'k--', label = 'DCM', color = 'green')
ax.plot(r['X2'], r['DIFF_CHECK_Q'], 'k--', label = 'Quaternion', color = 'red')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Differences")
ax.set_title("Elemant differences")
legend = ax.legend(loc = 2)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_ed_biubiubiu.png', dpi = 100)

r = np.load("DCM_Quaternion_rotator_check3.npz")
fig, ax = plt.subplots()
ax.plot(r['X2'], r['DIFF_CHECK_DCM'], 'k--', label = 'DCM', color = 'green')
ax.plot(r['X2'], r['DIFF_CHECK_Q'], 'k--', label = 'Quaternion', color = 'red')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Differences")
ax.set_title("Angle differences")
legend = ax.legend(loc = 2)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_ad_biubiubiu.png', dpi = 100)
del r

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