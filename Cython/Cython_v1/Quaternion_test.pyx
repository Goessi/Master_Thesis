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
import csv
print("-------test on addition operation-------")
A = Quaternion(-1.0, 2.0, 3.0, 4.0)
B = Quaternion(5.0, 6.0, 7.0, 8.0)
C = A + B
print(A)
print(B)
print(C)

print("-------test on subtraction operation-------")
E = Quaternion(1.3, 2.0, 3.0, 4.0)
F = Quaternion(-2.0, 6.0, 7.0, 8.0)
G = E - F
print(E)
print(F)
print(G)

print("-------test on dot production-------")
E = Quaternion(1.1, 2.0, 3.0, 4.0)
F = Quaternion(-2.0, 6.0, 7.0, 8.0)
G = E.dot(F)
print(E)
print(F)
print(G)

print("-------test on multiply-------")
E = Quaternion(1.0, 2.0, 3.0, 4.0)
F = Quaternion(-2.0, 6.0, 7.0, 8.0)
G = E * F
print(E)
print(F)
print(G)

print("-------test on scalar multiply-------")
E = Quaternion(1.0, 2.0, 3.0, 4.0)
F = 100
G = E.scalar_mul(F)
print(E)
print(F)
print(G)

print("-------test on norm-------")
E = Quaternion(1.0, 2.0, 3.0, 4.0)
print(E.norm())

print("-------test on conjugate-------")
E = Quaternion(1.0, 2.0, 3.0, 4.0)
print(E.conj())

print("-------test on Quaternion to DCM-------")
E = Quaternion(0.0, 1.0, 1.0, 1.0)
D = E.toDCM()
print(D)

print("-------test on DCM to Quaternion-------")
print(f.DCMtoQuaternion(D))
print("-------Quaternion-----------------------")
(X1, X2, MEANDIFF, MINDIFF, MAXDIFF, TIME) = f.Quaternion_rotation_precision(30, 10, np.pi, np.pi, np.pi)
np.savez("Quaternion_rotation_precision.npz", X1 = X1, X2 = X2, MEANDIFF = MEANDIFF, MINDIFF = MINDIFF, MAXDIFF = MAXDIFF, TIME = TIME)
fig, ax = plt.subplots()
ax.loglog(X1, MEANDIFF, 'k--', label = 'Mean', color = 'red')
ax.loglog(X1, MINDIFF, 'k:', label = 'MIN', color = 'green')
ax.loglog(X1, MAXDIFF, 'k', label = 'MAX', color = 'blue')
ax.set_xlabel("1/Num_Full_Rotations, [1/degrees]")
ax.set_ylabel("Differences")
ax.set_title("Mean, Max, Min differences(Quaternion), loglog_plot")
legend = ax.legend(loc = 2, bbox_to_anchor = (1.05, 1.0),borderaxespad = 0.)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('Quaternion_rotation_Precision_loglog.png', dpi = 100)

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

#fig, ax = plt.subplots()
#ax.semilogy(X1, TIME, 'o--', label = 'time', color = 'red')
#ax.set_xlabel("Times, 1/n")
#ax.set_ylabel("Time in seconds")
#ax.set_title("Time in seconds(Quaternion)_semilogy")
#legend = ax.legend(loc = 2, bbox_to_anchor = (1.05, 1.0),borderaxespad = 0.)
#legend.get_frame().set_facecolor('C0')
#plt.show()
#fig.set_size_inches(18.5, 10.5)
#fig.savefig('Time_Quaternion2.png', dpi = 100)

print("-------DCM-----------------------")
(X1, X2, MEANDIFF, MINDIFF, MAXDIFF, TIME) = f.DCM_rotation_precision(30, 10, np.pi, np.pi, np.pi)
np.savez("DCM_rotation_precision.npz", X1 = X1, X2 = X2, MEANDIFF = MEANDIFF, MINDIFF = MINDIFF, MAXDIFF = MAXDIFF, TIME = TIME)
fig, ax = plt.subplots()
ax.loglog(X1, MEANDIFF, 'k--', label = 'Mean', color = 'red')
ax.loglog(X1, MINDIFF, 'k:', label = 'MIN', color = 'green')
ax.loglog(X1, MAXDIFF, 'k', label = 'MAX', color = 'blue')
ax.set_xlabel("1/Num_Full_Rotations, [1/degrees]")
ax.set_ylabel("Differences")
ax.set_title("Mean, Max, Min differences(DCM), loglog_plot")
legend = ax.legend(loc = 2, bbox_to_anchor = (1.05, 1.0), borderaxespad = 0.)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_rotation_Precision_loglog.png', dpi = 100)

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

#fig, ax = plt.subplots()
#ax.semilogy(X1, TIME, 'o--', label = 'time', color = 'green')
#ax.set_xlabel("Times, 1/n")
#ax.set_ylabel("Time in seconds")
#ax.set_title("Time in seconds(DCM)_semilogy")
#legend = ax.legend(loc = 2, bbox_to_anchor = (1.05, 1.0),borderaxespad = 0.)
#legend.get_frame().set_facecolor('C0')
#plt.show()
#fig.set_size_inches(18.5, 10.5)
#fig.savefig('Time_DCM2.png', dpi = 100)

(X1, X2, MEANDIFF_Q, MINDIFF_Q, MAXDIFF_Q, MEANDIFF_DCM, MINDIFF_DCM, MAXDIFF_DCM) = f.Quaternion_DCM_rotation_precision(30, 10, np.pi, np.pi, np.pi)
np.savez("Quaternion_DCM_rotation_precision.npz", X1 = X1, X2 = X2, MEANDIFF_Q = MEANDIFF_Q, MINDIFF_Q = MINDIFF_Q, MAXDIFF_Q = MAXDIFF_Q, MEANDIFF_DCM = MEANDIFF_DCM, MINDIFF_DCM = MINDIFF_DCM, MAXDIFF_DCM = MAXDIFF_DCM)
fig, ax = plt.subplots()
ax.loglog(X1, MEANDIFF_Q, 'k--', label = 'Quaternion', color = 'red')
ax.loglog(X1, MEANDIFF_DCM, 'k--', label = 'DCM', color = 'green')
ax.set_xlabel("1/Num_Full_Rotations, [1/degrees]")
ax.set_ylabel("Differences")
ax.set_title("Mean Differences of DCM and Quaternion, loglog_plot")
legend = ax.legend(loc = 2)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_Precision_loglog.png', dpi = 100)

fig, ax = plt.subplots()
ax.semilogy(X1, MEANDIFF_Q, 'k--', label = 'Quaternion', color = 'red')
ax.semilogy(X1, MEANDIFF_DCM, 'k--', label = 'DCM', color = 'green')
ax.set_xlabel("1/Num_Full_Rotations, [1/degrees]")
ax.set_ylabel("Differences")
ax.set_title("Mean Differences of DCM and Quaternion, semilog_plot")
legend = ax.legend(loc = 2)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_Precision_semilog.png', dpi = 100)

print("--------test on rotator differences-------")
(X1, X2, DIAGONAL_CHECK_Q, OFF_DIAGONAL_CHECK_Q, DIAGONAL_CHECK_DCM, OFF_DIAGONAL_CHECK_DCM) = f.DCM_Quaternion_rotator_check0(30, 10, np.pi / 2, np.pi / 2, np.pi / 2)
np.savez("DCM_Quaternion_rotator_check0.npz", X1 = X1, X2 = X2, DIAGONAL_CHECK_Q = DIAGONAL_CHECK_Q, OFF_DIAGONAL_CHECK_Q = OFF_DIAGONAL_CHECK_Q, DIAGONAL_CHECK_DCM = DIAGONAL_CHECK_DCM, OFF_DIAGONAL_CHECK_DCM = OFF_DIAGONAL_CHECK_DCM)

fig, ax = plt.subplots()
ax.plot(X2, DIAGONAL_CHECK_Q, 'k--', label = 'Quaternion', color = 'red')
ax.plot(X2, DIAGONAL_CHECK_DCM, 'k--', label = 'DCM', color = 'green')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Diagonal_Differences")
ax.set_title("DCM and Quaternion in diagonal check")
legend = ax.legend(loc = 1)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_diagonal_check.png', dpi = 100)

fig, ax = plt.subplots()
ax.plot(X2, OFF_DIAGONAL_CHECK_Q, 'k--', label = 'Quaternion', color = 'red')
ax.plot(X2, OFF_DIAGONAL_CHECK_DCM, 'k--', label = 'DCM', color = 'green')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Off_Diagonal_Differences")
ax.set_title("DCM and Quaternion in off diagonal check")
legend = ax.legend(loc = 2)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_off_diagonal_check.png', dpi = 100)

(X1, X2, ORTHONORMALITY_COL_Q, ORTHONORMALITY_ROW_Q, ORTHONORMALITY_COL_DCM, ORTHONORMALITY_ROW_DCM) = f.DCM_Quaternion_rotator_check1(30, 10, np.pi / 2, np.pi / 2, np.pi / 2)
np.savez("DCM_Quaternion_rotator_check1.npz", X1 = X1, X2 = X2, ORTHONORMALITY_COL_Q = ORTHONORMALITY_COL_Q, ORTHONORMALITY_ROW_Q = ORTHONORMALITY_ROW_Q, ORTHONORMALITY_COL_DCM = ORTHONORMALITY_COL_DCM, ORTHONORMALITY_ROW_DCM = ORTHONORMALITY_ROW_DCM)

fig, ax = plt.subplots()
ax.plot(X2, ORTHONORMALITY_ROW_Q, 'k--', label = 'Quaternion', color = 'red')
ax.plot(X2, ORTHONORMALITY_ROW_DCM, 'k--', label = 'DCM', color = 'green')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Orthonormality_check_row")
ax.set_title("Orthonormality_check_row")
legend = ax.legend(loc = 2)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_orthonormality_row_check.png', dpi = 100)

fig, ax = plt.subplots()
ax.plot(X2, ORTHONORMALITY_COL_Q, 'k--', label = 'Quaternion', color = 'red')
ax.plot(X2, ORTHONORMALITY_COL_DCM, 'k--', label = 'DCM', color = 'green')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Orthonormality_check_col")
ax.set_title("Orthonormality_check_col")
legend = ax.legend(loc = 2)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_orthonormality_col_check.png', dpi = 100)


(X1, X2, DIFF_CHECK_Q, DIFF_CHECK_DCM) = f.DCM_Quaternion_rotator_check2(30, 10, np.pi / 2, np.pi / 2, np.pi / 2)
np.savez("DCM_Quaternion_rotator_check2.npz", X1 = X1, X2 = X2, DIFF_CHECK_Q = DIFF_CHECK_Q, DIFF_CHECK_DCM = DIFF_CHECK_DCM)

fig, ax = plt.subplots()
ax.plot(X2, DIFF_CHECK_Q, 'k--', label = 'Quaternion', color = 'red')
ax.plot(X2, DIFF_CHECK_DCM, 'k--', label = 'DCM', color = 'green')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Differences")
ax.set_title("DCM and Quaternion rotator check, elemant differences")
legend = ax.legend(loc = 2, bbox_to_anchor = (1.05, 1.0),borderaxespad = 0.)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_ed.png', dpi = 100)

(X1, X2, DIFF_CHECK_Q, DIFF_CHECK_DCM) = f.DCM_Quaternion_rotator_check3(30, 10, np.pi / 2, np.pi / 2, np.pi / 2)
np.savez("DCM_Quaternion_rotator_check3.npz", X1 = X1, X2 = X2, DIFF_CHECK_Q = DIFF_CHECK_Q, DIFF_CHECK_DCM = DIFF_CHECK_DCM)

fig, ax = plt.subplots()
ax.plot(X2, DIFF_CHECK_Q, 'k--', label = 'Quaternion', color = 'red')
ax.plot(X2, DIFF_CHECK_DCM, 'k--', label = 'DCM', color = 'green')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Differences")
ax.set_title("DCM and Quaternion rotator check, angle differences")
legend = ax.legend(loc = 2, bbox_to_anchor = (1.05, 1.0),borderaxespad = 0.)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_ad.png', dpi = 100)

print("-----------------------------------file store test---------------------------------------------------")
r = np.load("Quaternion_rotation_precision.npz")

fig, ax = plt.subplots()
ax.loglog(r['X1'], r['MEANDIFF'], 'k--', label = 'Mean', color = 'red')
ax.loglog(r['X1'], r['MINDIFF'], 'k:', label = 'MIN', color = 'green')
ax.loglog(r['X1'], r['MAXDIFF'], 'k', label = 'MAX', color = 'blue')
ax.set_xlabel("1/Num_Full_Rotations, [1/degrees]")
ax.set_ylabel("Differences")
ax.set_title("Mean, Max, Min differences(Quaternion), loglog_plot")
legend = ax.legend(loc = 2, bbox_to_anchor = (1.05, 1.0),borderaxespad = 0.)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('Quaternion_rotation_Precision_loglog_biubiubiu.png', dpi = 100)

fig, ax = plt.subplots()
ax.plot(r['X2'], r['TIME'], 'o--', label = 'time', color = 'red')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Time in seconds")
ax.set_title("Time in seconds(Quaternion)_plot")
legend = ax.legend(loc = 2, bbox_to_anchor = (1.05, 1.0),borderaxespad = 0.)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('Time_Quaternion1_biubiubiu.png', dpi = 100)
del r

r = np.load("DCM_rotation_precision.npz")

fig, ax = plt.subplots()
ax.loglog(r['X1'], r['MEANDIFF'], 'k--', label = 'Mean', color = 'red')
ax.loglog(r['X1'], r['MINDIFF'], 'k:', label = 'MIN', color = 'green')
ax.loglog(r['X1'], r['MAXDIFF'], 'k', label = 'MAX', color = 'blue')
ax.set_xlabel("1/Num_Full_Rotations, [1/degrees]")
ax.set_ylabel("Differences")
ax.set_title("Mean, Max, Min differences(DCM), loglog_plot")
legend = ax.legend(loc = 2, bbox_to_anchor = (1.05, 1.0), borderaxespad = 0.)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_rotation_Precision_loglog_biubiubiu.png', dpi = 100)

fig, ax = plt.subplots()
ax.plot(r['X2'], r['TIME'], 'o--', label = 'time', color = 'green')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Time in seconds")
ax.set_title("Time in seconds(DCM)_plot")
legend = ax.legend(loc = 2, bbox_to_anchor = (1.05, 1.0),borderaxespad = 0.)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('Time_DCM1_biubiubiu.png', dpi = 100)
del r

r = np.load("Quaternion_DCM_rotation_precision.npz")

fig, ax = plt.subplots()
ax.loglog(r['X1'], r['MEANDIFF_Q'], 'k--', label = 'Quaternion', color = 'red')
ax.loglog(r['X1'], r['MEANDIFF_DCM'], 'k--', label = 'DCM', color = 'green')
ax.set_xlabel("1/Num_Full_Rotations, [1/degrees]")
ax.set_ylabel("Differences")
ax.set_title("Mean Differences of DCM and Quaternion, loglog_plot")
legend = ax.legend(loc = 2)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_Precision_loglog_biubiubiu.png', dpi = 100)

fig, ax = plt.subplots()
ax.semilogy(r['X1'], r['MEANDIFF_Q'], 'k--', label = 'Quaternion', color = 'red')
ax.semilogy(r['X1'], r['MEANDIFF_DCM'], 'k--', label = 'DCM', color = 'green')
ax.set_xlabel("1/Num_Full_Rotations, [1/degrees]")
ax.set_ylabel("Differences")
ax.set_title("Mean Differences of DCM and Quaternion, semilog_plot")
legend = ax.legend(loc = 2)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_Precision_semilog_biubiubiu.png', dpi = 100)
del r

r = np.load("DCM_Quaternion_rotator_check0.npz")

fig, ax = plt.subplots()
ax.plot(r['X2'], r['DIAGONAL_CHECK_Q'], 'k--', label = 'Quaternion', color = 'red')
ax.plot(r['X2'], r['DIAGONAL_CHECK_DCM'], 'k--', label = 'DCM', color = 'green')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Diagonal_Differences")
ax.set_title("DCM and Quaternion in diagonal check")
legend = ax.legend(loc = 1)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_diagonal_check_biubiubiu.png', dpi = 100)

fig, ax = plt.subplots()
ax.plot(r['X2'], r['OFF_DIAGONAL_CHECK_Q'], 'k--', label = 'Quaternion', color = 'red')
ax.plot(r['X2'], r['OFF_DIAGONAL_CHECK_DCM'], 'k--', label = 'DCM', color = 'green')
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
ax.plot(r['X2'], r['ORTHONORMALITY_ROW_Q'], 'k--', label = 'Quaternion', color = 'red')
ax.plot(r['X2'], r['ORTHONORMALITY_ROW_DCM'], 'k--', label = 'DCM', color = 'green')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Orthonormality_check_row")
ax.set_title("Orthonormality_check_row")
legend = ax.legend(loc = 2)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_orthonormality_row_check_biubiubiu.png', dpi = 100)

fig, ax = plt.subplots()
ax.plot(r['X2'], r['ORTHONORMALITY_COL_Q'], 'k--', label = 'Quaternion', color = 'red')
ax.plot(r['X2'], r['ORTHONORMALITY_COL_DCM'], 'k--', label = 'DCM', color = 'green')
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
ax.plot(r['X2'], r['DIFF_CHECK_Q'], 'k--', label = 'Quaternion', color = 'red')
ax.plot(r['X2'], r['DIFF_CHECK_DCM'], 'k--', label = 'DCM', color = 'green')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Differences")
ax.set_title("DCM and Quaternion rotator check, elemant differences")
legend = ax.legend(loc = 2, bbox_to_anchor = (1.05, 1.0),borderaxespad = 0.)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_ed_biubiubiu.png', dpi = 100)

r = np.load("DCM_Quaternion_rotator_check3.npz")
fig, ax = plt.subplots()
ax.plot(r['X2'], r['DIFF_CHECK_Q'], 'k--', label = 'Quaternion', color = 'red')
ax.plot(r['X2'], r['DIFF_CHECK_DCM'], 'k--', label = 'DCM', color = 'green')
ax.set_xlabel("Number of Full Rotations")
ax.set_ylabel("Differences")
ax.set_title("DCM and Quaternion rotator check, angle differences")
legend = ax.legend(loc = 2, bbox_to_anchor = (1.05, 1.0),borderaxespad = 0.)
legend.get_frame().set_facecolor('C0')
plt.show()
fig.set_size_inches(18.5, 10.5)
fig.savefig('DCM_Quaternion_ad_biubiubiu.png', dpi = 100)
del r