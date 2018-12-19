# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 09:56:28 2018

test on quaternion

@author: JingQIN
"""
import math
import random
import matplotlib.pyplot as plt
import numpy as np
from Quaternion import Quaternion

print("-------test on addition operation-------")
A = Quaternion(-1.0,2.0,3.0,4.0)
B = Quaternion(5.0,6.0,7.0,8.0)
C = A + B
print(A)
print(B)
print(C)
print("-------test on subtraction operation-------")
E = Quaternion(1.0,2.0,3.0,4.0)
F = Quaternion(-2.0,6.0,7.0,8.0)
G = E - F
print(E)
print(F)
print(G)
print("-------test on dot production-------")
E = Quaternion(1.0,2.0,3.0,4.0)
F = Quaternion(-2.0,6.0,7.0,8.0)
G = E.dot(F)
print(E)
print(F)
print(G)
print("-------test on multiply-------")
E = Quaternion(1.0,2.0,3.0,4.0)
F = Quaternion(-2.0,6.0,7.0,8.0)
G = E*F
print(E)
print(F)
print(G)
print("-------test on scalar multiply-------")
E = Quaternion(1.0,2.0,3.0,4.0)
F = 100
G = E.scalar_mul(F)
print(E)
print(F)
print(G)
print("-------test on norm-------")
E = Quaternion(1.0,2.0,3.0,4.0)
print(E.norm())
print("-------test on conjugate-------")
E = Quaternion(1.0,2.0,3.0,4.0)
print(E.conj())
#print("-------test on rotation pricision------")
#
#X1 = []
#X2 = []
#MEANDIFF = []
#MINDIFF = []
#MAXDIFF = []
#for N in range (1,30):
#    meandiff=0.
#    maxdiff=0.
#    mindiff=1e9
#    for r in range(0,10000):
#        y = random.random()
#        p = Quaternion(0.0, 0.0, y, 0.)
#        p_zero = Quaternion(0.0, 0.0, y, 0.)
#        steps= 180*N
#        for i in range(0,steps):
#            p = p.rotator(np.pi/steps,[1.0, 0.0, 0.0])
#        for i in range(0,steps):
#            p = p.rotator(np.pi/steps,[0.0, 1.0, 0.0])
#        for i in range(0,steps):
#            p = p.rotator(np.pi/steps,[0.0, 0.0, 1.0])               
#        x=(p_zero - p).norm()
#        meandiff+=x
#        mindiff=min(mindiff,x)
#        maxdiff=max(maxdiff,x)
#    X1.append(1.0/N)
#    X2.append(N)
#    MEANDIFF.append(meandiff/100)
#    MINDIFF.append(mindiff)
#    MAXDIFF.append(maxdiff)
#
#fig, ax = plt.subplots()
#ax.semilogy(X1, MEANDIFF, 'k--', label='Mean', color = 'red')
#ax.semilogy(X1, MINDIFF, 'k:', label='MIN', color = 'green')
#ax.semilogy(X1, MAXDIFF, 'k', label='MAX', color = 'blue')
#ax.set_xlabel("Times")
#ax.set_ylabel("Differences")
#legend = ax.legend(loc=2, bbox_to_anchor=(1.05,1.0),borderaxespad = 0.)
#legend.get_frame().set_facecolor('C0')
#plt.show()
#fig.set_size_inches(18.5, 10.5)
#fig.savefig('test2.png', dpi=100)
print("-------test on toDCM-------")
E = Quaternion(0.0,1.0,1.0,1.0)
D = E.toDCM()
print(D)
print("-------test on DCM to Quaternion-------")
def DCMtoQuaternion(aListOfLists):
    q0s = 0.5*math.sqrt(aListOfLists[0][0] + aListOfLists[1][1] + aListOfLists[2][2] + 1)
    q0x = 0.5*math.sqrt(aListOfLists[0][0] - aListOfLists[1][1] - aListOfLists[2][2] + 1)
    q0y = 0.5*math.sqrt(-aListOfLists[0][0] + aListOfLists[1][1] - aListOfLists[2][2] + 1)
    q0z = 0.5*math.sqrt(-aListOfLists[0][0] - aListOfLists[1][1] + aListOfLists[2][2] + 1)
    max_q = max([q0s, q0x, q0y, q0z])
    
    if max_q == q0s:
        qx = (aListOfLists[2][1] - aListOfLists[1][2])/(4*q0s)
        qy = (aListOfLists[0][2] - aListOfLists[2][0])/(4*q0s)
        qz = (aListOfLists[1][0] - aListOfLists[0][1])/(4*q0s)
        return Quaternion(q0s, qx, qy, qz)
    elif max_q == q0x:
        qs = (aListOfLists[2][1] - aListOfLists[1][2])/(4*q0x)
        qy = (aListOfLists[1][0] + aListOfLists[0][1])/(4*q0x)
        qz = (aListOfLists[0][2] + aListOfLists[2][0])/(4*q0x)
        return Quaternion(qs, q0x, qy, qz)
    elif max_q == q0y:
        qs = (aListOfLists[0][2] - aListOfLists[2][0])/(4*q0y)
        qx = (aListOfLists[1][0] + aListOfLists[0][1])/(4*q0y)
        qz = (aListOfLists[2][1] + aListOfLists[1][2])/(4*q0y)
        return Quaternion(qs, qx, q0y, qz)
    elif max_q == q0z:
        qs = (aListOfLists[1][0] - aListOfLists[0][1])/(4*q0z)
        qx = (aListOfLists[0][2] + aListOfLists[2][0])/(4*q0z)
        qy = (aListOfLists[2][1] + aListOfLists[1][2])/(4*q0z)
        return Quaternion(qs, qx, qy, q0z)
print(DCMtoQuaternion(D))



