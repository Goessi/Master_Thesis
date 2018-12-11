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
print("-------test on norm-------")
E = Quaternion(1.0,2.0,3.0,4.0)
print(E.norm())
print("-------test on conjugate-------")
E = Quaternion(1.0,2.0,3.0,4.0)
print(E.conj())
print("-------test on rotation pricision------")
num = []
error = []
n = 10
for m in range(100):
    n = n*n
    num.append(n)

N = 100
q1 = Quaternion.rotator(np.pi/(200),[1.0, 0.0, 0.0])
print(q1)
q2 = Quaternion.rotator(np.pi/(200),[0.0, 1.0, 0.0])
print(q2)
q3 = Quaternion.rotator(np.pi/(200),[0.0, 0.0, 1.0])
print(q3)
p = [1.0, 0.0, 0.0]
p_o = [1.0, 0.0, 0.0]
for i in range(100):
    p = q3.mul_vector(p)
#    print(p)
#    print('-----------')
#    p = p.norm_q()
#    print(p)
    p = Quaternion(0.0, p[0], p[1], p[2])
    p = p.norm_q()
    p = [p.b, p.c, p.d]

for j in range(100):
    p = q1.mul_vector(p)
    p = Quaternion(0.0, p[0], p[1], p[2])
    p = p.norm_q()
    p = [p.b, p.c, p.d]
#    p = p.norm_q()

for k in range(100):
    p = q2.mul_vector(p)
    p = Quaternion(0.0, p[0], p[1], p[2])
    p = p.norm_q()
    p = [p.b, p.c, p.d]
#    p = p.norm_q()
err = pow(p_o[0]-p[0], 2) + pow(p_o[1]-p[1], 2) + pow(p_o[2]-p[2], 2)
err = math.sqrt(err)
error.append(err)







