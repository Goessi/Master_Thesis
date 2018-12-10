# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 09:56:28 2018

test on quaternion

@author: JingQIN
"""
import math
import random
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
x = random.random()
y = random.random()
z = random.random()
p = Quaternion(0.0, x, y, z)
p = p.norm_quaternion()
p_zero = p
print(p)
print(p_zero)
print("---------------------")
p = p.rotation(360, [1.0, 0.0, 0.0])
p = p.rotation(359, [0.0, 1.0, 0.0])
p = p.rotation(360, [0.0, 0.0, 1.0])
print(p)
print(p_zero)
print(p_zero - p)
