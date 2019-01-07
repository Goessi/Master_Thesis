# -*- coding: utf-8 -*-
"""
Created on Mon Jan  7 18:53:44 2019

@author: JingQIN
"""
import math
import random
import functions as f
import Quaternion as q
def extracAngleQuaternion(aQuaternion, x, y, z):
    '''
    extract angles from a DCM and calculate angle differences 
    
    arguments:
        aQuaternion: a Quaternion
    '''
    q0 = aQuaternion.a
    q1 = aQuaternion.b
    q2 = aQuaternion.c
    q3 = aQuaternion.d
    
    alpha = math.atan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (pow(q1, 2) + pow(q2, 2)))
    theta = math.asin(2 * (q0 * q2 - q3 * q1))
    phi = math.atan2(2 * (q0 * q3 + q2 * q1), 1 - 2 * (pow(q2, 2) + pow(q3, 2)))
    diff = abs((alpha - x) + (theta - y) + (phi - z)) / 3
    return diff

x = random.random()
y = random.random()
z = random.random()

m1 = f.computeDCM(x, [1, 0, 0])
m2 = f.computeDCM(y, [0, 1, 0])
m3 = f.computeDCM(z, [0, 0, 1])

DCM = m3 @ m2 @ m1
Q = f.DCMtoQuaternion(DCM)

(x1, y1, z1) = extracAngleQuaternion(Q)
print(abs((x1 - x) + (y1 - y) + (z1 - z)) / 3)