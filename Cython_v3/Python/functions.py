# -*- coding: utf-8 -*-
"""
Created on Sat Dec 29 20:15:07 2018

@author: JingQIN
"""
import math
import random
import numpy as np
import time
from Quaternion import Quaternion
import timeit

def computeDCM(theta, vectors):
    '''
    compute standard DCM, a list of lists, 3*3 matrix
    
    arguments:
        theta: rotation angle, radians
        vectors: coordinates of f, [f1, f2 ,f3]        
    '''
    DCM = [[0] * 3 for i in range(3)]
    f1 = vectors[0]
    f2 = vectors[1]
    f3 = vectors[2]
    cosTheta = math.cos(theta)
    sinTheta = math.sin(theta)
    oneCosT = 1 - cosTheta
    DCM[0][0] = cosTheta + pow(f1, 2) * oneCosT
    DCM[0][1] = f1 * f2 * oneCosT + f3 * sinTheta
    DCM[0][2] = f1 * f3 * oneCosT - f2 * sinTheta
    DCM[1][0] = f1 * f2 * oneCosT - f3 * sinTheta
    DCM[1][1] = cosTheta + pow(f2, 2) * oneCosT
    DCM[1][2] = f2 * f3 * oneCosT + f1 * sinTheta
    DCM[2][0] = f1 * f3 * oneCosT + f2 * sinTheta
    DCM[2][1] = f2 * f3 * oneCosT - f1 * sinTheta
    DCM[2][2] = cosTheta + pow(f3, 2) * oneCosT

    return np.array(DCM)

def Quaternion_rotation_precision(N, R, x_theta, y_theta, z_theta):
    '''
    calculate precision for Quaternion rotation
    
    arguments:
        N: how many steps for each round, steps is N in loop
        R: how many rounds you want to compute
        x_theta: total rotation angle around x-axis, radian
        y_theta: total rotation angle around y-axis, radian
        z_theta: total rotation angle around z-axis, radian
    
    '''
    time1 = time.clock()
    X1 = []
    X2 = []
    MEANDIFF = []
    MINDIFF = []
    MAXDIFF = []
    TIME = []
    CONST = (x_theta + y_theta + z_theta) * 180 / np.pi
    TOTALTIME = []
    
    for n in range (1,N):
        meandiff = 0.0
        maxdiff = 0.0
        mindiff = 1e9
        time_step = 0
        for r in range(0,R):
            x = random.random()
            y = random.random()
            z = random.random()
            p = Quaternion(0.0, x, y, z)
            print(p)
            p_zero = Quaternion(0.0, x, y, z)
            steps = n
            aTime = timeit.default_timer()
            for i in range(0, steps):
                p = p.rotator(x_theta / steps, [1.0, 0.0, 0.0])
            for j in range(0, steps):
                p = p.rotator(y_theta / steps, [0.0, 1.0, 0.0])
            for k in range(0, steps):
                p = p.rotator(z_theta / steps, [0.0, 0.0, 1.0])  
            bTime = timeit.default_timer()
            time_step += (bTime - aTime) / (3 * n)
            diff = (p_zero - p).norm()
            meandiff += diff
            mindiff = min(mindiff, diff)
            maxdiff = max(maxdiff, diff)        
        TIME.append(time_step / R)
        CONST1 = n * R
        X1.append(1.0 / (CONST1 * CONST))
        X2.append(CONST1)
        MEANDIFF.append(meandiff / R)
        MINDIFF.append(mindiff)
        MAXDIFF.append(maxdiff)
    TOTALTIME.append(time.clock() - time1)
    return X1, X2, MEANDIFF, MINDIFF, MAXDIFF, TIME, TOTALTIME


def DCM_rotation_precision(N, R, x_theta, y_theta, z_theta):
    '''
    calculate precision for DCM rotation
    
    arguments:
        N: how many steps for each round, steps is N in loop
        R: how many rounds you want to compute
        x_theta: total rotation angle around x-axis, radian
        y_theta: total rotation angle around y-axis, radian
        z_theta: total rotation angle around z-axis, radian
    
    '''
    time1 = time.clock()
    X1 = []
    X2 = []
    MEANDIFF = []
    MINDIFF = []
    MAXDIFF = []
    TIME = []
    CONST = (x_theta + y_theta + z_theta) * 180 / np.pi
    TOTALTIME = []
    
    for n in range (1,N):
        meandiff = 0.0
        maxdiff = 0.0
        mindiff = 1e9
        time_step = 0
        for r in range(0, R):
            x = random.random()
            y = random.random()
            z = random.random()
            v = [x, y, z]
            v_zero = [x, y, z]
            print(v)
            steps = n
            m1 = computeDCM(x_theta / steps, [1, 0, 0])
            m2 = computeDCM(y_theta / steps, [0, 1, 0])
            m3 = computeDCM(z_theta / steps, [0, 0, 1])
            aTime = timeit.default_timer()
            for i in range(0, steps):
                v = np.dot(m1, v)
            for j in range(0, steps):
                v = np.dot(m2, v)
            for k in range(0, steps):
                v = np.dot(m3, v)   
            bTime = timeit.default_timer()
            time_step += (bTime - aTime) / (3 * n)
            diff = math.sqrt(pow(v[0] - v_zero[0], 2) + pow(v[1] - v_zero[1], 2) + pow(v[2] - v_zero[2], 2))
            meandiff += diff
            mindiff = min(mindiff, diff)
            maxdiff = max(maxdiff, diff)        
        TIME.append(time_step / R)
        CONST1 = n * R
        X1.append(1.0 / (CONST1 * CONST))
        X2.append(CONST1)
        MEANDIFF.append(meandiff / R)
        MINDIFF.append(mindiff)
        MAXDIFF.append(maxdiff)
    TOTALTIME.append(time.clock() - time1)
    return X1, X2, MEANDIFF, MINDIFF, MAXDIFF, TIME, TOTALTIME



