# -*- coding: utf-8 -*-
"""
Created on Sat Dec 29 20:15:07 2018

@author: JingQIN
"""
import math
import random
import numpy as np
from Quaternion_Decimal import Quaternion
import timeit
import decimal

def computeDCM(theta, vectors):
    '''
    compute standard DCM, a list of lists, 3*3 matrix
    
    arguments:
        theta: rotation angle, radians
        vectors: coordinates of f, [f1, f2 ,f3]        
    '''
    DCM = [[0] * 3 for i in range(3)]
    f1 = decimal.Decimal(vectors[0])
    f2 = decimal.Decimal(vectors[1])
    f3 = decimal.Decimal(vectors[2])
    cosTheta = decimal.Decimal(math.cos(theta))
    sinTheta = decimal.Decimal(math.sin(theta))
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

def computeDCM_angle(alpha, beta, gamma):
    '''
    compute standard DCM, a list of lists, 3*3 matrix
    
    arguments:
        theta: rotation angle
        beta: rotation angle
        gamma: rotation angle
    '''
    DCM = [[0] * 3 for i in range(3)]
    sa = decimal.Decimal(math.sin(alpha))
    sb = decimal.Decimal(math.sin(beta))
    sg = decimal.Decimal(math.sin(gamma))
    
    ca = decimal.Decimal(math.cos(alpha))
    cb = decimal.Decimal(math.cos(beta))
    cg = decimal.Decimal(math.cos(gamma))
    
    DCM[0][0] = cb * cg
    DCM[0][1] = cb * sg
    DCM[0][2] = -sb
    DCM[1][0] = -sg * ca + sb * cg * sa
    DCM[1][1] = cg * ca + sb * sg * sa
    DCM[1][2] = cb * sa
    DCM[2][0] = sg * sa + sb * cg * ca
    DCM[2][1] = -cg * sa + sb * sg * ca
    DCM[2][2] = cb * ca

    return np.array(DCM)

def DCMtoQuaternion(aListOfLists):
    '''
    transfer from DCM to Quaternion, DCM is a list of lists
    
    argument:
        aListOfLists: a DCM    
    '''
    q0s = decimal.Decimal(0.5 * math.sqrt(aListOfLists[0][0] + aListOfLists[1][1] + aListOfLists[2][2] + 1))
    q0x = decimal.Decimal(0.5 * math.sqrt(aListOfLists[0][0] - aListOfLists[1][1] - aListOfLists[2][2] + 1))
    q0y = decimal.Decimal(0.5 * math.sqrt(-aListOfLists[0][0] + aListOfLists[1][1] - aListOfLists[2][2] + 1))
    q0z = decimal.Decimal(0.5 * math.sqrt(-aListOfLists[0][0] - aListOfLists[1][1] + aListOfLists[2][2] + 1))
    max_q = max([q0s, q0x, q0y, q0z])
    
    if max_q == q0s:
        qx = decimal.Decimal(aListOfLists[1][2] - aListOfLists[2][1]) / (4 * q0s) 
        qy = decimal.Decimal(aListOfLists[2][0] - aListOfLists[0][2]) / (4 * q0s)
        qz = decimal.Decimal(aListOfLists[0][1] - aListOfLists[1][0]) / (4 * q0s)
        return Quaternion(q0s, qx, qy, qz)
    elif max_q == q0x:
        qs = decimal.Decimal(aListOfLists[1][2] - aListOfLists[2][1]) / (4 * q0x)
        qy = decimal.Decimal(aListOfLists[0][1] + aListOfLists[1][0]) / (4 * q0x)
        qz = decimal.Decimal(aListOfLists[2][0] + aListOfLists[0][2]) / (4 * q0x)
        return Quaternion(qs, q0x, qy, qz)
    elif max_q == q0y:
        qs = decimal.Decimal(aListOfLists[2][0] - aListOfLists[0][2]) / (4 * q0y)
        qx = decimal.Decimal(aListOfLists[0][1] + aListOfLists[1][0]) / (4 * q0y)
        qz = decimal.Decimal(aListOfLists[1][2] + aListOfLists[2][1]) / (4 * q0y)
        return Quaternion(qs, qx, q0y, qz)
    elif max_q == q0z:
        qs = decimal.Decimal(aListOfLists[0][1] - aListOfLists[1][0]) / (4 * q0z)
        qx = decimal.Decimal(aListOfLists[2][0] + aListOfLists[0][2]) / (4 * q0z)
        qy = decimal.Decimal(aListOfLists[1][2] + aListOfLists[2][1]) / (4 * q0z)
        return Quaternion(qs, qx, qy, q0z)

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
    X1 = []
    X2 = []
    MEANDIFF = []
    MINDIFF = []
    MAXDIFF = []
    TIME = []
    CONST = (x_theta + y_theta + z_theta) * 180 / np.pi
    
    for n in range (1,N):
        meandiff = 0.0
        maxdiff = 0.0
        mindiff = 1e9
        time = 0
        for r in range(0,R):
            x = decimal.Decimal(random.random())
            y = decimal.Decimal(random.random())
            z = decimal.Decimal(random.random())
            zero = decimal.Decimal(0.0)
            p = Quaternion(zero, x, y, z)
            print(p)
            p_zero = Quaternion(zero, x, y, z)
            steps = n
            aTime = timeit.default_timer()
            for i in range(0, steps):
                p = p.rotator(x_theta / steps, [1.0, 0.0, 0.0])
            for j in range(0, steps):
                p = p.rotator(y_theta / steps, [0.0, 1.0, 0.0])
            for k in range(0, steps):
                p = p.rotator(z_theta / steps, [0.0, 0.0, 1.0])  
            bTime = timeit.default_timer()
            time += (bTime - aTime) / (3 * n)
            diff = (p_zero - p).norm()
            meandiff += diff
            mindiff = min(mindiff, diff)
            maxdiff = max(maxdiff, diff)        
        TIME.append(time / R)
        CONST1 = n * R
        X1.append(1.0 / (CONST1 * CONST))
        X2.append(CONST1)
        MEANDIFF.append(meandiff / R)
        MINDIFF.append(mindiff)
        MAXDIFF.append(maxdiff)
    return X1, X2, MEANDIFF, MINDIFF, MAXDIFF, TIME


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
    X1 = []
    X2 = []
    MEANDIFF = []
    MINDIFF = []
    MAXDIFF = []
    TIME = []
    CONST = (x_theta + y_theta + z_theta) * 180 / np.pi
    
    for n in range (1,N):
        meandiff = 0.0
        maxdiff = 0.0
        mindiff = 1e9
        time = 0
        for r in range(0, R):
            x = decimal.Decimal(random.random())
            y = decimal.Decimal(random.random())
            z = decimal.Decimal(random.random())
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
            time += (bTime - aTime) / (3 * n)
            diff = math.sqrt(pow(v[0] - v_zero[0], 2) + pow(v[1] - v_zero[1], 2) + pow(v[2] - v_zero[2], 2))
            meandiff += diff
            mindiff = min(mindiff, diff)
            maxdiff = max(maxdiff, diff)        
        TIME.append(time / R)
        CONST1 = n * R
        X1.append(1.0 / (CONST1 * CONST))
        X2.append(CONST1)
        MEANDIFF.append(meandiff / R)
        MINDIFF.append(mindiff)
        MAXDIFF.append(maxdiff)
    return X1, X2, MEANDIFF, MINDIFF, MAXDIFF, TIME

def Quaternion_DCM_rotation_precision(N, R, x_theta, y_theta, z_theta):
    '''
    calculate precision for Quaternion and DCM rotation, using same initial vector
    
    arguments:
        N: how many steps for each round, steps is N in loop
        R: how many rounds you want to compute
        x_theta: double total rotation angle around x-axis, radian
        y_theta: double total rotation angle around y-axis, radian
        z_theta: double total rotation angle around z-axis, radian
    
    '''
    X1 = []
    X2 = []
    MEANDIFF_Q = []
    MINDIFF_Q = []
    MAXDIFF_Q = []
    MEANDIFF_DCM = []
    MINDIFF_DCM = []
    MAXDIFF_DCM = []
    CONST = (x_theta + y_theta + z_theta) * 180 / np.pi
    
    for n in range (1,N):
        meandiff_Q = 0.0
        maxdiff_Q = 0.0
        mindiff_Q = 1e9
        meandiff_DCM = 0.0
        maxdiff_DCM = 0.0
        mindiff_DCM = 1e9
        for r in range(0,R):
            x = decimal.Decimal(random.random())
            y = decimal.Decimal(random.random())
            z = decimal.Decimal(random.random())
            zero = decimal.Decimal(0.0)
            p = Quaternion(zero, x, y, z)
            print(p)
            p_zero = Quaternion(zero, x, y, z)
            v = [x, y, z]
            v_zero = [x, y, z]
            steps = n
            m1 = computeDCM(x_theta / steps, [1, 0, 0])
            m2 = computeDCM(y_theta / steps, [0, 1, 0])
            m3 = computeDCM(z_theta / steps, [0, 0, 1])       
            for i in range(0, steps):
                p = p.rotator(x_theta / steps, [1.0, 0.0, 0.0])
                v = np.dot(m1, v)
            for i in range(0, steps):
                p = p.rotator(y_theta / steps, [0.0, 1.0, 0.0])
                v = np.dot(m2, v)
            for i in range(0, steps):
                p = p.rotator(z_theta / steps, [0.0, 0.0, 1.0])
                v = np.dot(m3, v)               
            x_Q = (p_zero - p).norm()
            meandiff_Q += x_Q
            mindiff_Q = min(mindiff_Q, x_Q)
            maxdiff_Q = max(maxdiff_Q, x_Q)
            x_DCM = math.sqrt(pow(v[0] - v_zero[0], 2) + pow(v[1] - v_zero[1], 2) + pow(v[2] - v_zero[2], 2))
            meandiff_DCM += x_DCM
            mindiff_DCM = min(mindiff_DCM, x_DCM)
            maxdiff_DCM = max(maxdiff_DCM, x_DCM)
        CONST1 = n * R
        X1.append(1.0 / (CONST1 * CONST))
        X2.append(CONST1)
        MEANDIFF_Q.append(meandiff_Q / R)
        MINDIFF_Q.append(mindiff_Q)
        MAXDIFF_Q.append(maxdiff_Q)
        MEANDIFF_DCM.append(meandiff_DCM / R)
        MINDIFF_DCM.append(mindiff_DCM)
        MAXDIFF_DCM.append(maxdiff_DCM)
    return X1, X2, MEANDIFF_Q, MINDIFF_Q, MAXDIFF_Q, MEANDIFF_DCM, MINDIFF_Q, MAXDIFF_Q



