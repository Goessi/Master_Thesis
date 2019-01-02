# -*- coding: utf-8 -*-
"""
Created on Sat Dec 29 20:15:07 2018

@author: JingQIN
"""
import math
import random
import numpy as np
from Quaternion import Quaternion
def Quaternion_rotation_precision(N, R, x_double_theta, y_double_theta, z_double_theta):
    '''
    calculate precision for Quaternion rotation
    
    arguments:
        N: how many steps for each round, steps is 180*N in loop
        R: how many rounds you want to compute
        x_double_theta: double total rotation angle around x-axis, radian
        y_double_theta: double total rotation angle around y-axis, radian
        z_double_theta: double total rotation angle around z-axis, radian
    
    '''
    X1 = []
    X2 = []
    MEANDIFF = []
    MINDIFF = []
    MAXDIFF = []
    for n in range (1,N):
        meandiff = 0.0
        maxdiff = 0.0
        mindiff = 1e9
        for r in range(0,R):
            x = random.random()
            y = random.random()
            z = random.random()
            p = Quaternion(0.0, x, y, z)
            p_zero = Quaternion(0.0, x, y, z)
            steps = 180 * n
            for i in range(0, steps):
                p = p.rotator(x_double_theta / steps,[1.0, 0.0, 0.0])
            for j in range(0, steps):
                p = p.rotator(y_double_theta / steps,[0.0, 1.0, 0.0])
            for k in range(0, steps):
                p = p.rotator(z_double_theta / steps,[0.0, 0.0, 1.0])               
            diff = (p_zero - p).norm()
            meandiff += diff
            mindiff = min(mindiff, diff)
            maxdiff = max(maxdiff, diff)
        X1.append(1.0 / n)
        X2.append(n)
        MEANDIFF.append(meandiff / R)
        MINDIFF.append(mindiff)
        MAXDIFF.append(maxdiff)
    return X1, X2, MEANDIFF, MINDIFF, MAXDIFF

def DCMtoQuaternion(aListOfLists):
    '''
    transfer from DCM to Quaternion, DCM is a list of lists
    
    argument:
        aListOfLists -- DCM    
    '''
    q0s = 0.5 * math.sqrt(aListOfLists[0][0] + aListOfLists[1][1] + aListOfLists[2][2] + 1)
    q0x = 0.5 * math.sqrt(aListOfLists[0][0] - aListOfLists[1][1] - aListOfLists[2][2] + 1)
    q0y = 0.5 * math.sqrt(-aListOfLists[0][0] + aListOfLists[1][1] - aListOfLists[2][2] + 1)
    q0z = 0.5 * math.sqrt(-aListOfLists[0][0] - aListOfLists[1][1] + aListOfLists[2][2] + 1)
    max_q = max([q0s, q0x, q0y, q0z])
    
    if max_q == q0s:
        qx = (aListOfLists[2][1] - aListOfLists[1][2]) / (4 * q0s)
        qy = (aListOfLists[0][2] - aListOfLists[2][0]) / (4 * q0s)
        qz = (aListOfLists[1][0] - aListOfLists[0][1]) / (4 * q0s)
        return Quaternion(q0s, qx, qy, qz)
    elif max_q == q0x:
        qs = (aListOfLists[2][1] - aListOfLists[1][2]) / (4 * q0x)
        qy = (aListOfLists[1][0] + aListOfLists[0][1]) / (4 * q0x)
        qz = (aListOfLists[0][2] + aListOfLists[2][0]) / (4 * q0x)
        return Quaternion(qs, q0x, qy, qz)
    elif max_q == q0y:
        qs = (aListOfLists[0][2] - aListOfLists[2][0]) / (4 * q0y)
        qx = (aListOfLists[1][0] + aListOfLists[0][1]) / (4 * q0y)
        qz = (aListOfLists[2][1] + aListOfLists[1][2]) / (4 * q0y)
        return Quaternion(qs, qx, q0y, qz)
    elif max_q == q0z:
        qs = (aListOfLists[1][0] - aListOfLists[0][1]) / (4 * q0z)
        qx = (aListOfLists[0][2] + aListOfLists[2][0]) / (4 * q0z)
        qy = (aListOfLists[2][1] + aListOfLists[1][2]) / (4 * q0z)
        return Quaternion(qs, qx, qy, q0z)
    
def computeDCM(theta, vectors):
    '''
    compute standard DCM, a list of lists, 3*3 matrix
    
    arguments:
        theta -- rotation angle, radians
        vectors -- coordinates of f, [f1, f2 ,f3]        
    '''
    DCM = [[0] * 3 for i in range(3)]
    f1 = vectors[0]
    f2 = vectors[1]
    f3 = vectors[2]
    cosTheta = math.cos(theta)
    sinTheta = math.sin(theta)
    DCM[0][0] = cosTheta + pow(f1, 2) * (1 - cosTheta)
    DCM[0][1] = f1 * f2 * (1 - cosTheta) + f3 * sinTheta
    DCM[0][2] = f1 * f3 * (1 - cosTheta) - f2 * sinTheta
    DCM[1][0] = f1 * f2 * (1 - cosTheta) - f3 * sinTheta
    DCM[1][1] = cosTheta + pow(f2, 2) * (1 - cosTheta)
    DCM[1][2] = f2 * f3 * (1 - cosTheta) + f1 * sinTheta
    DCM[2][0] = f1 * f3 * (1 - cosTheta) + f2 * sinTheta
    DCM[2][1] = f2 * f3 * (1 - cosTheta) - f1 * sinTheta
    DCM[2][2] = cosTheta + pow(f3, 2) * (1 - cosTheta)
    return np.array(DCM)

def DCM_rotation_precision(N, R, x_double_theta, y_double_theta, z_double_theta):
    '''
    calculate precision for DCM rotation
    
    arguments:
        N: how many steps for each round, steps is 180*N in loop
        R: how many rounds you want to compute
        x_double_theta: double total rotation angle around x-axis, radian
        y_double_theta: double total rotation angle around y-axis, radian
        z_double_theta: double total rotation angle around z-axis, radian
    
    '''
    X1 = []
    X2 = []
    MEANDIFF = []
    MINDIFF = []
    MAXDIFF = []
    for n in range (1,N):
        meandiff = 0.0
        maxdiff = 0.0
        mindiff = 1e9
        for r in range(0, R):
            x = random.random()
            y = random.random()
            z = random.random()
            v = [x, y, z]
            v_zero = [x, y, z]
            steps= 180 * n
            m1 = computeDCM(x_double_theta/ steps, [1, 0, 0])
            m2 = computeDCM(y_double_theta/ steps, [0, 1, 0])
            m3 = computeDCM(z_double_theta/ steps, [0, 0, 1])
            for i in range(0, steps):
                v = np.dot(m1, v)
            for j in range(0, steps):
                v = np.dot(m2, v)
            for k in range(0, steps):
                v = np.dot(m3, v)              
            diff = math.sqrt(pow(v[0] - v_zero[0], 2) + pow(v[1] - v_zero[1], 2) + pow(v[2] - v_zero[2], 2))
            meandiff += diff
            mindiff = min(mindiff, diff)
            maxdiff = max(maxdiff, diff)
        X1.append(1.0 / n)
        X2.append(n)
        MEANDIFF.append(meandiff / R)
        MINDIFF.append(mindiff)
        MAXDIFF.append(maxdiff)
    return X1, X2, MEANDIFF, MINDIFF, MAXDIFF

def Quaternion_DCM_precision(N, R, x_double_theta, y_double_theta, z_double_theta):
    '''
    calculate precision for Quaternion and DCM rotation, using same initial vector
    
    arguments:
        N: how many steps for each round, steps is 180*N in loop
        R: how many rounds you want to compute
        x_double_theta: double total rotation angle around x-axis, radian
        y_double_theta: double total rotation angle around y-axis, radian
        z_double_theta: double total rotation angle around z-axis, radian
    
    '''
    X1 = []
    X2 = []
    MEANDIFF_Q = []
    MINDIFF_Q = []
    MAXDIFF_Q = []
    MEANDIFF_DCM = []
    MINDIFF_DCM = []
    MAXDIFF_DCM = []
    
    for n in range (1,N):
        meandiff_Q = 0.0
        maxdiff_Q = 0.0
        mindiff_Q = 1e9
        meandiff_DCM = 0.0
        maxdiff_DCM = 0.0
        mindiff_DCM = 1e9
        for r in range(0,R):
            x = random.random()
            y = random.random()
            z = random.random()
            p = Quaternion(0.0, x, y, z)
            p_zero = Quaternion(0.0, x, y, z)
            v = [x, y, z]
            v_zero = [x, y, z]
            steps = 180 * n
            m1 = computeDCM(x_double_theta / steps, [1, 0, 0])
            m2 = computeDCM(y_double_theta / steps, [0, 1, 0])
            m3 = computeDCM(z_double_theta / steps, [0, 0, 1])       
            for i in range(0, steps):
                p = p.rotator(x_double_theta / steps,[1.0, 0.0, 0.0])
                v = np.dot(m1, v)
            for i in range(0, steps):
                p = p.rotator(y_double_theta / steps,[0.0, 1.0, 0.0])
                v = np.dot(m2, v)
            for i in range(0, steps):
                p = p.rotator(z_double_theta / steps,[0.0, 0.0, 1.0])
                v = np.dot(m3, v)               
            x_Q = (p_zero - p).norm()
            meandiff_Q += x_Q
            mindiff_Q = min(mindiff_Q, x_Q)
            maxdiff_Q = max(maxdiff_Q, x_Q)
            x_DCM = math.sqrt(pow(v[0] - v_zero[0], 2) + pow(v[1] - v_zero[1], 2) + pow(v[2] - v_zero[2], 2))
            meandiff_DCM += x_DCM
            mindiff_DCM = min(mindiff_DCM, x_DCM)
            maxdiff_DCM = max(maxdiff_DCM, x_DCM)
        X1.append(1.0 / n)
        X2.append(n)
        MEANDIFF_Q.append(meandiff_Q / R)
        MINDIFF_Q.append(mindiff_Q)
        MAXDIFF_Q.append(maxdiff_Q)
        MEANDIFF_DCM.append(meandiff_DCM / R)
        MINDIFF_DCM.append(mindiff_DCM)
        MAXDIFF_DCM.append(maxdiff_DCM)
    return X1, X2, MEANDIFF_Q, MINDIFF_Q, MAXDIFF_Q, MEANDIFF_DCM, MINDIFF_Q, MAXDIFF_Q

def DCM_diagonal_check(aDCM):
    '''
    calculate differences  between aDCM disgonal and ones
    '''
    diff = math.sqrt(pow(aDCM[0][0] - 1, 2) + pow(aDCM[1][1] - 1, 2) + pow(aDCM[2][2] - 1, 2))
    return diff

def DCM_off_diagonal_check(aDCM):
    '''
    calculate differences between aDCM off-diagonal and zeros
    '''
    diff = math.sqrt(pow(aDCM[0][1], 2) + pow(aDCM[0][2], 2) + pow(aDCM[1][0], 2) + pow(aDCM[1][2], 2) + pow(aDCM[2][0], 2) + pow(aDCM[2][1], 2))
    return diff

def DCM_orthonormality_check(aDCM, direction):
    '''
    calculate differences between columns dot products or rows dot products
    
    arguments:
        aDCM -- a 3*3 rotation matrix
        direction -- define if columns or rows in calculation, 0 is columns product, 1 is rows product'
    '''
    if direction == 0:
        dot_product01 = [aDCM[0][0]*aDCM[0][1], aDCM[1][0]*aDCM[1][1], aDCM[2][0]*aDCM[2][1]]
        dot_product12 = [aDCM[0][1]*aDCM[0][2], aDCM[1][1]*aDCM[1][2], aDCM[2][1]*aDCM[2][2]]
        dot_product02 = [aDCM[0][0]*aDCM[0][2], aDCM[1][0]*aDCM[1][2], aDCM[2][0]*aDCM[2][2]]
    if direction == 1:
        dot_product01 = [aDCM[0][0]*aDCM[1][0], aDCM[0][1]*aDCM[1][1], aDCM[0][2]*aDCM[1][2]]
        dot_product12 = [aDCM[1][0]*aDCM[2][0], aDCM[1][1]*aDCM[2][1], aDCM[1][2]*aDCM[2][2]]
        dot_product02 = [aDCM[0][0]*aDCM[2][0], aDCM[0][1]*aDCM[2][1], aDCM[0][2]*aDCM[2][2]]
    
    L = math.sqrt(pow(dot_product01[0], 2) + pow(dot_product01[1], 2) + pow(dot_product01[2], 2)) + math.sqrt(pow(dot_product12[0], 2) + pow(dot_product12[1], 2) + pow(dot_product12[2], 2)) + math.sqrt(pow(dot_product02[0], 2) + pow(dot_product02[1], 2) + pow(dot_product02[2], 2))
    return L

def DCM_Quaternion_rotator_check(N, R, x_double_theta, y_double_theta, z_double_theta):
    '''
    calculate how many differences from rotation operator after rotation
    
    arguments:
        N: how many steps for each round, steps is 180*N in loop
        R: how many rounds you want to compute
        x_double_theta: double total rotation angle around x-axis, radian
        y_double_theta: double total rotation angle around y-axis, radian
        z_double_theta: double total rotation angle around z-axis, radian    
    '''
    X1 = []
    X2 = []
    DIAGONAL_CHECK_Q = []
    OFF_DIAGONAL_CHECK_Q = []
    ORTHONORMALITY_COL_Q = []
    ORTHONORMALITY_ROW_Q = []
    DIAGONAL_CHECK_DCM = []
    OFF_DIAGONAL_CHECK_DCM = []
    ORTHONORMALITY_COL_DCM = []
    ORTHONORMALITY_ROW_DCM = []
    
    for n in range (1,N):
        diagonal_check_q = 0
        off_diagonal_check_q = 0
        orthonormality_col_q = 0
        orthonormality_row_q = 0
        diagonal_check_dcm = 0
        off_diagonal_check_dcm = 0
        orthonormality_col_dcm = 0
        orthonormality_row_dcm = 0
        for r in range(0,R):
            DCM = np.array([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
            DCM_zero = np.array([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
            Q = DCMtoQuaternion(DCM)
            Q_zero = DCMtoQuaternion(DCM)
            steps = 180 * n
            m1 = computeDCM(x_double_theta / steps, [1, 0, 0])
            m2 = computeDCM(y_double_theta / steps, [0, 1, 0])
            m3 = computeDCM(z_double_theta / steps, [0, 0, 1])
            q1 = Quaternion(math.cos(x_double_theta / (2*steps)), math.sin(x_double_theta / (2*steps)), 0, 0).norm_q()
            q2 = Quaternion(math.cos(x_double_theta / (2*steps)), 0, math.sin(x_double_theta / (2*steps)), 0).norm_q()
            q3 = Quaternion(math.cos(x_double_theta / (2*steps)), 0, 0, math.sin(x_double_theta / (2*steps))).norm_q()
            for i in range(0, steps):
                DCM = m1 * DCM
                Q = (q1 * Q).norm_q()
            for j in range(0, steps):
                DCM = m2 * DCM
                Q = (q2 * Q).norm_q()
            for k in range(0, steps):
                DCM = m3 * DCM
                Q = (q3 * Q).norm_q()
            Q = Q.toDCM()
            diagonal_check_q += DCM_diagonal_check(Q)
            off_diagonal_check_q += DCM_off_diagonal_check(Q)
            orthonormality_col_q += DCM_orthonormality_check(Q, 0)
            orthonormality_row_q += DCM_orthonormality_check(Q, 1)
            diagonal_check_dcm += DCM_diagonal_check(DCM)
            off_diagonal_check_dcm += DCM_off_diagonal_check(DCM)
            orthonormality_col_dcm += DCM_orthonormality_check(DCM, 0)
            orthonormality_row_dcm += DCM_orthonormality_check(DCM, 1)
        X1.append(1.0 / n)
        X2.append(n)
        DIAGONAL_CHECK_Q.append(diagonal_check_q / R)
        OFF_DIAGONAL_CHECK_Q.append(off_diagonal_check_q / R)
        ORTHONORMALITY_COL_Q.append(orthonormality_col_q / R)
        ORTHONORMALITY_ROW_Q.append(orthonormality_row_q / R)
        DIAGONAL_CHECK_DCM.append(diagonal_check_dcm / R)
        OFF_DIAGONAL_CHECK_DCM.append(off_diagonal_check_dcm / R)
        ORTHONORMALITY_COL_DCM.append(orthonormality_col_dcm / R)
        ORTHONORMALITY_ROW_DCM.append(orthonormality_row_dcm / R)
    return X1, X2, DIAGONAL_CHECK_Q, OFF_DIAGONAL_CHECK_Q, ORTHONORMALITY_COL_Q, ORTHONORMALITY_ROW_Q, DIAGONAL_CHECK_DCM, OFF_DIAGONAL_CHECK_DCM, ORTHONORMALITY_COL_DCM, ORTHONORMALITY_ROW_DCM