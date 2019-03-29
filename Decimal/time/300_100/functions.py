# -*- coding: utf-8 -*-
"""
Created on Sat Dec 29 20:15:07 2018

@author: JingQIN
"""
import math
import random
import numpy as np
from Quaternion import Quaternion
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

def computeDCM_angle(alpha, beta, gamma):
    '''
    compute standard DCM, a list of lists, 3*3 matrix
    
    arguments:
        theta: rotation angle
        beta: rotation angle
        gamma: rotation angle
    '''
    DCM = [[0] * 3 for i in range(3)]
    sa = math.sin(alpha)
    sb = math.sin(beta)
    sg = math.sin(gamma)
    
    ca = math.cos(alpha)
    cb = math.cos(beta)
    cg = math.cos(gamma)
    
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

def DCM_diagonal_check(aDCM, bDCM):
    '''
    calculate differences  between aDCM disgonal and ones, mean of euclidean distance 
    '''
    diff = (aDCM[0][0] + aDCM[1][1] + aDCM[2][2] - bDCM[0][0] - bDCM[1][1] - bDCM[2][2]) / 3
    return diff

def DCM_off_diagonal_check(aDCM, bDCM):
    '''
    calculate differences between aDCM off-diagonal and zeros
    '''
    diff = (aDCM[0][1] + aDCM[0][2] + aDCM[1][0] + aDCM[1][2] + aDCM[2][0] + aDCM[2][1]- (bDCM[0][1] + bDCM[0][2] + bDCM[1][0] + bDCM[1][2] + bDCM[2][0] + bDCM[2][1])) / 6
    return diff

def DCM_Quaternion_rotator_check0(N, R, x_theta, y_theta, z_theta):
    '''
    calculate how many differences from rotation operator after rotation
    
    arguments:
        N: how many steps for each round, steps is N in loop
        R: how many rounds you want to compute
        x_theta: total rotation angle around x-axis, radian
        y_theta: total rotation angle around y-axis, radian
        z_theta: total rotation angle around z-axis, radian    
    '''
    X1 = []
    X2 = []
    DIAGONAL_CHECK_Q = []
    OFF_DIAGONAL_CHECK_Q = []
    DIAGONAL_CHECK_DCM = []
    OFF_DIAGONAL_CHECK_DCM = []
    CONST = (x_theta + y_theta + z_theta) * 180 / np.pi
    
    for n in range (1, N):
        diagonal_check_q = 0
        off_diagonal_check_q = 0
        diagonal_check_dcm = 0
        off_diagonal_check_dcm = 0
        
        for r in range(0, R):
            x = random.random()
            y = random.random()
            z = random.random()
            
            DCM = computeDCM_angle(x, y ,z)
            DCM_zero = computeDCM_angle(x, y ,z)
            print(DCM)
            Q = DCMtoQuaternion(DCM)
            Q_zero = DCMtoQuaternion(DCM_zero)
            print(Q)
            steps = n
            m1 = computeDCM(x_theta / steps, [1, 0, 0])
            m2 = computeDCM(y_theta / steps, [0, 1, 0])
            m3 = computeDCM(z_theta / steps, [0, 0, 1])
            m4 = computeDCM(-y_theta / steps, [0, 1, 0])
            zero = decimal.Decimal(0.0)
            q1 = Quaternion(decimal.Decimal(math.cos(x_theta / (2 * steps))), decimal.Decimal(math.sin(x_theta / (2 * steps))), zero, zero).norm_q()
            q2 = Quaternion(decimal.Decimal(math.cos(x_theta / (2 * steps))), zero, decimal.Decimal(math.sin(x_theta / (2 * steps))), zero).norm_q()
            q3 = Quaternion(decimal.Decimal(math.cos(x_theta / (2 * steps))), zero, zero, decimal.Decimal(math.sin(x_theta / (2 * steps)))).norm_q()
            q4 = Quaternion(decimal.Decimal(math.cos(-y_theta / (2 * steps))), zero, decimal.Decimal(math.sin(-y_theta / (2 * steps))), zero).norm_q()
            for i in range(0, steps):
                DCM = DCM @ m1
                Q = (q1 * Q).norm_q()
            for j in range(0, steps):
                DCM = DCM @ m2
                Q = (q2 * Q).norm_q()
            for k in range(0, steps):
                DCM = DCM @ m3
                Q = (q3 * Q).norm_q()
            for k in range(0, steps):
                DCM = DCM @ m4
                Q = (q4 * Q).norm_q()
            Q = Q.toDCM()
            Q_zero = Q_zero.toDCM()
            diagonal_check_q += DCM_diagonal_check(Q, Q_zero)
            off_diagonal_check_q += DCM_off_diagonal_check(Q, Q_zero)
            diagonal_check_dcm += DCM_diagonal_check(DCM, DCM_zero)
            off_diagonal_check_dcm += DCM_off_diagonal_check(DCM, DCM_zero)
        CONST1 = n * R
        X1.append(1.0 / (CONST1 * CONST))
        X2.append(CONST1)
        DIAGONAL_CHECK_Q.append(diagonal_check_q / R)
        OFF_DIAGONAL_CHECK_Q.append(off_diagonal_check_q / R)
        DIAGONAL_CHECK_DCM.append(diagonal_check_dcm / R)
        OFF_DIAGONAL_CHECK_DCM.append(off_diagonal_check_dcm / R)
    return X1, X2, DIAGONAL_CHECK_Q, OFF_DIAGONAL_CHECK_Q, DIAGONAL_CHECK_DCM, OFF_DIAGONAL_CHECK_DCM

def DCM_orthonormality_check(aDCM, direction):
    '''
    calculate differences between columns dot products or rows dot products
    
    arguments:
        aDCM: a 3*3 rotation matrix
        direction: define if columns or rows in calculation, 0 is columns product, 1 is rows product'
    '''
    if direction == 0:
        dot_product01 = aDCM[0][0] * aDCM[0][1] + aDCM[1][0] * aDCM[1][1] + aDCM[2][0] * aDCM[2][1]
        dot_product12 = aDCM[0][1] * aDCM[0][2] + aDCM[1][1] * aDCM[1][2] + aDCM[2][1] * aDCM[2][2]
        dot_product02 = aDCM[0][0] * aDCM[0][2] + aDCM[1][0] * aDCM[1][2] + aDCM[2][0] * aDCM[2][2]
    if direction == 1:
        dot_product01 = aDCM[0][0] * aDCM[1][0] + aDCM[0][1] * aDCM[1][1] + aDCM[0][2] * aDCM[1][2]
        dot_product12 = aDCM[1][0] * aDCM[2][0] + aDCM[1][1] * aDCM[2][1] + aDCM[1][2] * aDCM[2][2]
        dot_product02 = aDCM[0][0] * aDCM[2][0] + aDCM[0][1] * aDCM[2][1] + aDCM[0][2] * aDCM[2][2]
        
    L = math.sqrt(pow(dot_product01, 2) + pow(dot_product12, 2) + pow(dot_product02, 2))
    return L

def DCM_Quaternion_rotator_check1(N, R, x_theta, y_theta, z_theta):
    '''
    calculate how many differences from rotation operator after rotation
    
    arguments:
        N: how many steps for each round, steps is N in loop
        R: how many rounds you want to compute
        x_theta: total rotation angle around x-axis, radian
        y_theta: total rotation angle around y-axis, radian
        z_theta: total rotation angle around z-axis, radian    
    '''
    X1 = []
    X2 = []
    ORTHONORMALITY_COL_Q = []
    ORTHONORMALITY_ROW_Q = []
    ORTHONORMALITY_COL_DCM = []
    ORTHONORMALITY_ROW_DCM = []
    CONST = (x_theta + y_theta + z_theta) * 180 / np.pi
    
    for n in range (1,N):
        orthonormality_col_q = 0
        orthonormality_row_q = 0
        orthonormality_col_dcm = 0
        orthonormality_row_dcm = 0
        for r in range(0,R):
            x = random.random()
            y = random.random()
            z = random.random()

            DCM = computeDCM_angle(x, y ,z)
            DCM_zero = computeDCM_angle(x, y ,z)
            print(DCM)
            Q = DCMtoQuaternion(DCM)
            Q_zero = DCMtoQuaternion(DCM_zero)
            print(Q)
            steps = n
            m1 = computeDCM(x_theta / steps, [1, 0, 0])
            m2 = computeDCM(y_theta / steps, [0, 1, 0])
            m3 = computeDCM(z_theta / steps, [0, 0, 1])
            m4 = computeDCM(-y_theta / steps, [0, 1, 0])
            zero = decimal.Decimal(0.0)
            q1 = Quaternion(decimal.Decimal(math.cos(x_theta / (2 * steps))), decimal.Decimal(math.sin(x_theta / (2 * steps))), zero, zero).norm_q()
            q2 = Quaternion(decimal.Decimal(math.cos(x_theta / (2 * steps))), zero, decimal.Decimal(math.sin(x_theta / (2 * steps))), zero).norm_q()
            q3 = Quaternion(decimal.Decimal(math.cos(x_theta / (2 * steps))), zero, zero, decimal.Decimal(math.sin(x_theta / (2 * steps)))).norm_q()
            q4 = Quaternion(decimal.Decimal(math.cos(-y_theta / (2 * steps))), zero, decimal.Decimal(math.sin(-y_theta / (2 * steps))), zero).norm_q()
            for i in range(0, steps):
                DCM = DCM @ m1
                Q = (q1 * Q).norm_q()
            for j in range(0, steps):
                DCM = DCM @ m2
                Q = (q2 * Q).norm_q()
            for k in range(0, steps):
                DCM = DCM @ m3
                Q = (q3 * Q).norm_q()
            for k in range(0, steps):
                DCM = DCM @ m4
                Q = (q4 * Q).norm_q()
            Q = Q.toDCM()
            Q_zero = Q_zero.toDCM()

            orthonormality_col_q += DCM_orthonormality_check(Q, 0)
            orthonormality_row_q += DCM_orthonormality_check(Q, 1)
            orthonormality_col_dcm += DCM_orthonormality_check(DCM, 0)
            orthonormality_row_dcm += DCM_orthonormality_check(DCM, 1)
        CONST1 = n * R
        X1.append(1.0 / (CONST1 * CONST))
        X2.append(CONST1)
        ORTHONORMALITY_COL_Q.append(orthonormality_col_q / R)
        ORTHONORMALITY_ROW_Q.append(orthonormality_row_q / R)
        ORTHONORMALITY_COL_DCM.append(orthonormality_col_dcm / R)
        ORTHONORMALITY_ROW_DCM.append(orthonormality_row_dcm / R)
    return X1, X2, ORTHONORMALITY_COL_Q, ORTHONORMALITY_ROW_Q, ORTHONORMALITY_COL_DCM, ORTHONORMALITY_ROW_DCM

def DCM_check (aDCM,bDCM): 
    '''
    calculate differences between each elements in 2 DCMs, mean of Manhattan Distance
    '''
    diff = abs(aDCM[0][0] - bDCM[0][0]) + abs(aDCM[0][1] - bDCM[0][1]) + abs(aDCM[0][2] - bDCM[0][2]) + \
    abs(aDCM[1][0] - bDCM[1][0]) + abs(aDCM[1][1] - bDCM[1][1]) + abs(aDCM[1][2] - bDCM[1][2]) + \
    abs(aDCM[2][0] - bDCM[2][0]) + abs(aDCM[2][1] - bDCM[2][1]) + abs(aDCM[2][2] - bDCM[2][2])
    return diff / 9.0

def DCM_Quaternion_rotator_check2(N, R, x_theta, y_theta, z_theta):
    '''
    calculate how many differences from rotation operator after rotation
    
    arguments:
        N: how many steps for each round, steps is N in loop
        R: how many rounds you want to compute
        x_theta: total rotation angle around x-axis, radian
        y_theta: total rotation angle around y-axis, radian
        z_theta: total rotation angle around z-axis, radian    
    '''
    X1 = []
    X2 = []
    DIFF_CHECK_Q = []
    DIFF_CHECK_DCM = []
    CONST = (x_theta + y_theta + z_theta) * 180 / np.pi

    
    for n in range (1, N):
        diff_chk_q = 0
        diff_chk_dcm = 0

        for r in range(0, R):
            x = random.random()
            y = random.random()
            z = random.random()

            DCM = computeDCM_angle(x, y, z)
            DCM_zero = DCM
            print(DCM)
            Q = DCMtoQuaternion(DCM)
            print(Q)
            steps = n
            m1 = computeDCM(x_theta / steps, [1, 0, 0])
            m2 = computeDCM(y_theta / steps, [0, 1, 0])
            m3 = computeDCM(z_theta / steps, [0, 0, 1])
            m4 = computeDCM(-y_theta / steps, [0, 1, 0])
            zero = decimal.Decimal(0.0)
            q1 = Quaternion(decimal.Decimal(math.cos(x_theta / (2 * steps))), decimal.Decimal(math.sin(x_theta / (2 * steps))), zero, zero).norm_q()
            q2 = Quaternion(decimal.Decimal(math.cos(y_theta / (2 * steps))), zero, decimal.Decimal(math.sin(y_theta / (2 * steps))), zero).norm_q()
            q3 = Quaternion(decimal.Decimal(math.cos(z_theta / (2 * steps))), zero, zero, decimal.Decimal(math.sin(z_theta / (2 * steps)))).norm_q()
            q4 = Quaternion(decimal.Decimal(math.cos(-y_theta / (2 * steps))), zero, decimal.Decimal(math.sin(-y_theta / (2 * steps))), zero).norm_q()
            for i in range(0, steps):
                DCM = DCM @ m1
                Q = (q1 * Q).norm_q()
            for j in range(0, steps):
                DCM = DCM @ m2
                Q = (q2 * Q).norm_q()
            for k in range(0, steps):
                DCM = DCM @ m3
                Q = (q3 * Q).norm_q()
            for k in range(0, steps):
                DCM = DCM @ m4
                Q = (q4 * Q).norm_q()
            Q = Q.toDCM();
            diff_chk_q += DCM_check(Q, DCM_zero)
            diff_chk_dcm += DCM_check(DCM, DCM_zero)
        CONST1 = n * R
        X1.append(1.0 / (CONST1 * CONST))
        X2.append(CONST1)
        DIFF_CHECK_Q.append(diff_chk_q / R)
        DIFF_CHECK_DCM.append(diff_chk_dcm / R)
        
    return X1, X2, DIFF_CHECK_Q, DIFF_CHECK_DCM

def extracAngleQuaternion(aQuaternion, x, y, z):
    '''
    extract angles from a Quaternion and calculate angle differences 
    
    arguments:
        aQuaternion: a Quaternion
        x: original theta
        y: original beta
        z: original phi
    '''
    q0 = aQuaternion.a
    q1 = aQuaternion.b
    q2 = aQuaternion.c
    q3 = aQuaternion.d
    
    alpha = math.atan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (pow(q1, 2) + pow(q2, 2)))
    beta = math.asin(2 * (q0 * q2 - q3 * q1))
    gamma = math.atan2(2 * (q0 * q3 + q2 * q1), 1 - 2 * (pow(q2, 2) + pow(q3, 2)))
#    diff = abs((alpha - x) + (beta - y) + (gamma - z)) / 3
    return alpha, beta, gamma

def extracAngleDCM(aDCM, x, y, z):
    '''
    extract angles from a DCM and calculate angle differences 
    
    arguments:
        aQuaternion: a Quaternion
        x: original theta
        y: original beta
        z: original gamma
    '''
    alpha = math.atan2(aDCM[1][2], aDCM[2][2])
    beta = math.asin(-aDCM[0][2])
    gamma = math.atan2(aDCM[0][1], aDCM[0][0])
    #diff = abs((x - alpha) + (y - beta) + (z - phi)) / 3
    return alpha, beta, gamma

def DCM_Quaternion_rotator_check3(N, R, x_theta, y_theta, z_theta):
    '''
    calculate how many differences in angles from rotation operator after rotation, radians
    
    arguments:
        N: how many steps for each round, steps is N in loop
        R: how many rounds you want to compute
        x_theta: total rotation angle around x-axis, radian
        y_theta: total rotation angle around y-axis, radian
        z_theta: total rotation angle around z-axis, radian    
    '''
    X1 = []
    X2 = []
    DIFF_CHECK_Q = []
    DIFF_CHECK_DCM = []
    CONST = (x_theta + y_theta + z_theta) * 180 / np.pi

    
    for n in range (1, N):
        diff_chk_q = 0
        diff_chk_dcm = 0

        for r in range(0, R):
            x = random.random()
            y = random.random()
            z = random.random()

            DCM = computeDCM_angle(x, y, z)
            print(DCM)
            Q = DCMtoQuaternion(DCM)
            print(Q)
            steps = n
            m1 = computeDCM(x_theta / steps, [1, 0, 0])
            m2 = computeDCM(y_theta / steps, [0, 1, 0])
            m3 = computeDCM(z_theta / steps, [0, 0, 1])
            m4 = computeDCM(-y_theta / steps, [0, 1, 0])
            zero = decimal.Decimal(0.0)
            q1 = Quaternion(decimal.Decimal(math.cos(x_theta / (2 * steps))), decimal.Decimal(math.sin(x_theta / (2 * steps))), zero, zero).norm_q()
            q2 = Quaternion(decimal.Decimal(math.cos(y_theta / (2 * steps))), zero, decimal.Decimal(math.sin(y_theta / (2 * steps))), zero).norm_q()
            q3 = Quaternion(decimal.Decimal(math.cos(z_theta / (2 * steps))), zero, zero, decimal.Decimal(math.sin(z_theta / (2 * steps)))).norm_q()
            q4 = Quaternion(decimal.Decimal(math.cos(-y_theta / (2 * steps))), zero, decimal.Decimal(math.sin(-y_theta / (2 * steps))), zero).norm_q()
            for i in range(0, steps):
                DCM = DCM @ m1
                Q = (q1 * Q).norm_q()
            for j in range(0, steps):
                DCM = DCM @ m2
                Q = (q2 * Q).norm_q()
            for k in range(0, steps):
                DCM = DCM @ m3
                Q = (q3 * Q).norm_q()
            for q in range(0, steps):
                DCM = DCM @ m4
                Q = (q4 * Q).norm_q()
#            diff_chk_q += extracAngleQuaternion(Q, x, y, z)
            (x0, y0, z0) = extracAngleQuaternion(Q, x, y, z)
            diff_chk_q += x - x0 + y - y0 + z - z0
#            diff_chk_dcm += extracAngleDCM(DCM, x, y, z)
            (x1, y1, z1) = extracAngleDCM(DCM, x, y, z)
            diff_chk_dcm += x - x1 + y - y1 + z - z1
        CONST1 = n * R
        X1.append(1.0 / (CONST1 * CONST))
        X2.append(CONST1)
        DIFF_CHECK_Q.append(diff_chk_q / R)
        DIFF_CHECK_DCM.append(diff_chk_dcm / R)
        
    return X1, X2, DIFF_CHECK_Q, DIFF_CHECK_DCM

