# -*- coding: utf-8 -*-
"""
Created on Sat Dec 29 20:15:07 2018

@author: JingQIN
"""
import random
from Quaternion import Quaternion
import timeit
import numpy as np

cdef extern from "math.h":
    double sqrt(double m)
    float cosf(float theta)
    float sinf(float theta)
from numpy cimport ndarray
cimport numpy as np
cimport cython
from libcpp.string cimport string
from libcpp.vector cimport vector


def computeDCM(float theta, vector[float] vectors):
    '''
    compute standard DCM, a list of lists, 3*3 matrix
    
    arguments:
        theta: rotation angle, radians
        vectors: coordinates of f, [f1, f2 ,f3]        
    '''
    DCM = [[0] * 3 for i in range(3)]
    cdef float f1 = vectors[0]
    cdef float f2 = vectors[1]
    cdef float f3 = vectors[2]
    cdef float cosTheta = cosf(theta)
    cdef float sinTheta = sinf(theta)
    cdef float oneCosT = 1 - cosTheta
    DCM[0][0] = cosTheta + pow(f1, 2) * oneCosT
    DCM[0][1] = f1 * f2 * oneCosT + f3 * sinTheta
    DCM[0][2] = f1 * f3 * oneCosT - f2 * sinTheta
    DCM[1][0] = f1 * f2 * oneCosT - f3 * sinTheta
    DCM[1][1] = cosTheta + pow(f2, 2) * oneCosT
    DCM[1][2] = f2 * f3 * oneCosT + f1 * sinTheta
    DCM[2][0] = f1 * f3 * oneCosT + f2 * sinTheta
    DCM[2][1] = f2 * f3 * oneCosT - f1 * sinTheta
    DCM[2][2] = cosTheta + pow(f3, 2) * oneCosT

    return ndarray(DCM)

def Quaternion_rotation_precision(int N, int R, float x_theta, float y_theta, float z_theta):
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
    cdef float CONST = (x_theta + y_theta + z_theta) * 180 / np.pi
    cdef int CONST1 = 0
    cdef float meandiff = 0.0
    cdef float maxdiff = 0.0
    cdef float mindiff = 0.0
    cdef float time = 0.0
    cdef float x = 0.0
    cdef float y = 0.0
    cdef float z = 0.0
    cdef int steps = 0
    cdef float diff = 0.0
    
    for n in range (1,N):
        meandiff = 0.0
        maxdiff = 0.0
        mindiff = 1e9
        time = 0
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


def DCM_rotation_precision(int N, int R, float x_theta, float y_theta, float z_theta):
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
    cdef float CONST = (x_theta + y_theta + z_theta) * 180 / np.pi
    cdef int CONST1 = 0
    cdef float meandiff = 0.0
    cdef float maxdiff = 0.0
    cdef float mindiff = 0.0
    cdef float time = 0.0
    cdef float x = 0.0
    cdef float y = 0.0
    cdef float z = 0.0
    v = []
    cdef int steps = 0
    cdef float diff = 0.0
    
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
            diff = sqrt(pow(v[0] - v_zero[0], 2) + pow(v[1] - v_zero[1], 2) + pow(v[2] - v_zero[2], 2))
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
