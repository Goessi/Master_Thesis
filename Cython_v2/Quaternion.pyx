# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 22:04:11 2018

first version of Master thesis codes, initial model

@author: JingQIN
"""

cdef extern from "math.h":
    double sqrt(double m)
    float cosf(float theta)
    float sinf(float theta)
from numpy cimport ndarray
cimport numpy as np
cimport cython
from libcpp.string cimport string
from libcpp.vector cimport vector


cdef class Quaternion(object):
    """class of Quaternion that do the simple operations
    
    Attributes:
        a -- a float parameter of real part
        b -- a float parameter of fundamental quaternion unit i
        c -- a float parameter of fundamental quaternion unit j
        d -- a float parameter of fundamental quaternion unit k
    """    
    cdef float a
    cdef float b
    cdef float c
    cdef float d
    def __init__(self, float a, float b, float c, float d):
        '''initial Quaternion class with 4 floats'''
        #assert type(a) == decimal.Decimal and type(b) == decimal.Decimal and type(c) == decimal.Decimal and type(d) == decimal.Decimal
        
        self.a = a
        self.b = b
        self.c = c
        self.d = d        
    
    def __add__(self, Quaternion other):
        '''compute Quaternion objects addition
        
        arguments:
            other -- another Quaternion object
        '''
        return Quaternion(self.a + other.a, self.b + other.b, self.c + other.c, self.d + other.d)
    
    def __sub__(self, Quaternion other):
        '''compute Quaternion objects subtraction
        
        arguments:
            other -- another Quaternion object
        '''
        return Quaternion(self.a - other.a, self.b - other.b, self.c - other.c, self.d - other.d)
        
    def __mul__(self, Quaternion other):
        '''compute Quaternion objects multiple
        
        arguments:
            other -- another Quaternion object
        '''
        cdef float a = self.a * other.a - self.b * other.b - self.c * other.c - self.d * other.d
        cdef float b = self.a * other.b + self.b * other.a + self.c * other.d - self.d * other.c
        cdef float c = self.a * other.c - self.b * other.d + self.c * other.a + self.d * other.b
        cdef float d = self.a * other.d + self.b * other.c - self.c * other.b + self.d * other.a
        return Quaternion(a, b, c, d)
    
    cpdef scalar_mul(self, float scalar):
        '''compute Quaternion objects multiple
        
        arguments:
            scalar -- a scalar
        '''
        cdef float a = self.a * scalar
        cdef float b = self.b * scalar
        cdef float c = self.c * scalar
        cdef float d = self.d * scalar
        
        return Quaternion(a, b, c, d)
    
    cpdef dot(self, Quaternion other):
        '''compute Quaternion objects dot production
        
        arguments:
            other -- another Quaternion object
        '''
        return Quaternion(self.a * other.a, self.b * other.b, self.c * other.c, self.d * other.d)
    
    cpdef norm(self):
        '''
        compute Quaternion object norm
        '''
        return sqrt(pow(self.a, 2) + pow(self.b, 2) + pow(self.c, 2) + pow(self.d, 2))
    
    cpdef norm_q(self):
        '''
        compute normalized Quaternion
        '''
        cdef float mynorm = self.norm()
        cdef Quaternion my_norm_q = Quaternion(self.a / mynorm, self.b / mynorm, self.c / mynorm, self.d / mynorm)
        return my_norm_q
    
    cpdef conj(self):
        '''
        compute Quaternion object complex conjugate
        '''
        a = self.a
        b = -self.b
        c = -self.c
        d = -self.d
        return Quaternion(a, b, c, d)
    
    cpdef rotator(self, float theta, vector[float] vectors):
        '''
        from angle and vectors, compute a quaternion
        
        arguments:
            theta -- rotation angle, radians
            vectors -- indicates rotation aixs, list, like [1, 0, 0]
        '''
        
        cdef float sum_v = sum([v * v for v in vectors])
        cdef float norm_v = sqrt(sum_v)
        vectors = [v / norm_v for v in vectors]
            
        cdef float a = cosf(theta / 2.)
        cdef float b = vectors[0] * sinf(theta / 2.)
        cdef float c = vectors[1] * sinf(theta / 2.)
        cdef float d = vectors[2] * sinf(theta / 2.)
        
        r = Quaternion(a, b, c, d)
        
        return r * self * r.conj()
    
    cpdef toDCM(self):
        '''
        compute a Quaternion object to a DCM, a list of lists
        specifically, a list of three 1*3 list, normalized 
        '''
        cdef float q0 = self.a
        cdef float q1 = self.b
        cdef float q2 = self.c
        cdef float q3 = self.d
        
        cdef float C11 = pow(q0, 2) + pow(q1, 2) - pow(q2, 2) - pow(q3, 2)
        cdef float C12 = 2 * (q1 * q2 + q0 * q3)
        cdef float C13 = 2 * (q1 * q3 - q0 * q2)        
        cdef float C21 = 2 * (q1 * q2 - q0 * q3)
        cdef float C22 = pow(q0, 2) - pow(q1, 2) + pow(q2, 2) - pow(q3, 2)
        cdef float C23 = 2 * (q2 * q3 + q0 * q1)        
        cdef float C31 = 2 * (q1 * q3 + q0 * q2)
        cdef float C32 = 2 * (q2 * q3 - q0 * q1)
        cdef float C33 = pow(q0, 2) - pow(q1, 2) - pow(q2, 2) + pow(q3, 2)
        cdef float C3_norm = sqrt(pow(C31, 2) + pow(C32, 2) + pow(C33, 2))
        cdef float C1_norm = sqrt(pow(C11, 2) + pow(C12, 2) + pow(C13, 2))
        cdef float C2_norm = sqrt(pow(C21, 2) + pow(C22, 2) + pow(C23, 2))
        DCM = [[C11 / C1_norm, C12 / C1_norm, C13 / C1_norm],[C21 / C2_norm, C22 / C2_norm, C23 / C2_norm],[C31 / C3_norm, C32 / C3_norm, C33 / C3_norm]]
        return ndarray(DCM)
        
        
    def __str__(self):
        ''' document printing'''
        parameters = {'':self.a, 'i':self.b, 'j':self.c, 'k':self.d}
        cdef int count = 0
        w = ''
        for k,v in parameters.items():
            if v != 0:
                if count == 0:
                    w = w + str(v) + k
                    count += 1
                else:
                    if v < 0:
                        w = w + str(v) + k
                    else:
                        w = w + '+' + str(v) + k                
        return w
    
