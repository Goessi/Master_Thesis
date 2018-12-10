# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 22:04:11 2018

first version of Master thesis codes, initial model

@author: JingQIN
"""

import math

class Quaternion(object):
    """class of Quaternion that do the simple operations
    
    Attributes:
        a: a float parameter of real part
        b: a float parameter of fundamental quaternion unit i
        c: a float parameter of fundamental quaternion unit j
        d: a float parameter of fundamental quaternion unit k
    """    
    def __init__(self, a, b, c, d):
        '''initial Quaternion class with 4 floats'''
        assert type(a) == float and type(b) == float and type(c) == float and type(d) == float
        
        self.a = a
        self.b = b
        self.c = c
        self.d = d
    
    def __add__(self, other):
        '''compute Quaternion objects addition
        
        arguments:
            other -- another Quaternion object
        '''
        return Quaternion(self.a + other.a, self.b + other.b, self.c + other.c, self.d + other.d)
    
    def __sub__(self, other):
        '''compute Quaternion objects subtraction
        
        arguments:
            other -- another Quaternion object
        '''
        return Quaternion(self.a - other.a, self.b - other.b, self.c - other.c, self.d - other.d)
    
    def __mul__(self, other):
        '''compute Quaternion objects multiple
        
        arguments:
            other -- another Quaternion object
        '''
        a = self.a*other.a - self.b*other.b - self.c*other.c - self.d*other.d
        b = self.b*other.a + self.a*other.b - self.d*other.c - self.c*other.d
        c = self.c*other.a + self.d*other.b + self.a*other.c - self.b*other.d
        d = self.d*other.a - self.c*other.b + self.b*other.c + self.a*other.d
        return Quaternion(a, b, c, d)
    
    def dot(self, other):
        '''compute Quaternion objects dot production
        
        arguments:
            other -- another Quaternion object
        '''
        return Quaternion(self.a * other.a, self.b * other.b, self.c * other.c, self.d * other.d)
    
    def norm(self):
        '''
        compute Quaternion object norm
        '''
        return math.sqrt(pow(self.a, 2) + pow(self.b, 2) + pow(self.c, 2) + pow(self.d, 2))
    
    def norm_quaternion(self):
        '''
        compute normalized Quaternion object
        '''
        return Quaternion(self.a/self.norm(), self.b/self.norm(), self.c/self.norm(), self.d/self.norm())
    
    def conj(self):
        '''
        compute Quaternion object complex conjugate
        '''
        a = self.a
        b = -self.b
        c = -self.c
        d = -self.d
        return Quaternion(a, b, c, d)
    
    def rotation(self, alpha, axis):
        '''compute Quaternion object rotation according to [x, y, z]
        
        arguements:
            alpha -- rotation angle around the x axis, degrees
            axis -- a list indicates which axi is used, [1, 0, 0] refers to x axis
        '''
        if alpha%360 == 0:
            return self
        else:
            for i in range(len(axis)):
                axis[i] = axis[i]*math.sin(alpha*math.pi/360)
            q = Quaternion(math.cos(alpha*math.pi/360), axis[0], axis[1], axis[2])
            return q*self*q.conj()
    
    def __str__(self):
        ''' document printing'''
        parameters = {'':self.a, 'i':self.b, 'j':self.c, 'k':self.d}
        print(parameters)
        count = 0
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
    
