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
    
#    def __mul__(self, other):
#        '''compute Quaternion objects multiple
#        
#        arguments:
#            other -- another Quaternion object
#        '''
#        a = self.a*other.a - self.b*other.b - self.c*other.c - self.d*other.d
#        b = self.b*other.a + self.a*other.b - self.d*other.c - self.c*other.d
#        c = self.c*other.a + self.d*other.b + self.a*other.c - self.b*other.d
#        d = self.d*other.a - self.c*other.b + self.b*other.c + self.a*other.d
#        return Quaternion(a, b, c, d)
        
    def __mul__(self, other):
        '''compute Quaternion objects multiple
        
        arguments:
            other -- another Quaternion object
        '''
        a = self.a*other.a - self.b*other.b - self.c*other.c - self.d*other.d
        b = self.a*other.b + self.b*other.a + self.c*other.d - self.d*other.c
        c = self.a*other.c - self.b*other.d + self.c*other.a + self.d*other.b
        d = self.a*other.d + self.b*other.c - self.c*other.b + self.d*other.a
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
    
    def norm_q(self):
        '''
        compute normalized Quaternion
        '''
        mynorm = self.norm()
        my_norm_q = Quaternion(self.a/mynorm, self.b/mynorm, self.c/mynorm, self.d/mynorm)
        return my_norm_q
    
    def conj(self):
        '''
        compute Quaternion object complex conjugate
        '''
        a = self.a
        b = -self.b
        c = -self.c
        d = -self.d
        return Quaternion(a, b, c, d)
    
    def rotator(self, theta, vectors):
        '''
        from angle and vectors, compute a quaternion
        
        arguments:
            theta -- rotation angle, radians
            vectors -- indicates rotation aixs, list, like [1, 0, 0]
        '''
        
        # normalize the vector to nearly 1
        sum_v = sum([ v*v for v in vectors])
        #if abs(sum_v - 1.0) > 0.00000001:
        norm_v = math.sqrt(sum_v)
        vectors = [v/norm_v for v in vectors]
            
        a = math.cos(theta/2.)
        b = vectors[0]*math.sin(theta/2.)
        c = vectors[1]*math.sin(theta/2.)
        d = vectors[2]*math.sin(theta/2.)
        
        r = Quaternion(a, b, c, d)
        
        return r*self*r.conj()
        
    def __str__(self):
        ''' document printing'''
        parameters = {'':self.a, 'i':self.b, 'j':self.c, 'k':self.d}
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
    
