# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 22:04:11 2018

first version of Master thesis codes, initial model

@author: JingQIN
"""

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
        return Quaternion(self.a+other.a, self.b+other.b, self.c+other.c, self.d+other.d)
    
    def __str__(self):
        ''' document printing'''
        return str(self.a) + ' + ' + str(self.b) + 'i + ' + str(self.c) + 'j + ' + str(self.d) + 'k'
    
