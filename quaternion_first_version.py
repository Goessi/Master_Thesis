# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 22:04:11 2018

first version of Master thesis codes, initial model

@author: JingQIN
"""

class Quaternion(object):
    """
    
    class of Quaternion that do the simple operations
    
    """    
    def __init__(self, a, b, c, d):
        assert type(a) == float and type(b) == float and type(c) == float and type(d) == float
        
        self.a = a
        self.b = b
        self.c = c
        self.d = d
    
    def __add__(self, other):
        return Quaternion(self.a+other.a, self.b+other.b, self.c+other.c, self.d+other.d)
    
    def __str__(self):
        return str(self.a) + ' + ' + str(self.b) + 'i + ' + str(self.c) + 'j + ' + str(self.d) + 'k'
    
A = Quaternion(1.0,2.0,3.0,4.0)
B = Quaternion(5.0,6.0,7.0,8.0)
C = A + B
print(A)
print(B)
print(C)