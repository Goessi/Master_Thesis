# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 09:56:28 2018

test on quaternion

@author: JingQIN
"""

from Quaternion import Quaternion

print("-------test on addition operation-------")
A = Quaternion(1.0,2.0,3.0,4.0)
B = Quaternion(5.0,6.0,7.0,8.0)
C = A + B
print(A)
print(B)
print(C)