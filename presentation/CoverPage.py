# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 15:46:27 2019

@author: JingQIN
"""

import matplotlib.pyplot as plt
import numpy as np
import math
import mpl_toolkits.mplot3d.axes3d as p3
import random

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
def vector_location(v, x_theta, y_theta, z_theta):
    '''
    calculate locaton for vector
    
    arguments:
        v[x,y,z]: rotated vector
        x_theta: total rotation angle around x-axis, radian
        y_theta: total rotation angle around y-axis, radian
        z_theta: total rotation angle around z-axis, radian
    
    '''
        
    v0 = v
    V = [v]
    steps = 180
    m1 = computeDCM(x_theta / steps, [1, 0, 0])
    m2 = computeDCM(y_theta / steps, [0, 1, 0])
    m3 = computeDCM(z_theta / steps, [0, 0, 1])
    for i in range(0, steps):
        v = np.dot(m1, v)
        V.append(v)
    for j in range(0, steps):
        v = np.dot(m2, v)
        V.append(v)
    for k in range(0, steps):
        v = np.dot(m3, v)
        V.append(v)
    V.append(v0)
    return V

def DCM_rotation_precision(x_theta, y_theta, z_theta):
    '''
    calculate precision for DCM rotation
    
    arguments:
        x_theta: total rotation angle around x-axis, radian
        y_theta: total rotation angle around y-axis, radian
        z_theta: total rotation angle around z-axis, radian
    
    '''
    X1 = []
    Y1 = []
    Z1 = []
    X2 = []
    Y2 = []
    Z2 = []
    X3 = []
    Y3 = []
    Z3 = []
    x = 0.5
    y = 0.5
    z = 0.5
    v = [x,y,z]
    X1.append(x)
    Y1.append(y)
    Z1.append(z)
    steps = 180
    m1 = computeDCM(x_theta / steps, [1, 0, 0])
    m2 = computeDCM(y_theta / steps, [0, 1, 0])
    m3 = computeDCM(z_theta / steps, [0, 0, 1])
    for i in range(0, steps):
        v = np.dot(m1, v)
        X1.append(v[0])
        Y1.append(v[1])
        Z1.append(v[2])
    for j in range(0, steps):
        v = np.dot(m2, v)
        X2.append(v[0])
        Y2.append(v[1])
        Z2.append(v[2])
    for k in range(0, steps):
        v = np.dot(m3, v)
        X3.append(v[0])
        Y3.append(v[1])
        Z3.append(v[2])
    X3.append(0.5)
    Y3.append(0.5)
    Z3.append(0.5)
    return X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3

x1, y1, z1, x2, y2, z2, x3, y3, z3 = DCM_rotation_precision(np.pi, np.pi, np.pi)
x1 = np.array(x1)
y1 = np.array(y1)
z1 = np.array(z1)
x2 = np.array(x2)
y2 = np.array(y2)
z2 = np.array(z2)
x3 = np.array(x3)
y3 = np.array(y3)
z3 = np.array(z3)

VX = vector_location([1,0,0],2*np.pi, 2*np.pi, 2*np.pi)
VY = vector_location([0,1,0],2*np.pi, 2*np.pi, 2*np.pi)
VZ = vector_location([0,0,1],2*np.pi, 2*np.pi, 2*np.pi)

fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')
ax1.set_facecolor('xkcd:black')
ax1.set_xlim(-1.1, 1.1)
ax1.set_ylim(-1.1, 1.1)
ax1.set_zlim(-1.1, 1.1)
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_zlabel("z")
plt.rcParams['animation.html'] = 'html5'
ax1.view_init(elev=-20., azim=235)
ax1.set_aspect('equal')
ax1.axis('off')
ax1.plot([0, 0], [0,2.5],[0,0],color='green')
ax1.plot([0,0],[2.5,2.4],[0,0.1],color='green')
ax1.plot([0,0],[2.5,2.4],[0,-0.1],color='green')

ax1.plot([0, 0], [0,0],[0,1.7],color='red')
ax1.plot([0,0.07],[0,0],[1.7,1.6],color='red')
ax1.plot([0,-0.07],[0,0],[1.7,1.6],color='red')

ax1.plot([0,2], [0,0],[0,0],color='blue')
ax1.plot([2,1.9],[0,0],[0,0.1],color='blue')
ax1.plot([2,1.9],[0,0],[0,-0.1],color='blue')
for j in range(0,len(VX),10):
    ax1.plot([0, VX[j][0]], [0,VY[j][0]],[0,VZ[j][0]],color='paleturquoise')
    ax1.plot([0, VX[j][1]], [0,VY[j][1]],[0,VZ[j][1]],color='palegreen')
    ax1.plot([0, VX[j][2]], [0,VY[j][2]],[0,VZ[j][2]],color='lightcoral')
for j in range(0, len(x1)):
    ax1.scatter(x1[:j], y1[:j], z1[:j], '.', color='steelblue',s=2)
for j in range(0, len(x2)):
    ax1.scatter(x2[:j], y2[:j], z2[:j], '.', color='lime',s=2)
for j in range(0, len(x3)):
    ax1.scatter(x3[:j], y3[:j], z3[:j], '.', color='maroon',s=2)
#ax2 = fig.add_subplot(121, projection='3d')
#ax2.set_xlim(-1.1, 1.1)
#ax2.set_ylim(-1.1, 1.1)
#ax2.set_zlim(-1.1, 1.1)
#ax2.set_xlabel("x")
#ax2.set_ylabel("y")
#ax2.set_zlabel("z")
#plt.rcParams['animation.html'] = 'html5'
#ax2.view_init(elev=-20., azim=235)
#ax2.set_aspect('equal')
#ax2.axis('off')
#
#ax2.plot([0, 0], [0,2.5],[0,0],color='green')
#ax2.plot([0,0],[2.5,2.4],[0,0.1],color='green')
#ax2.plot([0,0],[2.5,2.4],[0,-0.1],color='green')
#
#ax2.plot([0, 0], [0,0],[0,1.7],color='red')
#ax2.plot([0,0.07],[0,0],[1.7,1.6],color='red')
#ax2.plot([0,-0.07],[0,0],[1.7,1.6],color='red')
#
#ax2.plot([0,2], [0,0],[0,0],color='blue')
#ax2.plot([2,1.9],[0,0],[0,0.1],color='blue')
#ax2.plot([2,1.9],[0,0],[0,-0.1],color='blue')
#for j in range(0, len(x1)):
#    ax2.scatter(x1[:j], y1[:j], z1[:j], '.', color='steelblue',s=2)
#for j in range(0, len(x2)):
#    ax2.scatter(x2[:j], y2[:j], z2[:j], '.', color='lime',s=2)
#for j in range(0, len(x3)):
#    ax2.scatter(x3[:j], y3[:j], z3[:j], '.', color='maroon',s=2)

plt.savefig("coverPage2.png", dpi=400)
