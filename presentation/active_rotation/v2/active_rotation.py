# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 12:44:58 2019

@author: JingQIN
"""

import math
import random
import numpy as np
from numpy import array
from matplotlib import pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
# presentation data 

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
def DCM_rotation_precision(x_theta, y_theta, z_theta):
    '''
    calculate precision for DCM rotation
    
    arguments:
        x_theta: total rotation angle around x-axis, radian
        y_theta: total rotation angle around y-axis, radian
        z_theta: total rotation angle around z-axis, radian
    
    '''
    X = []
    Y = []
    Z = []
    x = 0.5
    y = 0.5
    z = 0.5
    v = [x, y, z]
    X.append(x)
    Y.append(y)
    Z.append(z)
    steps = 180
    m1 = computeDCM(x_theta / steps, [1, 0, 0])
    m2 = computeDCM(y_theta / steps, [0, 1, 0])
    m3 = computeDCM(z_theta / steps, [0, 0, 1])
    for i in range(0, steps):
        v = np.dot(m1, v)
        X.append(v[0])
        Y.append(v[1])
        Z.append(v[2])
    for j in range(0, steps):
        v = np.dot(m2, v)
        X.append(v[0])
        Y.append(v[1])
        Z.append(v[2])
    for k in range(0, steps):
        v = np.dot(m3, v)
        X.append(v[0])
        Y.append(v[1])
        Z.append(v[2])
    X.append(0.5)
    Y.append(0.5)
    Z.append(0.5)
    return X, Y, Z

def plot_axis(ax):
    ax.plot([0, 0], [0,2.5],[0,0],color='green')
    ax.plot([0,0],[2.5,2.4],[0,0.1],color='green')
    ax.plot([0,0],[2.5,2.4],[0,-0.1],color='green')
    
    
    ax.plot([0, 0], [0,0],[0,1.7],color='red')
    ax.plot([0,0.07],[0,0],[1.7,1.6],color='red')
    ax.plot([0,-0.07],[0,0],[1.7,1.6],color='red')
    
    ax.plot([0,2], [0,0],[0,0],color='blue')
    ax.plot([2,1.9],[0,0],[0,0.1],color='blue')
    ax.plot([2,1.9],[0,0],[0,-0.1],color='blue')
    
x, y, z = DCM_rotation_precision(np.pi, np.pi, np.pi)
x = array(x)
y = array(y)
z = array(z)
fig = plt.figure()
ax = p3.Axes3D(fig)

ax.set_xlim(-1.1, 1.1)
ax.set_ylim(-1.1, 1.1)
ax.set_zlim(-1.1, 1.1)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
plt.rcParams['animation.html'] = 'html5'
ax.view_init(elev=-20., azim=235)
ax.set_aspect('equal')
for j in range(0,183):
#    ax.plot([-2,2], [0,0],[0,0],color='steelblue')
#    ax.plot([0, 0], [-2,2],[0,0],color='lime')
#    ax.plot([0, 0], [0,0],[-2,2],color='maroon')
    plot_axis(ax)
    ax.plot(x[:j], y[:j], z[:j], 'o-', color='steelblue')
    plt.savefig("movie%d.png" %j)
for i in range(1,181):
#    ax.plot([-2,2], [0,0],[0,0],color='steelblue')
#    ax.plot([0, 0], [-2,2],[0,0],color='lime')
#    ax.plot([0, 0], [0,0],[-2,2],color='maroon')
    plot_axis(ax)
    ax.plot(x[:j+i], y[:j+i], z[:j+i], 'o-', color='lime')
    plt.savefig("movie%d.png" %(j+i))
for k in range(1,180):
#    ax.plot([-2,2], [0,0],[0,0],color='steelblue')
#    ax.plot([0, 0], [-2,2],[0,0],color='lime')
#    ax.plot([0, 0], [0,0],[-2,2],color='maroon')
    plot_axis(ax)
    ax.plot(x[:j+i+k], y[:j+i+k], z[:j+i+k], 'o-', color='maroon')
    plt.savefig("movie%d.png" %(j+i+k))
