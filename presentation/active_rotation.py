# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 16:46:03 2019

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
    X = []
    Y = []
    Z = []
    x = random.random()
    y = random.random()
    z = random.random()
    v = [x, y, z]
    X.append(x)
    Y.append(y)
    Z.append(z)
    steps = 10
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
    X.append(x)
    Y.append(y)
    Z.append(z)
    return X, Y, Z

x, y, z = DCM_rotation_precision(30, 10, np.pi, np.pi, np.pi)
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

lns = []
lns.append(ax.plot(x[:0], y[:0], z[:0], 'o-', color='steelblue'))
for i in range(10):
    lns.append(ax.plot(x[:i+1], y[:i+1], z[:i+1], 'o-', color='steelblue'))
for i in range(10):
    lns.append(ax.plot(x[:i+11], y[:i+11], z[:i+11], 'o-', color='lime'))
for i in range(10):
    lns.append(ax.plot(x[:i+21], y[:i+21], z[:i+21], 'o-', color='gold'))
lns.append(ax.plot(x[:31], y[:31], z[:31], 'o-', color='steelblue'))

line_ani = animation.ArtistAnimation(fig, lns, interval=1000, blit=True)

line_ani.save('active_rotation.gif',writer='imagemagick',fps=1000/100)


