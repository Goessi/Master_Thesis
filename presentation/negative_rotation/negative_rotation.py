# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 10:13:27 2019

@author: JingQIN
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import math
import mpl_toolkits.mplot3d.axes3d as p3

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

VX = vector_location([1,0,0],np.pi, np.pi, np.pi)
VY = vector_location([0,1,0],np.pi, np.pi, np.pi)
VZ = vector_location([0,0,1],np.pi, np.pi, np.pi)


    
    

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
for j in range(len(VX)):
    ax.plot([0, VX[j][0]], [0,VY[j][0]],[0,VZ[j][0]],color='steelblue')
    ax.plot([0, VX[j][1]], [0,VY[j][1]],[0,VZ[j][1]],color='lime')
    ax.plot([0, VX[j][2]], [0,VY[j][2]],[0,VZ[j][2]],color='maroon')
    plt.savefig("movie%d.png" %j)
    del ax.lines[0]
    del ax.lines[0]
    del ax.lines[0]

#line_ani = animation.ArtistAnimation(fig, lns, interval=1000, blit=True)
#
#line_ani.save('passive_rotation.gif',writer='imagemagick',fps=1000/100)