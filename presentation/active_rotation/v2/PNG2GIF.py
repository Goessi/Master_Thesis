# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 12:08:48 2019

@author: JingQIN
"""

import os
import imageio
images = []
filenames = []

for i in range(541):
    name = ("movie%d.png" %i)
    filenames.append(name)

for filename in filenames:
    images.append(imageio.imread(filename))

imageio.mimsave('movie.gif', images, duration=0.02)