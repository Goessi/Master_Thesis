# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 09:56:28 2018

test on quaternion

@author: JingQIN
"""

import matplotlib.pyplot as plt
from Quaternion import Quaternion
import functions as f
import numpy as np
import csv

print("-------Quaternion-----------------------")
X1, X2, MEANDIFF, MINDIFF, MAXDIFF, TIME, TOTALTIME = f.Quaternion_rotation_precision(3000, 1000, np.pi, np.pi, np.pi)
np.savez("Quaternion_rotation_precision.npz", X1 = X1, X2 = X2, MEANDIFF = MEANDIFF, MINDIFF = MINDIFF, MAXDIFF = MAXDIFF, TIME = TIME, TOTALTIME = TOTALTIME)


print("-------DCM-----------------------")
X1, X2, MEANDIFF, MINDIFF, MAXDIFF, TIME, TOTALTIME = f.DCM_rotation_precision(3000, 1000, np.pi, np.pi, np.pi)
np.savez("DCM_rotation_precision.npz", X1 = X1, X2 = X2, MEANDIFF = MEANDIFF, MINDIFF = MINDIFF, MAXDIFF = MAXDIFF, TIME = TIME, TOTALTIME = TOTALTIME)

