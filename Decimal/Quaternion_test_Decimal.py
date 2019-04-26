# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 09:56:28 2018

test on quaternion

@author: JingQIN
"""

import matplotlib.pyplot as plt
import functions_Decimal as f
import numpy as np


print("-------Quaternion-----------------------")
(X1, X2, MEANDIFF, MINDIFF, MAXDIFF, TIME) = f.Quaternion_rotation_precision(3000, 1000, np.pi, np.pi, np.pi)
np.savez("Quaternion_rotation_precision.npz", X1 = X1, X2 = X2, MEANDIFF = MEANDIFF, MINDIFF = MINDIFF, MAXDIFF = MAXDIFF, TIME = TIME)
print("-------DCM-----------------------")
(X1, X2, MEANDIFF, MINDIFF, MAXDIFF, TIME) = f.DCM_rotation_precision(3000, 1000, np.pi, np.pi, np.pi)
np.savez("DCM_rotation_precision.npz", X1 = X1, X2 = X2, MEANDIFF = MEANDIFF, MINDIFF = MINDIFF, MAXDIFF = MAXDIFF, TIME = TIME)

(X1, X2, MEANDIFF_Q, MINDIFF_Q, MAXDIFF_Q, MEANDIFF_DCM, MINDIFF_DCM, MAXDIFF_DCM) = f.Quaternion_DCM_rotation_precision(3000, 1000, np.pi, np.pi, np.pi)
np.savez("Quaternion_DCM_rotation_precision.npz", X1 = X1, X2 = X2, MEANDIFF_Q = MEANDIFF_Q, MINDIFF_Q = MINDIFF_Q, MAXDIFF_Q = MAXDIFF_Q, MEANDIFF_DCM = MEANDIFF_DCM, MINDIFF_DCM = MINDIFF_DCM, MAXDIFF_DCM = MAXDIFF_DCM)



