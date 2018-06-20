# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 16:37:54 2018

@author: TPflumm
"""
#==============================================================================
#TEST TRIANGLE Package
#==============================================================================
import triangle
import triangle.plot
import numpy as np
import matplotlib.pyplot as plt
R = 10
theta = np.linspace(0, -2*np.pi, 36)[:-1]
pts = np.vstack((R*np.cos(theta), R*np.sin(theta))).T
A = dict(vertices=pts)
B = triangle.triangulate(A, 'qa1')
triangle.plot.compare(plt, A, B)
plt.show()

