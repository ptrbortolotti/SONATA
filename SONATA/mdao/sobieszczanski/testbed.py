# -*- coding: utf-8 -*-
from __future__ import print_function

from openmdao.api import IndepVarComp, Component, Problem, Group
from openmdao.api import ScipyOptimizer
from openmdao.api import ExecComp
from openmdao.api import view_tree, view_connections

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
plt.close("all")

# Make data.
x1 = np.arange(0, 2.1, 0.1)
x2 = np.arange(-1, 1.1, 0.1)
x1, x2 = np.meshgrid(x1, x2)
f = (x1-1)**2 + x2**2

a = float(2)
g = x1-(x2**2)/a

# Plot the surface.
# You can force all the contours to be the same color.
plt.figure()
CS = plt.contour(x1, x2, f, colors='k')
plt.clabel(CS, fontsize=9, inline=1)
CS = plt.contour(x1, x2, g,20, colors='r')
plt.clabel(CS, fontsize=9, inline=1)
plt.title('Single color - negative contours dashed')
plt.show()

