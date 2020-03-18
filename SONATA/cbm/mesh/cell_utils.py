# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 14:10:22 2018

@author: TPflumm
"""

# Third party modules
import numpy as np

# First party modules
from SONATA.cbm.topo.utils import calc_angle_between


def calc_cell_angles(cell):
    corners = []
    for node in cell.nodes:
        corners.append(node.coordinates)
    corners = np.asarray(corners)
    temp = []
    for i in range(0, corners.shape[0]):
        if i == corners.shape[0] - 1:  # last point
            v1 = corners[i - 1] - corners[i]
            v2 = corners[0] - corners[i]
        else:
            v1 = corners[i - 1] - corners[i]
            v2 = corners[i + 1] - corners[i]
        temp.append(calc_angle_between(v1, v2))
    return np.array(temp)
