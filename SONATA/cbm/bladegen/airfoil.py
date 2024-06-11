# -*- coding: utf-8 -*-
"""
Created on Mon Sep 04 10:45:43 2017

@author: TPflumm
"""
# Third party modules
from OCC.Core.BRepBuilderAPI import (BRepBuilderAPI_MakeEdge,
                                     BRepBuilderAPI_MakeWire,)
from OCC.Core.Geom import Geom_BezierCurve
from OCC.Core.GeomAPI import GeomAPI_Interpolate
# PythonOCC Libraries
from OCC.Core.gp import gp_Pnt, gp_Vec

# First party modules
# SONATA modules:
from SONATA.cbm.fileIO.readinput import UIUCAirfoil
from SONATA.cbm.topo.utils import (TColgp_HArray1OfPnt_from_nparray,
                                   point_list_to_TColgp_Array1OfPnt,)


class Airfoil(object):
    class_counter = 0

    def __init__(self, name=None, trailing_edge_style=0):
        self.id = self.__class__.class_counter
        self.__class__.class_counter += 1