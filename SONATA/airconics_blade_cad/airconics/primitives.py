# -*- coding: utf-8 -*-
"""
Airfoil primitive class definitions

Created on Fri Dec  4 13:31:35 2015

@author: pchambers
"""
from OCC.Core.GC import GC_MakeSegment
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_OX, gp_OY

# import AirCONICStools as act
from pkg_resources import resource_string, resource_exists
import numpy as np


import sys
sys.path.insert(0,'/Users/rfeil/work/6_SONATA/SONATA')  # import sys path to import 'SONATA' & 'job' modules
import SONATA.blade_cad.airconics.AirCONICStools as act

# Classes
# -----------------------------------------------------------------------------
class Airfoil(object):
    """Class for defining a range of spline-fitted airfoil curves

        Parameters
        ----------

        Attributes
        ----------
        points : array of scalar, shape (N, 3)
            The x-y-z coordinates of points on the airfoils surface

        Curve - OCC.Geom.Handle_Geom_BsplineCurve
            The generated airfoil spline

        """
    def __init__(self,
                 inpPoints):

        # Initialise the x and z airfoil surface points as None
        self._points = inpPoints
        self.Curve = self._fitAirfoiltoPoints()

    @property
    def points(self):
        return self._points

    @points.setter
    def points(self, newpoints):
        # Updating this value recalls the _fitAirfoiltoPoints function. Does
        # not transform the airfoil
        self._points = newpoints

    def _fitAirfoiltoPoints(self):
        """ Fits an OCC curve to current airfoil.points.

        airfoil.points should be N by 2, the array of x, z points on the
        airfoil's surface.

        Returns
        -------
        Curve : OCC.Geom.Geom_BSplineCurve
            the generated spline
        """
        points3d = np.column_stack([self._points[:, 0], self._points[:, 1], self._points[:, 2]])
        Curve = act.points_to_bspline(points3d)
        return Curve
