# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 14:28:42 2017

@author: TPflumm
"""
from OCC.Geom2dAPI import Geom2dAPI_PointsToBSpline

from SONATA.topo.BSplineLst_utils import get_BSplineLst_Pnt2d, intersect_BSplineLst_with_BSpline
from SONATA.topo.utils import point2d_list_to_TColgp_Array1OfPnt2d

class Web(object):
    """
    Does nothing more than demonstrate syntax.
    
    This is an example of how a Pythonic human-readable docstring can
    get parsed by doxypypy and marked up with Doxygen commands as a
    regular input filter to Doxygen.
    
    Args:
    arg1: A positional argument.
    arg2: Another positional argument.
    
    Kwargs:
    kwarg: A keyword argument.
    
    Returns:
    A string holding the result.
    
    Raises:
    ZeroDivisionError, AssertionError, & ValueError.
    
    Examples:
    >>> myfunction(2, 3)
    '5 - 0, whatever.'
    >>> myfunction(5, 0, 'oops.')
    Traceback (most recent call last):
    ...
    ZeroDivisionError: integer division or modulo by zero
    >>> myfunction(4, 1, 'got it.')
    '5 - 4, got it.'
    >>> myfunction(23.5, 23, 'oh well.')
    Traceback (most recent call last):
    ...
    AssertionError
    >>> myfunction(5, 50, 'too big.')
    Traceback (most recent call last):
    ...
    ValueError
    """

   
    def __init__(self, ID, Pos1, Pos2, Segment0_BSplineLst, Segment0_Boundary_BSplineLst):
        self.ID = ID
        self.Pos1 = Pos1
        self.Pos2 = Pos2
        self.Segment0_BSplineLst = Segment0_BSplineLst
        self.Segment0_Boundary_BSplineLst = Segment0_Boundary_BSplineLst
        
        self.Pos1_Pnt2d = get_BSplineLst_Pnt2d(self.Segment0_BSplineLst,self.Pos1,0,1)
        self.Pos2_Pnt2d = get_BSplineLst_Pnt2d(self.Segment0_BSplineLst,self.Pos2,0,1)
        self.BSpline_Line = Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([self.Pos1_Pnt2d,self.Pos2_Pnt2d])).Curve().GetObject()
        self.BSplineLst = [self.BSpline_Line]
        
        self.wr_nodes =[] 
        self.wl_nodes = []
        
        #Intersect Segment0_Boundary_BSplineLst with self.BSpline
        [self.IntPnts,self.IntPnts_Pnt2d] = intersect_BSplineLst_with_BSpline(self.Segment0_Boundary_BSplineLst,self.BSpline_Line)
        
        