#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 10:13:43 2019

@author: gu32kij
"""
# Third party modules
from OCC.Core.GeomAdaptor import GeomAdaptor_Curve
from OCC.Core.GeomAPI import geomapi_To3d
from OCC.Core.gp import gp_Pnt, gp_Vec


def bsplinelst_to3d(BSplineLst2d, pln):
    """ Transforms the 2d BSplineLst to a 3d BSplineLst on the plane pln
    
    Parameters
    ----------
    BSplineLst : [Geom2d_BSplineCurve]
        List of a Geom2d_BplineCurves 
        
    pln : gp_Pln 
        
    Returns
    -------
    BSplineLst3d : [Geom_BSplineCurve]
        list of 3d Geom_BSplinecurves
    
    """
    curves = [geomapi_To3d(s, pln) for s in BSplineLst2d]
    return [GeomAdaptor_Curve(c).BSpline() for c in curves]


def pnt_to3d(pnt2d):
    """
    transforms a 2d point to th 3d point on the xy plane (z=0)
    
    """
    return gp_Pnt(pnt2d.X(), pnt2d.Y(), 0)


def vec_to3d(vec2d):
    """
    transforms a 2d vec to th 3d point on the xy plane (z=0)
    
    """
    return gp_Vec(vec2d.X(), vec2d.Y(), 0)
