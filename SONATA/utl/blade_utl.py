#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 10:20:51 2019

@author: gu32kij
"""
import numpy as np
from scipy.interpolate import interp1d

#PythonOCC Libraries
from OCC.TopoDS import TopoDS_Wire,  TopoDS_Vertex
from OCC.BRepOffsetAPI import BRepOffsetAPI_ThruSections



def interp_loads(loads, grid_loc):
    """
    Interpolates the loads at the given radial station (grid location)
    
    Parameters
    ----------
    loads : dict
        dictionary of the following keys and values, (default=None)
        for detailed information see the VABSConfig documentation or the 
        VABS user manual
        F : nparray([[grid, F1, F2, F3]]) 
        M : nparray([[grid, M1, M2, M3]]) 
        f : nparray([[grid, f1, f2, f2]])
        df : nparray([[grid, f1', f2', f3']])
        dm :  nparray([[grid, m1', m2', m3']])
        ddf : nparray([[grid, f1'', f2'', f3'']])
        ddm : nparray([[grid, m1'', m2'', m3'']])
        
    grid_loc : float
        location of interpolation
        
    Returns
    ----------
    sectional_load : dict
    dictionary of the following keys and values, (default=None)
        for detailed information see the VABSConfig documentation or the 
        VABS user manual
        F : nparray([F1, F2, F3]) 
        M : nparray([M1, M2, M3]) 
        f : nparray([f1, f2, f2])
        df : nparray([f1', f2', f3'])
        dm :  nparray([m1', m2', m3'])
        ddf : nparray([f1'', f2'', f3''])
        ddm : nparray([m1'', m2'', m3''])
        
    """
    
    d = {}
    for k,item in loads.items():
        fit = interp1d(item[:,0], item[:,1:], axis=0)
        d[k] = fit(grid_loc)
    return d

def interp_airfoil_position(airfoil_position, airfoils, grid_loc):
    """
    TBD:Docstring
    
    """
    
    #TBD: Extrapolation
    if grid_loc in airfoil_position[0]:
        afname = airfoil_position[1][airfoil_position[0].index(grid_loc)]
        return next((x for x in airfoils if x.name == afname), None)

    #find closest value:
    min_idx = np.argmin([abs(x-grid_loc) for x in airfoil_position[0]])
    min_val = airfoil_position[0][min_idx]
    
    if grid_loc > min_val:
        iv_idx = (min_idx,min_idx+1)
    else:
        iv_idx = (min_idx-1, min_idx)
    
    iv_val = tuple(airfoil_position[0][iv_idx[0]:iv_idx[1]+1])
    iv_af = tuple(airfoil_position[1][iv_idx[0]:iv_idx[1]+1])
    k = (grid_loc-iv_val[0]) / (iv_val[1]-iv_val[0])
    
    #select af from airfoils
    af1 = next((x for x in airfoils if x.name == iv_af[0]), None)
    af2 = next((x for x in airfoils if x.name == iv_af[1]), None)
    
    if af1 == af2:
        return af1
    
    #return transformed airfoil
    return af1.transformed(af2,k,200)


def make_loft(elements, ruled=False, tolerance=1e-6, continuity=4, check_compatibility=True):
    """
    A set of sections that are used to generate a surface with the 
        BRepOffsetAPI_ThruSections function from OCC
        
    Parameters
    ----------
    elements : list of OCC.TopoDS_Wire or TopoDS_Vertex
        A set of sections that are used to generate a surface with the 
        BRepOffsetAPI_ThruSections function from OCC
        
    ruled : bool
    tolerance : float
    continuity : int
    check_compatibility : bool
    
    Returns:
    ----------
    loft : TopoDS_Shape
        surface of the ThruSections Loft
    """
    
    sections = BRepOffsetAPI_ThruSections(False, ruled, tolerance)
    for i in elements:
        if isinstance(i, TopoDS_Wire):
            sections.AddWire(i)
        elif isinstance(i, TopoDS_Vertex):
            sections.AddVertex(i)
        else:
            raise TypeError('elements is a list of TopoDS_Wire or TopoDS_Vertex, found a %s fool' % i.__class__)

    sections.CheckCompatibility(check_compatibility)
    sections.SetContinuity(continuity)
    sections.Build()
    loft = sections.Shape()
    
    return loft
    