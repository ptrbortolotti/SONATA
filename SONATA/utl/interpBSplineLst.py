#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 15:22:29 2019

@author: gu32kij
"""
import numpy as np
from scipy.interpolate import interp1d

from OCC.Geom import Geom_Plane
from OCC.gp import gp_Dir, gp_Pnt, gp_Trsf, gp_Vec, gp_Ax2, gp_Ax3

from SONATA.utl.blade_utl import check_uniformity 
from SONATA.cbm.topo.BSplineLst_utils import intersect_BSplineLst_with_plane

class interpBSplineLst(object):
    """
    
    
    """
    
    def __init__(self, BSplineLst, xgrid, xvalues):
        """
        Parameters
        --------
        BSplineLst : 
        xgrid : array
        xvalues : array
        """
        
        self.BSplineLst = BSplineLst
        if check_uniformity(xgrid, xvalues) == False:
                print('WARNING:\t The reference axis is not uniformly defined along x!')                
        self._f_xint = interp1d(xgrid, xvalues, bounds_error=False, fill_value='extrapolate')

    
    def interpolate(self, grid):
        """
        interpolates the BSplineLst at a certain radial station by intersecting
        with a plane.
        
        Parameters
        --------
        grid : float
            nondimensional x coordinate in the blade ref. frame
        
        Retruns
        --------
        res, resCoords : (array,array)
            tuple of arrays. res contains the coordinates, while resCoords 
            contains the BSpline List parameter coordinates i and u.
            
        """
        
        if not isinstance(grid,np.ndarray):
            grid = np.atleast_1d(grid)
        res = []
        resCoords = []
        for gp in grid:
            plane = Geom_Plane(gp_Pnt(float(self._f_xint(gp)),0,0), gp_Dir(1,0,0))
            IntCoords, IntPnts = intersect_BSplineLst_with_plane(self.BSplineLst, plane)
            if len(IntPnts) > 1:
                print('WARNING:\t More than one intersection point was found!')
        
            try:
                res.append(IntPnts[0].Coord())
                resCoords.append(IntCoords[0])
            except:
                print('No Intersection Found')
                raise Exception

        return (np.asarray(res),np.asarray(resCoords))
    
    
    def interpolate_curvature(self, grid, ax2=None):
        """
        interpolates the BSplineLst at a certain radial station and returns the 
        curvature kappa2 and kappa3 in the local Ax2 frame if the option is 
        passed.   
        
        
        Parameters
        --------
        grid : float
            nondimensional x coordinate in the blade ref. frame
        ax2 : gp_Ax2, optional
            right handed coordinate system. 
        
        Returns
        --------
        (k2,k3) : tuple
            tuple of moment strain measures (curvature) about the given 
            coorindate systems (x2 and x3) axis.
                    
        """        
        p = gp_Pnt()
        v1 = gp_Vec()
        v2 = gp_Vec()
        
        para = self.interpolate(grid)[1]
        i = int(para[0][0])
        u = para[0][1]
        s = self.BSplineLst[i]
        s.D2(u,p,v1,v2)
        
        #transform curvature vector to local Ax2
        if ax2:
            trsf = gp_Trsf()
            trsf.SetTransformation(gp_Ax3(ax2))
            k_vec = v2.Transformed(trsf)
        else:
            k_vec = v2
            
        #swap the axis of the curvature vector, return moment strain measures. 
        k2 = k_vec.Coord()[2] 
        k3 = k_vec.Coord()[1]
        return (k2,k3)
                
    
if __name__ == '__main__':
    pass