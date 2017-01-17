# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 14:28:42 2017

@author: TPflumm
"""
from OCC.Geom2dAPI import Geom2dAPI_PointsToBSpline

from BSplineLst_utils import get_BSplineLst_Pnt2d, intersect_BSplineLst_with_BSpline
from utils import point2d_list_to_TColgp_Array1OfPnt2d

class Web(object):
   
    def __init__(self, ID, Pos1, Pos2, Segment0_BSplineLst, Segment0_Boundary_BSplineLst):
        self.ID = ID
        self.Pos1 = Pos1
        self.Pos2 = Pos2
        self.Segment0_BSplineLst = Segment0_BSplineLst
        self.Segment0_Boundary_BSplineLst = Segment0_Boundary_BSplineLst
        
        self.Pos1_Pnt2d = get_BSplineLst_Pnt2d(self.Segment0_BSplineLst,self.Pos1,0,1)
        self.Pos2_Pnt2d = get_BSplineLst_Pnt2d(self.Segment0_BSplineLst,self.Pos2,0,1)
        self.BSpline_Line = Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([self.Pos1_Pnt2d,self.Pos2_Pnt2d])).Curve().GetObject()
        
        #Intersect Segment0_Boundary_BSplineLst with self.BSpline
        [self.IntPnts,self.IntPnts_Pnt2d] = intersect_BSplineLst_with_BSpline(self.Segment0_Boundary_BSplineLst,self.BSpline_Line)
                