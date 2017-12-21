# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 14:28:42 2017

@author: TPflumm
"""
from OCC.Geom2dAPI import Geom2dAPI_PointsToBSpline

from SONATA.topo.BSplineLst_utils import get_BSplineLst_Pnt2d, intersect_BSplineLst_with_BSpline
from SONATA.topo.utils import point2d_list_to_TColgp_Array1OfPnt2d

class Web(object):
      #Build Webs:
        #TODO: CHECK IF WEB DEFINITION INTERSECT EACH OTHER
        #TODO: SORT WEBS BY POS1 VALUES:
    def __init__(self, ID, Pos1, Pos2, Segment0):
        self.ID = ID
        self.Pos1 = Pos1
        self.Pos2 = Pos2
        
        lid = Segment0.LayerLst[-1].ID
        self.Pos1_Pnt2d = Segment0.get_Pnt2d(lid,Pos1)
        self.Pos2_Pnt2d = Segment0.get_Pnt2d(lid,Pos2)
        self.BSplineLst = [Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([self.Pos1_Pnt2d,self.Pos2_Pnt2d])).Curve().GetObject()]

        self.wr_nodes =[] 
        self.wl_nodes = []
        
        #self.BackBSplineLst = trim_BSplineLst(BSplineLst, S1, S2, start, end):
        #Intersect Segment0_Boundary_BSplineLst with self.BSpline
        #[self.IntPnts,self.IntPnts_Pnt2d] = intersect_BSplineLst_with_BSpline(self.Segment0_Boundary_BSplineLst,self.BSpline_Line)
