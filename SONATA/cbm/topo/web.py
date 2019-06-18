# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 14:28:42 2017

@author: TPflumm
"""
from OCC.Geom2dAPI import Geom2dAPI_PointsToBSpline
from OCC.gp import gp_Pnt, gp_Pnt2d

from SONATA.cbm.topo.BSplineLst_utils import get_BSplineLst_Pnt2d, intersect_BSplineLst_with_BSpline
from SONATA.cbm.topo.utils import point2d_list_to_TColgp_Array1OfPnt2d
from SONATA.cbm.topo.para_Geom2d_BsplineCurve import ParaLst_from_BSplineLst, BSplineLst_from_ParaLst


class Web(object):
      #Build Webs:
        #TODO: CHECK IF WEB DEFINITION INTERSECT EACH OTHER
        #TODO: SORT WEBS BY POS1 VALUES:
    def __init__(self, ID, Pos1, Pos2, SegmentLst):
        self.ID = ID
        self.Pos1 = Pos1
        self.Pos2 = Pos2
        
        self.Pos1_Pnt2d = SegmentLst[0].get_Pnt2d(self.Pos1, SegmentLst)
        self.Pos2_Pnt2d = SegmentLst[0].get_Pnt2d(self.Pos2, SegmentLst)
        self.BSplineLst = [Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([self.Pos1_Pnt2d,self.Pos2_Pnt2d])).Curve()]

        self.wr_nodes =[] 
        self.wl_nodes = []
        self.wr_cells =[] 
        self.wl_cells = []
        
        #self.BackBSplineLst = trim_BSplineLst(BSplineLst, S1, S2, start, end):
        #Intersect Segment0_Boundary_BSplineLst with self.BSpline
        #[self.IntPnts,self.IntPnts_Pnt2d] = intersect_BSplineLst_with_BSpline(self.Segment0_Boundary_BSplineLst,self.BSpline_Line)
        
    def __getstate__(self):
        """Return state values to be pickled."""
        self.Para_BSplineLst = ParaLst_from_BSplineLst(self.BSplineLst)
        self.Pos1_coordinates = (self.Pos1_Pnt2d.X(),self.Pos1_Pnt2d.Y())
        self.Pos2_coordinates = (self.Pos2_Pnt2d.X(),self.Pos2_Pnt2d.Y())
        return (self.ID, self.Pos1, self.Pos2, self.Pos1_coordinates, self.Pos2_coordinates, self.Para_BSplineLst, self.wr_nodes, self.wl_nodes)   
    
    
    def __setstate__(self, state):
        """Restore state from the unpickled state values."""
        (self.ID, self.Pos1, self.Pos2, self.Pos1_coordinates, self.Pos2_coordinates, self.Para_BSplineLst, self.wr_nodes, self.wl_nodes)  = state
        self.BSplineLst = BSplineLst_from_ParaLst(self.Para_BSplineLst)
        self.Pos1_Pnt2d = gp_Pnt2d(self.Pos1_coordinates[0],self.Pos1_coordinates[1])
        self.Pos2_Pnt2d = gp_Pnt2d(self.Pos2_coordinates[0],self.Pos2_coordinates[1])
        
        