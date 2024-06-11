# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 14:28:42 2017

@author: TPflumm
"""
from OCC.Core.Geom2dAPI import Geom2dAPI_PointsToBSpline
from OCC.Core.gp import gp_Pnt2d, gp_Vec2d

from SONATA.cbm.topo.BSplineLst_utils import find_BSplineLst_coordinate
from SONATA.cbm.topo.utils import point2d_list_to_TColgp_Array1OfPnt2d
from SONATA.cbm.topo.para_Geom2d_BsplineCurve import ParaLst_from_BSplineLst, BSplineLst_from_ParaLst

from OCC.Core.Geom2dConvert import geom2dconvert_CurveToBSplineCurve
from OCC.Core.Geom2d import Geom2d_BezierCurve


#
# from OCC.Display.SimpleGui import init_display
# display, start_display, add_menu, add_function_to_menu = init_display()
# from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge2d

class Web(object):
      #Build Webs:
        #TODO: CHECK IF WEB DEFINITION INTERSECT EACH OTHER
        #TODO: SORT WEBS BY POS1 VALUES:
    def __init__(self, ID, Pos1, Pos2, curvature=0, SegmentLst=None):
        """

        ID:         web id (?)
        Pos1:       Pos of start point of web
        Pos2:       Pos of end point of web
        curvature:  value that defines the chordwise displacement of the ellipsis, i.e. the curvature of a web
        SegmentLst: List of segments

        """

        self.ID = ID
        self.Pos1 = Pos1
        self.Pos2 = Pos2
        
        if curvature is None:
            self.curvature = 0
        else: self.curvature = curvature
        
        self.Pos1_Pnt2d = SegmentLst[0].get_Pnt2d(self.Pos1, SegmentLst)    # p1 - start point at arc
        self.Pos2_Pnt2d = SegmentLst[0].get_Pnt2d(self.Pos2, SegmentLst)    # p2 - end point at arc

        BSplineLst_straight = [Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([self.Pos1_Pnt2d, self.Pos2_Pnt2d])).Curve()]


        if self.curvature == 0:  # if no curvature is assigned or curvature is assigned to be zero
            self.BSplineLst = BSplineLst_straight
        else:
            self.BSplineLst = bend_web_BSplineLst(BSplineLst_straight, self.Pos1_Pnt2d, self.Pos2_Pnt2d, self.curvature)
            
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
        

def bend_web_BSplineLst(BSplineLst_straight, Pos1_Pnt2d, Pos2_Pnt2d, curvature=0.0):
    """
    

    Parameters
    ----------
    BSplineLst_straight : TYPE
        DESCRIPTION.
    Pos1_Pnt2d : TYPE
        DESCRIPTION.
    Pos2_Pnt2d : TYPE
        DESCRIPTION.
    curvature : TYPE, optional
        DESCRIPTION. The default is 0.

    Returns
    -------
    web_BSplineLst : TYPE
        DESCRIPTION.

    """

            
    # Determine mid point between start and end arc position and translate by given curvature offset
    # len = get_BSplineLst_length(self.BSplineLst)
    [i_mid, u_mid] = find_BSplineLst_coordinate(BSplineLst_straight, 0.5, 0, 1)  # get midpoint of web
    s_mid = BSplineLst_straight[i_mid]  # allocate spline location of start point
    p_mid = gp_Pnt2d()          # init p - Mid point coordinates
    v_mid = gp_Vec2d()          # init v - normal vector
    s_mid.D1(u_mid, p_mid, v_mid)
    n_mid_pt = v_mid.GetNormal()
    n_mid_pt.Normalize()
    n_mid_pt.Multiply(curvature)
    p_curvature = p_mid.Translated(n_mid_pt)    # p_curvature - offset point in the middle of the straight line that connects the start and end arc positions

    # Option 1 - connect points directly via BSpline
    # self.BSplineLst = [Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([self.Pos1_Pnt2d, p_curvature, self.Pos2_Pnt2d])).Curve()]

    # ---------- S ------------
    #            x
    #            x
    #            x
    #            x
    #            x
    #            x
    # ---------- E -----------


    # Option 2 - Determine Bezier curve
    bezier_mount_pts_offset = 0.05  # determine vertical offset from start & end bezier mounting points; greater than 0 in order to account for Suction- and pressure innerside curvatures of innermost plies in segment 0
    # translate START arc position by curvature
    [i_start_arc_offset, u_start_arc_offset] = find_BSplineLst_coordinate(BSplineLst_straight, bezier_mount_pts_offset, 0, 1)
    s_start_arc_offset = BSplineLst_straight[i_start_arc_offset]    # allocate spline location of start point
    p_start_arc_offset = gp_Pnt2d()          # init p - Start arc offset point coordinates
    v_start_arc_offset = gp_Vec2d()          # init v - normal vector
    s_start_arc_offset.D1(u_start_arc_offset, p_start_arc_offset, v_start_arc_offset)
    n_start_arc_offset = v_start_arc_offset.GetNormal()
    n_start_arc_offset.Normalize()
    n_start_arc_offset.Multiply(curvature)
    p_start_bezier = p_start_arc_offset.Translated(n_start_arc_offset)    # p_start_bezier - offset point
    # translate END arc position by curvature
    [i_end_arc_offset, u_end_arc_offset] = find_BSplineLst_coordinate(BSplineLst_straight, 1-bezier_mount_pts_offset, 0, 1)
    s_end_arc_offset = BSplineLst_straight[i_end_arc_offset]    # allocate spline location of start point
    p_end_arc_offset = gp_Pnt2d()          # init p - Start arc offset point coordinates
    v_end_arc_offset = gp_Vec2d()          # init v - normal vector
    s_end_arc_offset.D1(u_end_arc_offset, p_end_arc_offset, v_end_arc_offset)
    n_end_arc_offset = v_end_arc_offset.GetNormal()
    n_end_arc_offset.Normalize()
    n_end_arc_offset.Multiply(curvature)
    p_end_bezier = p_end_arc_offset.Translated(n_end_arc_offset)    # p_end_bezier - offset point


    Bezier_PntList = [Pos1_Pnt2d, p_start_bezier, p_curvature, p_end_bezier, Pos2_Pnt2d]
    # Bezier_PntList = [self.Pos1_Pnt2d, p_start_bezier, p_curvature, self.Pos2_Pnt2d]
    Bezier = Geom2d_BezierCurve(point2d_list_to_TColgp_Array1OfPnt2d(Bezier_PntList))
    

    # spline = BRepBuilderAPI_MakeEdge2d(self.BSplineLst)
    # spline.Build()
    # display.DisplayShape(spline.Shape(), update=True)


    #     --------- x  S     -----------    S  - Start arc position at arc
    # Ps        x      |                    Ps - p_start_bezier mounting point: orthogonal offset from S by curvature value; additional offset along web direction by bezier_mount_pts_offset
    #        x         |
    #      x           |
    #     x            |
    #     x            |
    # P1  x <----------o                    P1 - p_curvature: orthogonal offset from web mid point by curvature value
    #     x            |
    #     x            |
    #      x           |
    #        x         |
    # Pe        x      |                    Pe - p_end_bezier mounting point: orthogonal offset from S by curvature value; additional offset along web direction by bezier_mount_pts_offset
    #     --------- x  E     -----------    E  - End arc position at arc  
    
    web_BSplineLst = [geom2dconvert_CurveToBSplineCurve(Bezier)]
    return web_BSplineLst