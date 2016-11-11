# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 13:03:08 2016

@author: TPflumm
"""

import math
import numpy as np   

from OCC.gp import gp_Pnt, gp_OX, gp_OY,gp_OZ,gp_Vec, gp_Trsf, gp_DZ, gp_Ax2, gp_Ax3, gp_Pnt2d, gp_Dir2d, gp_Ax2d, gp_Dir, gp_XYZ, gp_Circ, gp_Ax1
from OCC.GC import GC_MakeArcOfCircle, GC_MakeSegment
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge,BRepBuilderAPI_MakeWire,BRepBuilderAPI_Transform
from OCC.TopoDS import topods, TopoDS_Edge

#from OCC.Display.SimpleGui import init_display
#display, start_display, add_menu, add_function_to_menu = init_display()

#Frequently Used Functions 
from core_geometry_utils import *
from core_operations_utils import *


def tab_definition(ref,tab_thickness,tab_length, tab_radius, tab_radius_angle):  
    
#    tab_thickness = 0.008  # times basic chord length
#    tab_length = 0.04      # times basic chord length
#    tab_radius = 0.02      # times basic chord length
#    tab_radius_angle = math.radians(-45)
#    ref = np.array([0.0, 0.75, 0.0])
               
    P1 = ref+(0,0,tab_thickness/2)
    P2 = ref+(0,-tab_length,tab_thickness/2)
    P3 = ref+(0,-tab_length,tab_thickness/2+tab_radius)
    
    
    #OCC - Environment
    ref_Ax = gp_Dir(0.0, -1.0, 0.0)
    ref_Pnt = np_array_to_gp_Pnt(ref) 
    P1_gpPnt = np_array_to_gp_Pnt(P1) 
    P2_gpPnt = np_array_to_gp_Pnt(P2) 
    P3_gpPnt = np_array_to_gp_Pnt(P3) 
    P4_gpPnt = P2_gpPnt.Rotated(gp_Ax1(P3_gpPnt,gp_Dir(1,0,0)),tab_radius_angle)
    
    
    E1 = BRepBuilderAPI_MakeEdge(ref_Pnt,P1_gpPnt)
    E2 = BRepBuilderAPI_MakeEdge(P1_gpPnt,P2_gpPnt)
    
    AX2 = gp_Ax2(P3_gpPnt,gp_Dir(1,0,0))
    tab_circle = gp_Circ(AX2,tab_radius)
    ArcOfCircle = GC_MakeArcOfCircle(tab_circle, P4_gpPnt, P2_gpPnt, True)
    ARC1 = BRepBuilderAPI_MakeEdge(ArcOfCircle.Value())
    
    # Create a wire out of the edges
    topWire = BRepBuilderAPI_MakeWire(E1.Edge(), E2.Edge(), ARC1.Edge())
    
    # Set up the mirror
    aTrsf = gp_Trsf()
    MirrorAxis= gp_OY()
    aTrsf.SetMirror(MirrorAxis)
    
    # Apply the mirror transformation
    aBRespTrsf = BRepBuilderAPI_Transform(topWire.Wire(), aTrsf)
    # Get the mirrored shape back out of the transformation and convert back to a wire,  A wire instead of a generic shape now
    MirroredWire = topods.Wire(aBRespTrsf.Shape())    
    
    #topWire.Add(aMirroredWire)
    # Combine the two constituent wires
    tab_wire = BRepBuilderAPI_MakeWire()
    tab_wire.Add(topWire.Wire())
    tab_wire.Add(MirroredWire)
    tab_wire = tab_wire.Wire()

    return tab_wire


##map(display.DisplayShape,[ref_Pnt,P1_gpPnt,P2_gpPnt,P3_gpPnt,P4_gpPnt])
##display.DisplayShape(E1.Shape())
##display.DisplayShape(E2.Shape())
##display.DisplayShape(ARC1.Shape())
##display.DisplayShape(topWire.Shape(),color='BLUE')
##display.DisplayShape(aBRespTrsf.Shape(),color='GREEN')
#display.DisplayShape(ref_Pnt)
#display.DisplayShape(tab_wire,color='WHITE')
#
#display.FitAll()
#start_display()
