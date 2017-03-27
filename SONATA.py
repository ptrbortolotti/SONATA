# -*- coding: utf-8 -*-
"""
THIS IS THE SONATA EXECUTION FILE!
@author: TPflumm
"""

#Basic PYTHON Modules:
import numpy as np       
import matplotlib.pyplot as plt
import pickle
from functools import partial

#PythonOCC Libraries
from OCC.Display.SimpleGui import init_display
from OCC.AIS import AIS_Shape
from OCC.gp import gp_Pnt2d,  gp_Trsf2d, gp_Vec2d, gp_Pnt
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.TColgp import TColgp_HArray1OfPnt2d, TColgp_Array1OfPnt2d
from OCC.Geom2d import Geom2d_TrimmedCurve, Geom2d_BezierCurve
from OCC.Geom2d import Handle_Geom2d_BSplineCurve_DownCast
from OCC.Geom2dAPI import Geom2dAPI_InterCurveCurve, Geom2dAPI_Interpolate, Geom2dAPI_PointsToBSpline
from OCC.Geom2dConvert import geom2dconvert_CurveToBSplineCurve
from OCC.GCPnts import GCPnts_AbscissaPoint, GCPnts_QuasiUniformDeflection, GCPnts_TangentialDeflection, GCPnts_QuasiUniformAbscissa
from OCC.Quantity import Quantity_Color
from OCC.Graphic3d import (Graphic3d_EF_PDF,
                           Graphic3d_EF_SVG,
                           Graphic3d_EF_TEX,
                           Graphic3d_EF_PostScript,
                           Graphic3d_EF_EnhPostScript)
from OCC.STEPControl import STEPControl_Writer, STEPControl_AsIs, STEPControl_GeometricCurveSet
from OCC.Interface import Interface_Static_SetCVal
from OCC.IFSelect import IFSelect_RetDone


#SONATA modules:
from SONATA.fileIO.readinput import section_config
from SONATA.fileIO.CADinput import import_2d_stp, import_3d_stp
 
from SONATA.topo.segment import Segment
from SONATA.topo.layer import Layer
from SONATA.topo.wire_utils import build_wire_from_BSplineLst
from SONATA.topo.web import Web
from SONATA.topo.weight import Weight
from SONATA.topo.utils import calc_DCT_angles, TColgp_HArray1OfPnt2d_from_nparray, discrete_stepsize, curvature_of_curve, getID, Pnt2dLst_to_npArray,unique_rows, P2Pdistance, point2d_list_to_TColgp_HArray1OfPnt2d, point2d_list_to_TColgp_Array1OfPnt2d
from SONATA.topo.BSplineLst_utils import find_BSplineLst_coordinate, get_BSpline_length, get_BSplineLst_length, \
                            get_BSplineLst_Pnt2d, discretize_BSplineLst, BSplineLst_from_dct, copy_BSpline, \
                            findPnt_on_2dcurve, set_BSplineLst_to_Origin, seg_boundary_from_dct, copy_BSplineLst,\
                            trim_BSplineLst, get_BSplineLst_D2, trim_BSplineLst_by_Pnt2d, intersect_BSplineLst_with_BSpline
                            

from SONATA.display.display_utils import export_to_JPEG, export_to_PNG, export_to_PDF, export_to_SVG, export_to_PS, export_to_EnhPS, export_to_TEX, \
                          export_to_BMP,export_to_TIFF, show_coordinate_system, display_SONATA_SegmentLst
###############################################################################
#                           M    A    I    N                                  #
###############################################################################

#==========================
#DISPLAY CONFIG:
display, start_display, add_menu, add_function_to_menu = init_display('wx')
display.Context.SetDeviationAngle(0.000001)       # 0.001 default. Be careful to scale it to the problem.
display.Context.SetDeviationCoefficient(0.000001) # 0.001 default. Be careful to scale it to the problem. 
#show_coordinate_system(display) #CREATE AXIS SYSTEM for Visualization


#==========================
#READ INPUT:
print "STATUS:\t Reading Crossection Configuration File"
filename = 'sec_config.input'
Configuration = section_config(filename)

#Normalize layer thickness to chord lenght
#for i, item in enumerate(Configuration.SEG_Layup):
#    item[:,2] =  item[:,2]/Configuration.SETUP_scale_factor

#==========================                  
#Initialize Segments and sort the according to ID 
SegmentLst = []   #List of Segment Objects
for i,item in enumerate(Configuration.SEG_ID):
    if item == 0:        
        if Configuration.SETUP_input_type == 0:   #0) Airfoil from UIUC Database  --- naca23012
            SegmentLst.append(Segment(item, Layup = Configuration.SEG_Layup[i], CoreMaterial = Configuration.SEG_CoreMaterial[i], scale_factor = Configuration.SETUP_scale_factor, Theta = Configuration.SETUP_Theta, OCC=False, airfoil = Configuration.SETUP_datasource))
        
        elif Configuration.SETUP_input_type == 1: #1) Geometry from .dat file --- AREA_R250.dat
            SegmentLst.append(Segment(item, Layup = Configuration.SEG_Layup[i], CoreMaterial = Configuration.SEG_CoreMaterial[i], scale_factor = Configuration.SETUP_scale_factor,  Theta = Configuration.SETUP_Theta, OCC=False, filename = Configuration.SETUP_datasource))
        
        elif Configuration.SETUP_input_type == 2: #2)2d .step or .iges  --- AREA_R230.stp
            BSplineLst = import_2d_stp(Configuration.SETUP_datasource, Configuration.SETUP_scale_factor,Configuration.SETUP_Theta)
            SegmentLst.append(Segment(item, Layup = Configuration.SEG_Layup[i], CoreMaterial = Configuration.SEG_CoreMaterial[i], Theta = Configuration.SETUP_Theta, OCC=True, Boundary = BSplineLst))
        
        elif Configuration.SETUP_input_type == 3: #3)3D .step or .iges and radial station of crosssection --- AREA_Blade.stp, R=250
            BSplineLst = import_3d_stp(Configuration.SETUP_datasource,Configuration.SETUP_radial_station,Configuration.SETUP_scale_factor,Configuration.SETUP_Theta)
            SegmentLst.append(Segment(item, Layup = Configuration.SEG_Layup[i], CoreMaterial = Configuration.SEG_CoreMaterial[i], Theta = Configuration.SETUP_Theta, OCC=True, Boundary = BSplineLst))  

        else:
            print 'ERROR: \t WRONG input_type'
 
    else:
        SegmentLst.append(Segment(item, Layup = Configuration.SEG_Layup[i], CoreMaterial = Configuration.SEG_CoreMaterial[i],Theta = Configuration.SETUP_Theta))
sorted(SegmentLst, key=getID)  

# ============================================================================= 
#               Build SEGMENT 0:
# =============================================================================
SegmentLst[0].build_wire()
SegmentLst[0].build_layers()
SegmentLst[0].determine_final_boundary()    #Determine Boundary from Segment 0:
    
#============================================================================= 
#               Build Webs:
# =============================================================================
#TODO: CHECK IF WEB DEFINITION INTERSECT EACH OTHER
#TODO: SORT WEBS BY POS1 VALUES:

#Create WEB Object
WebLst = []
if Configuration.SETUP_NbOfWebs > 0:
    for i in range(0,Configuration.SETUP_NbOfWebs):
        print 'STATUS: \t Building Web %s' %(i+1)
        WebLst.append(Web(Configuration.WEB_ID[i],Configuration.WEB_Pos1[i],Configuration.WEB_Pos2[i],SegmentLst[0].BSplineLst, SegmentLst[0].final_Boundary_BSplineLst))
    sorted(SegmentLst, key=getID)  
    
# ============================================================================= 
#               Build remaining SEGMENTS 
# =============================================================================
if Configuration.SETUP_NbOfWebs > 0:
    for i,seg in enumerate(SegmentLst[1:],start=1):
        seg.build_segment_boundary_from_WebLst(WebLst,SegmentLst[0].final_Boundary_BSplineLst)
        seg.build_layers()


#============================================================================= 
#               Balance Weight
# =============================================================================
if Configuration.SETUP_BalanceWeight == True:
    print 'STATUS: \t Building Balance Weight'   
    BW = Weight(0,Configuration.BW_XPos,Configuration.BW_YPos,Configuration.BW_Diameter,Configuration.BW_MatID)

    
# ============================================================================= 
#               DISPLAY and EXPORT to STEP_AP203 ::
# =============================================================================

def display_custome_shape(shape,linewidth,transparency,RGB):
    s = shape
    ais_shp = AIS_Shape(s)
    ais_shp.SetWidth(linewidth)
    ais_shp.SetTransparency(transparency)
    ais_shp.SetColor(Quantity_Color(RGB[0], RGB[1], RGB[2], 0))
    ais_context = display.GetContext().GetObject()
    ais_context.Display(ais_shp.GetHandle())
    return None

# initialize the STEP exporter
step_writer = STEPControl_Writer()
Interface_Static_SetCVal("write.step.schema", "AP203")


# transfer shapes and display them in the viewer
#display.DisplayShape(SegmentLst[0].wire, color="BLACK")
display_custome_shape(SegmentLst[0].wire,2,0,[0,0,0])



display.DisplayShape(SegmentLst[0].BSplineLst[0].StartPoint())

#step_writer.Transfer(SegmentLst[0].wire, STEPControl_AsIs)
for i,seg in enumerate(SegmentLst):
    display_custome_shape(seg.wire,2,0,[0,0,0])
    k = 0
    for j,layer in enumerate(seg.LayerLst):
        [R,G,B,T] =  plt.cm.jet(k*50)
        
        if i==0:
            #None
            display.DisplayColoredShape(layer.wire, Quantity_Color(R, G, B, 0),update=True)
            #display.DisplayShape(layer.wire, color="BLACK")
        elif i==1:
            #None
            display.DisplayColoredShape(layer.wire, Quantity_Color(R, G, B, 0),update=True)
            #display.DisplayShape(layer.wire, color="BLACK")
        else:
            None
            display.DisplayColoredShape(layer.wire, Quantity_Color(R, G, B, 0),update=True)
        #step_writer.Transfer(layer.wire, STEPControl_AsIs)
        #Boundary_BSplineLst = layer.Boundary_BSplineLst
        #Wire = build_wire_from_BSplineLst(Boundary_BSplineLst)    
        #display.DisplayColoredShape(Wire, Quantity_Color(R, G, B, 0))
        k = k+1;
        if k>5:
            k = 0
        #item.get_pnt2d(0,)

#display.DisplayShape(BW.Curve, color="BLACK")        
#step_writer.Transfer(BW.Curve, STEPControl_AsIs)        
        
#status = step_writer.Write("SONATA.stp")    
#assert(status == IFSelect_RetDone)

# ============================================================================= 
#               SAVE THE The SegmentLst as pickle:
# =============================================================================
output_filename = filename.replace('.input', '.pkl')

with open(output_filename, 'wb') as output:
    pickle.dump(SegmentLst, output, protocol=pickle.HIGHEST_PROTOCOL)

#======================================================================
display.set_bg_gradient_color(20,6,111,200,200,200)
show_coordinate_system(display,5)

add_menu('screencapture')
add_function_to_menu('screencapture','export to PDF', partial(export_to_PDF,display))
add_function_to_menu('screencapture','export to SVG', partial(export_to_SVG,display))
add_function_to_menu('screencapture','export to PS', partial(export_to_PS,display))
add_function_to_menu('screencapture','export to EnhPS', partial(export_to_EnhPS,display))
add_function_to_menu('screencapture','export to TEX', partial(export_to_TEX,display))
add_function_to_menu('screencapture','export to BMP', partial(export_to_BMP,display))
add_function_to_menu('screencapture', 'export to PNG', partial(export_to_PNG,display))
add_function_to_menu('screencapture', 'export to JPEG', partial(export_to_JPEG,display))
add_function_to_menu('screencapture', 'export to TIFF', partial(export_to_TIFF,display))

display.View_Top()
display.FitAll()
start_display()
