# -*- coding: utf-8 -*-
"""
THIS IS THE SONATA EXECUTION FILE!
@author: TPflumm
"""

#Basic Libraries:
import math
import numpy as np       
import matplotlib.pyplot as plt
import shapely.geometry as shp
from scipy.optimize import leastsq

#Third Party Libaries: OCC Libraries
from OCC.Display.SimpleGui import init_display

#Own Modules:
from readinput import section_config 
from display import show_coordinate_system
from segment import Segment
from layer import Layer


#Basic Libraries:
import numpy as np

#PythonOCC Libraries
from OCC.gp import gp_Pnt2d,  gp_Trsf2d, gp_Vec2d
from OCC.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.TColgp import TColgp_HArray1OfPnt2d, TColgp_Array1OfPnt2d
from OCC.Geom2d import Geom2d_TrimmedCurve, Geom2d_BezierCurve
from OCC.Geom2d import Handle_Geom2d_BSplineCurve_DownCast
from OCC.Geom2dAPI import Geom2dAPI_InterCurveCurve, Geom2dAPI_Interpolate, Geom2dAPI_PointsToBSpline
from OCC.Geom2dConvert import geom2dconvert_CurveToBSplineCurve
from OCC.GCPnts import GCPnts_AbscissaPoint, GCPnts_QuasiUniformDeflection, GCPnts_TangentialDeflection
from OCC.Geom2dAPI import Geom2dAPI_Interpolate
from OCC.Quantity import Quantity_Color
from OCC.Graphic3d import (Graphic3d_EF_PDF,
                           Graphic3d_EF_SVG,
                           Graphic3d_EF_TEX,
                           Graphic3d_EF_PostScript,
                           Graphic3d_EF_EnhPostScript)

#Own Libraries:
from utils import calc_DCT_angles, TColgp_HArray1OfPnt2d_from_nparray, discrete_stepsize, curvature_of_curve
from BSplineLst_utils import find_BSplineLst_coordinate, get_BSpline_length, get_BSplineLst_length, \
                            get_BSplineLst_Pnt2d, discretize_BSplineLst, BSplineLst_from_dct, copy_BSpline, \
                            findPnt_on_2dcurve, set_BSplineLst_to_Origin, seg_boundary_from_dct, copy_BSplineLst,\
                            trim_BSplineLst, get_BSplineLst_D2
from wire_utils import build_wire_from_BSplineLst
from layer import Layer
from utils import getID, Pnt2dLst_to_npArray,unique_rows, P2Pdistance, point2d_list_to_TColgp_HArray1OfPnt2d, point2d_list_to_TColgp_Array1OfPnt2d
from offset import shp_parallel_offset
from readinput import UIUCAirfoil2d, AirfoilDat2d
from cutoff import cutoff_layer




def export_to_PDF(event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    f.Export('SONATA_export.pdf', Graphic3d_EF_PDF)
    display.set_bg_gradient_color(20,6,111,200,200,200)
    
def export_to_SVG(event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    f.Export('SONATA_export.svg', Graphic3d_EF_SVG)
    display.set_bg_gradient_color(20,6,111,200,200,200)

def export_to_PS(event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    f.Export('SONATA_export.ps', Graphic3d_EF_PostScript)
    display.set_bg_gradient_color(20,6,111,200,200,200)

def export_to_EnhPS(event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    f.Export('SONATA_export_enh.ps', Graphic3d_EF_EnhPostScript)
    display.set_bg_gradient_color(20,6,111,200,200,200)

def export_to_TEX(event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    f.Export('SONATA_export.tex', Graphic3d_EF_TEX)
    display.set_bg_gradient_color(20,6,111,200,200,200)

def print_xy_click(shp, *kwargs):
    for shape in shp:
        print("Shape selected: ", shape)
    print(kwargs)


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
print "STATUS:\t Reading Input"
filename = 'sec_config.input'
Configuration = section_config(filename)

#Normalize layer thickness to chord lenght
for i, item in enumerate(Configuration.SEG_Layup):
    item[:,2] =  item[:,2]/Configuration.SETUP_chord

#==========================                  
#Initialize Segments and sort the according to ID 
SegmentLst = []   #List of Segment Objects
for i,item in enumerate(Configuration.SEG_ID):
    if item == 0:
        #SegmentLst.append(Segment(item, Layup = Configuration.SEG_Layup[i], CoreMaterial = Configuration.SEG_CoreMaterial[i], OCC=False, airfoil = 'naca23012'))
        SegmentLst.append(Segment(item, Layup = Configuration.SEG_Layup[i], CoreMaterial = Configuration.SEG_CoreMaterial[i], OCC=False, filename = 'AREA_R250.dat'))
    else:
        SegmentLst.append(Segment(item, Layup = Configuration.SEG_Layup[i], CoreMaterial = Configuration.SEG_CoreMaterial[i]))
sorted(SegmentLst, key=getID)
        


# ============================================================================= 
#               Build SEGMENT 0:
# =============================================================================
SegmentLst[0].build_wire()
SegmentLst[0].build_layers()


        
# ============================================================================= 
#               Build Webs:
# =============================================================================

WEB_Pos1_Pnt2d = get_BSplineLst_Pnt2d(SegmentLst[0].BSplineLst,Configuration.WEB_Pos1[0],0,1)
WEB_Pos2_Pnt2d = get_BSplineLst_Pnt2d(SegmentLst[0].BSplineLst,Configuration.WEB_Pos2[0],0,1)
WEB_BSpline = Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([WEB_Pos1_Pnt2d,WEB_Pos2_Pnt2d])).Curve().GetObject()

#Determine Boundary from Segment 0:
new_Boundary_BSplineLst = []
new_Boundary_BSplineLst += trim_BSplineLst(SegmentLst[0].LayerLst[-1].Boundary_BSplineLst, 0, SegmentLst[0].LayerLst[-1].S1, 0, 1)  #start und ende der lage
new_Boundary_BSplineLst += copy_BSplineLst(SegmentLst[0].LayerLst[-1].BSplineLst)
new_Boundary_BSplineLst += trim_BSplineLst(SegmentLst[0].LayerLst[-1].Boundary_BSplineLst, SegmentLst[0].LayerLst[-1].S2, 1, 0, 1)  #start und ende der lage
new_Boundary_BSplineLst = set_BSplineLst_to_Origin(new_Boundary_BSplineLst)


tolerance=1e-10
IntPnts = []
IntPnts_Pnt2d = []
for i,item in enumerate(new_Boundary_BSplineLst): 
    intersection = Geom2dAPI_InterCurveCurve(WEB_BSpline.GetHandle(), item.GetHandle(),tolerance)
    for j in range(1,intersection.NbPoints()+1):
            IntPnt = intersection.Point(j)
            IntPnts_Pnt2d.append(IntPnt)
            u = findPnt_on_2dcurve(IntPnt,item.GetHandle())
            IntPnts.append([i,u])
    
            
P1 = IntPnts_Pnt2d[0]
P2 = IntPnts_Pnt2d[1]
Boundary_WEB_BSplineLst = [Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([P1,P2])).Curve().GetObject()]
BSplineLst = copy_BSplineLst(new_Boundary_BSplineLst)

# =============================================================================
#CREATE SEGMENT 1 BOUNDARY
trimmed_Boundary1= []
para1 =  IntPnts[0]
para2 =  IntPnts[1]
for i,item in enumerate(BSplineLst):       
     First = item.FirstParameter() 
     Last =  item.LastParameter() 
     if para1[0] == i and para2[0] != i:
         BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
         BSplineCopy.Segment(para1[1],Last)
         trimmed_Boundary1.append(BSplineCopy)
         
     elif (para1[0] != i and para2[0] != i) and (para1[0] < i and para2[0] > i):
         BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
         BSplineCopy.Segment(First,Last)
         trimmed_Boundary1.append(BSplineCopy)
         
     elif para1[0] == i and para2[0] == i:
         BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
         BSplineCopy.Segment(para1[1],para2[1])
         trimmed_Boundary1.append(BSplineCopy)
         break

     elif para1[0] != i and para2[0] == i:
         BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
         BSplineCopy.Segment(First,para2[1])
         trimmed_Boundary1.append(BSplineCopy)
         break

Boundary1_WEB_BSplineLst = [Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([P2,P1])).Curve().GetObject()]
Boundary1_BSplineLst = trimmed_Boundary1 + Boundary1_WEB_BSplineLst                                  
Boundary1_BSplineLst = set_BSplineLst_to_Origin(Boundary1_BSplineLst)                          
# =============================================================================
#CREATE SEGMENT 2 BOUNDARY
trimmed_Boundary2 = []
front_BSplineLst = []
rear_BSplineLst = []
para1 =  IntPnts[1]
para2 =  IntPnts[0]

for i,item in enumerate(BSplineLst):
    First = item.FirstParameter() 
    Last =  item.LastParameter()
    if para1[0] == i and para2[0] != i: #Okay
         BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
         BSplineCopy.Segment(para1[1],Last)
         rear_BSplineLst.append(BSplineCopy)
         
    elif (para1[0] != i and para2[0] != i) and (para1[0] > i and para2[0] > i):
         BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
         BSplineCopy.Segment(First,Last)
         front_BSplineLst.append(BSplineCopy)
         
    elif (para1[0] != i and para2[0] != i) and (para1[0] < i and para2[0] < i):
         BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
         BSplineCopy.Segment(First,Last)
         rear_BSplineLst.append(BSplineCopy)
         
    elif para1[0] == i and para2[0] == i: #Okay
         BSplineCopy1 = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
         BSplineCopy1.Segment(para1[1],Last)
         trimmed_Boundary2.append(BSplineCopy1)
         BSplineCopy2 = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
         BSplineCopy2.Segment(First,para2[1])
         trimmed_Boundary2.append(BSplineCopy2)
         break

    elif para1[0] != i and para2[0] == i: #Okay
         BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
         BSplineCopy.Segment(First,para2[1])
         front_BSplineLst.append(BSplineCopy)
         #break

if front_BSplineLst and rear_BSplineLst:
    trimmed_Boundary2 =  rear_BSplineLst + front_BSplineLst

elif front_BSplineLst and not rear_BSplineLst:
    trimmed_Boundary2 = front_BSplineLst
    
elif not front_BSplineLst and rear_BSplineLst:
    trimmed_Boundary2 = rear_BSplineLst  
 
  
Boundary2_WEB_BSplineLst = [Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([P1,P2])).Curve().GetObject()] 
Boundary2_BSplineLst = trimmed_Boundary2 + Boundary2_WEB_BSplineLst                 
Boundary2_BSplineLst = set_BSplineLst_to_Origin(Boundary2_BSplineLst)


SegmentLst[1].BSplineLst = Boundary1_BSplineLst
SegmentLst[2].BSplineLst = Boundary2_BSplineLst


SegmentLst[1].build_wire()
SegmentLst[2].build_wire()

# ============================================================================= 
#               Build SEGMENT 1 & 2:
# =============================================================================
SegmentLst[1].build_layers()
SegmentLst[2].build_layers()





# ============================================================================= 
#               DISPLAY :
# =============================================================================


#DISPLAY LAYERS


for i,seg in enumerate(SegmentLst):
    display.DisplayShape(seg.wire, color="BLACK")
    k = 0
    for j,layer in enumerate(seg.LayerLst):
        [R,G,B,T] =  plt.cm.jet(k*50)
        display.DisplayColoredShape(layer.wire, Quantity_Color(R, G, B, 0))
        Boundary_BSplineLst = layer.Boundary_BSplineLst
        Wire = build_wire_from_BSplineLst(Boundary_BSplineLst)    
        #display.DisplayColoredShape(Wire, Quantity_Color(R, G, B, 0))
        k = k+1;
        if k>5:
            k = 0
        #item.get_pnt2d(0,)







    


#======================================================================
#PLOT
f = display.View.View().GetObject()

display.register_select_callback(print_xy_click)
display.set_bg_gradient_color(20,6,111,200,200,200)
add_menu('screencapture')
add_function_to_menu('screencapture', export_to_PDF)
add_function_to_menu('screencapture', export_to_SVG)
add_function_to_menu('screencapture', export_to_PS)
add_function_to_menu('screencapture', export_to_EnhPS)
add_function_to_menu('screencapture', export_to_TEX)
display.View_Top()
display.FitAll()

  
start_display()