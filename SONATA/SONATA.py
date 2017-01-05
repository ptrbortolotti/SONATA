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
from OCC.Geom2dAPI import Geom2dAPI_InterCurveCurve, Geom2dAPI_Interpolate, Geom2dAPI_PointsToBSpline
from OCC.Geom2dConvert import geom2dconvert_CurveToBSplineCurve
from OCC.GCPnts import GCPnts_AbscissaPoint, GCPnts_QuasiUniformDeflection, GCPnts_TangentialDeflection
from OCC.Geom2dAPI import Geom2dAPI_Interpolate
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
display.Context.SetDeviationAngle(0.0000002)       # 0.001 default. Be careful to scale it to the problem.
display.Context.SetDeviationCoefficient(0.0000002) # 0.001 default. Be careful to scale it to the problem. 
#show_coordinate_system(display) #CREATE AXIS SYSTEM for Visualization


#==========================
#READ INPUT:
print "STATUS:\t Read Input"
filename = 'sec_config.input'
Configuration = section_config(filename)
Configuration.SEG_Layup[0][:,2] =  Configuration.SEG_Layup[0][:,2]/Configuration.SETUP_chord
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
        
#==========================
#Build Segment 0
Segment0 = SegmentLst[0].copy()
test = Segment0.Projection


Layer1 = Layer(1,SegmentLst[0].BSplineLst, Segment0.Layup[0][0], Segment0.Layup[0][1],Segment0.Layup[0][2],Segment0.Layup[0][3],Segment0.Layup[0][4],cutoff_style= 1.24, join_style=1, name = 'Skin Layer 1')
#Layer2 = Layer(0002,Boundary2)
print Layer1.cutoff_style
#=========================================================================

Segment0 = SegmentLst[0].copy()
Segment0.build_wire()
LayerLst = []
Layup = Segment0.Layup



for i in range(1,5):
    print "i:",i
    
    #get_boundary_layer
    if i == 1:
        Boundary_BSplineLst = Segment0.BSplineLst
    
    else:
        layup_projection = Segment0.Projection[i-2]
        Boundary_BSplineLst = []
        for j,item in enumerate(layup_projection):
            Boundary_to_trim = copy_BSplineLst(LayerLst[int(item[2])-1].Boundary_BSplineLst)
 
            if item[2]==i:
                Boundary_BSplineLst += LayerLst[int(item[2])-1].BSplineLst
            else:
                Boundary_BSplineLst += trim_BSplineLst(Boundary_to_trim, item[0], item[1], 0, 1)
    
            #print int(item[2])-1
            #Wire = build_wire_from_BSplineLst(Trimmed)          
            #display.DisplayShape(Wire, color="YELLOW")
              
   
    #CREATE LAYER Object
    tmp_Layer = Layer(i,Boundary_BSplineLst, Segment0.Layup[i-1][0], Segment0.Layup[i-1][1],Segment0.Layup[i-1][2],Segment0.Layup[i-1][3],Segment0.Layup[i-1][4],cutoff_style= 2, join_style=1, name = 'test')   
    tmp_Layer.build_layer() 
    LayerLst.append(tmp_Layer)     
    OffsetBSplineLst = tmp_Layer.BSplineLst
    S1 = tmp_Layer.S1
    S2 = tmp_Layer.S2
    OffsetWire = build_wire_from_BSplineLst(OffsetBSplineLst)
    display.DisplayShape(OffsetWire, color="BLUE")
    
    
    
    
    
    
    
#check if OffsetBSplineLst is closed:
Offset_StartPnt = get_BSplineLst_Pnt2d(OffsetBSplineLst,S1,S1,S2)
Offset_EndPnt = get_BSplineLst_Pnt2d(OffsetBSplineLst,S2,S1,S2)
if not Offset_StartPnt.IsEqual(Offset_EndPnt, 1e-9):
    Trimm1 = SegmentLst[0].copy()
    Trimm2 = SegmentLst[0].copy()
    if S1>0:
        Trimm1 = Trimm1.trim(0,S1,0,1)
    else: 
        Trimm1 = []
    if S2<1:
        Trimm2 = Trimm2.trim(S2,1,0,1)
    else:
        Trimm2 = []

    newList =  Trimm1 + OffsetBSplineLst + Trimm2
    

else:
    newList = OffsetBSplineLst






#plt.figure(1)
#plt.plot(*npArray.T, color='green', marker='.')
#plt.plot(*offlinepts.T, color='blue', marker='.')
#plt.axis('equal')
#plt.show()
        
OffsetWire = build_wire_from_BSplineLst(OffsetBSplineLst)
newBoundary = build_wire_from_BSplineLst(newList)

display.DisplayShape(Segment0.wire)
#display.DisplayShape(Trimmed_Wire2, color="GREEN")
#display.DisplayShape(OffsetWire, color="BLUE")
#display.DisplayShape(newBoundary, color="CYAN", update=True)



##SECOND LAYER TEST
#copynewList = copy_BSplineLst(newList)
#
#S1 = 0.3
#S2 = 0.7
#Trimmed_BSplineLst = trim_BSplineLst(copynewList,S1,S2,0,1)
#Trimmed_Wire3 = build_wire_from_BSplineLst(Trimmed_BSplineLst)
#npArray = discretize_BSplineLst(Trimmed_BSplineLst,1e-05)   
#offlinepts = shp_parallel_offset(npArray,0.005)
#OffsetBSplineLst = BSplineLst_from_dct(offlinepts)
#OffsetBSplineLst = trim_BSplineLst(OffsetBSplineLst, S1+0.003, S2-0.003, S1, S2)
#
#
#StartPnt1 = get_BSplineLst_Pnt2d(Trimmed_BSplineLst,S1,S1,S2)
#EndPnt1 = get_BSplineLst_Pnt2d(Trimmed_BSplineLst,S2,S1,S2)
#StartPnt2 = get_BSplineLst_Pnt2d(OffsetBSplineLst,S1,S1,S2)
#EndPnt2 = get_BSplineLst_Pnt2d(OffsetBSplineLst,S2,S1,S2)
#
## the first 
#array = TColgp_Array1OfPnt2d(1, 2)
#array.SetValue(1, StartPnt1)
#array.SetValue(2, StartPnt2)
#tmp_bspline1 = Geom2dAPI_PointsToBSpline(array).Curve().GetObject()
#
## the second 
#array = TColgp_Array1OfPnt2d(1, 2)
#array.SetValue(1, EndPnt2)
#array.SetValue(2, EndPnt1)
#tmp_bspline2 = Geom2dAPI_PointsToBSpline(array).Curve().GetObject()
#
#OffsetBSplineLst.insert(0,tmp_bspline1)
#OffsetBSplineLst.append(tmp_bspline2)
#
#OffsetWire2 = build_wire_from_BSplineLst(OffsetBSplineLst)
#
##plt.figure(2)
##plt.plot(*npArray.T, color='green', marker='.')
##plt.plot(*offlinepts.T, color='blue', marker='.')
##plt.axis('equal')
##plt.show()
#        
#
#
#
#
#
##display.DisplayShape(OffsetBSplineLst[0],  color="ORANGE")
##display.DisplayShape(OffsetBSplineLst[1],  color="YELLOW")
#
#
##display.DisplayShape(Trimmed_Wire3, color="GREEN")
#display.DisplayShape(OffsetWire2, color="GREEN")

#CHECK OFFSET BSPlineLst intersects original bounding BSplineLst:
#OrgBSplineLst = Segment0.BSplineLst
#for i,OffsetItem in enumerate(OffsetBSplineLst):
#    for j,OrgItem in enumerate(OrgBSplineLst):
#        Inter = Geom2dAPI_InterCurveCurve(OffsetItem.GetHandle(), OrgItem.GetHandle(), 1.0e-10)
#        if Inter.NbPoints() == 1:
#            IntPnt = Inter.Point(1)
#            display.DisplayShape(IntPnt)      
#            u = findPnt_on_2dcurve(IntPnt,OffsetItem)
#            BSplineCopy = copy_BSpline(OffsetItem)
#            last = BSplineCopy.LastParameter()
#            BSplineCopy.Segment(u,last)
#            OffsetBSplineLst[j] = BSplineCopy
#            print u





#pnt2d = get_BSplineLst_Pnt2d(Segment0.BSplineLst,0.399,0,1)
#display.DisplayShape(pnt2d, color = 'CYAN')

#pnt2d = get_BSplineLst_Pnt2d(Trimmed_BSplineLst,0.995,S1,S2)
#display.DisplayShape(pnt2d, color = 'BLUE')







#Pnt = Segment0.trim(0.4,0.6)    

#print get_BSplineLst_length(tmp_BSplineLst)
#print find_BSplineLst_coordinate(tmp_BSplineLst,1)
#Pnt =  get_BSplineLst_Pnt(tmp_BSplineLst,0.9995)


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