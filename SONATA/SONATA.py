# -*- coding: utf-8 -*-
"""
THIS IS THE SONATA EXECUTION FILE!
@author: TPflumm
"""

#Basic Libraries:
import math
import numpy as np       
import matplotlib.pyplot as plt

#Third Party Libaries: OCC Libraries
from OCC.Display.SimpleGui import init_display

#Own Modules:
from readinput import section_config 
from display import show_coordinate_system
from utils import *
from segment import Segment
from layer import Layer




#Basic Libraries:
import numpy as np
import copy       

#PythonOCC Libraries
from OCC.gp import gp_Pnt2d,  gp_Trsf2d, gp_Vec2d
from OCC.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.Geom2d import Geom2d_TrimmedCurve
from OCC.GCPnts import GCPnts_AbscissaPoint
from OCC.Geom2dAPI import Geom2dAPI_Interpolate


#Own Libraries:
from utils import calc_DCT_angles, TColgp_HArray1OfPnt2d_from_nparray, discrete_stepsize, curvature_of_curve
from BSplineLst_utils import  find_BSplineLst_coordinate, get_BSpline_length, get_BSplineLst_length, get_BSplineLst_Pnt2d
from wire_utils import build_wire_from_BSplineLst
from layer import Layer
from utils import getID


###############################################################################
#                           M    A    I    N                                  #
###############################################################################


#==========================
#DISPLAY CONFIG:
display, start_display, add_menu, add_function_to_menu = init_display('wx')
display.Context.SetDeviationAngle(0.000001)       # 0.001 default. Be careful to scale it to the problem.
display.Context.SetDeviationCoefficient(0.000001) # 0.001 default. Be careful to scale it to the problem. 
show_coordinate_system(display) #CREATE AXIS SYSTEM for Visualization


#==========================
#READ INPUT:
print "STATUS:\t Read Input"
filename = 'sec_config.input'
Configuration = section_config(filename)

#==========================                  
#Initialize Segments and sort the according to ID 
SegmentLst = []   #List of Segment Objects
for i,item in enumerate(Configuration.SEG_ID):
    if item == 0:
        SegmentLst.append(Segment(item, Layup = Configuration.SEG_Layup[i], CoreMaterial = Configuration.SEG_CoreMaterial[i], OCC=False, airfoil = 'naca23012'))
    else:
        SegmentLst.append(Segment(item, Layup = Configuration.SEG_Layup[i], CoreMaterial = Configuration.SEG_CoreMaterial[i]))

sorted(SegmentLst, key=getID)
        
#==========================
#Build Segment 0
tmp_Segment = SegmentLst[0]
Layer1 = Layer(0001,SegmentLst[0].BSplineLst, tmp_Segment.Layup[0][0], tmp_Segment.Layup[0][1],tmp_Segment.Layup[0][2],tmp_Segment.Layup[0][3],tmp_Segment.Layup[0][4])
Layer1.trim_to_coords()


#Trim Layer 1
#Discretize Layer 1
#Offset Layer 1
#Rebuild Layer 1






#=========================================================================
S1 = 0.3
S2 = 0.7

Segment0 = SegmentLst[0]
Segment0.build_wire()

Trimmed_BSplineLst = Segment0.trim(S1,S2)
Trimmed_Wire2 = build_wire_from_BSplineLst(Trimmed_BSplineLst)

display.DisplayShape(Segment0.wire)
display.DisplayShape(Trimmed_Wire2, color="GREEN")



#Discretize_Layer:

    
Trimmed_BSplineLst_length = get_BSplineLst_length(Trimmed_BSplineLst)    
for i,item in enumerate(Trimmed_BSplineLst):
    BSpline_length = get_BSpline_length(item)

    
BSplineLst = Segment0.BSplineLst 



#Discretize Trimmed_BSplineLst
BSplineLst = Trimmed_BSplineLst
tmp_spline = BSplineLst[0]
first = tmp_spline.FirstParameter()
last = tmp_spline.LastParameter()

#the_get_BSplineLst_Point function must get the start and end arguments of the Spline!
pnt2d = get_BSplineLst_Pnt2d(BSplineLst,0.99)
display.DisplayShape(pnt2d, color = 'BLUE')



S1 = 0.3
S2 = 0.7 

S = 0
accuracy = 100
idx_old = 0
while S <= 1:
    pnt2d = get_BSplineLst_Pnt2d(BSplineLst,S)
    #grab vertices!    
    [idx,U] = find_BSplineLst_coordinate(BSplineLst,S)
    if idx > idx_old:
        BSpline = BSplineLst[idx_old]
        pnt2d = BSpline.EndPoint()
        display.DisplayShape(pnt2d, color = 'RED')
               
    else:  
        BSpline = BSplineLst[idx]
        BSplineLst_length = get_BSplineLst_length(BSplineLst)
        #print S, curvature_of_curve(BSpline,U)
        #display.DisplayShape(pnt2d)
        S += BSplineLst_length/accuracy * discrete_stepsize(curvature_of_curve(BSpline,U))

    idx_old = idx
#grab last vertice    
pnt2d = get_BSplineLst_Pnt2d(BSplineLst,1)    
display.DisplayShape(pnt2d, color = 'RED')










#Pnt = Segment0.trim(0.4,0.6)    

#print get_BSplineLst_length(tmp_BSplineLst)
#print find_BSplineLst_coordinate(tmp_BSplineLst,1)
#Pnt =  get_BSplineLst_Pnt(tmp_BSplineLst,0.9995)


#======================================================================
#PLOT



display.set_bg_gradient_color(20,6,111,200,200,200)
display.View_Top()
display.FitAll()
start_display()