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
from utils import *
from segment import Segment
from layer import Layer


#Basic Libraries:
import numpy as np

#PythonOCC Libraries
from OCC.gp import gp_Pnt2d,  gp_Trsf2d, gp_Vec2d
from OCC.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.Geom2d import Geom2d_TrimmedCurve 
from OCC.Geom2dAPI import Geom2dAPI_InterCurveCurve
from OCC.GCPnts import GCPnts_AbscissaPoint
from OCC.Geom2dAPI import Geom2dAPI_Interpolate


#Own Libraries:
from utils import calc_DCT_angles, TColgp_HArray1OfPnt2d_from_nparray, discrete_stepsize, curvature_of_curve
from BSplineLst_utils import find_BSplineLst_coordinate, get_BSpline_length, get_BSplineLst_length, \
                            get_BSplineLst_Pnt2d, discretize_BSplineLst, BSplineLst_from_dct, copy_BSpline, \
                            findPnt_on_2dcurve, set_BSplineLst_to_Origin
from wire_utils import build_wire_from_BSplineLst
from layer import Layer
from utils import getID, Pnt2dLst_to_npArray
from offset import shp_parallel_offset


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
#tmp_Segment = SegmentLst[0]
#Layer1 = Layer(0001,SegmentLst[0].BSplineLst, tmp_Segment.Layup[0][0], tmp_Segment.Layup[0][1],tmp_Segment.Layup[0][2],tmp_Segment.Layup[0][3],tmp_Segment.Layup[0][4])
#Layer1.trim_to_coords()

#Trim Layer 1
#Discretize Layer 1
#Offset Layer 1
#Rebuild Layer 1

#=========================================================================
Segment0 = SegmentLst[0].copy()
Segment0.build_wire()


S1 = 0.2
S2 = 0.7


#OrgBSplineLst = Segment0.BSplineLst
#npArray = discretize_BSplineLst(OrgBSplineLst,0,1,200)
#offlinepts = shp_parallel_offset(npArray,0.03)
#OffsetBSplineLst = BSplineLst_from_dct(offlinepts)

Trimmed_BSplineLst = Segment0.trim(S1,S2,0,1)
Trimmed_Wire2 = build_wire_from_BSplineLst(Trimmed_BSplineLst)
npArray = discretize_BSplineLst(Trimmed_BSplineLst,S1,S2,200)
offlinepts = shp_parallel_offset(npArray,0.01)
OffsetBSplineLst = BSplineLst_from_dct(offlinepts)



#display.DisplayShape(OffsetBSplineLst[0],  color="ORANGE")
#display.DisplayShape(OffsetBSplineLst[1],  color="YELLOW")

plt.plot(*npArray.T, color='black', marker='.')
plt.plot(*offlinepts.T, color='red', marker='.')
plt.axis('equal')
plt.show()
        
OffsetWire = build_wire_from_BSplineLst(OffsetBSplineLst)

display.DisplayShape(Segment0.wire)
display.DisplayShape(Trimmed_Wire2, color="GREEN")
display.DisplayShape(OffsetWire, color="BLUE")


#Pnt = get_BSplineLst_Pnt2d(OrgBSplineLst,1,0,1)
#display.DisplayShape(Pnt, color="ORANGE")


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



display.set_bg_gradient_color(20,6,111,200,200,200)
display.View_Top()
display.FitAll()
start_display()