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
from utils import calc_DCT_angles, TColgp_HArray1OfPnt2d_from_nparray
from BSplineLst_utils import  find_BSplineLst_coordinate, get_BSpline_length, get_BSplineLst_length
from wire_utils import build_wire_from_BSplineLst
from layer import Layer

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
section = section_config()
section.read_config(filename)  

#==========================
#Build Segment 0
print "STATUS:\t Build Segment 0"

Segment0 = Segment()
Segment0.BSplineLst_from_airfoil_database(section.SETUP_Airfoil, 150)
Segment0.build_wire()

S1 = 0.4
S2 = 0.6
Trimmed_Wire1 = Segment0.trim_SEGwire(S1,S2)

S1 = 0.3
S2 = 0.9995
Trimmed_BSplineLst = Segment0.trim(S1,S2)
Trimmed_Wire2 = build_wire_from_BSplineLst(Trimmed_BSplineLst)

#Segment0.build_wire() 
display.DisplayShape(Segment0.wire)
display.DisplayShape(Trimmed_Wire2, color="GREEN")
display.DisplayShape(Trimmed_Wire1, color="ORANGE")

Layer0 = Layer(Segment0.BSplineLst,section.SEG_Layup[1][1][0],section.SEG_Layup[1][1][1],section.SEG_Layup[1][1][2],section.SEG_Layup[1][1][3],section.SEG_Layup[1][1][4])

#Discretize_BSplineLst:

def radius_of_curve(curve,u):
    p = gp_Pnt2d()
    v1 = gp_Vec2d()
    v2 = gp_Vec2d()
    curve.D2(u,p,v1,v2)
    eq1 = v1.Crossed(v2)      
    eq2 = v1.Magnitude()
    curvature = abs(eq1)/abs(eq2**3)
    radius = 1/curvature
    return radius
    
def curvature_of_curve(curve,u):
    p = gp_Pnt2d()
    v1 = gp_Vec2d()
    v2 = gp_Vec2d()
    curve.D2(u,p,v1,v2)
    eq1 = v1.Crossed(v2)      
    eq2 = v1.Magnitude()
    curvature = abs(eq1)/abs(eq2**3)
    return curvature 

def discrete_stepsize(kappa):
    beta = 0.8
    stepsize = 1-beta*math.tanh(kappa)
    return stepsize
    
    
def find_BSpline_coordinate(BSpline,s):
    # Be careful, s stands for the lenght coordinate of a single BSpline, while S represents the Global Coordinate!
    BSpline_length = get_BSpline_length(BSpline)
    tolerance=1e-10
    Adaptor = Geom2dAdaptor_Curve(BSpline.GetHandle())
    if s <= BSpline_length:
             tmp = GCPnts_AbscissaPoint(tolerance, Adaptor, s, 0)
             U = tmp.Parameter()
    return U             
             

    
Trimmed_BSplineLst_length = get_BSplineLst_length(Trimmed_BSplineLst)    
 
for i,item in enumerate(Trimmed_BSplineLst):
    BSpline_length = get_BSpline_length(item)
    
    












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