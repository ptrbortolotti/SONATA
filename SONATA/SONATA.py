# -*- coding: utf-8 -*-
"""
THIS IS THE SONATA EXECUTION FILE!
@author: TPflumm
"""

#-------------------------------
#          H E A D E R
#-------------------------------

#Basic Libraries
import math
import numpy as np       
import matplotlib.pyplot as plt

from OCC.Display.SimpleGui import init_display
from OCC.gp import gp_Pnt2d
from OCC.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.GCPnts import GCPnts_AbscissaPoint

from readinput import section_config 
from display import show_coordinate_system
from utils import *
from segment import Segment
from layer import Layer



###############################################################################
#                           M    A    I    N                                  #
###############################################################################


#======================================================================
#DISPLAY CONFIG:
display, start_display, add_menu, add_function_to_menu = init_display('wx')
display.Context.SetDeviationAngle(0.000001)       # 0.001 default. Be careful to scale it to the problem.
display.Context.SetDeviationCoefficient(0.000001) # 0.001 default. Be careful to scale it to the problem. 
show_coordinate_system(display) #CREATE AXIS SYSTEM for Visualization

#======================================================================
#READ INPUT:
print "STATUS:\t Read Input"
filename = 'sec_config.input'
section = section_config()
section.read_config(filename)  

#======================================================================
#Build Segment 0
print "STATUS:\t Build Segment 0"

Segment0 = Segment()
Segment0.BSplineLst_from_airfoil_database(section.SETUP_Airfoil)
#Segment0.seg_boundary_from_file('naca671215')
Segment0.build_wire()   



tmp_BSplineLst = Segment0.BSplineLst
    
print get_BSplineLst_length(tmp_BSplineLst)
print find_BSplineLst_coordinate(tmp_BSplineLst,1)
Pnt =  get_BSplineLst_Pnt(tmp_BSplineLst,0.9995)


#======================================================================
#PLOT


display.DisplayShape(Segment0.wire)
display.DisplayShape(Pnt)

display.set_bg_gradient_color(20,6,111,200,200,200)
display.View_Top()
display.FitAll()
start_display()