# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 13:36:27 2016

@author: TPflumm
"""

#-------------------------------
#          H E A D E R
#-------------------------------

#Basic Libraries
import math
import numpy as np       

from read_inputfile import * 
from core_geometry_utils import *
from core_operations_utils import *


from OCC.Geom import Geom_OffsetCurve, Geom_BezierCurve, Geom_Plane, Geom_TrimmedCurve, Geom_Curve 
from OCC.BRepCheck import BRepCheck_Wire
#-------------------------------
#          FUNCTIONS
#-------------------------------


#======================================================
#       MAIN
#======================================================
#START    
    
#======================================================================
#READ INPUT  
   
filename = 'sec_config.input'
section = section_config()
section.read_config(filename)  

#======================================================================
#Build Segment 0
print "STATUS:\t Build Segment 0"
#GET Boundary_DCT for Segment 0
section.SEG_Boundary_DCT.append(UIUCAirfoil2d(section.SETUP_Airfoil))   

ID = 0


#def build_segment(section,ID)

#===================
#read input

#Check if either a Discrete of a Topo Boundary is defined continoue with build process
if (0 < len(section.SEG_Boundary_DCT)) or (0 < len(section.SEG_Boundary_OCC)) :
    print "INFO:\t Boundary defined for segment 0"
    if not (0 < len(section.SEG_Boundary_OCC)):
        print 'INFO:\t The segment boundary is only defined by discrete values. Applying interpolation algorithm to obtain OCC definition'
        #Interpolation Algorithm to Obtain OCC definition from discrete boundary values. 
        #TBD: Analyse discrete curve, check for edges and periodicity to create a OCC definition!
        tmp_data = section.SEG_Boundary_DCT[0]
        tmp_harray = TColgp_HArray1OfPnt2d_from_nparray(tmp_data)                                   #Put Data into Harray1OfPnt2d
        tmp_interpolation = Geom2dAPI_Interpolate(tmp_harray.GetHandle(), False, 0.001)             #Interpolate datapoints to bspline
        tmp_interpolation.Perform()                                               
        tmp_bspline = tmp_interpolation.Curve().GetObject()
        
        #Convert to 3D SPACE onto x,y plane
        P = gp_Pnt(0,0,0)
        V = gp_Dir(gp_Vec(0,0,1))
        Plane = Geom_Plane(P, V)
        
        tmp_edge = BRepBuilderAPI_MakeEdge(tmp_bspline.GetHandle(),Plane.GetHandle())
        #display.DisplayShape(edge.Shape(), update=True, color='WHITE')
        tmp_wire = BRepBuilderAPI_MakeWire(tmp_edge.Edge())        
        
        
        section.SEG_Boundary_OCC.append(tmp_wire)   
    else: 
        print 'INFO:\t The segment boundary is defined in OCC definition' 
        
else:
    print "ERROR:\t No boundary defined for segment 0"


#===================
#read check if seg_boundary is closed!
#if not close it
test = BRepCheck_Wire(tmp_wire.Wire())
print test.Closed()     #Returns 0 for NoError
                        #Returns 27 for 



#======================================================================
if section.SETUP_NbOfWebs>0:
    print "STATUS:\t Continue  Building Segments:"
    
else:
    print "STATUS:\t Fill Segment 0 with Core"

#======================================================================
#Balance Weight
if section.SETUP_BalanceWeight == True:
    print "STATUS:\t Create Balanace Weight"
    print "STATUS:\t Trim Geometry to Balance Weight"

#======================================================================
#Create MESH
    


#======================================================================
#EXPORT
print "STATUS:\t Exporting PATRAN file Format"
print "STATUS:\t Exporting .STP and .IGS file Format"

#END    