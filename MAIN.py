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
import matplotlib.pyplot as plt

from read_inputfile import * 
from core_geometry_utils import *
from core_operations_utils import *

from OCC.gp import gp_Pnt, gp_Vec,  gp_Pln, gp_Dir, gp_Trsf, gp_Ax1, gp_OX, gp_Ax3, gp_Ax2, gp_Circ, gp_OY
from OCC.gp import gp_Pnt2d, gp_Vec2d, gp_XY, gp_Lin2d, gp_Dir2d
from OCC.TColgp import TColgp_Array1OfPnt, TColgp_Array1OfPnt2d, TColgp_HArray1OfPnt2d
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge,BRepBuilderAPI_MakeWire 
from OCC.Geom import Geom_Plane
from OCC.Geom2dAPI import Geom2dAPI_PointsToBSpline, Geom2dAPI_Interpolate, Geom2dAPI_InterCurveCurve 
from OCC.BRepCheck import BRepCheck_Wire
from OCC.Display.SimpleGui import init_display
from OCC.BRepGProp import brepgprop_LinearProperties
from OCC.GProp import GProp_GProps
#-------------------------------
#          FUNCTIONS
#-------------------------------


#======================================================
#       MAIN
#======================================================
#START    
###############################################################################
display, start_display, add_menu, add_function_to_menu = init_display()
display.Context.SetDeviationAngle(0.000001)      # 0.001 default
display.Context.SetDeviationCoefficient(0.000001) # 0.001 default




#======================================================================
#READ INPUT  
   
filename = 'sec_config.input'
section = section_config()
section.read_config(filename)  

#======================================================================
#Build Segment 0
print "STATUS:\t Build Segment 0"
#GET Boundary_DCT for Segment 0
section.SEG_Boundary_DCT.append(UIUCAirfoil2d(section.SETUP_Airfoil).T)   
ID = 0

#def build_segment(section,ID)
#===================
#read input 
DCT_data = section.SEG_Boundary_DCT[0]

#Neccesary functions:
def calc_angle_between(v1, v2):
    return np.degrees(math.atan2(np.linalg.norm(np.cross(v1,v2)), np.dot(v1,v2)))

def calc_DCT_angles(DCT_data):
    temp = []
    for i in range(0,DCT_data.shape[0]):
        if i == 0: #first point
            v1 = DCT_data[i-2]-DCT_data[i] 
            v2 = DCT_data[i+1]-DCT_data[i]
        elif i == DCT_data.shape[0]-1: #last point
            v1 = DCT_data[i-1]-DCT_data[i] 
            v2 = DCT_data[1]-DCT_data[i]
        else:
            v1 = DCT_data[i-1]-DCT_data[i]
            v2 = DCT_data[i+1]-DCT_data[i]
        temp.append(calc_angle_between(v1,v2))
    return np.array(temp)


def display_points_of_array(array):
    for j in range(array.Lower(),array.Upper()+1):
        p = array.Value(j)
        display.DisplayShape(p, update=False)            


#Check if either a Discrete of a Topo Boundary is defined continue with build process
if (0 < len(section.SEG_Boundary_DCT)) or (0 < len(section.SEG_Boundary_OCC)) :
    print "INFO:\t Boundary defined for segment 0"
    if not (0 < len(section.SEG_Boundary_OCC)):
        print 'INFO:\t The segment boundary is only defined by discrete values. Applying interpolation algorithm to obtain OCC definition'
        #Interpolation Algorithm to Obtain OCC definition from discrete boundary values. 
        
        #Check if DCT_Definition is closed, if not: close it
        if not np.array_equal(DCT_data[0],DCT_data[-1]):
            print 'INFO:\t Closing open discrete boundary definition'
            DCT_data = np.concatenate((DCT_data,DCT_data[0:1,:]),axis=0)
            section.SEG_Boundary_DCT[0] = DCT_data #update it to section definition            
            
        #Find corners and edges of data
        DCT_angles = calc_DCT_angles(DCT_data)
        
        #Split DCT_data into steady segments 
        min_degree = 140   #allowed angle in discrete representation before starting to split
        
        corners = []
        for i in range(0,DCT_angles.shape[0]):
            if DCT_angles[i] < min_degree: 
                corners.append(i)
        NbCorners = np.size(corners)
        
        
        DCT_Segments = []        
        if NbCorners == 0:
            DCT_Segments[0] = DCT_data
        
        if NbCorners > 0:
            for i in range(0,NbCorners-1):
                print i,corners[i]
                DCT_Segments.append(DCT_data[corners[i]:corners[i+1]+1])
                
                
        #plt.plot(DCT_data[:,0],DCT_data[:,1], color='black', marker='.')
        #for i in range(0,len(DCT_Segments)):
        #    plt.plot(DCT_Segments[i][:,0],DCT_Segments[i][:,1],marker='o')
        
                
        #If no corner is found, curve is periodic!
        #If one corner is found, only one segment exists with that is non periodic
        #If two corner are found, two segments exists that are non periodic 
                        #TBD: Analyse discrete curve, check for edges and periodicity to create a OCC definition!
        
        #FOR EACH discretized Segment create OCC definition by interpolation        
        
        wire = 	BRepBuilderAPI_MakeWire()
        for i,item in enumerate(DCT_Segments):
            print i
            data = item.T
            print data
            tmp_harray = TColgp_HArray1OfPnt2d_from_nparray(data)
            
            if NbCorners == 0:
                tmp_interpolation = Geom2dAPI_Interpolate(tmp_harray.GetHandle(), True, 0.0000001)             #Interpolate datapoints to bspline
            else:     
                tmp_interpolation = Geom2dAPI_Interpolate(tmp_harray.GetHandle(), False, 0.0000001)             #Interpolate datapoints to bspline
                 
            tmp_interpolation.Perform()                                               
            tmp_bspline = tmp_interpolation.Curve().GetObject()
            
            #Convert to 3D SPACE onto x,y plane
            P = gp_Pnt(0,0,0)
            V = gp_Dir(gp_Vec(0,0,1))
            Plane = Geom_Plane(P, V)
            tmp_edge = BRepBuilderAPI_MakeEdge(tmp_bspline.GetHandle(),Plane.GetHandle())
            wire.Add(tmp_edge.Edge())
            
        section.SEG_Boundary_OCC.append(wire.Wire())
        #display.DisplayShape(wire.Wire(), color='YELLOW')

    else: 
        print 'INFO:\t The segment boundary is defined in OCC definition' 
        
else:
    print "ERROR:\t No boundary defined for segment 0"
    
    
#======================================================================
#      C R E A T E   S E G.   C O O R D S.
#======================================================================
tmp_wire = section.SEG_Boundary_OCC[0]
display.DisplayShape(wire.Wire(), color='WHITE')
tmp_shape = wire.Shape()
prop = GProp_GProps()
test = brepgprop_LinearProperties(tmp_shape, prop)

BRepAdaptor_Curve2d 

wire.Shape()

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


#======================================================================
#PLOT




#plt.axis('equal')
#plt.show()


#display.FitAll()
#start_display()
#END    