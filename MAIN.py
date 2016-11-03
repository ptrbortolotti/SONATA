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
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge,BRepBuilderAPI_MakeWire,BRepBuilderAPI_MakeFace
from OCC.BRepOffsetAPI import BRepOffsetAPI_MakeOffset
from OCC.BRepLib import BRepLib_MakeFace 
from OCC.BRep import BRep_Tool_Pnt
from OCC.Geom import Geom_Plane, Geom_TrimmedCurve, Geom_RectangularTrimmedSurface
from OCC.GeomAPI import GeomAPI_IntCS, GeomAPI_ExtremaCurveSurface
from OCC.Geom2dAPI import Geom2dAPI_PointsToBSpline, Geom2dAPI_Interpolate, Geom2dAPI_InterCurveCurve 
from OCC.BRepCheck import BRepCheck_Wire
from OCC.Display.SimpleGui import init_display
from OCC.BRepGProp import BRepGProp_EdgeTool,brepgprop
from OCC.GProp import GProp_GProps
from OCC.BRepAdaptor import BRepAdaptor_Curve2d, BRepAdaptor_Curve, 	BRepAdaptor_CompCurve
from OCC.GCPnts import GCPnts_AbscissaPoint
from OCC.IntCurvesFace import IntCurvesFace_Intersector

from OCCUtils.Construct import (make_closed_polygon, make_n_sided,
                                make_vertex, make_face, trim_wire,make_offset)
from OCCUtils.Topology import WireExplorer, Topo
from OCCUtils.base import GlobalProperties, BaseObject
from OCCUtils.types_lut import ShapeToTopology


#-------------------------------
#          FUNCTIONS
#-------------------------------


#======================================================
#       MAIN
#======================================================
#START    
###############################################################################
display, start_display, add_menu, add_function_to_menu = init_display('wx')
display.Context.SetDeviationAngle(0.000001)      # 0.001 default
display.Context.SetDeviationCoefficient(0.000001) # 0.001 default

#CREATE AXIS SYSTEM for Visualization
COSY = gp_Ax3()	
O  = gp_Pnt(0., 0., 0.)
p1 = gp_Pnt(0.3,0.,0.)
p2 = gp_Pnt(0.,0.1,0.)
p3 = gp_Pnt(0.,0.,0.1)

h1 = BRepBuilderAPI_MakeEdge(O,p1).Shape()
h2 = BRepBuilderAPI_MakeEdge(O,p2).Shape()
h3 = BRepBuilderAPI_MakeEdge(O,p3).Shape()

display.DisplayShape(O,color='BLACK')
display.DisplayShape(h1,color='RED')
display.DisplayShape(h2,color='GREEN')
display.DisplayShape(h3,color='BLUE')


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
                #print i,corners[i]
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
            #print i
            data = item.T
            #print data
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
#      DETERMINE  S E G.  with  ORIGIN:
#====================================================================== 
 
def set_BoundaryWire_to_Origin(TopoDS_wire):
    '''
    The Origin is determined by the most right Intersection Point of the X-Axis with the segment boundary 
    '''   
    Face = BRepLib_MakeFace(gp_Pln(gp_Pnt(0,0,0), gp_Dir(0,1,0)),-2,2,-2,2).Face()
    tolerance=1e-10
    intPnts = []
    idx = 0
    for edg in WireExplorer(TopoDS_wire).ordered_edges():       #Iterate over Wire! 
        Adaptor = BRepAdaptor_Curve(edg)
        First = Adaptor.FirstParameter()
        Last = Adaptor.LastParameter()
        Adaptor3d_HCurve = Adaptor.Trim(First,Last,tolerance)
        intersection = IntCurvesFace_Intersector(Face, tolerance)
        intersection.Perform(Adaptor3d_HCurve,First,Last)
        if intersection.IsDone():
            for j in range (1,intersection.NbPnt()+1):
                XValue = intersection.Pnt(j).X()
                W = intersection.WParameter(j)
                intPnts.append([idx,W,XValue])
        idx += 1   
    
            
    #Determine Origin as point                 
    IntPntsarray = np.asarray(intPnts)  #idx,W,XValue
    OriEdgePnt = IntPntsarray[np.argmax(IntPntsarray[:,2]),:]                    
     
    
    #Reorder Sequence of edges of wire
    Owire =  BRepBuilderAPI_MakeWire()
    idx = 0                 
    for edg in WireExplorer(TopoDS_wire).ordered_edges():       #Iterate over Wire!        
        if idx  == OriEdgePnt[0]:
            Adaptor = BRepAdaptor_Curve(edg)
            First = Adaptor.FirstParameter() 
            Last =  Adaptor.LastParameter()
            BSplineCurve1 = Adaptor.BSpline().GetObject()
            BSplineCurve1.Segment(OriEdgePnt[1],Last)
            BSplineCurve2 = Adaptor.BSpline().GetObject()
            BSplineCurve2.Segment(First,OriEdgePnt[1])
            tmp_edge1 = BRepBuilderAPI_MakeEdge(BSplineCurve1.GetHandle())
            tmp_edge2 = BRepBuilderAPI_MakeEdge(BSplineCurve2.GetHandle())
            Owire.Add(tmp_edge1.Edge())
            
        elif idx > OriEdgePnt[0]:
            Owire.Add(edg)
        
        else:
            None
        idx += 1
      
        
    idx = 0
    for edg in WireExplorer(TopoDS_wire).ordered_edges():
        if idx < OriEdgePnt[0]:
            Owire.Add(edg)
        else:
            None
        idx += 1
           
    Owire.Add(tmp_edge2.Edge())       
    Owire.Build()

    return Owire.Wire()
    
#======================================================================
#      DETERMINE  S E G.   C O O R D S.
#======================================================================
'''  Determin Seg. Coords. -> create Layer -> Determine new. Seg. Coords, ...
TBD: think of data structure to store information about layup, with all wires, lists of edges, and so on. 
 TBD:     
    

neccesarry functions:   - find_innermost_boundary
                        - gp_Pnt2d P= get_wire_point(TopoDS_Wire, std_real s)
                        - std_real s = get_wire_coord(TopoDS_Wire, gp_Pnt2d P)
                        - virtual_void = get_wire_D2 (TopoDS_Wire,std_real s, gp_Pnt2d P, gp_Vec2d V1, gp_Vec2d V2)
                        - std_real L = get_wire_length(TopoDS_Wire)

                        - TopoDS_Wire = trim_wire_to_interval(TopoDS_Wire, real S1,real S2)
                        - discretize            
                         
'''


def get_wire_length(TopoDS_wire):   #std_real L = get_wire_length(TopoDS_Wire)
    AdaptorComp = BRepAdaptor_CompCurve(TopoDS_wire, True)
    length = AdaptorComp.LastParameter()-AdaptorComp.FirstParameter()
    return length
  
def get_wire_Pnt2d(TopoDS_wire, S):     #gp_Pnt2d P= get_wire_point(TopoDS_Wire, std_real s)
    length = get_wire_length(TopoDS_wire)
    AdaptorComp = BRepAdaptor_CompCurve(TopoDS_wire, True)
    P = AdaptorComp.Value(S*length)
    return gp_Pnt2d(P.X(),P.Y())

def get_wire_Pnt(TopoDS_wire, S):     #gp_Pnt2d P= get_wire_point(TopoDS_Wire, std_real s)
    length = get_wire_length(TopoDS_wire)
    AdaptorComp = BRepAdaptor_CompCurve(TopoDS_wire, True)
    P = AdaptorComp.Value(S*length)
    return P 
    
    
def find_coordinate_on_ordered_edges(TopoDS_wire,S):
    WireLength = get_wire_length(TopoDS_wire)      
    CummLength = 0
    tolerance=1e-10
    idx = 0
    for edg in WireExplorer(TopoDS_wire).ordered_edges():       #Iterate over Wire! 
         Adaptor = BRepAdaptor_Curve(edg)
         CummLength += GCPnts_AbscissaPoint().Length(Adaptor, tolerance)
         if S*WireLength <= CummLength:
             U0 = 0
             dist = GCPnts_AbscissaPoint().Length(Adaptor, tolerance)-(CummLength-S*WireLength)
             Test = GCPnts_AbscissaPoint(tolerance, Adaptor, dist, U0)
             U = Test.Parameter()
             break
         idx += 1
    return [idx,U]  #Return index of edge and parameter U on edge!


def trim_wire_to_interval(TopoDS_wire, S1, S2):
    twire =  BRepBuilderAPI_MakeWire()
    para1 =  find_coordinate_on_ordered_edges(TopoDS_wire, S1)
    para2 =  find_coordinate_on_ordered_edges(TopoDS_wire, S2)
    idx = 0
    for edg in WireExplorer(TopoDS_wire).ordered_edges():       #Iterate over Wire! 
         Adaptor = BRepAdaptor_Curve(edg)
         First = Adaptor.FirstParameter() 
         Last =  Adaptor.LastParameter()
         #CummLength += GCPnts_AbscissaPoint().Length(Adaptor, tolerance)
         if para1[0] == idx and para2[0] != idx:
             BSplineCurve = Adaptor.BSpline().GetObject()
             BSplineCurve.Segment(para1[1],Last,)
             tmp_edge = BRepBuilderAPI_MakeEdge(BSplineCurve.GetHandle())
             #display.DisplayShape(tmp_edge.Edge(), color='CYAN') 
             twire.Add(tmp_edge.Edge())
             
         elif (para1[0] != idx and para2[0] != idx) and (para1[0] < idx and para2[0] > idx):
             BSplineCurve = Adaptor.BSpline().GetObject()
             tmp_edge = BRepBuilderAPI_MakeEdge(BSplineCurve.GetHandle())
             twire.Add(tmp_edge.Edge()) 
             
         elif para1[0] == idx and para2[0] == idx:
             BSplineCurve = Adaptor.BSpline().GetObject()
             BSplineCurve.Segment(para1[1],para2[1])
             tmp_edge = BRepBuilderAPI_MakeEdge(BSplineCurve.GetHandle())
             twire.Add(tmp_edge.Edge())
             break
    
         elif para1[0] != idx and para2[0] == idx:
             BSplineCurve = Adaptor.BSpline().GetObject()
             BSplineCurve.Segment(First,para2[1])
             tmp_edge = BRepBuilderAPI_MakeEdge(BSplineCurve.GetHandle())
             twire.Add(tmp_edge.Edge())
             break
         idx += 1
    return twire.Wire()

    
#======================================================================
#      MAIN
#======================================================================
    
 
tmp_wire = set_BoundaryWire_to_Origin(section.SEG_Boundary_OCC[0])
Owire = tmp_wire
AdaptorComp = BRepAdaptor_CompCurve(tmp_wire, True)


   
          
#testwire = BRepBuilderAPI_MakeWire(tmp_edge1.Edge()).Wire()
display.DisplayShape(get_wire_Pnt2d(Owire,0.00051))
display.DisplayShape(get_wire_Pnt2d(Owire,1))             
display.DisplayShape(Owire, color='BLUE') 
S1 = 0.4
S2 = 0.6  



pnt2d_S1 = get_wire_Pnt2d(tmp_wire,S1)
pnt2d_S2 = get_wire_Pnt2d(tmp_wire,S2)
twire = trim_wire_to_interval(tmp_wire,S1,S2)
display.DisplayShape(pnt2d_S1)   
display.DisplayShape(pnt2d_S2) 




    
#display.DisplayShape(tmp_wire, color='WHITE')
#display.DisplayShape(twire, color='ORANGE') 
#display.DisplayShape(Owire.Wire(), color='RED') 


   



         








 

##======================================================================
#if section.SETUP_NbOfWebs>0:
#    print "STATUS:\t Continue  Building Segments:"
#    
#else:
#    print "STATUS:\t Fill Segment 0 with Core"
#
##======================================================================
##Balance Weight
#if section.SETUP_BalanceWeight == True:
#    print "STATUS:\t Create Balanace Weight"
#    print "STATUS:\t Trim Geometry to Balance Weight"
#
##======================================================================
##Create MESH
#    
#
#
##======================================================================
##EXPORT
#print "STATUS:\t Exporting PATRAN file Format"
#print "STATUS:\t Exporting .STP and .IGS file Format"


#======================================================================
#PLOT




#plt.axis('equal')
#plt.show()

display.set_bg_gradient_color(20,6,111,200,200,200)
display.View_Top()
display.FitAll()
start_display()
#END    