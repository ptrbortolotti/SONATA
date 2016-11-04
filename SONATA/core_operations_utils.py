# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 14:10:57 2016

@author: TPflumm
"""
#-------------------------------
#          H E A D E R
#-------------------------------


import urllib2                                       # help in opening URL
import numpy as np           

#Python OCC Libraries
from OCC.gp import gp_Pnt, gp_Vec, gp_Pnt2d, gp_Pln, gp_Dir, gp_Trsf,gp_Ax1, gp_OX, gp_Ax3
from OCC.Geom2dAPI import Geom2dAPI_PointsToBSpline, Geom2dAPI_Interpolate 
from OCC.GeomAPI import geomapi
from OCC.TColgp import TColgp_Array1OfPnt, TColgp_Array1OfPnt2d, TColgp_HArray1OfPnt2d, TColgp_HArray1OfPnt 
from OCC.BRepBuilderAPI import (BRepBuilderAPI_MakeEdge,
                                BRepBuilderAPI_MakeVertex,
                                BRepBuilderAPI_MakeWire,
                                BRepBuilderAPI_Transform,
                                BRepBuilderAPI_MakeFace)
from OCC.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCC.Display.SimpleGui import init_display
from OCC.TopoDS import topods, TopoDS_Edge, TopoDS_Compound
from OCC.Display.WebGl import x3dom_renderer
from OCC.Display.SimpleGui import init_display

#from OCCUtils.Construct import make_loft
#-------------------------------
#          FUNCTIONS
#-------------------------------
                                 


def UIUCAirfoil(name):
    foil_dat_url = 'http://www.ae.illinois.edu/m-selig/ads/coord_seligFmt/%s.dat' % name
    f = urllib2.urlopen(foil_dat_url)
    temp_x = []
    temp_y = []
    temp_z = []
    for line in f.readlines()[1:]:                                                                  # The first line contains info only
        line = line.lstrip().rstrip().replace('    ', ' ').replace('   ', ' ').replace('  ', ' ')   # do some cleanup on the data (mostly dealing with spaces)
        data = line.split(' ')   
        temp_y.append(float(data[0]))                                                               # data[0] = x coord.
        temp_z.append(float(data[1]))                                                               # data[1] = y coord.
        temp_x.append(float(0))                                                                     # data[2] = z coord 
    return np.array([temp_x,temp_y,temp_z])                                                         # return AirfoilCoordinate as np.arrray

def UIUCAirfoil2d(name):
    foil_dat_url = 'http://www.ae.illinois.edu/m-selig/ads/coord_seligFmt/%s.dat' % name
    f = urllib2.urlopen(foil_dat_url)
    temp_x = []
    temp_y = []
    for line in f.readlines()[1:]:                                                                  # The first line contains info only
        line = line.lstrip().rstrip().replace('    ', ' ').replace('   ', ' ').replace('  ', ' ')   # do some cleanup on the data (mostly dealing with spaces)
        data = line.split(' ')   
        temp_x.append(float(data[0]))                                                               # data[0] = x coord.
        temp_y.append(float(data[1]))                                                               # data[1] = y coord.
    return np.array([temp_x,temp_y])    

def AirfoilDat(name):
    string = "%s.dat" %(name)
    f = open(string)
    temp_x = []
    temp_y = []
    temp_z = []
    for line in f.readlines()[1:]:                                                                  # The first line contains info only
        line = line.lstrip().rstrip().replace('    ', ' ').replace('   ', ' ').replace('  ', ' ')   # do some cleanup on the data (mostly dealing with spaces)
        data = line.split(' ')   
        temp_y.append(float(data[0]))                                                               # data[0] = x coord.
        temp_z.append(float(data[1]))                                                               # data[1] = y coord.
        temp_x.append(float(0))                                                                     # data[2] = z coord 
    return np.array([temp_x,temp_y,temp_z])                                                         # return AirfoilCoordinate as np.arrray    
    
def AirfoilDat2d(name):
    string = "%s.dat" %(name)
    f = open(string)
    temp_x = []
    temp_y = []
    for line in f.readlines()[1:]:                                                                  # The first line contains info only
        line = line.lstrip().rstrip().replace('    ', ' ').replace('   ', ' ').replace('  ', ' ')   # do some cleanup on the data (mostly dealing with spaces)
        data = line.split(' ')   
        temp_x.append(float(data[0]))                                                             
        temp_y.append(float(data[1]))                                                                                                       
    return np.array([temp_x,temp_y])                                                                # return AirfoilCoordinate as np.arrray    

def np_array_to_gp_Pnt(array):              
    Pnt = gp_Pnt(array[0], array[1], array[2])
    return Pnt


def TColgp_HArray1OfPnt_from_nparray(data):
    #Create TColgp_HArray1OfPnt from np.array															
    # data[0] = x coord,  data[1] = y coord,  data[2] = z coord 
    harray = TColgp_HArray1OfPnt(1, np.ma.size(data,1))
    for index in range(0, np.ma.size(data,1)):
        i = index+1
        harray.SetValue(i, gp_Pnt(float(data[0,index]),float(data[1,index]),float(data[2,index])))
    return harray
    
def TColgp_Array1OfPnt2d_from_nparray(data):
    #Create TColgp_HArray1OfPnt from np.array															
    # data[0] = x coord,  data[1] = y coord, 
    array = TColgp_Array1OfPnt2d(1, np.ma.size(data,1))
    for index in range(0, np.ma.size(data,1)):
        i = index+1
        array.SetValue(i, gp_Pnt2d(float(data[0,index]),float(data[1,index])))
    return array  
    
def TColgp_HArray1OfPnt2d_from_nparray(data):
    #Create TColgp_HArray1OfPnt from np.array															
    # data[0] = x coord,  data[1] = y coord, 
    harray = TColgp_HArray1OfPnt2d(1, np.ma.size(data,1))
    for index in range(0, np.ma.size(data,1)):
        i = index+1
        harray.SetValue(i, gp_Pnt2d(float(data[0,index]),float(data[1,index])))
    return harray  


def point_list_to_TColgp_Array1OfPnt(li):
    pts = TColgp_Array1OfPnt(0, len(li) - 1)
    for n, i in enumerate(li):
        pts.SetValue(n, i)
    return pts

def point2d_list_to_TColgp_Array1OfPnt2d(li):
    return _Tcol_dim_1(li, TColgp_Array1OfPnt2d)

def point2d_list_to_TColgp_HArray1OfPnt2d(li):
    return _Tcol_dim_1(li, TColgp_HArray1OfPnt2d)
    
def point_list_to_TColgp_Array1OfPnt(li):
    return  _Tcol_dim_1(li, TColgp_Array1OfPnt)
    
def point_list_to_TColgp_HArray1OfPnt(li):
    return  _Tcol_dim_1(li, TColgp_HArray1OfPnt)      

def _Tcol_dim_1(li, _type):
    pts = _type(0, len(li)-1)
    for n, i in enumerate(li):
        pts.SetValue(n, i)
    return pts

def create_list_of_2Dpoints(data):                                             #CREATE OCC array of Points
    pts_2d = []
    for index,value in np.ndenumerate(data[0,:]):
        pts_2d.append(gp_Pnt2d(float(data[0,index]),float(data[1,index])))
    return pts_2d

def create_list_of_points(data):                                             #CREATE OCC array of Points
    pts = []
    for index,value in np.ndenumerate(data[0,:]):
        pts.append(gp_Pnt(float(data[0,index]),float(data[1,index]),float(data[2,index])))
    return pts

def create_array_of_points(data):                                             #CREATE OCC array of Points
    pts = []
    for index,value in np.ndenumerate(data[0,:]):
        
        pts.append(gp_Pnt(float(data[0,index]),float(data[1,index]),float(data[2,index])))
    return pts

   
def make_wire(*args):
    # if we get an iterable, than add all edges to wire builder
    if isinstance(args[0], list) or isinstance(args[0], tuple):
        wire = BRepBuilderAPI_MakeWire()
        for i in args[0]:
            wire.Add(i)
        wire.Build()
        return wire.Wire()
    wire = BRepBuilderAPI_MakeWire(*args)
    return wire.Wire()

def make_edge(*args):
    edge = BRepBuilderAPI_MakeEdge(*args)
    result = edge.Edge()
    return result

def rotate_TopoDS_wire(wire,ax1,angle):    #for top
    aTrsf = gp_Trsf()    
    aTrsf.SetRotation(ax1,angle)
    aBRespTrsf = BRepBuilderAPI_Transform(wire, aTrsf)
    RotatedWire = aBRespTrsf.Shape()    
    rotWire = topods.Wire(RotatedWire)
    return rotWire
    
def translate_TopoDS_wire(wire,gp_Pnt1,gp_Pnt2):
    aTrsf = gp_Trsf()
    aTrsf.SetTranslation(gp_Pnt1,gp_Pnt2)
    aBRespTrsf = BRepBuilderAPI_Transform(wire, aTrsf)
    atraslatedShape = aBRespTrsf.Shape()
    translateWire = topods.Wire(atraslatedShape)
    return translateWire

def scale_TopoDS_wire(wire,gp_Pnt1,factor):
    aTrsf = gp_Trsf()
    aTrsf.SetScale(gp_Pnt1,factor)
    aBRespTrsf = BRepBuilderAPI_Transform(wire, aTrsf)
    aScaledShape = aBRespTrsf.Shape()
    scaledWire = topods.Wire(aScaledShape)
    return scaledWire
 

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    
    
 
  


def gp_Pnt2d_to_npArray(Ptn2d):
        vector = np.array([Ptn2d.X(),Ptn2d.Y()])
        return vector

def gp_Vec2d_to_npArray(Vec2d):
        vector = np.array([Vec2d.X(),Vec2d.Y()])
        return vector
        
def np_GetNormal2d(Vec2d):
        vector = np.array([-Vec2d[1],Vec2d[0]])
        return vector
        
def npArray_to_gp_Pnt2d(vector):
        Pnt2d = gp_Pnt2d(vector[0],vector[1])
        return Pnt2d     
