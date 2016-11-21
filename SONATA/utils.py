## 2016 Tobias Pflumm (tobias.pflumm@tum.de)
## Note that this is adapted from pythonocc-Utils.
''' This Module describes all helpful utility functionalities of the code'''

#Basic Libraries:
import math 
import numpy as np 
                          
#Python OCC Libraries
from OCC.gp import gp_Pnt, gp_Pnt2d, gp_Trsf, gp_Vec2d
from OCC.TColgp import (TColgp_Array1OfPnt, TColgp_Array1OfPnt2d, TColgp_HArray1OfPnt2d, 
                        TColgp_HArray1OfPnt )
from OCC.BRepBuilderAPI import (BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire,
                                BRepBuilderAPI_Transform )
from OCC.TopoDS import topods 

#Own Modules:

    
#######################UTILITIE FUNCTIONS######################################    
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
    '''Returns the angle in radians between vectors 'v1' and 'v2'::'''
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

#Neccesary functions:
def calc_angle_between(v1, v2):
    '''Returns the angle in degree between vectors 'v1' and 'v2'''
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
    min_step = 0.1
    stretch = 4
    stepsize = (1-(1-min_step)*math.tanh(kappa/stretch))
    return stepsize
        

def Pnt2dLst_to_npArray(Pnt2dLst):
    lst_tmp = []
    for i,item in enumerate(Pnt2dLst):
        lst_tmp.append([item.X(),item.Y()])
    return np.asarray(lst_tmp)
    
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
        
def getID(custom):
    return custom.ID        