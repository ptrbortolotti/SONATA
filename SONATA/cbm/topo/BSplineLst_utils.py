#Basic Libraries:
import math
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

#PythonOCC Libraries
from OCC.gp import gp_Pnt2d, gp_Dir2d, gp_Vec2d
from OCC.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.GCPnts import GCPnts_AbscissaPoint, GCPnts_QuasiUniformDeflection,GCPnts_UniformDeflection, GCPnts_TangentialDeflection, GCPnts_QuasiUniformAbscissa
from OCC.Geom2dAPI import Geom2dAPI_Interpolate, Geom2dAPI_InterCurveCurve, Geom2dAPI_ProjectPointOnCurve, Geom2dAPI_PointsToBSpline
from OCC.Geom2d import Handle_Geom2d_BSplineCurve_DownCast, Geom2d_Line
from OCC.GeomAPI import GeomAPI_Interpolate
#Own Libraries:
from SONATA.cbm.topo.utils import calc_DCT_angles, TColgp_HArray1OfPnt2d_from_nparray, Pnt2dLst_to_npArray, \
                    discrete_stepsize, curvature_of_curve, isclose, unique_rows, \
                    P2Pdistance, PolygonArea, TColgp_Array1OfPnt2d_from_nparray, TColgp_HArray1OfPnt_from_nparray

from SONATA.cbm.topo.para_Geom2d_BsplineCurve import Para_Geom2d_BSplineCurve

###############################################################################
# BSpline and BSplineLst Utilities
###############################################################################

def corners_of_BSplineLst(BSplineLst):
    corners = [] 
    for item in BSplineLst:
        corners.append(item.EndPoint())
        
    corners.pop(-1)
    return corners #gp_Pnt2d Lst
    

def ProjectPointOnBSplineLst(BSplineLst,Pnt2d,tolerance_distance=100):
    p2 = []
    for idx,item in enumerate(BSplineLst):
        projection2 = Geom2dAPI_ProjectPointOnCurve(Pnt2d,item.GetHandle())
        for j in range(1,projection2.NbPoints()+1):
            if projection2.Distance(j)<=tolerance_distance:
                p2.append([projection2.Point(j), idx, projection2.Parameter(j),projection2.Distance(j)])
            else: None   
    
    if p2:
        p2 = np.asarray(p2)
        min_index = p2[:,3].argmin() 
        return p2[min_index,:]
    else:
        return []

def distance_on_BSplineLst(BSplineLst,para1,para2):
#    para1 = findPnt_on_BSplineLst(P1,BSplineLst)
#    para2 = findPnt_on_BSplineLst(P2,BSplineLst)
    
#    if closed:
#        print para1,para2
#        para2 = [len(BSplineLst)-1, BSplineLst[-1].LastParameter()]
#        print para1,para2
    tol=1e-7
    #print para1,para2
    
    Distance = 0
    for i,item in enumerate(BSplineLst):
        Adaptor = Geom2dAdaptor_Curve(item.GetHandle())
        First = item.FirstParameter() 
        Last =  item.LastParameter() 
        if para1[0] == i and para2[0] != i:
            Length = GCPnts_AbscissaPoint().Length(Adaptor, para1[1],Last, tol)
            Distance += Length
             
        elif (para1[0] != i and para2[0] != i) and (para1[0] < i and para2[0] > i):
            Length = GCPnts_AbscissaPoint().Length(Adaptor, First, Last, tol)
            Distance += Length
             
        elif para1[0] == i and para2[0] == i:
            Length = GCPnts_AbscissaPoint().Length(Adaptor, para1[1],para2[1], tol)
            Distance += Length
            break
    
    #NOTE: Somehow the execption if para2 is close to Last doesn't work properly!!!!!!! Similar to the trimm function     
        elif para1[0] != i and para2[0] == i:
            if isclose(para2[1],Last):
                 Length = GCPnts_AbscissaPoint().Length(Adaptor, First,Last, tol)
                 Distance += Length
            else:
                 Length = GCPnts_AbscissaPoint().Length(Adaptor, First,para2[1], tol)
                 Distance += Length
            break

    return Distance		


def isPnt_on_2dcurve(Pnt2d,Curve2d,tolerance=1e-6):
    projection = Geom2dAPI_ProjectPointOnCurve(Pnt2d,Curve2d)
    Trigger = False
    for i in range(1,projection.NbPoints()+1):
        if projection.Distance(i) <= tolerance:
            Trigger = True
        else:
            None           
    return Trigger


def isPnt_on_BSplineLst(Pnt2d,BSplineLst,tolerance=1e-6):
    Bool = False
    for item in BSplineLst:
        if isPnt_on_2dcurve(Pnt2d,item.GetHandle(),tolerance):
            Bool = True
            break
    return Bool

def findPnt_on_2dcurve(Pnt2d,Curve2d):
    projection = Geom2dAPI_ProjectPointOnCurve(Pnt2d,Curve2d)
    u = projection.LowerDistanceParameter()
    return u
    
def findPnt_on_BSplineLst(Pnt2d,BSplineLst):
    for i,item in enumerate(BSplineLst):
        if isPnt_on_2dcurve(Pnt2d,item.GetHandle()):
            u = findPnt_on_2dcurve(Pnt2d,item.GetHandle())
            coordinates = [i,u]
            break
        else:
            coordinates = None
    return coordinates

    
def get_BSpline_length(BSpline):
    tolerance=1e-10
    first = BSpline.FirstParameter()
    last = BSpline.LastParameter()
    Adaptor = Geom2dAdaptor_Curve(BSpline.GetHandle())
    Length = GCPnts_AbscissaPoint().Length(Adaptor, first, last, tolerance)
    return Length

    
def copy_BSpline(BSpline):
    BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(BSpline.Copy()).GetObject()
    return BSplineCopy
  
    
def copy_BSplineLst(BSplineLst):
    copied_BSplineLst = []
    for i,item in enumerate(BSplineLst):
        copied_BSplineLst.append(copy_BSpline(item))
    return copied_BSplineLst 
    
    
def get_BSplineLst_length(BSplineLst):
    #Return the cummulated length of the BSplineLst
    CummLength = 0
    for i,item in enumerate(BSplineLst):
         CummLength += get_BSpline_length(item)
    return CummLength  


def equidistant_Points_on_BSplineLst(BSplineLst,minLen):
    ''' minLen = the minimum distance between points. Should be determined by the segment0.BSplineLstLenght/Resolution
    '''
    Pnt2dLst = []
    for item in BSplineLst:
        Adaptor = Geom2dAdaptor_Curve(item.GetHandle())
        length = get_BSpline_length(item)
        NbPoints = int(length//minLen)+2  
        discretization = GCPnts_QuasiUniformAbscissa(Adaptor,NbPoints)
        
        for j in range(1, NbPoints+1):
                para = discretization.Parameter(j)
                Pnt = gp_Pnt2d()
                item.D0(para,Pnt)
                Pnt2dLst.append(Pnt)

    return Pnt2dLst

    
def find_BSpline_coordinate(BSpline,s):
    # Be careful, s stands for the lenght coordinate of a single BSpline, while S represents the Global Coordinate!
    BSpline_length = get_BSpline_length(BSpline)
    tolerance=1e-10
    Adaptor = Geom2dAdaptor_Curve(BSpline.GetHandle())
    if s <= BSpline_length:
             tmp = GCPnts_AbscissaPoint(tolerance, Adaptor, s, 0)
             u = tmp.Parameter()
    return u                            

def find_BSplineLst_coordinate(BSplineLst, S, start, end):
    '''The function find_BSplineLst_coordinate returns the list index of 
    the bspline and its parameter where the coordinate is located.
    
    Args:
        S: (float) Coordinate to be found 
        start: (float) start of the interval where the BSplineLst is defined
        end: (float) end of the interval where the BSplineLst is defined
        
    Returns: 
        [i,U]: #Return [index of the bspline,parameter U on the bspline]
        
    '''
    BSplineLstLength = get_BSplineLst_length(BSplineLst)
    x = BSplineLstLength*(S-start)/(end-start)
    CummLength = 0
    tolerance=1e-8
    #print 'S: '+ str(S)
    for i,item in enumerate(BSplineLst):
        first = item.FirstParameter()
        last = item.LastParameter()
        Adaptor = Geom2dAdaptor_Curve(item.GetHandle())
        CummLength += get_BSpline_length(item)
        if x < CummLength or isclose(x,CummLength):
             dist = x - (CummLength - GCPnts_AbscissaPoint().Length(Adaptor, first, last, tolerance))
             tmp = GCPnts_AbscissaPoint(tolerance, Adaptor, dist, first)
             U = tmp.Parameter()
             break
    #print 'idx: '+ str(i) + '     U: '+ str(U) 
    return [i,U]

def find_BSplineLst_pos(BSplineLst,para):
    idx = para[0]
    U = para[1]
    tol=1e-7
    BSplineLst_length = get_BSplineLst_length(BSplineLst)
    CummLength = 0
    for x in range(0, idx):
        #print 'x',x
        CummLength += get_BSpline_length(BSplineLst[x])    

    Adaptor = Geom2dAdaptor_Curve(BSplineLst[idx].GetHandle())
    First = BSplineLst[idx].FirstParameter() 
    partial_length = GCPnts_AbscissaPoint().Length(Adaptor, First, U, tol)     
    S = (CummLength + partial_length)/float(BSplineLst_length)
    return S


def reverse_BSplineLst(BSplineLst):
    BSplineLst = list(reversed(BSplineLst))
    for i, item in enumerate(BSplineLst):
        item.Reverse()
    return BSplineLst

  
def get_BSplineLst_Pnt2d(BSplineLst,S,start,end):
    if start>end: #transfrom the interval if start>end to a basic interval [0..]
        D = 1-start
        S = (S+D)%1
        end = (end+D)%1
        start = (start+D)%1 
        
    P = gp_Pnt2d()
    [idx,U] = find_BSplineLst_coordinate(BSplineLst,S, start, end)
    BSplineLst[idx].D0(U,P)
    return P

def get_BSplineLst_D2(BSplineLst,S, start, end):
    if start>end: #transfrom the interval if start>end to a basic interval [0..]
        D = 1-start
        S = (S+D)%1
        end = (end+D)%1
        start = (start+D)%1 
    
    P = gp_Pnt2d()
    V1 = gp_Vec2d()
    V2 = gp_Vec2d()
    [idx,U] = find_BSplineLst_coordinate(BSplineLst,S, start, end)
    BSplineLst[idx].D2(U,P,V1,V2)
    return P,V1,V2

    
def intersect_BSplineLst_with_BSpline(BSplineLst,BSpline):
    tolerance=1e-10
    IntPnts = []
    IntPnts_Pnt2d = []
    for i,item in enumerate(BSplineLst): 
        intersection = Geom2dAPI_InterCurveCurve(BSpline.GetHandle(), item.GetHandle(),tolerance)
        for j in range(1,intersection.NbPoints()+1):
                IntPnt = intersection.Point(j)
                IntPnts_Pnt2d.append(IntPnt)
                u = findPnt_on_2dcurve(IntPnt,item.GetHandle())
                IntPnts.append([i,u])
    return IntPnts, IntPnts_Pnt2d

    

def BSplineLst_Orientation(BSplineLst, NbPoints=31):   

    #FASTER BUT NOT EQUIDISTANT APPROACH!    
    Pnt2dLst = []   
    for item in BSplineLst:
        First = item.FirstParameter()
        Pnt2dLst.append(item.StartPoint())
        Last = item.LastParameter()
        inc = (Last-First)/NbPoints
        para = First+inc
        while para<=Last:
            p = gp_Pnt2d()
            item.D0(para,p)
            Pnt2dLst.append(p)
            para += inc
    b = Pnt2dLst_to_npArray(Pnt2dLst)   
        
#   Place Equidistant Points on BSplineLst:                             #THIS APPROACH WAS SUPER SLOW!
#    a = np.linspace(0.0, 1.0, num=NbPoints, endpoint=True)
#    Pnt2dLst = []
#    for s in a:
#        Pnt2dLst.append(get_BSplineLst_Pnt2d(BSplineLst,s, 0, 1))
#    b = Pnt2dLst_to_npArray(Pnt2dLst)
    
    #Calculate of Polygon.
    area = PolygonArea(b)
    #print area
    if area > 0: #counter-clockwise=True
        return True
    elif area < 0: #clockwise=False
        return False
    else: return None
    
    
def trim_BSplineLst(BSplineLst, S1, S2, start, end):
    ''' the trim_BSplineLst function trims the BSplineLst that is defined on 
    the interval between start and end to the interval defined by S1 and S2.
    
    Args:
        BSplineLst: (list of geom2d_bsplines) is the object to be modified.
        S1: (float) start of the interval to which the object is trimmed.
        S2: (float) end of the interval to which the object is trimmed.
        start: (float) start of the interval where the BSplineLst is defined
        end: (float) end of the interval where the BSplineLst is defined
            
    Returns: 
        trimmed_BSplineLst: (list of geom2d_bsplines) is the modified object.
            
    Example: a BSplineLst is defined between 0.0 and 1.0 and shall be trimmed 
            to a interval S1=0.3 and S2=0.7
    '''
    if start>end: #transfrom the interval if start>end to a basic interval [0..]
        D = 1-start
        S1 = (S1+D)%1
        S2 = (S2+D)%1
        end = (end+D)%1
        start = (start+D)%1 
    
    if S1 > S2:
        trimmed_BSplineLst = []
        front_BSplineLst = []
        rear_BSplineLst = []
        para1 =  find_BSplineLst_coordinate(BSplineLst, S1, start, end)
        para2 =  find_BSplineLst_coordinate(BSplineLst, S2, start, end)
        for i,item in enumerate(BSplineLst):
            First = item.FirstParameter() 
            Last =  item.LastParameter()
            if para1[0] == i and para2[0] != i: #Okay
                 BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
                 if isclose(para1[1],Last):
                     pass 
                 elif isclose(para1[1],First):
                     BSplineCopy.Segment(First,Last)
                     trimmed_BSplineLst.append(BSplineCopy)
                 else:                    
                     BSplineCopy.Segment(para1[1],Last)
                     trimmed_BSplineLst.append(BSplineCopy)
                 
            elif (para1[0] != i and para2[0] != i) and (para1[0] > i and para2[0] > i):
                 BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
                 BSplineCopy.Segment(First,Last)
                 front_BSplineLst.append(BSplineCopy)
                 
            elif (para1[0] != i and para2[0] != i) and (para1[0] < i and para2[0] < i):
                 BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
                 BSplineCopy.Segment(First,Last)
                 rear_BSplineLst.append(BSplineCopy)
                 
            elif para1[0] == i and para2[0] == i: #Okay
                 BSplineCopy1 = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
                 BSplineCopy1.Segment(para1[1],Last)
                 trimmed_BSplineLst.append(BSplineCopy1)
                 BSplineCopy2 = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
                 BSplineCopy2.Segment(First,para2[1])
                 trimmed_BSplineLst.append(BSplineCopy2)
                 break
        
            elif para1[0] != i and para2[0] == i: #Okay
                 BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
                 BSplineCopy.Segment(First,para2[1])
                 front_BSplineLst.append(BSplineCopy)
                 #break
        
        if front_BSplineLst and rear_BSplineLst:
            trimmed_BSplineLst =  rear_BSplineLst + front_BSplineLst
        
        elif front_BSplineLst and not rear_BSplineLst:
            trimmed_BSplineLst = front_BSplineLst
            
        elif not front_BSplineLst and rear_BSplineLst:
            trimmed_BSplineLst = rear_BSplineLst 
             

    elif S1 == start and S2 == end:
        trimmed_BSplineLst = copy_BSplineLst(BSplineLst)
        
    elif S2 == start or S1 == S2:
        trimmed_BSplineLst = []
    
    elif S2 > S1:
        trimmed_BSplineLst = []
        #check if para1 and para2 is out of bounds!
        para1 =  find_BSplineLst_coordinate(BSplineLst, S1, start, end)
        para2 =  find_BSplineLst_coordinate(BSplineLst, S2, start, end)
        for i,item in enumerate(BSplineLst):       
             First = item.FirstParameter() 
             Last =  item.LastParameter() 
             if para1[0] == i and para2[0] != i:
                 BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
                 if isclose(para1[1],Last):
                     pass 
                 elif isclose(para1[1],First):
                     BSplineCopy.Segment(First,Last)
                     trimmed_BSplineLst.append(BSplineCopy)
                 else: 
                     BSplineCopy.Segment(para1[1],Last)
                     trimmed_BSplineLst.append(BSplineCopy)
                 
             elif (para1[0] != i and para2[0] != i) and (para1[0] < i and para2[0] > i):
                 BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
                 BSplineCopy.Segment(First,Last)
                 trimmed_BSplineLst.append(BSplineCopy)
                 
             elif para1[0] == i and para2[0] == i:
                 BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
                 BSplineCopy.Segment(para1[1],para2[1])
                 trimmed_BSplineLst.append(BSplineCopy)
                 break
        
             elif para1[0] != i and para2[0] == i:
                 BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
                 if isclose(para2[1],Last):
                     trimmed_BSplineLst.append(BSplineCopy)
                 else:
                     BSplineCopy.Segment(First,para2[1])
                     trimmed_BSplineLst.append(BSplineCopy)
                 break
             

    return trimmed_BSplineLst		


def trim_BSplineLst_by_coordinates(BSplineLst, para1,para2):
    trimmed_BSplineLst = []
#    para1 =  find_BSplineLst_coordinate(BSplineLst, S1, start, end)
#    para2 =  find_BSplineLst_coordinate(BSplineLst, S2, start, end)
    for i,item in enumerate(BSplineLst):       
         First = item.FirstParameter() 
         Last =  item.LastParameter() 
         if para1[0] == i and para2[0] != i:
             BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
             BSplineCopy.Segment(para1[1],Last)
             trimmed_BSplineLst.append(BSplineCopy)
             
         elif (para1[0] != i and para2[0] != i) and (para1[0] < i and para2[0] > i):
             BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
             BSplineCopy.Segment(First,Last)
             trimmed_BSplineLst.append(BSplineCopy)
             
         elif para1[0] == i and para2[0] == i:
             BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
             BSplineCopy.Segment(para1[1],para2[1])
             trimmed_BSplineLst.append(BSplineCopy)
             break
    
    #NOTE: Somehow the execption if para2 is close to Last doesn't work properly!!!!!!!     
         elif para1[0] != i and para2[0] == i:
             BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
             if isclose(para2[1],Last):
                 trimmed_BSplineLst.append(BSplineCopy)
             else:
                 BSplineCopy.Segment(First,para2[1])
                 trimmed_BSplineLst.append(BSplineCopy)
             break
         
    #             elif para1[0] != i and para2[0] == i:
    #                 BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
    #                 BSplineCopy.Segment(First,para2[1])
    #                 trimmed_BSplineLst.append(BSplineCopy)
    #                 break
    #                
                
    return trimmed_BSplineLst	


    
def trim_BSplineLst_by_Pnt2d(BSplineLst,Pos1_Pnt2d,Pos2_Pnt2d):
    trimmed_BSplineLst= []
    front_BSplineLst = []
    rear_BSplineLst =[]
    para1 = findPnt_on_BSplineLst(Pos1_Pnt2d,BSplineLst)
    para2 = findPnt_on_BSplineLst(Pos2_Pnt2d,BSplineLst)
    
    #if para2>para1
    if (para2[0] > para1[0]) or ((para2[0] == para1[0]) and (para2[1]>para1[1])):
        for i,item in enumerate(BSplineLst):  
            First = item.FirstParameter() 
            Last =  item.LastParameter() 
            if para1[0] == i and para2[0] != i:
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
                BSplineCopy.Segment(para1[1],Last)
                trimmed_BSplineLst.append(BSplineCopy)
                 
            elif (para1[0] != i and para2[0] != i) and (para1[0] < i and para2[0] > i):
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
                BSplineCopy.Segment(First,Last)
                trimmed_BSplineLst.append(BSplineCopy)
                 
            elif para1[0] == i and para2[0] == i:
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
                BSplineCopy.Segment(para1[1],para2[1])
                trimmed_BSplineLst.append(BSplineCopy)
                break
        
            elif para1[0] != i and para2[0] == i:
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
                BSplineCopy.Segment(First,para2[1])
                trimmed_BSplineLst.append(BSplineCopy)
                break
    
    #if para1>para2
    else:
        #para1, para2 = para2, para1
        for i,item in enumerate(BSplineLst):
            First = item.FirstParameter() 
            Last =  item.LastParameter()
            if para1[0] == i and para2[0] != i: #Okay
                 BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
                 BSplineCopy.Segment(para1[1],Last)
                 rear_BSplineLst.append(BSplineCopy)
                 
            elif (para1[0] != i and para2[0] != i) and (para1[0] > i and para2[0] > i):
                 BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
                 BSplineCopy.Segment(First,Last)
                 front_BSplineLst.append(BSplineCopy)
                 
            elif (para1[0] != i and para2[0] != i) and (para1[0] < i and para2[0] < i):
                 BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
                 BSplineCopy.Segment(First,Last)
                 rear_BSplineLst.append(BSplineCopy)
                 
            elif para1[0] == i and para2[0] == i: #Okay
                 BSplineCopy1 = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
                 BSplineCopy1.Segment(para1[1],Last)
                 trimmed_BSplineLst.append(BSplineCopy1)
                 BSplineCopy2 = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
                 BSplineCopy2.Segment(First,para2[1])
                 trimmed_BSplineLst.append(BSplineCopy2)
                 break
        
            elif para1[0] != i and para2[0] == i: #Okay
                 BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
                 BSplineCopy.Segment(First,para2[1])
                 front_BSplineLst.append(BSplineCopy)
                 #break

        if front_BSplineLst and rear_BSplineLst:
            trimmed_BSplineLst =  rear_BSplineLst + front_BSplineLst
        
        elif front_BSplineLst and not rear_BSplineLst:
            trimmed_BSplineLst = front_BSplineLst
            
        elif not front_BSplineLst and rear_BSplineLst:
            trimmed_BSplineLst = rear_BSplineLst

    return trimmed_BSplineLst



    
def seg_boundary_from_dct(DCT_data,angular_deflection = 30):
    #OUTDATED!!!
    #Check if DCT_Definition is closed, if not: close it
    if not np.allclose(DCT_data[0],DCT_data[-1]):
        print('INFO:\t Closing open discrete boundary definition')
        DCT_data = np.concatenate((DCT_data,DCT_data[0:1,:]),axis=0)
    
    #Find corners and edges of data
    DCT_angles = calc_DCT_angles(DCT_data)
    print(DCT_data.shape, DCT_angles.shape)
    corners = []
    for i in range(0,DCT_angles.shape[0]):
        if DCT_angles[i] < (180 - angular_deflection) or DCT_angles[i] > (180 + angular_deflection): 
            corners.append(i)
    NbCorners = np.size(corners)
    
    #Segment by Corners    
    print('NbCorners = ', NbCorners)
    DCT_Segments = []
    if NbCorners == 0:
        DCT_Segments.append(DCT_data)
    if NbCorners > 0:
        DCT_Segments.append(DCT_data[0:corners[0]+1])
        for i in range(0,NbCorners-1):
            print('i, corners = ', i,corners[i])
            DCT_Segments.append(DCT_data[corners[i]:corners[i+1]+1])
        DCT_Segments.append(DCT_data[corners[-1]:-1])
    
    list_of_bsplines = []
    for i,item in enumerate(DCT_Segments):
        data = item.T
               
        tmp_harray = TColgp_HArray1OfPnt2d_from_nparray(data)
        if NbCorners == 0:
            #tmp_interpolation = Geom2dAPI_Interpolate(tmp_harray.GetHandle(), True, 1e-6)             #Interpolate datapoints to bspline
            tmp_interpolation = Geom2dAPI_Interpolate(tmp_harray.GetHandle(), False, 1e-6)             #Interpolate datapoints to bspline
        else:     
            tmp_interpolation = Geom2dAPI_Interpolate(tmp_harray.GetHandle(), False, 1e-6)             #Interpolate datapoints to bspline
        tmp_interpolation.Perform()                                               
        tmp_bspline = tmp_interpolation.Curve().GetObject()
        list_of_bsplines.append(tmp_bspline)
        
    return list_of_bsplines


def discretize_BSplineLst(BSplineLst, Deflection = 2e-4, AngularDeflection=0.02,CurvatureDeflection=0.06,MinimumOfPoints=50):
    Pnt2dLst = []
    for i,item in enumerate(BSplineLst):
        Adaptor = Geom2dAdaptor_Curve(item.GetHandle())
        discretization = GCPnts_QuasiUniformDeflection(Adaptor,Deflection, 1)	#GeomAbs_Shape Continuity: 1=C0, 2=G1, 3=C1, 3=G2,... 
        #discretization = GCPnts_TangentialDeflection(Adaptor,AngularDeflection,CurvatureDeflection,MinimumOfPoints)
        NbPoints = discretization.NbPoints()
        for j in range(1, NbPoints+1):
            Pnt = discretization.Value(j)
            Pnt2dLst.append(Pnt)
    
    a = Pnt2dLst_to_npArray(Pnt2dLst)
    b = np.around(a,10) # Evenly round to the given number of decimals. 
    
    #check if it is closed:
    if np.array_equal(b[0],b[-1]):
        closed = True
    else: closed = False
    
    npArray = unique_rows(b) # Remove possible doubles! 
    
    if closed:
        npArray = np.vstack((npArray,b[0]))
    else: None
    #print len(npArray)
    return npArray
             
    
def BSplineLst_from_dct(DCT_data, angular_deflection=15, closed=False, tol_interp = 1e-6, twoD = True):
    #Close DCT_date if closed == True
    if closed and not np.allclose(DCT_data[0],DCT_data[-1]):
        print('INFO:\t Closing open discrete definition')
        DCT_data = np.concatenate((DCT_data,DCT_data[0:1,:]),axis=0)   
    
    
    #TODP: check all datapoints and merge if allclose is true 
    
    #Find corners and edges of data
    #print(DCT_data)
    DCT_angles = calc_DCT_angles(DCT_data)
    #print(DCT_angles)
    corners = []
    for i in range(0,DCT_angles.shape[0]):
        if DCT_angles[i] < (180 - angular_deflection) or DCT_angles[i] > (180 + angular_deflection): 
            corners.append(i)
    NbCorners = np.size(corners)
    
    
    #Segmenting data according to corners    
    #====================================
    DCT_Segments = [] 
    
    #Closed:       
    if np.allclose(DCT_data[0],DCT_data[-1]): 
        closed = True
        if NbCorners == 0:
            DCT_Segments.append(DCT_data)
       
        elif NbCorners > 0:

            if corners[0] == 0 and corners[-1] == len(DCT_data)-1:
                for i in range(0,NbCorners-2):                
                   DCT_Segments.append(DCT_data[corners[i]:corners[i+1]+1])   
                DCT_Segments.append(DCT_data[corners[-2]:corners[-1]+1])         #last segment
            
            #regular case:
            else:
                DCT_Segments.append(DCT_data[:corners[0]+1])        #first               
                if NbCorners > 1:                                   #intermediate segments:
                   for i in range(0,NbCorners-1):                
                        DCT_Segments.append(DCT_data[corners[i]:corners[i+1]+1])             
                DCT_Segments.append(DCT_data[corners[-1]:])         #last 
                DCT_Segments[0] = np.concatenate((DCT_Segments[-1][:-1],DCT_Segments[0]))
                del DCT_Segments[-1]

    else:#Open:
        closed = False
        if NbCorners == 0:
            DCT_Segments.append(DCT_data)
        
        if NbCorners > 0:
            DCT_Segments.append(DCT_data[:corners[0]+1])        #first               
            if NbCorners > 1:                                   #intermediate segments:
               for i in range(0,NbCorners-1):                
                    DCT_Segments.append(DCT_data[corners[i]:corners[i+1]+1])             
            DCT_Segments.append(DCT_data[corners[-1]:])         #last segment

    
    list_of_bsplines = []
    for i,item in enumerate(DCT_Segments):
        data = item.T

        if twoD:
            tmp_harray = TColgp_HArray1OfPnt2d_from_nparray(data)
        else:
            tmp_harray = TColgp_HArray1OfPnt_from_nparray(data)
        
        try: 
            if twoD:      
                tmp_interpolation = Geom2dAPI_Interpolate(tmp_harray.GetHandle(), False, tol_interp)             #Interpolate datapoints to bspline
            else:
                tmp_interpolation = GeomAPI_Interpolate(tmp_harray.GetHandle(), False, tol_interp)
        except: 
            print('ERROR: BSplineLst_from_dct did not perform the Interpolation with tol_interp', tol_interp)#RAISED ERROR
#            plt.figure()
#            plt.clf()
#            plt.axis('equal')  
#            plt.plot(*item.T, marker='.')
#            for i, it in enumerate(item):
#                plt.annotate(i, (it[0],it[1]), color='black')   
#            plt.show()                   
                   
        tmp_interpolation.Perform()                              
        tmp_bspline = tmp_interpolation.Curve().GetObject()
        #print(tmp_bspline.Degree(), tmp_bspline.NbKnots())
        list_of_bsplines.append(tmp_bspline)
        
    return list_of_bsplines
    
    
    
def set_BSplineLst_to_Origin(BSplineLst, Theta=0):
    '''
    The Origin is determined by the most right Intersection Point of the X-Axis with the segment boundary, the axis is rotated backwards by Theta (in degree). 
    '''
    
    Theta=Theta/float(180)*np.pi
    x = math.cos(Theta)
    y = -math.sin(Theta)
    OPnt = gp_Pnt2d(0,0)
    dir2d = gp_Dir2d(x,y)
    axis = Geom2d_Line(OPnt,dir2d)   #x-axis
    v0 = gp_Vec2d(dir2d)
    tolerance=1e-10
    IntPnts = []

    for i,item in enumerate(BSplineLst): 
        First = item.FirstParameter()
        Last = item.LastParameter()
        intersection = Geom2dAPI_InterCurveCurve(axis.GetHandle(), item.GetHandle(),tolerance)
        for j in range(1,intersection.NbPoints()+1):
                IntPnt = intersection.Point(j)
                v1 = gp_Vec2d(OPnt,IntPnt) 
                vres = v0.Dot(v1)
                #distance = IntPnt.Distance(OPnt)
                u = findPnt_on_2dcurve(IntPnt,item.GetHandle())
                IntPnts.append([i,u,vres])
                
    #Determine Origin as point                 
    IntPntsarray = unique_rows(np.asarray(IntPnts)) #idx,W,XValue
    OriEdgePnt = IntPntsarray[np.argmax(IntPntsarray[:,2]),:]    
    #print(IntPntsarray,OriEdgePnt)   
             
    #Reorder Sequence of BSplines of BSplinesLst
    OBSplineLst =  []
    CorrectOrigin = False
              
    for i,item in enumerate(BSplineLst):     
        if i  == OriEdgePnt[0]:
            First = item.FirstParameter() 
            Last =  item.LastParameter()
            
            if isclose(OriEdgePnt[1],First) == True:
                OBSplineLst.append(item)
                CorrectOrigin = True
            
            elif isclose(OriEdgePnt[1],Last) == True:
                CorrectOrigin = False
                BSplineCurve2 = item
            
            else: 
                CorrectOrigin = False
                BSplineCurve1 = copy_BSpline(item)
                BSplineCurve1.Segment(OriEdgePnt[1],Last)
                BSplineCurve2 = copy_BSpline(item)
                BSplineCurve2.Segment(First,OriEdgePnt[1])
                OBSplineLst.append(BSplineCurve1)
   
        elif i > OriEdgePnt[0]:
            OBSplineLst.append(item)
        else:
            None

    for i,item in enumerate(BSplineLst):
        if i < OriEdgePnt[0]:
            OBSplineLst.append(item)
        else:
            None
            
    if CorrectOrigin == False:
        OBSplineLst.append(BSplineCurve2)    
    else: None

    return OBSplineLst




if __name__ == '__main__':
    pass
    

