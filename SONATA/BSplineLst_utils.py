#Basic Libraries:
import numpy as np

#PythonOCC Libraries
from OCC.gp import gp_Pnt2d
from OCC.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.GCPnts import GCPnts_AbscissaPoint
from OCC.Geom2dAPI import Geom2dAPI_Interpolate
from OCC.Geom2d import Handle_Geom2d_BSplineCurve_DownCast

#Own Libraries:
from utils import calc_DCT_angles, TColgp_HArray1OfPnt2d_from_nparray, Pnt2dLst_to_npArray, discrete_stepsize, curvature_of_curve


###############################################################################
# BSpline and BSplineLst Utilities
###############################################################################


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

    
def find_BSpline_coordinate(BSpline,s):
    # Be careful, s stands for the lenght coordinate of a single BSpline, while S represents the Global Coordinate!
    BSpline_length = get_BSpline_length(BSpline)
    tolerance=1e-10
    Adaptor = Geom2dAdaptor_Curve(BSpline.GetHandle())
    if s <= BSpline_length:
             tmp = GCPnts_AbscissaPoint(tolerance, Adaptor, s, 0)
             u = tmp.Parameter()
    return u                            

def find_BSplineLst_coordinate(BSplineLst,S, start, end):
    BSplineLstLength = get_BSplineLst_length(BSplineLst)
    x = BSplineLstLength*(S-start)/(end-start)
    CummLength = 0
    tolerance=1e-10
    #print 'S: '+ str(S)
    for i,item in enumerate(BSplineLst):
        first = item.FirstParameter()
        last = item.LastParameter()
        #print 'first:  ' , str(item.FirstParameter())
        #print 'last:  '  , str(item.LastParameter())
        Adaptor = Geom2dAdaptor_Curve(item.GetHandle())
        CummLength += get_BSpline_length(item)
        if x <= CummLength:
             dist = x - (CummLength - GCPnts_AbscissaPoint().Length(Adaptor, first, last, tolerance))
             tmp = GCPnts_AbscissaPoint(tolerance, Adaptor, dist, first)
             U = tmp.Parameter()
             break
    #print 'idx: '+ str(i) + '     U: '+ str(U) 
    return [i,U]  #Return index of edge and parameter U on edge!
 
  
def get_BSplineLst_Pnt2d(BSplineLst,S, start, end):
    P = gp_Pnt2d()
    [idx,U] = find_BSplineLst_coordinate(BSplineLst,S, start, end)
    BSplineLst[idx].D0(U,P)
    return P
	
	
def trim_BSplineLst(BSplineLst, S1, S2, start, end):
    trimmed_BSplineLst = []
    para1 =  find_BSplineLst_coordinate(BSplineLst, S1, start, end)
    para2 =  find_BSplineLst_coordinate(BSplineLst, S2, start, end)
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
    return trimmed_BSplineLst		

    
def seg_boundary_from_dct(DCT_data,min_degree):
    #Check if DCT_Definition is closed, if not: close it
    if not np.array_equal(DCT_data[0],DCT_data[-1]):
        print 'INFO:\t Closing open discrete boundary definition'
        DCT_data = np.concatenate((DCT_data,DCT_data[0:1,:]),axis=0)
        
    #Find corners and edges of data
    DCT_angles = calc_DCT_angles(DCT_data)
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
                  
    list_of_bsplines = []
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
        list_of_bsplines.append(tmp_bspline)
        
    return list_of_bsplines
    


def discretize_BSplineLst(BSplineLst,start,end,accuracy=100):
    #Discretize_Layer:
    Pnt2dLst = []
    S = start
    idx_old = 0
    while S <= end:
        pnt2d = get_BSplineLst_Pnt2d(BSplineLst,S,start,end)
        Pnt2dLst.append(pnt2d)
        
        #grab vertices!    
        [idx,U] = find_BSplineLst_coordinate(BSplineLst,S,start,end)
        if idx > idx_old:
            BSpline = BSplineLst[idx_old]
            pnt2d = BSpline.EndPoint()
            Pnt2dLst.append(pnt2d)
                   
        else:  
            BSpline = BSplineLst[idx]
            BSplineLst_length = get_BSplineLst_length(BSplineLst)
            S += BSplineLst_length/accuracy * discrete_stepsize(curvature_of_curve(BSpline,U))
        idx_old = idx
        
    #grab last vertice    
    pnt2d = get_BSplineLst_Pnt2d(BSplineLst,end, start,end)
    Pnt2dLst.append(pnt2d)
    npArray = Pnt2dLst_to_npArray(Pnt2dLst)
    return npArray