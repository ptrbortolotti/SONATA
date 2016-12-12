#Basic Libraries:
import numpy as np
from scipy.optimize import leastsq

#PythonOCC Libraries
from OCC.gp import gp_Pnt2d, gp_Dir2d
from OCC.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.GCPnts import GCPnts_AbscissaPoint, GCPnts_QuasiUniformDeflection
from OCC.Geom2dAPI import Geom2dAPI_Interpolate, Geom2dAPI_InterCurveCurve
from OCC.Geom2d import Handle_Geom2d_BSplineCurve_DownCast, Geom2d_Line

#Own Libraries:
from utils import calc_DCT_angles, TColgp_HArray1OfPnt2d_from_nparray, Pnt2dLst_to_npArray, \
                    discrete_stepsize, curvature_of_curve, isclose, unique_rows, \
                    P2Pdistance


###############################################################################
# BSpline and BSplineLst Utilities
###############################################################################
def findPnt_on_2dcurve(Pnt,curve,u0=0):
    def Pnt_Distance(u,Pnt):
        u = float(u[0])
        p = gp_Pnt2d()
        curve.D0(u,p)
        error = p.Distance(Pnt)
        return error
    
    y = leastsq(Pnt_Distance,u0,args=(Pnt))
    #error = Pnt_Distance(y,Pnt)
    #print('u:',float(y[0]),'Error:', error)
    u = float(y[0])
    return u

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
        Adaptor = Geom2dAdaptor_Curve(item.GetHandle())
        CummLength += get_BSpline_length(item)
        if x < CummLength or isclose(x,CummLength):
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
    if S1 > S2:
        trimmed_BSplineLst = []
        para1 =  find_BSplineLst_coordinate(BSplineLst, S1, start, end)
        para2 =  find_BSplineLst_coordinate(BSplineLst, S2, start, end)
        for i,item in enumerate(BSplineLst):
            First = item.FirstParameter() 
            Last =  item.LastParameter()
            if para1[0] == i and para2[0] != i: #Okay
                 BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
                 BSplineCopy.Segment(para1[1],Last)
                 trimmed_BSplineLst.append(BSplineCopy)
                 
            elif (para1[0] != i and para2[0] != i) and (para1[0] > i and para2[0] > i):
                 BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
                 BSplineCopy.Segment(First,Last)
                 trimmed_BSplineLst.append(BSplineCopy)
                 
            elif (para1[0] != i and para2[0] != i) and (para1[0] < i and para2[0] < i):
                 BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy()).GetObject()
                 BSplineCopy.Segment(First,Last)
                 trimmed_BSplineLst.append(BSplineCopy)
                 
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
                 trimmed_BSplineLst.append(BSplineCopy)
                 break
             

    elif S1 == start and S2 == end:
        trimmed_BSplineLst = copy_BSplineLst(BSplineLst)
    
    elif S2 > S1:
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

    
def seg_boundary_from_dct(DCT_data,angular_deflection = 30):
    #Check if DCT_Definition is closed, if not: close it
    if not np.array_equal(DCT_data[0],DCT_data[-1]):
        print 'INFO:\t Closing open discrete boundary definition'
        DCT_data = np.concatenate((DCT_data,DCT_data[0:1,:]),axis=0)
        
    #Find corners and edges of data
    DCT_angles = calc_DCT_angles(DCT_data)
    corners = []
    for i in range(0,DCT_angles.shape[0]):
        if DCT_angles[i] < (180 - angular_deflection) or DCT_angles[i] > (180 + angular_deflection): 
            corners.append(i)
    NbCorners = np.size(corners)
    
    #Segment by Corners    
    DCT_Segments = []
    if NbCorners == 0:
        DCT_Segments.append(DCT_data)
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
            #tmp_interpolation = Geom2dAPI_Interpolate(tmp_harray.GetHandle(), True, 0.0000001)             #Interpolate datapoints to bspline
            tmp_interpolation = Geom2dAPI_Interpolate(tmp_harray.GetHandle(), False, 0.0000001)             #Interpolate datapoints to bspline
        else:     
            tmp_interpolation = Geom2dAPI_Interpolate(tmp_harray.GetHandle(), False, 0.0000001)             #Interpolate datapoints to bspline
        tmp_interpolation.Perform()                                               
        tmp_bspline = tmp_interpolation.Curve().GetObject()
        list_of_bsplines.append(tmp_bspline)
        
    return list_of_bsplines

def discretize_BSplineLst(BSplineLst,Deflection=0.00001):
    Pnt2dLst = []
    for i,item in enumerate(BSplineLst):
        Adaptor = Geom2dAdaptor_Curve(item.GetHandle())
        discretization = GCPnts_QuasiUniformDeflection(Adaptor,Deflection,4)	#GeomAbs_Shape Continuity: 1=C0, 2=G1, 3=C1, 3=G2,... 
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

    return npArray
    
         
    
def BSplineLst_from_dct(DCT_data,angular_deflection=30):
           
    #Find corners and edges of data
    DCT_angles = calc_DCT_angles(DCT_data)
    corners = []
    for i in range(0,DCT_angles.shape[0]):
        if DCT_angles[i] < (180 - angular_deflection) or DCT_angles[i] > (180 + angular_deflection): 
            corners.append(i)
    NbCorners = np.size(corners)
    
    
    #Segmenting data according to corners    
    #====================================
    DCT_Segments = [] 
    
    #Closed:       
    if np.array_equal(DCT_data[0],DCT_data[-1]): 
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
        tmp_harray = TColgp_HArray1OfPnt2d_from_nparray(data)
        if closed == True and NbCorners == 0:
            tmp_interpolation = Geom2dAPI_Interpolate(tmp_harray.GetHandle(), False, 1e-06)             #Interpolate datapoints to bspline   
        else:
            tmp_interpolation = Geom2dAPI_Interpolate(tmp_harray.GetHandle(), False, 1e-06)             #Interpolate datapoints to bspline
            
        tmp_interpolation.Perform()                              
        tmp_bspline = tmp_interpolation.Curve().GetObject()
        list_of_bsplines.append(tmp_bspline)
        
    return list_of_bsplines
    
    
    
    
def set_BSplineLst_to_Origin(BSplineLst):
    '''
    The Origin is determined by the most right Intersection Point of the X-Axis with the segment boundary 
    '''   
    xaxis = Geom2d_Line(gp_Pnt2d(0,0),gp_Dir2d(1,0))   #x-axis
    tolerance=1e-10
    IntPnts = []
   
    for i,item in enumerate(BSplineLst): 
        First = item.FirstParameter()
        Last = item.LastParameter()
        intersection = Geom2dAPI_InterCurveCurve(xaxis.GetHandle(), item.GetHandle(),tolerance)
        for j in range(1,intersection.NbPoints()+1):
                IntPnt = intersection.Point(j)
                XValue = IntPnt.X()
                u = findPnt_on_2dcurve(IntPnt,item)
                IntPnts.append([i,u,XValue])
                
    #Determine Origin as point                 
    IntPntsarray = np.asarray(IntPnts)  #idx,W,XValue
    OriEdgePnt = IntPntsarray[np.argmax(IntPntsarray[:,2]),:]                    
     
    #Reorder Sequence of BSplines of BSplinesLst
    OBSplineLst =  []
    CorrectOrigin = False
              
    for i,item in enumerate(BSplineLst):     
        if i  == OriEdgePnt[0]:
            First = item.FirstParameter() 
            Last =  item.LastParameter()
            if isclose(OriEdgePnt[1],First) == False:
                CorrectOrigin = False
                BSplineCurve1 = copy_BSpline(item)
                BSplineCurve1.Segment(OriEdgePnt[1],Last)
                BSplineCurve2 = copy_BSpline(item)
                BSplineCurve2.Segment(First,OriEdgePnt[1])
                OBSplineLst.append(BSplineCurve1)
            else:
                OBSplineLst.append(item)
                CorrectOrigin = True
                
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