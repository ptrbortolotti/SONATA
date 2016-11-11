#Basic Libraries:
import numpy as np

#PythonOCC Libraries
from OCC.gp import gp_Pnt2d
from OCC.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.GCPnts import GCPnts_AbscissaPoint
from OCC.Geom2dAPI import Geom2dAPI_Interpolate
from OCC.Geom2d import Handle_Geom2d_BSplineCurve_DownCast

#Own Libraries:
from utils import calc_DCT_angles, TColgp_HArray1OfPnt2d_from_nparray


###############################################################################
# BSpline and BSplineLst Utilities
###############################################################################


def get_BSpline_length(BSpline):
    tolerance=1e-10
    Adaptor = Geom2dAdaptor_Curve(BSpline.GetHandle())
    Length = GCPnts_AbscissaPoint().Length(Adaptor, tolerance)
    return Length

	
def get_BSplineLst_length(BSplineLst):
    #Return the cummulated length of the BSplineLst
    CummLength = 0
    for i,item in enumerate(BSplineLst):
         CummLength += get_BSpline_length(item)
    return CummLength  
    

def find_BSplineLst_coordinate(BSplineLst,S):
    BSplineLstLength = get_BSplineLst_length(BSplineLst)      
    CummLength = 0
    tolerance=1e-10
    for i,item in enumerate(BSplineLst):
         Adaptor = Geom2dAdaptor_Curve(item.GetHandle())
         CummLength += get_BSpline_length(item)
         if S*BSplineLstLength <= CummLength:
             dist = GCPnts_AbscissaPoint().Length(Adaptor, tolerance)-(CummLength-S*BSplineLstLength)
             tmp = GCPnts_AbscissaPoint(tolerance, Adaptor, dist, 0)
             U = tmp.Parameter()
             break
    return [i,U]  #Return index of edge and parameter U on edge!
 
   
def get_BSplineLst_Pnt2d(BSplineLst,S):
    P = gp_Pnt2d()
    [idx,U] = find_BSplineLst_coordinate(BSplineLst,S)
    BSplineLst[idx].D0(U,P)
    return P
	
	
def trim_BSplineLst(BSplineLst, S1, S2):
    trimmed_BSplineLst = []
    para1 =  find_BSplineLst_coordinate(BSplineLst, S1)
    para2 =  find_BSplineLst_coordinate(BSplineLst, S2)
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