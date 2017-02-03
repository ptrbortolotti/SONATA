#Basic Libraries:
import numpy as np       
import matplotlib.pyplot as plt

#PythonOCC Libraries
from OCC.gp import gp_Pnt2d, gp_Pnt, gp_Pln, gp_Dir, gp_Vec, gp_Trsf
from OCC.GCPnts import GCPnts_AbscissaPoint
from OCC.BRepAdaptor import BRepAdaptor_CompCurve, BRepAdaptor_Curve
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeEdge, BRepBuilderAPI_Transform
from OCC.BRepLib import BRepLib_MakeFace
from OCC.BRepTools import BRepTools_WireExplorer
from OCC.IntCurvesFace import IntCurvesFace_Intersector
from OCC.Geom import Geom_Plane
from OCC.TopTools import TopTools_ListIteratorOfListOfShape, TopTools_ListOfShape
from OCC.TopoDS import topods

#Own Libraries:
from explorer import WireExplorer
from utils import Pnt2dLst_to_npArray, PolygonArea

###############################################################################
# Wire Utilities
###############################################################################

def Unique_EdgeLst(EdgeLst):
    Preci = 1e-3
    BSplineLst3D = []
    NonUniques = []
    for edg in EdgeLst:
        Adaptor = BRepAdaptor_Curve(edg)
        BSplineLst3D.append(Adaptor.BSpline().GetObject())
    
    for i,bspline in enumerate(BSplineLst3D):
        for j,item in enumerate(BSplineLst3D):
            if bspline.IsEqual(item.GetHandle(),Preci) and i!=j:
                NonUniques.append(j)
                
    return NonUniques
        
def Wire_Orientation(wire, NbPoints=31):
    #Place Equidistant Points on BSplineLst:
    #VERY SLOW
    a = np.linspace(0.0, 1.0, num=NbPoints, endpoint=True)
    Pnt2dLst = []
    for s in a:
        Pnt2d = get_wire_Pnt2d(wire,s)
        Pnt2dLst.append(Pnt2d)
    b = Pnt2dLst_to_npArray(Pnt2dLst)
    
    
    plt.figure(3)
    plt.clf()         
    plt.plot(*b.T, color='black', marker='.')
    plt.axis('equal')  
    plt.show()   
    
    
    
    #Calculate of Polygon.
    area = PolygonArea(b)
    #print area
    if area > 0: #counter-clockwise=True
        return True
    elif area < 0: #clockwise=False
        return False
    else: return None

      

def NbEdges_in_wire(Wire):
    ex = BRepTools_WireExplorer(Wire)
    counter = 0
    while ex.More():
        counter +=1
        ex.Next()
    return counter


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
    

def rotate_wire(wire,ax1,angle):    #for top
    aTrsf = gp_Trsf()    
    aTrsf.SetRotation(ax1,angle)
    aBRespTrsf = BRepBuilderAPI_Transform(wire, aTrsf)
    RotatedWire = aBRespTrsf.Shape()    
    rotWire = topods.Wire(RotatedWire)
    return rotWire


def translate_wire(wire,gp_Pnt1,gp_Pnt2):
    aTrsf = gp_Trsf()
    aTrsf.SetTranslation(gp_Pnt1,gp_Pnt2)
    aBRespTrsf = BRepBuilderAPI_Transform(wire, aTrsf)
    atraslatedShape = aBRespTrsf.Shape()
    translateWire = topods.Wire(atraslatedShape)
    return translateWire

def scale_wire(wire,gp_Pnt1,factor):
    aTrsf = gp_Trsf()
    aTrsf.SetScale(gp_Pnt1,factor)
    aBRespTrsf = BRepBuilderAPI_Transform(wire, aTrsf)
    aScaledShape = aBRespTrsf.Shape()
    scaledWire = topods.Wire(aScaledShape)
    return scaledWire

	
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


def trim_wire(TopoDS_wire, S1, S2):
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

    
#def build_wire_from_BSplineLst(BSplineLst): #Builds TopoDS_Wire from connecting BSplineSegments and returns it        
#    tmp_wire = 	BRepBuilderAPI_MakeWire()
#    for i,item in enumerate(BSplineLst):
#        P = gp_Pnt(0,0,0)
#        V = gp_Dir(gp_Vec(0,0,1))
#        Plane = Geom_Plane(P, V)
#        tmp_edge = BRepBuilderAPI_MakeEdge(item.GetHandle(),Plane.GetHandle())
#        tmp_wire.Add(tmp_edge.Edge())
#    return tmp_wire.Wire()


def build_wire_from_BSplineLst(BSplineLst): #Builds TopoDS_Wire from connecting BSplineSegments and returns it        
    P = gp_Pnt(0,0,0)
    V = gp_Dir(gp_Vec(0,0,1))
    Plane = Geom_Plane(P, V)          
    TopTools_EdgeLst = TopTools_ListOfShape()
    for i,item in enumerate(BSplineLst):
        tmp_edge = BRepBuilderAPI_MakeEdge(item.GetHandle(),Plane.GetHandle())
        #print tmp_edge.IsDone()
        TopTools_EdgeLst.Append(tmp_edge.Edge())
        
    wire = BRepBuilderAPI_MakeWire()
    wire.Add(TopTools_EdgeLst)
    wire.Build()
    return wire.Wire()


    
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
    

      