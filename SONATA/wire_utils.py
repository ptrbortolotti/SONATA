#Basic Libraries:
import numpy as np       

#PythonOCC Libraries
from OCC.gp import gp_Pnt2d, gp_Pnt, gp_Pln, gp_Dir, gp_Vec
from OCC.GCPnts import GCPnts_AbscissaPoint
from OCC.BRepAdaptor import BRepAdaptor_CompCurve, BRepAdaptor_Curve
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeEdge
from OCC.BRepLib import BRepLib_MakeFace
from OCC.IntCurvesFace import IntCurvesFace_Intersector
from OCC.Geom import Geom_Plane

#Own Libraries:
from explorer import WireExplorer

###############################################################################
# Wire Utilities
###############################################################################
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

    
def build_wire_from_BSplineLst(BSplineLst): #Builds TopoDS_Wire from connecting BSplineSegments and returns it        
    tmp_wire = 	BRepBuilderAPI_MakeWire()
    for i,item in enumerate(BSplineLst):
        P = gp_Pnt(0,0,0)
        V = gp_Dir(gp_Vec(0,0,1))
        Plane = Geom_Plane(P, V)
        tmp_edge = BRepBuilderAPI_MakeEdge(item.GetHandle(),Plane.GetHandle())
        tmp_wire.Add(tmp_edge.Edge())
    return tmp_wire.Wire()

    
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
    

      