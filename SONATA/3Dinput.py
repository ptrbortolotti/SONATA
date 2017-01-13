import sys
import os

from OCC.gp import gp_Pnt2d, gp_Pnt, gp_Pln, gp_Dir, gp_Vec
from OCC.GC import GC_MakePlane, GC_MakeEllipse
from OCC.Geom import Geom_Plane, Geom_RectangularTrimmedSurface
from OCC.STEPControl import STEPControl_Reader
from OCC.IGESControl import IGESControl_Reader
from OCC.StlAPI import StlAPI_Reader
from OCC.TopoDS import TopoDS_Shape
from OCC.BRepTools import breptools_Read
from OCC.BRep import BRep_Builder
from OCC.IFSelect import IFSelect_RetDone, IFSelect_ItemsByEntity
from OCC.Display.SimpleGui import init_display
from OCC.BRepAlgoAPI import BRepAlgoAPI_Section
from OCC.BRepAlgo import BRepAlgo_Section
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace
from OCC.TopExp import TopExp_Explorer
from OCC.TopTools import TopTools_ListIteratorOfListOfShape, TopTools_ListOfShape
import OCC.TopoDS as TopoDS
from OCC.BRepAdaptor import  BRepAdaptor_Curve2d, BRepAdaptor_Curve

from OCC.TopOpeBRep import TopOpeBRep_ShapeIntersector
from wire_utils import build_wire_from_BSplineLst

def load_stp(filename):
    step_reader = STEPControl_Reader()
    status = step_reader.ReadFile(filename)
     
    if status == IFSelect_RetDone:  # check status
        failsonly = False
        step_reader.PrintCheckLoad(failsonly, IFSelect_ItemsByEntity)
        step_reader.PrintCheckTransfer(failsonly, IFSelect_ItemsByEntity)
     
        ok = step_reader.TransferRoot(1)
        _nbs = step_reader.NbShapes()
        aResShape = step_reader.Shape(1)
    else:
        print("Error: can't read file.")
        sys.exit(0)

    return aResShape

def load_igs(filename):
    #READ IGS:
    iges_reader = IGESControl_Reader()
    status = iges_reader.ReadFile(filename)
    
    if status == IFSelect_RetDone:  # check status
        failsonly = False
        iges_reader.PrintCheckLoad(failsonly, IFSelect_ItemsByEntity)
        iges_reader.PrintCheckTransfer(failsonly, IFSelect_ItemsByEntity)
        ok = iges_reader.TransferRoots()
        aResShape = iges_reader.Shape(1)
    else:
        print("Error: can't read file.")
        sys.exit(0)

def load_stl(filename):
    #READ STL:
    stl_reader = StlAPI_Reader()
    stl_shp = TopoDS_Shape()
    stl_reader.Read(stl_shp, filename)
    return stl_shp    
    
def load_brep(filename):
    brep_shp = TopoDS_Shape()
    builder = BRep_Builder()
    breptools_Read(brep_shp, filename, builder)
    return brep_shp
 
 
def load_3D(filename):        
    name, file_extension = os.path.splitext(filename)
    if file_extension == '.stp':
        aResShape = load_stp(filename)   
    elif file_extension == '.igs':
         aResShape = load_igs(filename)
    elif file_extension == '.stl':
         aResShape = load_stl(filename)
    return aResShape


def Section_Shape_to_2dBSplineLst(TopoDS_Shape, Pnt, Dir):
    
    #Define the Point and direction of the slicing plane
    #Dir = gp_Dir(1., 0., 0.); Pnt = gp_Pnt(500,0,0)
    Pln = gp_Pln(Pnt, Dir)
    face = BRepBuilderAPI_MakeFace(Pln).Shape()
    TopoDS_Face = TopoDS.topods().Face(face)
    
    # Computes Shape/Plane intersection
    section = BRepAlgoAPI_Section(TopoDS_Shape, face)

    ex = TopExp_Explorer(section.Shape(),6,7) #Search for Edges(6), Exclude Vertices(7)
    
    results = []
    while ex.More():
        edge = TopoDS.topods().Edge(ex.Current())
        results.append(edge)
        #print "is null?", bool(edge.IsNull())
        ex.Next()
   
    BSplineLst = []
    for edg in results:
        #print "null now?", bool(edg.IsNull())
        Adaptor = BRepAdaptor_Curve2d(edg,TopoDS_Face)
        BSplineLst.append(Adaptor.BSpline().GetObject())
    
    return BSplineLst 
    
    
    
    
    
# =============================================================================
#            IF __MAIN__
# =============================================================================  
if __name__ == '__main__': 
    display, start_display, add_menu, add_function_to_menu = init_display('wx')
    display.Context.SetDeviationAngle(0.000001)       # 0.001 default. Be careful to scale it to the problem.
    display.Context.SetDeviationCoefficient(0.00001) # 0.001 default. Be careful to scale it to the problem. 
    display.set_bg_gradient_color(20,6,111,200,200,200)
    
    filename = 'AREA_Blatt_L.stl'
    aResShape = load_3D(filename) 
    display.DisplayShape(aResShape, color=None, transparency=0.7, update=True)
    
    # Define the direction
    D = gp_Dir(1., 0., 0.) 
    P = gp_Pnt(500,0,0)

    Pln = gp_Pln(P, D)
    face = BRepBuilderAPI_MakeFace(Pln,).Shape()
    display_face = BRepBuilderAPI_MakeFace(Pln,- 100., 100., -120., 120).Shape()
    TopoDS_Face = TopoDS.topods().Face(face)
    
    # Computes Shape/Plane intersection
    section = BRepAlgoAPI_Section(aResShape, face, False)       
    section.Approximation(True)
    section.Build()
    
    ex = TopExp_Explorer(section.Shape(),6,7) #Search for Edges(6), Exclude Vertices(7)
    results = []
    NbEdges = 0 
    while ex.More():
        edge = TopoDS.topods().Edge(ex.Current())
        results.append(edge)
        NbEdges += 1
        ex.Next()
    
    BSplineLst = []
    for edg in results:
        Adaptor = BRepAdaptor_Curve(edg,TopoDS_Face)
        BSplineLst.append(Adaptor.BSpline().GetObject())
       
    TopTools_EdgeLst = TopTools_ListOfShape()
    tmp_wire = 	BRepBuilderAPI_MakeWire()
    for i,item in enumerate(BSplineLst):
        Plane = Geom_Plane(P, D)
        tmp_edge = BRepBuilderAPI_MakeEdge(item.GetHandle())
        #print tmp_edge.IsDone()
        display.DisplayShape(tmp_edge.Edge())
        TopTools_EdgeLst.Append(tmp_edge.Edge())
        
    #wire = BRepBuilderAPI_MakeWire()
    #wire.Add(TopTools_EdgeLst)
    #wire.Build()
    #display.DisplayShape(wire.Wire())
    #wire = wire.Wire()
    
    #print "WIRE CLOSED: " + str(wire.Closed())
    #transformedwire = build_wire_from_BSplineLst(BSplineLst)     
        
        

# =============================================================================
#            DISPLAY
# =============================================================================


    
    #display.DisplayShape(gp_Pnt(500,0,0))
    #display.DisplayShape(wire)
    #display.DisplayShape(display_face, color="BLUE", transparency=0.8, update=True)
    #display.DisplayShape(transformedwire)
    #display.DisplayShape(section.Shape(),color="BLACK")
    #display.DisplayShape(section2.Shape(),color="GREEN")
    start_display()