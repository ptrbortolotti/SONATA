import sys
import os
import math

from OCC.gp import gp_Pnt2d, gp_Pnt, gp_Pln, gp_Dir, gp_Vec, gp_Trsf, gp_Ax3,gp_Ax1
from OCC.GC import GC_MakePlane, GC_MakeEllipse
from OCC.Geom import Geom_Plane, Geom_RectangularTrimmedSurface
from OCC.Geom2dAPI import Geom2dAPI_PointsToBSpline
from utils import point2d_list_to_TColgp_Array1OfPnt2d

from OCC.STEPControl import STEPControl_Reader
from OCC.IGESControl import IGESControl_Reader
from OCC.Graphic3d import (Graphic3d_EF_PDF,
                           Graphic3d_EF_SVG,
                           Graphic3d_EF_TEX,
                           Graphic3d_EF_PostScript,
                           Graphic3d_EF_EnhPostScript)
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
from OCC.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCC.TopOpeBRep import TopOpeBRep_ShapeIntersector


from wire_utils import build_wire_from_BSplineLst,rotate_wire
from BSplineLst_utils import get_BSpline_length, get_BSplineLst_length, \
                            find_BSplineLst_coordinate, get_BSplineLst_Pnt2d, \
                            trim_BSplineLst, seg_boundary_from_dct, set_BSplineLst_to_Origin, \
                            copy_BSplineLst, trim_BSplineLst_by_Pnt2d, copy_BSpline, \
                            reverse_BSplineLst

def export_to_PDF(event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('capture_pdf%s.pdf' % i):
        i += 1
    f.Export('capture_pdf%s.pdf' % i, Graphic3d_EF_PDF)
    print "EXPORT: \t Screencapture exported to capture_pdf%s.pdf" % i
    display.set_bg_gradient_color(20,6,111,200,200,200)
    
def export_to_SVG(event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('capture_svg%s.svg' % i):
        i += 1
    f.Export('capture_svg%s.svg' % i, Graphic3d_EF_SVG)
    print "EXPORT: \t Screencapture exported to capture_svg%s.svg" % i
    display.set_bg_gradient_color(20,6,111,200,200,200)
    
def export_to_PS(event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('capture_ps%s.ps' % i):
        i += 1
    f.Export('capture_ps%s.ps' % i, Graphic3d_EF_PostScript)
    print "EXPORT: \t Screencapture exported to capture_ps%s.ps" % i
    display.set_bg_gradient_color(20,6,111,200,200,200)

def export_to_EnhPS(event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('capture_Enh_ps%s.ps' % i):
        i += 1
    f.Export('capture_Enh_ps%s.ps' % i, Graphic3d_EF_EnhPostScript)
    print "EXPORT: \t Screencapture exported to capture_Enh_ps%s.ps" % i
    display.set_bg_gradient_color(20,6,111,200,200,200)
    
def export_to_TEX(event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('capture_tex%s.tex' % i):
        i += 1
    f.Export('capture_tex%s.tex' % i, Graphic3d_EF_TEX)
    print "EXPORT: \t Screencapture exported to capture_tex%s.tex" % i
    display.set_bg_gradient_color(20,6,111,200,200,200)
    
def export_to_BMP(event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('capture_bmp%s.bmp' % i):
        i += 1
    display.View.Dump('capture_bmp%s.bmp' % i)
    print "EXPORT: \t Screencapture exported to capture_bmp%s.bmp" % i
    display.set_bg_gradient_color(20,6,111,200,200,200)

def export_to_PNG(event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('capture_png%s.png' % i):
        i += 1
    display.View.Dump('capture_png%s.png' % i)
    print "EXPORT: \t Screencapture exported to capture_png%s.bmp" % i
    display.set_bg_gradient_color(20,6,111,200,200,200)

def export_to_JPEG(event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('capture_jpeg%s.jpeg' % i):
        i += 1
    display.View.Dump('capture_jpeg%s.jpeg' % i)
    print "EXPORT: \t Screencapture exported to capture_jpeg%s.jpeg" % i
    display.set_bg_gradient_color(20,6,111,200,200,200)

def export_to_TIFF(event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('capture_tiff%s.tiff' % i):
        i += 1
    display.View.Dump('capture_tiff%s.tiff' % i)
    print "EXPORT: \t Screencapture exported to capture_tiff%s.tiff" % i
    display.set_bg_gradient_color(20,6,111,200,200,200)



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


def stp2d_to_2dBSplineLst(TopoDS_Shape):
    # Define the direction
    Dir = gp_Dir(0., 0., 1.) 
    Pnt = gp_Pnt(0,0,0)

    Pln = gp_Pln(Pnt, Dir)
    face = BRepBuilderAPI_MakeFace(Pln).Shape()
    TopoDS_Face = TopoDS.topods().Face(face)

    ex = TopExp_Explorer(TopoDS_Shape,6,7) #Search for Edges(6), Exclude Vertices(7)
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
        #print Adaptor.GetType()
        
        if Adaptor.GetType() == 0: # gp_Lin2d 	
            P1 = Adaptor.Value(Adaptor.FirstParameter())
            P2 = Adaptor.Value(Adaptor.LastParameter())
            BSpline = Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([P1,P2])).Curve().GetObject()   
            BSplineLst.append(BSpline)
            
        elif Adaptor.GetType() == 6: #BSplineCurve
            BSplineLst.append(Adaptor.BSpline().GetObject())
    
    
    #Order BSplineLst Head-to-Tail
    rel_tol=1e-07
    NbSplines = len(BSplineLst)
    for i, item in enumerate(BSplineLst[:NbSplines-1]):
        EP1 = item.EndPoint()
        SP2 = BSplineLst[i+1].StartPoint()
        if not EP1.IsEqual(SP2,rel_tol):
            BSplineLst[i+1].Reverse()
        else:
            None
    
    #SET_Origin
    BSplineLst = set_BSplineLst_to_Origin(BSplineLst)

    return BSplineLst 
    

def stp3d_to_2dBSplineLst(TopoDS_Shape, R):
    #Define the Point and direction of the slicing plane
    Dir = gp_Dir(1., 0., 0.) 
    Pnt = gp_Pnt(R,0,0)
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
        display.DisplayShape(edge)
        ex.Next()
   
    BSplineLst = []
    for edg in results:
        #print "null now?", bool(edg.IsNull())
        Adaptor = BRepAdaptor_Curve2d(edg,TopoDS_Face)
        if Adaptor.GetType() == 0: # gp_Lin2d 	
            P1 = Adaptor.Value(Adaptor.FirstParameter())
            P2 = Adaptor.Value(Adaptor.LastParameter())
            BSpline = Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([P1,P2])).Curve().GetObject()   
            BSplineLst.append(BSpline)
            
        elif Adaptor.GetType() == 6: #BSplineCurve
            BSplineLst.append(Adaptor.BSpline().GetObject())
    

    return BSplineLst    
    
    
    
# =============================================================================
#            IF __MAIN__
# =============================================================================  
if __name__ == '__main__': 
    display, start_display, add_menu, add_function_to_menu = init_display('wx')
    display.Context.SetDeviationAngle(0.0001)       # 0.001 default. Be careful to scale it to the problem.
    display.Context.SetDeviationCoefficient(0.0001) # 0.001 default. Be careful to scale it to the problem. 
    display.set_bg_gradient_color(20,6,111,200,200,200)
    
    '''CREATE AXIS SYSTEM for Visualization'''
    O  = gp_Pnt(0., 0., 0.)
    p1 = gp_Pnt(10,0.,0.)
    p2 = gp_Pnt(0.,10.0,0.)
    p3 = gp_Pnt(0.,0.,10.0)
    
    h1 = BRepBuilderAPI_MakeEdge(O,p1).Shape()
    h2 = BRepBuilderAPI_MakeEdge(O,p2).Shape()
    h3 = BRepBuilderAPI_MakeEdge(O,p3).Shape()

    display.DisplayShape(O,color='BLACK')
    display.DisplayShape(h1,color='RED')
    display.DisplayShape(h2,color='GREEN')
    display.DisplayShape(h3,color='BLUE')
       
    
    filename = 'AREA_Blatt_L.stp'
    aResShape = load_3D(filename) 
    display.DisplayShape(aResShape, color=None, transparency=0.7, update=True)
    
    #=====================
    #IMPORT 2D Step File
    #=====================
#    BSplineLst = stp2d_to_2dBSplineLst(aResShape)        
#    BSplineLst = reverse_BSplineLst(BSplineLst)
#   
#    P = get_BSplineLst_Pnt2d(BSplineLst,0.0, 0, 1)
#    display.DisplayShape(P,color='BLACK') 
#    P = get_BSplineLst_Pnt2d(BSplineLst,0.005, 0, 1)
#    display.DisplayShape(P,color='WHITE')    
#
#    transformedwire = build_wire_from_BSplineLst(BSplineLst)
#    display.DisplayShape(transformedwire,color='BLUE')  
 

    #=====================
    #IMPORT 3D Step File and Slice it a Radial Station
    #=====================
   
    R = 1223
    BSplineLst = stp3d_to_2dBSplineLst(aResShape,R)

        
    #for i,item in enumerate(BSplineLst):
        #display.DisplayShape(item)
    
    #Transform BSplineLst to SONATA Coordinates
    wire = build_wire_from_BSplineLst(BSplineLst)
    rotwire1 = rotate_wire(wire,gp_Ax1(gp_Pnt(0,0,0),gp_Dir(0,0,1)),math.pi/2)
    rotwire2 = rotate_wire(rotwire1,gp_Ax1(gp_Pnt(0,0,0),gp_Dir(0,1,0)),math.pi)
    print "WIRE CLOSED: " + str(rotwire2.Closed())
    display.DisplayShape(rotwire2)
    
    Dir = gp_Dir(1., 0., 0.) 
    Pnt = gp_Pnt(R,0,0)
    Pln = gp_Pln(Pnt, Dir)
    face = BRepBuilderAPI_MakeFace(Pln,).Shape()
    display_face = BRepBuilderAPI_MakeFace(Pln,- 100., 100., -120., 120).Shape()
    display.DisplayShape(display_face, color="BLACK", transparency=0.8, update=True)


#TODO: Wire_to_BSplineLst
#TODO: def import_3d_stp(filename,R)
#TODO: def import_2d_stp(filename)


#
#    Pln = gp_Pln(P, D)
#    face = BRepBuilderAPI_MakeFace(Pln,).Shape()
#    display_face = BRepBuilderAPI_MakeFace(Pln,- 100., 100., -120., 120).Shape()
#    TopoDS_Face = TopoDS.topods().Face(face)
#    
#    # Computes Shape/Plane intersection
#    section = BRepAlgoAPI_Section(aResShape, face, False)       
#    section.Approximation(True)
#    section.Build()
#    
#    ex = TopExp_Explorer(section.Shape(),6,7) #Search for Edges(6), Exclude Vertices(7)
#    results = []
#    NbEdges = 0 
#    while ex.More():
#        edge = TopoDS.topods().Edge(ex.Current())
#        results.append(edge)
#        NbEdges += 1
#        ex.Next()
#    
#    BSplineLst = []
#    for edg in results:
#        Adaptor = BRepAdaptor_Curve(edg,TopoDS_Face)
#        BSplineLst.append(Adaptor.BSpline().GetObject())
#       
#    TopTools_EdgeLst = TopTools_ListOfShape()
#    tmp_wire = 	BRepBuilderAPI_MakeWire()
#    for i,item in enumerate(BSplineLst):
#        Plane = Geom_Plane(P, D)
#        tmp_edge = BRepBuilderAPI_MakeEdge(item.GetHandle())
#        #print tmp_edge.IsDone()
#        #display.DisplayShape(tmp_edge.Edge())
#        TopTools_EdgeLst.Append(tmp_edge.Edge())
#        
#    wire = BRepBuilderAPI_MakeWire()
#    wire.Add(TopTools_EdgeLst)
#    wire.Build()
#    display.DisplayShape(wire.Wire())
#    #wire = wire.Wire()
#    
#    #print "WIRE CLOSED: " + str(wire.Closed())
    #transformedwire = build_wire_from_BSplineLst(BSplineLst)         

# =============================================================================
#            DISPLAY
# =============================================================================

    f = display.View.View().GetObject()
    
    #display.DisplayShape(gp_Pnt(500,0,0))
    #display.DisplayShape(wire)
    #display.DisplayShape(display_face, color="BLUE", transparency=0.8, update=True)
    #display.DisplayShape(transformedwire)
    #display.DisplayShape(section.Shape(),color="BLACK")
    #display.DisplayShape(section2.Shape(),color="GREEN")
    #display.register_select_callback(print_xy_click)
    display.set_bg_gradient_color(20,6,111,200,200,200)
    add_menu('screencapture')
    add_function_to_menu('screencapture', export_to_PDF)
    add_function_to_menu('screencapture', export_to_SVG)
    add_function_to_menu('screencapture', export_to_PS)
    add_function_to_menu('screencapture', export_to_EnhPS)
    add_function_to_menu('screencapture', export_to_TEX)
    add_function_to_menu('screencapture', export_to_BMP)
    add_function_to_menu('screencapture', export_to_PNG)
    add_function_to_menu('screencapture', export_to_JPEG)
    add_function_to_menu('screencapture', export_to_TIFF)
    
    display.View_Top()
    #display.View_Iso()
    display.FitAll()
    start_display()