import sys
import os
import math
import numpy as np
import matplotlib.pyplot as plt

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
from OCC.Quantity import Quantity_Color
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
from OCC.GCPnts import GCPnts_AbscissaPoint, GCPnts_QuasiUniformDeflection, GCPnts_TangentialDeflection, GCPnts_QuasiUniformAbscissa
from OCC.Geom2dAPI import Geom2dAPI_Interpolate
from OCC.BRepTools import BRepTools_WireExplorer
from OCC.TopoDS import topods
from OCC.ShapeAnalysis import ShapeAnalysis_Wire, ShapeAnalysis_WireOrder
from OCC.ShapeFix import ShapeFix_Wire

from wire_utils import build_wire_from_BSplineLst,rotate_wire, translate_wire, NbEdges_in_wire, Unique_EdgeLst, Wire_Orientation
from explorer import WireExplorer
from utils import point2d_list_to_TColgp_HArray1OfPnt2d, Pnt2dLst_to_npArray, PolygonArea
from BSplineLst_utils import get_BSpline_length, get_BSplineLst_length, \
                            find_BSplineLst_coordinate, get_BSplineLst_Pnt2d, \
                            trim_BSplineLst, seg_boundary_from_dct, set_BSplineLst_to_Origin, \
                            copy_BSplineLst, trim_BSplineLst_by_Pnt2d, copy_BSpline, \
                            reverse_BSplineLst, BSplineLst_Orientation

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
    return aResShape

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



def EdgeLst_to_Wire(EdgeLst):  
    '''FROM UNORDERED LIST OF EDGES TO CONNECTING WIRE'''
    
     # filter out those entities that share the same TShape
     # but do *not* share the same orientation
#    filter_orientation_seq = []
#    for i in EdgeLst:
#        _present = False
#        for j in filter_orientation_seq:
#            if i.IsSame(j):
#                _present = True
#                break
#        if _present is False:
#            filter_orientation_seq.append(i)
#    EdgeLst = filter_orientation_seq
    
 
    pool = EdgeLst
    wire = BRepBuilderAPI_MakeWire()
    idx_pool = []
    for i in range(0,len(pool)):
        idx_pool.append(i)
        
    remaining_idx_pool = list(idx_pool)
    order_sheme = [0]
    
    wire = BRepBuilderAPI_MakeWire()
    tmp_wire = BRepBuilderAPI_MakeWire()
    
    while len(order_sheme) != len(pool):
        #Build Wire according to order_sheme
        for i,item in enumerate(order_sheme):
            wire.Add(pool[item])
         
        #Find the next element that fits to the temporary wire
        remaining_idx_pool = list(idx_pool)
        for idx in order_sheme:
            remaining_idx_pool.remove(idx) 
        
        tmp_wire = wire
        for idx in remaining_idx_pool: #BRepBuilderAPI_DisconnectedWire: The last edge which you attempted to add was not connected to the wire.
            tmp_wire.Add(pool[idx])
            if tmp_wire.Error() == 2:
                None
            else:
                order_sheme.append(idx)
                break
            tmp_wire = wire    
 
    final_wire = BRepBuilderAPI_MakeWire()
    #Build final_Wire according to order_sheme
    for i,item in enumerate(order_sheme):
            final_wire.Add(pool[item])

    #wire = BRepBuilderAPI_MakeWire()
    #wire.Add(TopTools_EdgeLst)
    final_wire.Build()
    
    #print NbEdges_in_wire(final_wire.Wire())
    
    if final_wire.Error() == 0:
        None
        #print 'No error occurred. The wire is correctly built.'
        
        
    return final_wire.Wire()



def wire_to_BSplineLst(wire):
    BSplineLst = []
    #Define the Point and direction of the slicing plane
    Dir = gp_Dir(0., 0., 1.) 
    Pnt = gp_Pnt(0,0,0)
    Pln = gp_Pln(Pnt, Dir)
    face = BRepBuilderAPI_MakeFace(Pln).Shape()
    TopoDS_Face = TopoDS.topods().Face(face)
    
    ex = BRepTools_WireExplorer(wire)
    while ex.More():
        edg = ex.Current()
        #print 
        ex.Next()
    
        #print "null now?", bool(edg.IsNull())
        Adaptor = BRepAdaptor_Curve2d(edg,TopoDS_Face)
        
        #print Adaptor.GetType()
        if Adaptor.GetType() == 0: # gp_Lin2d 	
            P1 = Adaptor.Value(Adaptor.FirstParameter())
            P2 = Adaptor.Value(Adaptor.LastParameter())
            BSpline = Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([P1,P2])).Curve().GetObject()   
            BSplineLst.append(BSpline)
        
        if Adaptor.GetType() == 1:
            NbPoints = 20        
            #print Adaptor.FirstParameter(),Adaptor.LastParameter()
            discretization = GCPnts_QuasiUniformAbscissa(Adaptor, NbPoints)  
            
            Pnt2dLst = []
            NbPoints = discretization.NbPoints()
            for j in range(1, NbPoints+1):
               para = discretization.Parameter(j)
               Pnt2dLst.append(Adaptor.Value(para))
               
            array = point2d_list_to_TColgp_Array1OfPnt2d(Pnt2dLst)                                               
            #tmp_interpolation = Geom2dAPI_Interpolate(harray.GetHandle(), False, 1e-04) 
            #tmp_interpolation.Perform()
            tmp_approximate = Geom2dAPI_PointsToBSpline(array)
            BSplineLst.append(tmp_approximate.Curve().GetObject())
            
        elif Adaptor.GetType() == 6: #BSplineCurve
            BSplineLst.append(Adaptor.BSpline().GetObject())
            
    return BSplineLst    


def order_BSplineLst_Head2Tail(BSplineLst,rel_tol=1e-07):    
    #Order BSplineLst Head-to-Tail
    NbSplines = len(BSplineLst)
    for i, item in enumerate(BSplineLst[:NbSplines-1]):
        EP1 = item.EndPoint()
        SP2 = BSplineLst[i+1].StartPoint()
        if not EP1.IsEqual(SP2,rel_tol):
            BSplineLst[i+1].Reverse()
        else:
            None
            
    return BSplineLst    



def stp2d_to_wire(TopoDS_Shape):
    ex = TopExp_Explorer(TopoDS_Shape,6,7) #Search for Edges(6), Exclude Vertices(7)
    EdgeLst = []
    counter = 0
    while ex.More():
        edge = TopoDS.topods().Edge(ex.Current())
        EdgeLst.append(edge)
        counter += 1
        #print "is null?", bool(edge.IsNull())
        #display.DisplayShape(edge)
        ex.Next()
    #print counter
    #NonUniques = Unique_EdgeLst(EdgeLst)
    wire = EdgeLst_to_Wire(EdgeLst)    
       
    #print wire.Closed()
    return wire
    
    

def stp3d_to_wire(TopoDS_Shape, R):
    #Define the Point and direction of the slicing plane
    Dir = gp_Dir(1., 0., 0.) 
    Pnt = gp_Pnt(R,0,0)
    Pln = gp_Pln(Pnt, Dir)
    face = BRepBuilderAPI_MakeFace(Pln).Shape()

    section = BRepAlgoAPI_Section(TopoDS_Shape, face)   # Computes Shape/Plane intersection

    ex = TopExp_Explorer(section.Shape(),6,7) #Search for Edges(6), Exclude Vertices(7)
    
    EdgeLst = []
    counter = 0
    #TopTools_EdgeLst = TopTools_ListOfShape()
    while ex.More():
        counter += 1
        edge = TopoDS.topods().Edge(ex.Current())
        EdgeLst.append(edge)
        #print "is null?", bool(edge.IsNull())
        #display.DisplayShape(edge)
        ex.Next()
        
    Wire = EdgeLst_to_Wire(EdgeLst)
    
    #DISPLAY SECTIONING PLANE
    #display_face = BRepBuilderAPI_MakeFace(Pln,- 100., 100., -120., 120).Shape()
    #display.DisplayShape(display_face, color="BLACK", transparency=0.8, update=True)
    return Wire
   
    

def import_2d_stp(filename):
    '''
    The 2D Shape must given in SONATA Coordinates!
    '''
    aResShape = load_3D(filename) 
    wire = stp2d_to_wire(aResShape)
    BSplineLst = wire_to_BSplineLst(wire)

    BSplineLst = order_BSplineLst_Head2Tail(BSplineLst)
    
    if BSplineLst_Orientation(BSplineLst) == False:
        BSplineLst = reverse_BSplineLst(BSplineLst)  
    
    #TODO: Rotate x-Axis to determine Origin
    BSplineLst = set_BSplineLst_to_Origin(BSplineLst) 
    BSplineLst = order_BSplineLst_Head2Tail(BSplineLst)

    print 'STATUS: \t IMPORT_3d_STP'
    print 'STATUS: \t CHECK ClosedWire: \t\t ', str(wire.Closed())
    print 'STATUS: \t CHECK Head2Tail: \t\t ', Check_BSplineLst_Head2Tail(BSplineLst)
    print 'STATUS: \t CHECK Counterclockwise: \t ', BSplineLst_Orientation(BSplineLst)
    return BSplineLst



def import_3d_stp(filename,R):
    '''
    The 3D Shape must given in Rotorblade Coordinates!
    '''
    print 'STATUS: \t IMPORT_3d_STP'
    aResShape = load_3D(filename) 
    wire = stp3d_to_wire(aResShape,R)   
    print 'STATUS: \t CHECK ClosedWire: \t\t ', str(wire.Closed())
    
    print 'STATUS: \t Transform wire to SONATA Crosssection Coordinates'
    #Transform wire to SONATA Crosssection Coordinates
    wire = translate_wire(wire,gp_Pnt(R,0,0),gp_Pnt(0,0,0))
    wire = rotate_wire(wire,gp_Ax1(gp_Pnt(0,0,0),gp_Dir(0,1,0)),math.pi/2)
    wire = rotate_wire(wire,gp_Ax1(gp_Pnt(0,0,0),gp_Dir(0,0,1)),math.pi/2)
    
    print 'STATUS: \t Transform wire to SONATA BSplineLst'
    BSplineLst = wire_to_BSplineLst(wire)
    
    print 'STATUS: \t Order BSplineLst Head2Tail'
    BSplineLst = order_BSplineLst_Head2Tail(BSplineLst)
    
    print 'STATUS: \t Check Orientation and Reverse it if Necessary'
    if BSplineLst_Orientation(BSplineLst,11) == False:                      #TODO: VERY SLOW CHANGE THE Get_Pnt2d from BSplineLst function to something else!!!!
        print 'Reversing BSplineLst'
        BSplineLst = reverse_BSplineLst(BSplineLst)  
    
    #TODO: Rotate x-Axis to determine Origin
    print 'STATUS: \t Set BSplineLst to Origin'
    BSplineLst = set_BSplineLst_to_Origin(BSplineLst) 
    BSplineLst = order_BSplineLst_Head2Tail(BSplineLst)
    

    print 'STATUS: \t CHECK Head2Tail: \t\t ', Check_BSplineLst_Head2Tail(BSplineLst)
    print 'STATUS: \t CHECK Counterclockwise: \t ', BSplineLst_Orientation(BSplineLst,11)
    
    return BSplineLst


def Check_BSplineLst_Head2Tail(BSplineLst):
    Head2Tail = True
    lin_tol=1e-07
    NbSplines = len(BSplineLst)
    for i, item in enumerate(BSplineLst[:NbSplines-1]):
        EP1 = item.EndPoint()
        SP2 = BSplineLst[i+1].StartPoint()
        if not EP1.IsEqual(SP2,lin_tol):
           #display.DisplayShape(item,color='RED')
           Head2Tail = False
        else:
            None
            #display.DisplayShape(item,color='BLUE')
    
    
    EP1 = BSplineLst[0].StartPoint()
    SP2 = BSplineLst[-1].EndPoint()
    if not EP1.IsEqual(SP2,lin_tol):
        #display.DisplayShape(BSplineLst[-1],color='RED')
        Head2Tail = False
    else:
        None
        #display.DisplayShape(BSplineLst[-1],color='BLUE')        
    
    return Head2Tail

    


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

    #=====================
    #IMPORT 2D Step File
    #=====================
    #BSplineLst = import_2d_stp('AREA_R230.stp')
    #BSplineLst = set_BSplineLst_to_Origin(BSplineLst)
    #for i,item in enumerate(BSplineLst):
       #display.DisplayShape(item) 

    #aResShape = load_3D('AREA_R230.stp') 
    #display.DisplayShape(aResShape)
    #transformedwire = build_wire_from_BSplineLst(BSplineLst)
    #display.DisplayShape(transformedwire,color='BLUE')  

    #=====================
    #IMPORT 3D Step File and Slice it a Radial Station
    #=====================
    BSplineLst = import_3d_stp('AREA_Blatt_L.stp',230)
    
#    for i,item in enumerate(BSplineLst):
#        [R,G,B,T] =  plt.cm.jet(i*10)
#        display.DisplayColoredShape(item, Quantity_Color(R, G, B, 0),update=True)
    
    
    #aResShape = load_3D('AREA_Blatt_L.stp') 
    #display.DisplayShape(aResShape, color=None, transparency=0.7, update=True)
    Wire = build_wire_from_BSplineLst(BSplineLst)
    display.DisplayShape(Wire,color="BLUE")


# DISPLAY
# =============================================================================

    f = display.View.View().GetObject()
    
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