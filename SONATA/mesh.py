import os
import numpy as np
import pickle
import matplotlib as plt

from OCC.AIS import AIS_Shape
from OCC.BRep import BRep_Builder, BRep_Tool
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeFace
from OCC.BRepAdaptor import BRepAdaptor_CompCurve
from OCC.gp import gp_Pnt2d, gp_Pnt
from OCC.GCPnts import GCPnts_QuasiUniformAbscissa
from OCC.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.Geom2dAPI import Geom2dAPI_ProjectPointOnCurve
from OCC.Display.SimpleGui import init_display
from OCC.Quantity import Quantity_Color
from OCC.Graphic3d import Graphic3d_EF_PDF, Graphic3d_EF_SVG, Graphic3d_EF_TEX, Graphic3d_EF_PostScript, Graphic3d_EF_EnhPostScript
from OCC.BRepLib import breplib_BuildCurves3d

from OCC.TopoDS import TopoDS_Compound, topods_Face, topods_Edge


from BSplineLst_utils import get_BSplineLst_length, get_BSpline_length, trim_BSplineLst
from wire_utils import build_wire_from_BSplineLst,get_wire_length
from utils import Pnt2dLst_to_npArray, unique_rows
from display_mesh import plot_mesh
#====================NOTES:==========================================================
'''
1) Distribute Equidistant Points on Layer_i
2) Find Cell where the Max.Error between Cell and the Layer is larger than allowed -> Split Cell at max deviation
3) Determine Area of Cells and Cell Angles -> Improve Cells 





'''
#===================LOAD CROSSSECTION==========================================
#LOAD .pkl data with SegmentLst
filename = 'sec_config.pkl'
with open(filename, 'rb') as handle:
    SegmentLst = pickle.load(handle)
    

#Build wires for each layer and segment
for seg in SegmentLst:
    seg.build_wire()
    for layer in seg.LayerLst:
        layer.build_wire()



class cell(object):
     
    def __init__(self):
        self.points
        self.faces
        self.neighbours
        self.Orientation
        self.MatID
        
    
class point(object):
     
    def __init__(self):
        self.coords
        self.cell_pointer
        self.face_pointer



def equidistant_Points_on_BSplineLst(BSplineLst, IC=False, **kwargs): #include corners
    ''' minLen = the minimum distance between points. Should be determined by the segment0.BSplineLstLenght/Resolution
        IC (Include Corners): True or False, Note: IC=True the the Points are not equidistant placed on BSplineLst!!!
        
        (BSplineLst, IC=False, either NbPoints or MinLen)
        (BSplineLst, IC=True,  MinLen)
        
    ''' 
    if IC==True:
        
        #KWARGS:
        if kwargs.get('minLen') !=  None:
            minLen = kwargs.get('minLen')
        
        Pnt2dLst = []
        for item in BSplineLst:
            Adaptor = Geom2dAdaptor_Curve(item.GetHandle())
            length = get_BSpline_length(item)
            NbPoints = int(length//minLen)+2  
            discretization = GCPnts_QuasiUniformAbscissa(Adaptor,NbPoints)
            
            for j in range(1, NbPoints):
                    para = discretization.Parameter(j)
                    Pnt = gp_Pnt2d()
                    item.D0(para,Pnt)
                    Pnt2dLst.append(Pnt)
            
        #add last point 
        para = discretization.Parameter(j+1)
        Pnt = gp_Pnt2d()
        item.D0(para,Pnt)
        Pnt2dLst.append(Pnt)
            

    else:
        wire = build_wire_from_BSplineLst(BSplineLst)
        wire_length = get_wire_length(wire)
        
        #KWARGS:
        if  kwargs.get('NbPoints') !=  None:
            NbPoints = kwargs.get('NbPoints')
            
        elif kwargs.get('minLen') !=  None and kwargs.get('NbPoints') == None:
            NbPoints = int(wire_length//kwargs.get('minLen'))+2  
                          
        AdaptorComp = BRepAdaptor_CompCurve(wire, True)
        discretization = GCPnts_QuasiUniformAbscissa(AdaptorComp,NbPoints)
        
        Pnt2dLst = []
        for j in range(1, NbPoints+1):
            para = discretization.Parameter(j)
            P = AdaptorComp.Value(para)
            Pnt2dLst.append(gp_Pnt2d(P.X(),P.Y()))
       
    return Pnt2dLst


#Point Projection on BSPLINELST 
#Pnt2d = a_Pnt2dLst[23]
#proPnt2dLst = []
#BSplineLst = Trimmed_BSplineLst 
#tol = 1e-3
#distance = (1+tol)*layer.thickness
##project_Pnt2d_on_BSplineLst(BSplineLst,Pnt2d,distance)
#for item in BSplineLst:
#    projection = Geom2dAPI_ProjectPointOnCurve(Pnt2d,item.GetHandle())
#    print projection.NbPoints()
#    for i in range(1,projection.NbPoints()+1):
#        print 'i:', i#
#        print projection.Distance(i)
#        if projection.Distance(i)<=distance:
#            proPnt2dLst.append(projection.Point(i))



#==============================================================================
Projection = SegmentLst[-1].Projection
Resolution = 300 # Nb of Points on Segment0
length = get_BSplineLst_length(SegmentLst[0].BSplineLst)
global_minLen = round(length/Resolution,5)
NbPoints = int(length//global_minLen)+2 
layer = SegmentLst[-1].LayerLst[-1]
BSplineLst = layer.BSplineLst
Boundary_BSplineLst = layer.Boundary_BSplineLst
S1 = layer.S1
S2 = layer.S2
Trimmed_BSplineLst = trim_BSplineLst(Boundary_BSplineLst, S1, S2, 0, 1)




a_Pnt2dLst = equidistant_Points_on_BSplineLst(BSplineLst, False, minLen = global_minLen)
a = Pnt2dLst_to_npArray(a_Pnt2dLst)
n = len(a_Pnt2dLst )

b_Pnt2dLst = equidistant_Points_on_BSplineLst(Trimmed_BSplineLst, False, NbPoints = n-2)
b = Pnt2dLst_to_npArray(b_Pnt2dLst)
Pnts = a_Pnt2dLst+b_Pnt2dLst

#CREATE mesh
nodes = np.vstack([a,b])
s_idx = len(a)
len_b = len(b)
elements = []
for i in range(0, len(b)-1):
    elements.append([i+3,i+2,s_idx+(i+1),s_idx+(i+2)])
data = np.zeros(len(elements))

plot_mesh(nodes,elements,data,True,False)



#DisplayMesh(nodex,elements)
builder = BRep_Builder()
comp = TopoDS_Compound()
WireBuilder = BRepBuilderAPI_MakeWire()
builder.MakeCompound(comp)
WireLst = []

for ele in elements:
    WireBuilder = BRepBuilderAPI_MakeWire()
    for i in range(0,len(ele)-1):
        node1 = ele[i]-1
        node2 = ele[i+1]-1
        Pnt1 = gp_Pnt(nodes[node1][0],nodes[node1][1],0)
        Pnt2 = gp_Pnt(nodes[node2][0],nodes[node2][1],0)
        me = BRepBuilderAPI_MakeEdge(Pnt1, Pnt2)
        if me.IsDone():
            WireBuilder.Add(me.Edge())
    
    node1 = ele[-1]-1
    node2 = ele[0]-1         
    Pnt1 = gp_Pnt(nodes[node1][0],nodes[node1][1],0)
    Pnt2 = gp_Pnt(nodes[node2][0],nodes[node2][1],0)
    me = BRepBuilderAPI_MakeEdge(Pnt1, Pnt2)
    if me.IsDone():
        WireBuilder.Add(me.Edge())         
        
    WireLst.append(WireBuilder.Wire())

FaceLst = [] 
for wire in WireLst:
    Face = BRepBuilderAPI_MakeFace(wire)
    FaceLst.append(Face.Face())






# HERE WILL BE THE MESHING ALGORITHM!



#====================DISPLAY===================================================
display, start_display, add_menu, add_function_to_menu = init_display('wx')
display.Context.SetDeviationAngle(0.000001)       # 0.001 default. Be careful to scale it to the problem.
display.Context.SetDeviationCoefficient(0.000001) # 0.001 default. Be careful to scale it to the problem. 
display.set_bg_gradient_color(20,6,111,200,200,200) 
    


for p in a_Pnt2dLst:
    display.DisplayShape(p,color='BLACK')

for p in b_Pnt2dLst:
    display.DisplayShape(p,color='GREEN')

#for p in proPnt2dLst:
#    display.DisplayShape(p,color='RED')

#for f in FaceLst:
#    breplib_BuildCurves3d(f)
#    display.DisplayShape(f)



# transfer shapes and display them in the viewer
display.DisplayShape(SegmentLst[0].wire, color="BLACK")
display.DisplayShape(SegmentLst[0].BSplineLst[0].StartPoint())

for i,seg in enumerate(SegmentLst):
    display.DisplayShape(seg.wire,color="BLACK")
    k = 0
    for j,layer in enumerate(seg.LayerLst):
        [R,G,B,T] =  plt.cm.jet(k*50)
        
        if i==0:
            display.DisplayColoredShape(layer.wire, Quantity_Color(R, G, B, 0),update=True)
            
        elif i==1:
            display.DisplayColoredShape(layer.wire, Quantity_Color(R, G, B, 0),update=True)

        else:
            display.DisplayColoredShape(layer.wire, Quantity_Color(R, G, B, 0),update=True)

        k = k+1;
        if k>5:
            k = 0








def display_custome_shape(shape,linewidth,transparency,RGB):
    s = shape
    ais_shp = AIS_Shape(s)
    ais_shp.SetWidth(linewidth)
    ais_shp.SetTransparency(transparency)
    ais_shp.SetColor(Quantity_Color(RGB[0], RGB[1], RGB[2], 0))
    ais_context = display.GetContext().GetObject()
    ais_context.Display(ais_shp.GetHandle())
    return None

for w in reversed(WireLst):
    display_custome_shape(w,1.6,0,[1,0,0])

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
    f.Export('capture_svg_%s.svg' % i, Graphic3d_EF_SVG)
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


'''CREATE AXIS SYSTEM for Visualization'''
O  = gp_Pnt(0., 0., 0.)
p1 = gp_Pnt(10.0,0.,0.)
p2 = gp_Pnt(0.,10.0,0.)
p3 = gp_Pnt(0.,0.,10.0)

h1 = BRepBuilderAPI_MakeEdge(O,p1).Shape()
h2 = BRepBuilderAPI_MakeEdge(O,p2).Shape()
h3 = BRepBuilderAPI_MakeEdge(O,p3).Shape()

display.DisplayShape(O,color='BLACK')
display.DisplayShape(h1,color='RED')
display.DisplayShape(h2,color='GREEN')
display.DisplayShape(h3,color='BLUE')
    
    
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
display.FitAll()
start_display()

