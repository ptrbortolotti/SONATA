#-------------------------------
#          H E A D E R
#-------------------------------

#Basic Libraries
import math
import numpy as np                                                              #fundamental package for scientific computing with Python (powerful N-Dimensional array objects)
import matplotlib.pyplot as plt                                                 #python 2D plotting library
from scipy.interpolate import interp1d
                                                                 # help in opening URL

#Python OCC Libraries
from OCC.gp import gp_Pnt, gp_Vec, gp_Pnt2d, gp_Pln, gp_Dir, gp_Trsf, gp_Ax1, gp_OX, gp_Ax3, gp_Ax2, gp_Circ, gp_OY
from OCC.Geom2dAPI import Geom2dAPI_PointsToBSpline, Geom2dAPI_Interpolate, Geom2dAPI_InterCurveCurve 
from OCC.GeomAPI import geomapi, GeomAPI_PointsToBSpline,GeomAPI_Interpolate
from OCC.TColgp import TColgp_Array1OfPnt, TColgp_Array1OfPnt2d, TColgp_HArray1OfPnt2d
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge,BRepBuilderAPI_MakeWire 
from OCC.BRepOffsetAPI import BRepOffsetAPI_MakeOffset 
from OCC.Geom import Geom_OffsetCurve, Geom_BezierCurve, Geom_Plane, Geom_TrimmedCurve, Geom_Curve 
from OCC.GeomAPI import GeomAPI_IntCS
from OCC.GC import GC_MakeArcOfCircle, GC_MakeSegment
from OCC.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCC.Display.SimpleGui import init_display
from OCC.TopoDS import topods, TopoDS_Edge, TopoDS_Compound, TopoDS_Vertex
from OCC.TopExp import TopExp_Explorer
from OCC.Display.WebGl import x3dom_renderer,threejs_renderer
from OCC.STEPControl import STEPControl_Writer, STEPControl_AsIs                #FOR STEP EXPORT
from OCC.Interface import Interface_Static_SetCVal
from OCC.IFSelect import IFSelect_RetDone


#Frequently Used Functions 
from core_geometry_utils import *
from core_operations_utils import *
from tab_definition import *
from Intersect_and_trim_utils import trim_2dcurve_selfintersectingLoop


#-------------------------------
#          FUNCTIONS
#-------------------------------

def display_points_of_array(array):
    for j in range(array.Lower(),array.Upper()+1):
        p = array.Value(j)
        display.DisplayShape(p, update=False)

#-------------------------------
#          M A I N
#-------------------------------
R = 1650                    #[mm]
R_0 = 120                   #[mm]
R_r = 150                   #[mm]
R_h = 250                   #[mm]

r_0 = R_0/float(R)          #[-]
r_r = R_r/float(R)          #[-]
r_h = R_h/float(R)          #[-]

airfoil = 'naca2410'  

###############################################################################
display, start_display, add_menu, add_function_to_menu = init_display(backend_str=None, size=(1920,1080))
# required for precise rendering of the circles
display.Context.SetDeviationAngle(0.00001)      # 0.001 default
display.Context.SetDeviationCoefficient(0.00001) # 0.001 default


#CREATE AXIS SYSTEM for Visualization
COSY = gp_Ax3()	
O  = gp_Pnt(0., 0., 0.)
p1 = gp_Pnt(1.,0.,0.)
p2 = gp_Pnt(0.,1.,0.)
p3 = gp_Pnt(0.,0.,1.)

h1 = BRepBuilderAPI_MakeEdge(O,p1).Shape()
h2 = BRepBuilderAPI_MakeEdge(O,p2).Shape()
h3 = BRepBuilderAPI_MakeEdge(O,p3).Shape()

display.DisplayShape(O,color='BLACK')
display.DisplayShape(h1,color='RED')
display.DisplayShape(h2,color='GREEN')
display.DisplayShape(h3,color='BLUE')

###############################################################################

#GET AIRFOIL POINTS AND PLACE THEM INTO AN NP.ARRAY called data   
data = UIUCAirfoil(airfoil)
#X = -(data[1]-0.25)             #Move Center of Airfoil to 1/4-chord with h2 pointing to the leading edge!
#data[1] = X

harray = TColgp_HArray1OfPnt_from_nparray(data)                    #Create TColgp_HArray1OfPnt from np.array
display_points_of_array(harray)

anInterpolation = GeomAPI_Interpolate(harray.GetHandle(), False, 0.001)
anInterpolation.Perform()
bspline_5 = anInterpolation.Curve()
edge_5 = BRepBuilderAPI_MakeEdge(bspline_5)
display.DisplayShape(edge_5.Shape(), update=True, color='WHITE')
aWire = BRepBuilderAPI_MakeWire(edge_5.Edge())


###################################
#Rotor Blade Skin
###################################

#dist = -0.003
#num_of_layers = 10
#
#geom_bspline = bspline_5
#for j in range(0,num_of_layers):
#    geom_bspline= Geom_OffsetCurve(geom_bspline, dist, gp_Dir(1., 0., 0.))
#    result = geom_bspline.IsCN(2)
#    display.DisplayShape(geom_bspline, color='BLUE', update=True)
#    geom_bspline = geom_bspline.GetHandle()


###################################
#C-SPAR Bezier Curve 
###################################

#P1 Definition:
x = 0.4
Plane = Geom_Plane(gp_Ax3(gp_Pnt(0, x, 0),gp_Dir(0,1.0,0)))
IntCS = GeomAPI_IntCS (bspline_5, Plane.GetHandle())
P1 = IntCS.Point(1)
display.DisplayColoredShape(P1, 'BLUE')

#P4 Definition:
x = 0.45
Plane = Geom_Plane(gp_Ax3(gp_Pnt(0, x, 0),gp_Dir(0,1.0,0)))
IntCS = GeomAPI_IntCS (bspline_5, Plane.GetHandle())
P4 = IntCS.Point(2)
display.DisplayColoredShape(P4, 'RED')

#P2 Definition:
alpha2 = 182*np.pi/180           #in rad
L2 = 0.4 
y = np.cos(alpha2)*L2  
z = np.sin(alpha2)*L2
P2_coord = gp_XYZ(0,y,z)
P2_coord.Add(P1.Coord())
P2 = gp_Pnt(P2_coord)
display.DisplayColoredShape(P2, 'YELLOW')

#P3 Definition:
alpha3 = 178*np.pi/180           #in rad
L3 = 0.5 
y = np.cos(alpha3)*L3  
z = np.sin(alpha3)*L3
P3_coord = gp_XYZ(0,y,z)
P3_coord.Add(P4.Coord())
P3 = gp_Pnt(P3_coord)
display.DisplayColoredShape(P3, 'GREEN')


array = TColgp_Array1OfPnt(1, 4)
array.SetValue(1, P1)
array.SetValue(2, P2)
array.SetValue(3, P3)
array.SetValue(4, P4)

curve = Geom_BezierCurve(array)
ME = BRepBuilderAPI_MakeEdge(curve.GetHandle())
GreenEdge = ME
display.DisplayColoredShape(GreenEdge.Edge(), 'BLACK')


###################################
#Rotor Blade Skin
###################################

dist = -0.004
num_of_layers = 8

offset_curves = []
geom_bspline = bspline_5
for j in range(0,num_of_layers):
    geom_bspline= Geom_OffsetCurve(geom_bspline, dist, gp_Dir(1., 0., 0.))
    result = geom_bspline.IsCN(2)
    display.DisplayShape(geom_bspline, color='BLUE', update=True)
    offset_curves.append(geom_bspline)     
    geom_bspline = geom_bspline.GetHandle()



Offset1 = offset_curves[0]
Inter = Geom2dAPI_InterCurveCurve(Offset1.GetHandle(), 1.0e-5)
print('Number of Intersection Points:', Inter.NbPoints())
#print('Number of tangetial Intersections:',Inter.NbSegments())
InterP1 = Inter.Point(1)
TrimmedWire = trim_2dcurve_selfintersectingLoop(InterP1)

display.DisplayShape(TrimmedWire.Wire(), color='BLACK')
#display.DisplayShape(Trim2, color='ORANGE')



###################################
#Test Offset Algorithms
###################################

testEdge = BRepBuilderAPI_MakeEdge(offset_curves[1].GetHandle())
testWire = BRepBuilderAPI_MakeWire(testEdge.Edge())

#offsetDistance = -0.004
#testoffset = make_offset(aWire.Wire(), offsetDistance, altitude=1, joinType=0)
#display.DisplayShape(testoffset,color='ORANGE')
    
###################################
#TRAILING EDGE TAB DEFINITION
###################################
#The tab is parallel to the end of the camber line, has a round end, 4-5% of chord length, thickness of ~0.8% of chord, radius of 2% of chord, 
'''Comment: This has to be completely functional for all kinds of different airfoil shapes!''' 
a = np.vstack((data[:,0],data[:,-1]))
b = np.vstack((data[:,1],data[:,-2]))
tab_pnt  = a.mean(axis=0)                       #tab_reference point
tab_pnt2 = b.mean(axis=0)
vec_0 = tab_pnt - tab_pnt2 
vec_1 = tab_pnt
tab_rotation_angle = -angle_between(vec_0,vec_1)       #angle between last profil points

tab_thickness = 0.008  # times basic chord lengtt,   has to be a multiple of layup and the same over total rotor blade!       
tab_length = 0.04      # times basic chord length
tab_radius = 0.03      # times basic chord length
tab_radius_angle = math.radians(-45)
#tab_pnt = np.array([0.0, -0.75, 0.0])

tab_template = tab_definition(tab_pnt,tab_thickness,tab_length, tab_radius, tab_radius_angle)
#display.DisplayShape(tab_wire,color='BLUE')

#Rotation
tab_gpPnt = np_array_to_gp_Pnt(tab_pnt)
ax1 = gp_Ax1(tab_gpPnt, gp_Dir(1., 0., 0.))
tab_wire  = rotate_TopoDS_wire(tab_template,ax1,tab_rotation_angle)
display.DisplayShape(tab_wire ,color='RED')





###################################
#Balance Weight
###################################



###########################
#       STP EXPORT
###########################
# initialize the STEP exporter
step_writer = STEPControl_Writer()
Interface_Static_SetCVal("write.step.schema", "AP203")

# transfer shapes and write file
step_writer.Transfer(aWire.Shape(), STEPControl_AsIs)
status = step_writer.Write("airfoil.stp")

assert(status == IFSelect_RetDone)







###########################
#       START DISPLAY
###########################
#
display.FitAll()
start_display()









