#-------------------------------
#          H E A D E R
#-------------------------------

#Basic Libraries
import math
import numpy as np                                                              #fundamental package for scientific computing with Python (powerful N-Dimensional array objects)
import matplotlib.pyplot as plt                                                 #python 2D plotting library
from scipy.interpolate import interp1d
from scipy.optimize import leastsq, broyden1, brentq, bisect, newton, fsolve,root                            # help in opening URL

#Python OCC Libraries
from OCC.gp import gp_Pnt, gp_Vec,  gp_Pln, gp_Dir, gp_Trsf, gp_Ax1, gp_OX, gp_Ax3, gp_Ax2, gp_Circ, gp_OY
from OCC.gp import gp_Pnt2d, gp_Vec2d, gp_XY, gp_Lin2d, gp_Dir2d
from OCC.gp import gp_Pnt, gp_Dir, gp_Vec
from OCC.GCE2d import  GCE2d_MakeSegment
from OCC.Geom2d import Geom2d_OffsetCurve, Geom2d_TrimmedCurve, Geom2d_Line, Geom2d_BezierCurve
from OCC.Geom2dConvert import Geom2dConvert_CompCurveToBSplineCurve
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





#======================================================
#       MAIN
#======================================================
if __name__ == '__main__':
        
    
    airfoil = 'naca23012'  
    
    ###############################################################################
    display, start_display, add_menu, add_function_to_menu = init_display()
    display.Context.SetDeviationAngle(0.000005)      # 0.001 default
    display.Context.SetDeviationCoefficient(0.000005) # 0.001 default
    
    #CREATE AXIS SYSTEM for Visualization
    COSY = gp_Ax3()	
    O  = gp_Pnt(0., 0., 0.)
    p1 = gp_Pnt(0.1,0.,0.)
    p2 = gp_Pnt(0.,0.1,0.)
    p3 = gp_Pnt(0.,0.,0.1)
    
    h1 = BRepBuilderAPI_MakeEdge(O,p1).Shape()
    h2 = BRepBuilderAPI_MakeEdge(O,p2).Shape()
    h3 = BRepBuilderAPI_MakeEdge(O,p3).Shape()
    
    display.DisplayShape(O,color='BLACK')
    display.DisplayShape(h1,color='RED')
    display.DisplayShape(h2,color='GREEN')
    display.DisplayShape(h3,color='BLUE')
 
      
    ###############################################################################
    #GET OUTER AIRFOIL DATA
    ###################################
    data = UIUCAirfoil2d(airfoil)                                                       #Get Airfoil Data from Database 
    harray = TColgp_HArray1OfPnt2d_from_nparray(data)                                   #Put Data into Harray1OfPnt2d
    #display_points_of_array(harray)                                                     #Display Points
    anInterpolation = Geom2dAPI_Interpolate(harray.GetHandle(), False, 0.001)           #Interpolate datapoints to bspline
    anInterpolation.Perform()                                               
    bspline = anInterpolation.Curve().GetObject()
    
    #Convert to 3D SPACE onto x,y plane
    P = gp_Pnt(0,0,0)
    V = gp_Dir(gp_Vec(0,0,1))
    Plane = Geom_Plane(P, V)
    
    edge = BRepBuilderAPI_MakeEdge(bspline.GetHandle(),Plane.GetHandle())
    #display.DisplayShape(edge.Shape(), update=True, color='WHITE')
    aWire = BRepBuilderAPI_MakeWire(edge.Edge())
    
    ###############################################################################
    #SPLIT AIRFOIL INTO SEGMENTS
    ###################################
    sep1 = 0.2
    sep2 = 0.9
    Dir = gp_Dir2d(gp_Vec2d(0,1))
    
    sepLine1 = Geom2d_Line(gp_Pnt2d(sep1,0),Dir)
    sepLine2 = Geom2d_Line(gp_Pnt2d(sep2,0),Dir)
    Inter1 = Geom2dAPI_InterCurveCurve(bspline.GetHandle(), sepLine1.GetHandle(), 1.0e-9)
    print('Number of Intersection Points:', Inter1.NbPoints())
    
    
    InterP1 = Inter1.Point(1)
    InterP2 = Inter1.Point(2)
    display.DisplayShape(InterP1,color='BLUE')
    display.DisplayShape(InterP2,color='RED')
    
    def findPnt_on_2dcurve(Pnt,curve,u0):
        def Pnt_Distance(u,Pnt):
            u = float(u[0])
            #print(u)
            p = gp_Pnt2d()
            curve.D0(u,p)
            error = p.Distance(Pnt)
            return error
        
        y = leastsq(Pnt_Distance,u0,args=(Pnt))
        error = Pnt_Distance(y,Pnt)
        print('u:',float(y[0]),'Error:', error)
        u = float(y[0])
        return u
        
    first = bspline.FirstParameter()
    last =  bspline.LastParameter()
        
    u1 = findPnt_on_2dcurve(InterP1,bspline,first+0.1)
    u2 = findPnt_on_2dcurve(InterP2,bspline,last-0.1)
    
    UP_curve = Geom2d_TrimmedCurve(bspline.GetHandle(), first, u1)
    LE_curve = Geom2d_TrimmedCurve(bspline.GetHandle(), u1, u2)
    LO_curve = Geom2d_TrimmedCurve(bspline.GetHandle(), u2, last)
    
    display.DisplayShape(UP_curve, color='CYAN')
    display.DisplayShape(LE_curve , color='GREEN')
    display.DisplayShape(LO_curve, color='BLUE')
    
    ###############################################################################
    #CREATE SKIN LAYER
    ###################################
    
    dist = -0.004
    num_of_layers = 3
    
    UP_offset_curves = []
    geom_bspline = UP_curve.GetHandle()
    for j in range(0,num_of_layers):
        geom_bspline= Geom2d_OffsetCurve(geom_bspline, dist)
        display.DisplayShape(geom_bspline, color='CYAN', update=True)
        UP_offset_curves.append(geom_bspline)     
        geom_bspline = geom_bspline.GetHandle()
    
    LO_offset_curves = []
    geom_bspline = LO_curve.GetHandle()
    for j in range(0,num_of_layers):
        geom_bspline= Geom2d_OffsetCurve(geom_bspline, dist)
        display.DisplayShape(geom_bspline, color='BLUE', update=True)
        LO_offset_curves.append(geom_bspline)     
        geom_bspline = geom_bspline.GetHandle()
    
    num_of_layers = 8
    LE_offset_curves = []
    geom_bspline = LE_curve.GetHandle()
    for j in range(0,num_of_layers):
        geom_bspline= Geom2d_OffsetCurve(geom_bspline, dist)
        display.DisplayShape(geom_bspline, color='GREEN', update=True)
        LE_offset_curves.append(geom_bspline)     
        geom_bspline = geom_bspline.GetHandle()
    
    
    
    ###############################################################################
    #REMOVE INTERSECTION ON TE
    ###################################
    
    UP_curve_1 = UP_offset_curves[0]
    LO_curve_1 = LO_offset_curves[0]
    TE_Inter = Geom2dAPI_InterCurveCurve(UP_curve_1.GetHandle() , LO_curve_1.GetHandle(), 1.0e-11)
    print('Number of Intersection Points:', TE_Inter.NbPoints())
    TE_InterP1 = TE_Inter.Point(1)
    display.DisplayShape(TE_InterP1,color='RED')
    
    UP_first = UP_curve.FirstParameter()
    UP_last =  UP_curve.LastParameter()
    
    LO_first = LO_curve.FirstParameter()
    LO_last =  LO_curve.LastParameter()
    
    UP_U = findPnt_on_2dcurve(TE_InterP1,UP_curve_1,first+0.1)
    LO_U = findPnt_on_2dcurve(TE_InterP1,LO_curve_1,last-0.1)
    
    UP_curve_1T = Geom2d_TrimmedCurve(UP_curve_1.GetHandle(), UP_last , UP_U )
    LE_curve_1T = Geom2d_TrimmedCurve(LO_curve_1.GetHandle(), LO_first , LO_U)
    
    display.DisplayShape(UP_curve_1T, color='RED')
    display.DisplayShape(LE_curve_1T , color='ORANGE')
    
    
    display.FitAll()
    start_display()


