#===================================================================================================================
#                               H   E   A   D   E   R                                
#===================================================================================================================

import numpy as np   
from scipy.optimize import leastsq
  
#==============GEOMETRY LIBs==================================================
from OCC.TColgp import TColgp_Array1OfPnt2d

from OCC.gp import gp_Pnt2d, gp_Vec2d, gp_XY, gp_Lin2d, gp_Dir2d
from OCC.gp import gp_Pnt, gp_Dir, gp_Vec
from OCC.GCE2d import  GCE2d_MakeSegment

from OCC.Geom2d import Geom2d_OffsetCurve, Geom2d_TrimmedCurve, Geom2d_Line, Geom2d_BezierCurve
from OCC.Geom2dAPI import Geom2dAPI_PointsToBSpline, Geom2dAPI_InterCurveCurve
from OCC.Geom import Geom_Plane

#==============BRep LIBs=====================================================
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire

#==============DISPLAY LIBs===================================================
#from OCC.Display.SimpleGui import init_display
#display, start_display, add_menu, add_function_to_menu = init_display()
#display.Context.SetDeviationAngle(0.00001)      # 0.001 default
#display.Context.SetDeviationCoefficient(0.00001) # 0.001 default

#===================================================================================================================
#                               B   O   D   Y                             
#===================================================================================================================


def display_points_of_array(array, COLOR):
    for j in range(array.Lower(),array.Upper()+1):
        p = array.Value(j)
        #display.DisplayShape(p, color=COLOR, update=False)

def gp_Pnt2d_to_npArray(Ptn2d):
        vector = np.array([Ptn2d.X(),Ptn2d.Y()])
        return vector

def gp_Vec2d_to_npArray(Vec2d):
        vector = np.array([Vec2d.X(),Vec2d.Y()])
        return vector
        
def np_GetNormal2d(Vec2d):
        vector = np.array([-Vec2d[1],Vec2d[0]])
        return vector
        
def npArray_to_gp_Pnt2d(vector):
        Pnt2d = gp_Pnt2d(vector[0],vector[1])
        return Pnt2d       

def trim_2dcurve_selfintersectingLoop(offset1,InterP1):
    #display.DisplayShape(InterP1, color='RED')       
    
    def Pnt_Distance(u,InterP1):
        u = float(u[0])
        #print(u)
        p = gp_Pnt2d()
        offset1.D0(u,p)
        error = p.Distance(InterP1)
        return error

    u0 = 0.2  #Starting Guess
    y = leastsq(Pnt_Distance,u0,args=(InterP1))
    error = Pnt_Distance(y,InterP1)
    print('u:',float(y[0]),'Error:', error)
    t1 = float(y[0])

    u0 = 0.9  #Starting Guess
    y = leastsq(Pnt_Distance,u0,args=(InterP1))
    error = Pnt_Distance(y,InterP1)
    print('u:',float(y[0]),'Error:', error)
    t2 = float(y[0])

    
    Trim1 = Geom2d_TrimmedCurve(offset1.GetHandle(), 0, t1)
    Trim2 = Geom2d_TrimmedCurve(offset1.GetHandle(), t2, 1)
    #display.DisplayShape(Trim1, color='CYAN')
    #display.DisplayShape(Trim2, color='CYAN')
     
    
    #Convert to 3D SPACE onto x,y plane
    P = gp_Pnt(0,0,0)
    V = gp_Dir(gp_Vec(0,0,1))
    Plane = Geom_Plane(P, V)
    
    aEdge1 = BRepBuilderAPI_MakeEdge(Trim1.GetHandle(),Plane.GetHandle())
    aEdge2 = BRepBuilderAPI_MakeEdge(Trim2.GetHandle(),Plane.GetHandle())
    aWire = BRepBuilderAPI_MakeWire(aEdge1.Edge(), aEdge2.Edge())
    #display.DisplayShape(aWire.Wire(), color='BLACK')
    
    return aWire


#==============================================================================
#                       M A I N
#==============================================================================
if __name__ == '__main__':
    

    #============================================
    #SUBDIVIDED CUBIC BEZIER SPLINE APPROXIMATION
    #============================================
    dist = -13
    
    P1 = gp_Pnt2d(0, 0)
    P2 = gp_Pnt2d(1, 1)
    P3 = gp_Pnt2d(2, 4)
    P4 = gp_Pnt2d(8, 9)
    P5 = gp_Pnt2d(9, 30)
    P6 = gp_Pnt2d(5, 10)
    P7 = gp_Pnt2d(-30, 36)
    
    array1 = TColgp_Array1OfPnt2d(1, 7)
    array1.SetValue(1, P1)
    array1.SetValue(2, P2)
    array1.SetValue(3, P3)
    array1.SetValue(4, P4)
    array1.SetValue(5, P5)
    array1.SetValue(6, P6)
    array1.SetValue(7, P7)
    
    
    curve1 = Geom2d_BezierCurve(array1)
    offset1 = Geom2d_OffsetCurve(curve1.GetHandle(), dist)              
    
    #display.DisplayShape(curve1, color='BLUE')
    #display.DisplayShape(offset1, color='YELLOW')
    display_points_of_array(array1,'GREEN')
    
    Inter = Geom2dAPI_InterCurveCurve(offset1.GetHandle(), 1.0e-5)
    print('Number of Intersection Points:', Inter.NbPoints())
    #print('Number of tangetial Intersections:',Inter.NbSegments())
    InterP1 = Inter.Point(1)
    
    
    
    TestWire = trim_2dcurve_selfintersectingLoop(offset1,InterP1)
    #display.DisplayShape(TestWire.Wire(), color='BLACK')
        
    
    #============================================
    # DISPLAY NORMAL VECTOR AND V2 
    #============================================    
    #nPoints = 100                                                             
    #dist = 1.2
    #
    #p = gp_Pnt2d()
    #v1 = gp_Vec2d()
    #v2 = gp_Vec2d()
    #s = np.linspace(0,1,nPoints)
    #    
    #error = np.zeros(nPoints)
    #for i,u in enumerate(s):  
    #        curve1.D2(u,p,v1,v2)
    #        normal = v1.GetNormal()
    #        Line = gp_Lin2d(p,gp_Dir2d(normal) )
    #        Normal_Segment = GCE2d_MakeSegment(Line, 0, dist)
    #        display.DisplayShape(Normal_Segment.Value(), color='WHITE')
    #        
    #        Line = gp_Lin2d(p,gp_Dir2d(v2) )
    #        Normal_Segment = GCE2d_MakeSegment(Line, 0, dist)
    #        
    #        scalarproduct = normal.Dot(v2)*dist
    #        if  scalarproduct > 0:  # CONCAVE
    #            display.DisplayShape(Normal_Segment.Value(), color='GREEN')
    #        
    #        else:   #CONVEX
    #            display.DisplayShape(Normal_Segment.Value(), color='RED')
                
    #display.FitAll()
    
    #start_display()