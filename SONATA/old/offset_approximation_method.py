from __future__ import print_function

import numpy as np   
import math
#from scipy.optimize import least_squares
from scipy.optimize import leastsq
  
from OCC.gp import gp_Pnt2d, gp_Vec2d, gp_XY, gp_Lin2d, gp_Dir2d
from OCC.GCE2d import  GCE2d_MakeSegment
from OCC.TColgp import TColgp_Array1OfPnt2d
from OCC.Geom2dAPI import Geom2dAPI_PointsToBSpline, Geom2dAPI_InterCurveCurve,Geom2dAPI_ProjectPointOnCurve
from OCC.Geom2d import Geom2d_OffsetCurve, Geom2d_TrimmedCurve, Geom2d_Line, Geom2d_BezierCurve
from OCC.Geom2dConvert import Geom2dConvert_BSplineCurveToBezierCurve


from OCC.Display.SimpleGui import init_display
display, start_display, add_menu, add_function_to_menu = init_display()
display.Context.SetDeviationAngle(0.00001)      # 0.001 default
display.Context.SetDeviationCoefficient(0.00005) # 0.001 default


def display_points_of_array(array, COLOR):
    for j in range(array.Lower(),array.Upper()+1):
        p = array.Value(j)
        display.DisplayShape(p, color=COLOR, update=False)

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


#--------------------------------------------------------
# Hoschek - Spline approximation of offset curves:
#--------------------------------------------------------


def BezierOffset(parameters,BezierApproximation,dist):
    Bez_NbPoles = BezierApproximation.NbPoles()
    Bez_Poles = TColgp_Array1OfPnt2d(1, Bez_NbPoles)
    BezierApproximation.Poles(Bez_Poles)
    V = Bez_Poles 
    
    lambda_t1 =  parameters[0]    
    lambda_t2 =  parameters[1]  
    
    Bez_NbPoles = BezierApproximation.NbPoles()
    Bez_Poles  = TColgp_Array1OfPnt2d(1, Bez_NbPoles)    
    BezierApproximation.Poles(Bez_Poles)    
        
    W = TColgp_Array1OfPnt2d(1, Bez_NbPoles)    
    p = gp_Pnt2d()
    t =  gp_Vec2d()
        
    #Set FIRST Boundary Point (W0):
    V0 = gp_Pnt2d_to_npArray(V.Value(1))
    first = BezierApproximation.FirstParameter()
    BezierApproximation.D1(first,p,t)
    t.Normalize()
    n = t.GetNormal()
    N = gp_Vec2d_to_npArray(n)
    W0 = V0+dist*N
    gp_W0 = npArray_to_gp_Pnt2d(W0)
    
    W.SetValue(1,gp_W0)
        
    #Set LAST Boundary Point (W0):
    V3 = gp_Pnt2d_to_npArray(V.Value(Bez_NbPoles))
    last = BezierApproximation.LastParameter()
    BezierApproximation.D1(last,p,t)
    t.Normalize()
    n = t.GetNormal()
    N = gp_Vec2d_to_npArray(n)
    W3 = V3+dist*N
    gp_W3 = npArray_to_gp_Pnt2d(W3)
        
    W.SetValue(4,gp_W3) 
    
    #--------------------------------------------------------
    #Set Initial Condition for W1 and W2
    #--------------------------------------------------------
    V1 = gp_Pnt2d_to_npArray(V.Value(2))
    V2= gp_Pnt2d_to_npArray(V.Value(Bez_NbPoles-1))
      
    W1 = W0+lambda_t1 *(V1-V0)
    W2 = W3+lambda_t2 *(V2-V3)
    
    gp_W1 = npArray_to_gp_Pnt2d(W1)
    gp_W2 = npArray_to_gp_Pnt2d(W2)
    
    W.SetValue(2,gp_W1)
    W.SetValue(3,gp_W2) 
    
    HoschekOffsetCurve = Geom2d_BezierCurve(W)
    
    #display.DisplayShape(HoschekOffsetCurve , color='BLUE')
    display_points_of_array(V,'RED')
    #display_points_of_array(W,'BLUE')
    return HoschekOffsetCurve


def BezierOffset_Error(parameters,X_Bezier,dist):
    Y_Bezier = BezierOffsetG2(parameters,X_Bezier,dist)
    #display.DisplayShape(Y_Bezier, color='CYAN')
    #--------------------------------------------------------
    #DEFINE MEASUREMENT OF ERROR 
    #--------------------------------------------------------
    
   #Spread equidistant points on Y_Bezier   
    n = X_Bezier.Degree()
    m = Y_Bezier.Degree()   
    k = 2*max(m,n)
    nPoints = k+1
    nPoints = 30                                                               #BUG: PROJECTION der equidistanten punkte von X auf Y. kann zu fehlern fueren, besonders wenn anzahl der Kontrollpunkte gering ist
    #print("equidistant points:",nPoints)

    p = gp_Pnt2d()
    v1 = gp_Vec2d()
    v2 = gp_Vec2d()
    s = np.linspace(0,1,nPoints)
        
    #EXACT OFFSET:
    Xd_Bezier = Geom2d_OffsetCurve(X_Bezier.GetHandle(), dist)              
    display.DisplayShape(Xd_Bezier, color='YELLOW')
    

    #OFFSET DISTANCE CORRECTION BY LEIHONG LI, for constant local area
    error = np.zeros(nPoints)
    for i,u in enumerate(s): 
        X_Bezier.D2(u,p,v1,v2)
        normal = v1.GetNormal()
        #display.DisplayShape(p,color='MAGENTA')
        radius = radius_of_curve(X_Bezier,u)         
        Projection = Geom2dAPI_ProjectPointOnCurve(p,Y_Bezier.GetHandle())
        nearestP = Projection.NearestPoint()
        dist_vec = gp_Vec2d (p, nearestP)
        Distance = Projection.LowerDistance()
        #print(Distance)
        #print(dist_vec.Magnitude())        
        scalarproduct = normal.Dot(dist_vec)
        #print("SCALARPRODUCT",scalarproduct)
        if scalarproduct < 0:
            Distance = -Projection.LowerDistance()
        
        print('Radius:',radius)
        VexVsCave = normal.Dot(v2)*dist
        if  VexVsCave > 0:  # CONCAVE: (Inside Curve)
            if radius > (2*abs(dist)):
                dist_corr = radius - math.sqrt(radius**2 - 2*radius*abs(dist))
                #dist_corr = dist
                print('Radius:',radius)  
                print ('h0:',dist,';\th:',dist_corr)
            else: 
                dist_corr =  2*dist  
        
        else:   #CONVEX: (Outside curve)
            if radius > (2*abs(dist)):
                dist_corr = -radius + math.sqrt(radius**2 + 2*radius*abs(dist))
                #dist_corr = dist
                print('Radius:',radius)  
                print ('h0:',dist,';\th:',dist_corr)
            else: 
                dist_corr = dist
                
        #print('Radius:',radius)    
        #print ('h0:',dist,';\th:',dist_corr)
        error[i] = Distance-dist_corr 
        #display.DisplayShape(p,color='ORANGE')
        
    #print(error)
    return error
    

def OffsetApprox(X_Bezier,dist):
    parameters = np.array([0.9,0.9,0.5,0.5])        #lambda_1,lambda_2,mu_1,mu_2
    x = leastsq(BezierOffset_Error, parameters, args=(X_Bezier,dist),full_output=1)
    Y_Bezier = BezierOffsetG2(x[0],X_Bezier, dist)
    error = BezierOffset_Error(x[0],X_Bezier, dist)
    print(error)
    return Y_Bezier
    
def radius_of_curve(curve,u):
    p = gp_Pnt2d()
    v1 = gp_Vec2d()
    v2 = gp_Vec2d()
    curve.D2(u,p,v1,v2)
    eq1 = v1.Crossed(v2)      
    eq2 = v1.Magnitude()
    curvature = abs(eq1)/abs(eq2**3)
    radius = 1/curvature
    return radius
    
def BezierOffsetG2(parameters,X_Bezier,dist):
    lambda_1 =  parameters[0]    
    lambda_2 =  parameters[1] 
    mu_1 =  parameters[2] 
    mu_2 =  parameters[3] 
    NbPoles = X_Bezier.NbPoles()
    Deg = X_Bezier.Degree()
    #print ("#INFO: Number of Poles of original Bezier Curve X(t):",NbPoles)
    #print ("#INFO: Degree of original Bezier Curve X(t):",Deg)    
    
    
    V =  TColgp_Array1OfPnt2d(1, NbPoles)    
    X_Bezier.Poles(V)
        
    W = TColgp_Array1OfPnt2d(1, 6) 
    p = gp_Pnt2d()
    t =  gp_Vec2d()
        
    #Set FIRST Boundary Point (W0):
    V0 = gp_Pnt2d_to_npArray(V.Value(1))
    first = X_Bezier.FirstParameter()
    X_Bezier.D1(first,p,t)
    t.Normalize()
    n = t.GetNormal()
    N = gp_Vec2d_to_npArray(n)
    W0 = V0+dist*N
    gp_W0 = npArray_to_gp_Pnt2d(W0)
    W.SetValue(1,gp_W0)
      
    #Set LAST Boundary Point (W5):
    VN = gp_Pnt2d_to_npArray(V.Value(NbPoles))
    last = X_Bezier.LastParameter()
    X_Bezier.D1(last,p,t)
    t.Normalize()
    n = t.GetNormal()
    N = gp_Vec2d_to_npArray(n)
    W5 = VN+dist*N
    gp_W5 = npArray_to_gp_Pnt2d(W5)  
    W.SetValue(6,gp_W5) 
   
    #--------------------------------------------------------
    #Set Initial Condition for W1, W2, W3 and W4
    #--------------------------------------------------------
    V1 = gp_Pnt2d_to_npArray(V.Value(2))
    V2 = gp_Pnt2d_to_npArray(V.Value(3))
    VN_1= gp_Pnt2d_to_npArray(V.Value(NbPoles-1))
    VN_2= gp_Pnt2d_to_npArray(V.Value(NbPoles-2))

    m = 5
    n = Deg 
    
    kappa_0 = 1 / radius_of_curve(X_Bezier,first)
    kappa_1 = 1 / radius_of_curve(X_Bezier,last )
    
    #print("kappa0,1:", kappa_0, kappa_1)    
    
    k0 = (m*(n-1)) / (n*(m-1)) * (1+kappa_0*dist)**(-1)
    k1 = (m*(n-1)) / (n*(m-1)) * (1+kappa_1*dist)**(-1)

    #print("k0,1:", k0, k1)      
    
    W1 = W0+lambda_1 *(V1-V0)
    W4 = W5+lambda_2 *(VN_1-VN)
    W2 = W0+lambda_1**2*k0*(V2-V1)+mu_1*(V1-V0)
    W3 = W5+lambda_2**2*k1*(VN_2-VN_1)+mu_2*(VN-VN_1)
    
    gp_W1 = npArray_to_gp_Pnt2d(W1)
    gp_W2 = npArray_to_gp_Pnt2d(W2)
    gp_W3 = npArray_to_gp_Pnt2d(W3)
    gp_W4 = npArray_to_gp_Pnt2d(W4)
      
    W.SetValue(2,gp_W1)
    W.SetValue(3,gp_W2) 
    W.SetValue(4,gp_W3)
    W.SetValue(5,gp_W4)     
    
    Y_Bezier = Geom2d_BezierCurve(W)
    
    #display.DisplayShape(Y_Bezier, color='CYAN')
    #display_points_of_array(V,'RED')
    #display_points_of_array(W,'BLUE')
    return Y_Bezier 




    
#============================================
#SUBDIVIDED CUBIC BEZIER SPLINE APPROXIMATION
#============================================


array0 = TColgp_Array1OfPnt2d(1, 4)
array0.SetValue(1, gp_Pnt2d(0, 0))
array0.SetValue(2, gp_Pnt2d(-3, 2))
array0.SetValue(3, gp_Pnt2d(-2, 3))
array0.SetValue(4, gp_Pnt2d(0, 3))

curve0 = Geom2d_BezierCurve(array0)

P1 = gp_Pnt2d(0, 0)
P2 = gp_Pnt2d(1, 1)
P3 = gp_Pnt2d(2, 4)
P4 = gp_Pnt2d(8, 9)
P5 = gp_Pnt2d(9, 16)
P6 = gp_Pnt2d(5, 25)
P7 = gp_Pnt2d(6, 36)

array1 = TColgp_Array1OfPnt2d(1, 7)
array1.SetValue(1, P1)
array1.SetValue(2, P2)
array1.SetValue(3, P3)
array1.SetValue(4, P4)
array1.SetValue(5, P5)
array1.SetValue(6, P6)
array1.SetValue(7, P7)

array2 = TColgp_Array1OfPnt2d(1, 4)
array2.SetValue(1, P4)
array2.SetValue(2, P5)
array2.SetValue(3, P6)
array2.SetValue(4, P7)

curve1 = Geom2d_BezierCurve(array1)
curve2 = Geom2d_BezierCurve(array2)

display.DisplayShape(curve1, color='BLUE')
#display.DisplayShape(curve2, color='CYAN')
display_points_of_array(array1,'GREEN')
#display_points_of_array(array2,'BLUE')






#============================================
#FUCTION GET NORMAL VECTOR:
#============================================    
nPoints = 100                                                             
dist = 1.2

p = gp_Pnt2d()
v1 = gp_Vec2d()
v2 = gp_Vec2d()
s = np.linspace(0,1,nPoints)
    
error = np.zeros(nPoints)
for i,u in enumerate(s):  
        curve1.D2(u,p,v1,v2)
        normal = v1.GetNormal()
        Line = gp_Lin2d(p,gp_Dir2d(normal) )
        Normal_Segment = GCE2d_MakeSegment(Line, 0, dist)
        display.DisplayShape(Normal_Segment.Value(), color='WHITE')
        
        Line = gp_Lin2d(p,gp_Dir2d(v2) )
        Normal_Segment = GCE2d_MakeSegment(Line, 0, dist)
        
        scalarproduct = normal.Dot(v2)*dist
        if  scalarproduct > 0:  # CONCAVE
            display.DisplayShape(Normal_Segment.Value(), color='GREEN')
        
        else:   #CONVEX
            display.DisplayShape(Normal_Segment.Value(), color='RED')
            







#=====================================================================
# OFFSET CALCULATION - MAIN
#=====================================================================


parameters = np.array([0.9,0.9,0.5,0.5])        #lambda_1,lambda_2,mu_1,mu_2
#BezierOffsetG2(parameters,curve1,dist)
#BezierOffset_Error(parameters,curve1,dist)
OffsetG2 = OffsetApprox(curve1,dist)
#OffsetG2_2 = OffsetApprox(OffsetG2,dist)
#OffsetG2_3 = OffsetApprox(OffsetG2_2,dist)
display.DisplayShape(OffsetG2, color='RED')
#display.DisplayShape(OffsetG2_2, color='ORANGE')
#display.DisplayShape(OffsetG2_3, color='GREEN')
#Xd_Bezier = Geom2d_OffsetCurve(curve1.GetHandle(), dist)
#display.DisplayShape(Xd_Bezier, color='YELLOW')


#Offset_1 = OffsetApprox(curve0, dist)
#Offset_2 = OffsetApprox(Offset_1,dist)
#Offset_3 = OffsetApprox(curve1,dist)
#Offset_4 = OffsetApprox(curve2,dist)
#Offset_5 = OffsetApprox(Offset_4,dist)
#Offset_6 = OffsetApprox(Offset_5,dist)

#display.EraseAll()
#display.DisplayShape(curve0, color='RED')
#display.DisplayShape(Offset_1, color='ORANGE')
#display.DisplayShape(Offset_2, color='YELLOW')
#display.DisplayShape(Offset_3, color='GREEN')
#display.DisplayShape(Offset_4, color='BLACK')
#display.DisplayShape(Offset_5, color='BLUE')
#display.DisplayShape(Offset_6, color='MAGENTA')


display.FitAll()

start_display()