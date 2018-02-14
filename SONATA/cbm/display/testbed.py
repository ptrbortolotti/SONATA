# -*- coding: utf-8 -*-
"""
Created on Tue Apr 04 14:49:03 2017

@author: TPflumm
"""

#Basic PYTHON Modules:
import numpy as np       
import matplotlib as plt
from datetime import datetime

from OCC.gp import gp_Pnt, gp_XOY,gp_Ax22d,gp_Pnt2d,gp_Dir,gp_Dir2d
from OCC.GeomAPI import GeomAPI_ProjectPointOnCurve
from OCC.Geom import Geom_Circle

from OCC.Geom2d import  Geom2d_Circle
from OCC.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.Geom2dAPI import Geom2dAPI_PointsToBSpline, Geom2dAPI_Interpolate

from OCC.GCPnts import GCPnts_QuasiUniformAbscissa

from OCC.Display.SimpleGui import init_display
from OCC.TColgp import (TColgp_Array1OfPnt, TColgp_Array1OfPnt2d, TColgp_HArray1OfPnt2d, 
                        TColgp_HArray1OfPnt )

def point2d_list_to_TColgp_Array1OfPnt2d(li):
    return _Tcol_dim_1(li, TColgp_Array1OfPnt2d)

def _Tcol_dim_1(li, _type):
    pts = _type(0, len(li)-1)
    for n, i in enumerate(li):
        pts.SetValue(n, i)
    return pts


display, start_display, add_menu, add_function_to_menu = init_display('wx')

point_to_project = gp_Pnt(1., 2., 0)
radius = 5.

# create a circle, centered at origin with a given radius
circle = Geom_Circle(gp_XOY(), radius)
circle2d = Geom2d_Circle(gp_Ax22d(gp_Pnt2d(0,0),gp_Dir2d(1,0),gp_Dir2d(0,1)) , radius)

Adaptor = Geom2dAdaptor_Curve(circle2d.GetHandle())
NbPoints = 100 
discretization = GCPnts_QuasiUniformAbscissa(Adaptor,NbPoints) 

Pnt2d_Lst = []
for j in range(1, discretization.NbPoints()):
        para = discretization.Parameter(j)
        P = gp_Pnt2d()
        circle2d.D0(para,P)
        Pnt2d_Lst.append(P)
        display.DisplayShape(P,color='WHITE')


BSpline = Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d(Pnt2d_Lst)).Curve().GetObject()   
display.DisplayShape(BSpline,color='BLACK')
display.DisplayShape(BSpline.StartPoint(),color='ORANGE')
BSpline2 = Geom2dAPI_Interpolate(point2d_list_to_TColgp_Array1OfPnt2d(Pnt2d_Lst)).Curve().GetObject()   

print(BSpline.StartPoint().X(), BSpline.StartPoint().Y())
print(BSpline.EndPoint().X(), BSpline.EndPoint().Y())
print(BSpline.IsClosed())



display.DisplayShape(circle)
display.DisplayShape(point_to_project, update=True)
display.DisplayMessage(point_to_project, "P")

# project the point P on the circle
projection = GeomAPI_ProjectPointOnCurve(point_to_project,circle.GetHandle())
# get the results of the projection
# the point
projected_point = projection.NearestPoint()
# the number of possible results
nb_results = projection.NbPoints()
print(("NbResults : %i" % nb_results))

pstring = "N : at Distance : %f" % projection.LowerDistance()
print(pstring)
#display.DisplayMessage(projected_point, pstring)

# thre maybe many different possible solutions
if nb_results > 0:
    for i in range(1, nb_results+1):
        Q = projection.Point(i)
        distance = projection.Distance(i)
        pstring = "Q%i: at Distance :%f" % (i, distance)
        print(pstring)
        display.DisplayShape(Q,color='Blue')
        #display.DisplayMessage(Q, pstring)


start_display()
