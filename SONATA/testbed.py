import numpy as np
import matplotlib.pyplot as plt



from OCC.gp import gp_Pnt2d
from OCC.Geom2dAPI import Geom2dAPI_Interpolate, Geom2dAPI_PointsToBSpline, Geom2dAPI_ProjectPointOnCurve
from OCC.TColgp import TColgp_HArray1OfPnt2d, TColgp_Array1OfPnt2d

from OCC.Display.SimpleGui import init_display
display, start_display, add_menu, add_function_to_menu = init_display('wx')
display.Context.SetDeviationAngle(0.0000002)       # 0.001 default. Be careful to scale it to the problem.
display.Context.SetDeviationCoefficient(0.0000002) # 0.001 default. Be careful to scale it to the problem. 


def findPnt_on_2dcurve(Pnt2d,Curve2d):
    projection = Geom2dAPI_ProjectPointOnCurve(Pnt2d,Curve2d)
    u = projection.LowerDistanceParameter()
    return u


if __name__ == '__main__': 
    
    # the first bspline
    array = TColgp_Array1OfPnt2d(1, 5)
    array.SetValue(1, gp_Pnt2d(0, 0))
    array.SetValue(2, gp_Pnt2d(1, 2))
    array.SetValue(3, gp_Pnt2d(2, 3))
    array.SetValue(4, gp_Pnt2d(4, 3))
    array.SetValue(5, gp_Pnt2d(5, 5))
    bspline_1 = Geom2dAPI_PointsToBSpline(array).Curve().GetObject()
    
    Curve2d = bspline_1
    P = gp_Pnt2d(2, 3)
    u = findPnt_on_2dcurve(P,bspline_1.GetHandle())
    print u


    
    for j in range(array.Lower(), array.Upper()+1):
        p = array.Value(j)
        display.DisplayShape(p, update=False)

    
        
    display.DisplayShape(bspline_1, update=False)
    
    
    display.View_Top()
    display.FitAll()
    start_display()