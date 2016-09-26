from StringIO import StringIO
import matplotlib.pyplot as plt
import numpy as np
import requests
import shapely.geometry as shp
import math, random



#Python OCC Libraries
from OCC.gp import gp_Pnt, gp_Vec,  gp_Pln, gp_Dir, gp_Trsf, gp_Ax1, gp_OX, gp_Ax3, gp_Ax2, gp_Circ, gp_OY
from OCC.gp import gp_Pnt2d, gp_Vec2d, gp_XY, gp_Lin2d, gp_Dir2d
from OCC.gp import gp_Pnt, gp_Dir, gp_Vec
from OCC.Geom2d import Geom2d_OffsetCurve, Geom2d_TrimmedCurve, Geom2d_Line, Geom2d_BezierCurve
from OCC.Geom2dConvert import Geom2dConvert_CompCurveToBSplineCurve
from OCC.Geom2dAPI import Geom2dAPI_PointsToBSpline, Geom2dAPI_Interpolate, Geom2dAPI_InterCurveCurve 
from OCC.GeomAPI import geomapi, GeomAPI_PointsToBSpline,GeomAPI_Interpolate
from OCC.TColgp import TColgp_Array1OfPnt, TColgp_Array1OfPnt2d, TColgp_HArray1OfPnt2d
from OCC.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.GCPnts import GCPnts_QuasiUniformAbscissa

from core_geometry_utils import *
from core_operations_utils import *


def equiPnt_curve2d(curve,NbPoints):
    Adaptor_curve2d = Geom2dAdaptor_Curve(curve.GetHandle())
    first = curve.FirstParameter()
    last = curve.LastParameter()
    equiPara = GCPnts_QuasiUniformAbscissa(Adaptor_curve2d, NbPoints, first,last)
    
    temp = []
    if equiPara.IsDone():
         for j in range(1,NbPoints+1):
             Pnt = curve.Value(equiPara.Parameter(j))
             temp.append([Pnt.X(), Pnt.Y()])
    return np.array(temp)

#TBD: use a maximum deviation do adapt equidistand point placement for refinement in high gradient regions.


def shp_parallel_offset(arrPts,dist):
    side = 'left'    
    afline = shp.LineString(arrPts)
    offline = afline.parallel_offset(dist,side,16,1,1.0)
    offlinepts = np.array(offline.coords)
    #plt.plot(*offlinepts.T, color='red', marker='.')
    return offlinepts



###############################################################################
#GET OUTER AIRFOIL DATA
###################################

if __name__ == '__main__':
      
    
    #def polygon_offset(curve): 
    airfoil = 'naca23012'  
    
    AFURL = 'ah79100c.dat'
    arrPnt = np.loadtxt(AFURL, skiprows=1)   
    
    data = UIUCAirfoil2d(airfoil)                                                       #Get Airfoil Data from Database 
    harray = TColgp_HArray1OfPnt2d_from_nparray(data)                                   #Put Data into Harray1OfPnt2d
    #display_points_of_array(harray)                                                     #Display Points
    anInterpolation = Geom2dAPI_Interpolate(harray.GetHandle(), False, 0.001)           #Interpolate datapoints to bspline
    anInterpolation.Perform()                                               
    bspline = anInterpolation.Curve().GetObject()
    
    #dist = -0.0025           
    arrPnt = equiPnt_curve2d(bspline,2000)
    
    #AFURL = 'ah79100c.dat'
    #arrPnt = np.loadtxt(AFURL, skiprows=1)   
    plt.plot(*arrPnt.T, color='black', marker='.')
    #offlinepts = shp_parallel_offset(arrPnt,dist)
    #plt.plot(*offlinepts.T, color='red', marker='.')
        
    
    
    dist = 0.0015
    num_of_layers = 41
    for j in range(0,num_of_layers):
        print j
        offlinepts = shp_parallel_offset(arrPnt,dist)
        plt.plot(*offlinepts.T, color='red')
        arrPnt = offlinepts
    
    #TBD: fill! => run loop until array  turns to zero!
    
               

#    # Create a Polygon from the nx2 array in `afpts`
#    afpoly = shp.Polygon(arrPnt)
#    
#    # Create offset airfoils, both inward and outward
#    poffafpoly = afpoly.buffer(0.031) # Outward offset
#    noffafpoly = afpoly.buffer(-0.002)  # Inward offset
#    noffafpoly2 = noffafpoly.buffer(-0.002)  # Inward offset
#    noffafpoly3 = noffafpoly2.buffer(-0.002)  # Inward offset
#    noffafpoly4 = noffafpoly3.buffer(-0.002) # Inward offset
#    noffafpoly5 = noffafpoly4.buffer(-0.002,32,1,1)  # Inward offset
#    
#    # Turn polygon points into numpy arrays for plotting
#    afpolypts = np.array(afpoly.exterior)
#    poffafpolypts = np.array(poffafpoly.exterior)
#    noffafpolypts = np.array(noffafpoly.exterior)
#    noffafpolypts2 = np.array(noffafpoly2.exterior)
#    noffafpolypts3 = np.array(noffafpoly3.exterior)
#    noffafpolypts4= np.array(noffafpoly4.exterior)
#    noffafpolypts5= np.array(noffafpoly5.exterior)
#    
#    # Plot points
#    plt.plot(*afpolypts.T, color='black', marker='.')
#    #plt.plot(*poffafpolypts.T, color='red')
#    plt.plot(*noffafpolypts.T, color='green', marker='.')
#    plt.plot(*noffafpolypts2.T, color='green', marker='.')
#    plt.plot(*noffafpolypts3.T, color='green', marker='.')
#    plt.plot(*noffafpolypts4.T, color='green', marker='.')
#    plt.plot(*noffafpolypts5.T, color='green', marker='.')
    plt.axis('equal')
    plt.show()


    


        