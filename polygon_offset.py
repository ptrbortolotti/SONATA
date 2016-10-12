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
from OCC.Display.SimpleGui import init_display

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


def unit_vector(vector):
    return vector / np.linalg.norm(vector)
    
        
def angle_between(v1, v2):
    return np.degrees( math.atan2(np.linalg.norm(np.cross(v1,v2)), np.dot(v1,v2)))

###############################################################################
#GET OUTER AIRFOIL DATA
###################################

def shp_discreteoffset(curve,dist,NbPoints):
    arrPnt = equiPnt_curve2d(curve,NbPoints)
    #print arrPnt
    offlinepts = shp_parallel_offset(arrPnt,dist)
    
    
    #TBD: if there is a edge! split array! and construct two or more curves by interpolation!    
    
    #print offlinepts
    harray = TColgp_HArray1OfPnt2d_from_nparray(np.transpose(offlinepts))
    #print harray    
    #print harray.Length()
    anInterpolation = Geom2dAPI_Interpolate(harray.GetHandle(), False, 0.000001)           #Interpolate datapoints to bspline
    anInterpolation.Perform()
    #print anInterpolation.IsDone()                                               
    bspline = anInterpolation.Curve().GetObject()
    return bspline



if __name__ == '__main__':
###############################################################################
    
    #def polygon_offset(curve): 
    airfoil = 'naca23012'  
    

    
    data = UIUCAirfoil2d(airfoil)                                                       #Get Airfoil Data from Database 
    harray = TColgp_HArray1OfPnt2d_from_nparray(data)                                   #Put Data into Harray1OfPnt2d
    #display_points_of_array(harray)                                                     #Display Points
    anInterpolation = Geom2dAPI_Interpolate(harray.GetHandle(), False, 0.001)           #Interpolate datapoints to bspline
    anInterpolation.Perform()                                               
    bspline = anInterpolation.Curve().GetObject()
    
    dist = 0.0013
    NbPoints = 2000

    
    
    arrPnt = equiPnt_curve2d(bspline,1000)
    plt.plot(*arrPnt.T, color='black', marker='.')
        
    
    test = shp_discreteoffset(bspline,dist,NbPoints)
    
    #AFURL = 'ah79100c.dat'
    #arrPnt = np.loadtxt(AFURL, skiprows=1)   
    
    dist = 0.008
    num_of_layers = 8
    for j in range(0,num_of_layers):
        print j
        offlinepts = shp_parallel_offset(arrPnt,dist)
        plt.plot(*offlinepts.T, color='red', marker='.')
        arrPnt = offlinepts
#    
    #TBD: fill! => run loop until array  turns to zero!
    
    #offlinepts = np.array([[0,0],[1,0],[2,1],[3,2],[5,1]])
    
    for i in range(1,offlinepts.shape[0]-1):
        v1 = offlinepts[i-1]-offlinepts[i]
        v2 = offlinepts[i+1]-offlinepts[i]
        print i, angle_between(v1,v2)
        if abs(angle_between(v1,v2)) < 150:
            print ("edge:", i, angle_between(v1,v2))
            plt.plot(*offlinepts[i].T, color='green',marker='x')
     
    #CHECK IF POLYLINE IS A CLOSED CONTOUR AND IF THAT IS THE CASE, CHECK THE ANGLE!
     
    
    #plt.plot(*offlinepts[0].T, color='BLUE',marker='>')
    #plt.plot(*offlinepts[-1].T, color='CYAN',marker='s')
    #vec1 = offlinepts[i]-offlinepts[0]         
        #vec2 = offlinepts[2]-offlinepts[1]
    

    
    #angle_between(vec1,vec2)
    #print angle_between(vec1,vec2)*180/np.pi
    
    
    
    
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
        



    


        