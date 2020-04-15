# Basic Libraries:
# Core Library modules
import math
import os

# Third party modules
import matplotlib.pyplot as plt
import numpy as np
from OCC.Core.GCPnts import (GCPnts_AbscissaPoint, GCPnts_QuasiUniformAbscissa,
                             GCPnts_QuasiUniformDeflection,
                             GCPnts_TangentialDeflection,
                             GCPnts_UniformDeflection,)
from OCC.Core.Geom import Geom_BSplineCurve, Geom_Curve
from OCC.Core.Geom2d import (Geom2d_BSplineCurve, Geom2d_Curve, Geom2d_Line,
                             Handle_Geom2d_BSplineCurve_DownCast,)
from OCC.Core.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.Core.Geom2dAPI import (Geom2dAPI_InterCurveCurve,
                                Geom2dAPI_Interpolate,
                                Geom2dAPI_PointsToBSpline,
                                Geom2dAPI_ProjectPointOnCurve,)
from OCC.Core.GeomAdaptor import GeomAdaptor_Curve
from OCC.Core.GeomAPI import (GeomAPI_IntCS, GeomAPI_Interpolate,
                              GeomAPI_ProjectPointOnCurve,)
# PythonOCC Libraries
from OCC.Core.gp import (gp_Dir, gp_Dir2d, gp_Pln, gp_Pnt,
                         gp_Pnt2d, gp_Vec, gp_Vec2d,)
from scipy.optimize import leastsq

# First party modules
from SONATA.cbm.topo.para_Geom2d_BsplineCurve import Para_Geom2d_BSplineCurve
# Own Libraries:
from SONATA.cbm.topo.utils import (P2Pdistance, Pnt2dLst_to_npArray,
                                   PolygonArea,
                                   TColgp_Array1OfPnt2d_from_nparray,
                                   TColgp_HArray1OfPnt2d_from_nparray,
                                   TColgp_HArray1OfPnt_from_nparray,
                                   calc_DCT_angles, curvature_of_curve,
                                   discrete_stepsize, fuse_rows, isclose,
                                   unique_rows,)

###############################################################################
# BSpline and BSplineLst Utilities
###############################################################################


def corners_of_BSplineLst(BSplineLst):
    corners = []
    for item in BSplineLst:
        corners.append(item.EndPoint())

    corners.pop(-1)
    return corners  # gp_Pnt2d Lst


def ProjectPointOnBSplineLst(BSplineLst, Pnt2d, tolerance_distance=100):
    """
    Projets a Pnt2d onto a list of BSplines with given tolerance_distance and 
    returns the a tuple of information of the closest Point
    
    Parameters
    ----------
    BSplineLst : list
        consecutive list of Geom2d_BSplineCurves 
    
    Pnt2d : gp_Pnt2d
        2d Point object of OCC.gp.gp_Pnt2d
    
    Returns
    ----------
    p2 : list 
     list of projection information [gp_Pnt2d, index of BSplineLst, curve parameter u, distance]
         
    
    """
    p2 = []
    for idx, item in enumerate(BSplineLst):
        projection2 = Geom2dAPI_ProjectPointOnCurve(Pnt2d, item)
        for j in range(1, projection2.NbPoints() + 1):
            if projection2.Distance(j) <= tolerance_distance:
                p2.append([projection2.Point(j), idx, projection2.Parameter(j), projection2.Distance(j)])
            else:
                None

    if p2:
        p2 = np.asarray(p2)
        min_index = p2[:, 3].argmin()
        return p2[min_index, :]
    else:
        return []


def distance_on_BSplineLst(BSplineLst, para1, para2):
    #    para1 = findPnt_on_BSplineLst(P1,BSplineLst)
    #    para2 = findPnt_on_BSplineLst(P2,BSplineLst)

    #    if closed:
    #        print para1,para2
    #        para2 = [len(BSplineLst)-1, BSplineLst[-1].LastParameter()]
    #        print para1,para2
    tol = 1e-7
    # print para1,para2

    Distance = 0
    for i, item in enumerate(BSplineLst):
        Adaptor = Geom2dAdaptor_Curve(item)
        First = item.FirstParameter()
        Last = item.LastParameter()
        if para1[0] == i and para2[0] != i:
            Length = GCPnts_AbscissaPoint().Length(Adaptor, para1[1], Last, tol)
            Distance += Length

        elif (para1[0] != i and para2[0] != i) and (para1[0] < i and para2[0] > i):
            Length = GCPnts_AbscissaPoint().Length(Adaptor, First, Last, tol)
            Distance += Length

        elif para1[0] == i and para2[0] == i:
            Length = GCPnts_AbscissaPoint().Length(Adaptor, para1[1], para2[1], tol)
            Distance += Length
            break

        # NOTE: Somehow the execption if para2 is close to Last doesn't work properly!!!!!!! Similar to the trimm function
        elif para1[0] != i and para2[0] == i:
            if isclose(para2[1], Last):
                Length = GCPnts_AbscissaPoint().Length(Adaptor, First, Last, tol)
                Distance += Length
            else:
                Length = GCPnts_AbscissaPoint().Length(Adaptor, First, para2[1], tol)
                Distance += Length
            break

    return Distance


# Obsolete
# def isPnt_on_2dcurve(Pnt2d,Curve2d,tolerance=1e-6):
#    projection = Geom2dAPI_ProjectPointOnCurve(Pnt2d,Curve2d)
#    Trigger = False
#    for i in range(1,projection.NbPoints()+1):
#        if projection.Distance(i) <= tolerance:
#            Trigger = True
#        else:
#            None
#    return Trigger


def isPnt_on_curve(Pnt, Curve, tolerance=1e-6):
    """
    checks if a point (either gp_Pnt or gp_Pnt2d) is on a curve (subclass of 
    Geom2d_Curve or Geom_Curve) by projecting the Point orthogonal on to the 
    curve and checing ifthe distance is below a tolerance tol.
    
    Parameter
    ---------
    Pnt : OCC.gp_Pnt or OCC.gp_Pnt2d
    Curve : OCC.Geom_Curve or OCC.Geom_Curve Handle
    tolerance : float, optional
        default tolerance 1e-6
        
    Return
    --------
    Trigger : bool    
    """
    if isinstance(Pnt, gp_Pnt2d) and issubclass(Curve.__class__, Geom2d_Curve):
        twoD = True
    elif isinstance(Pnt, gp_Pnt) and issubclass(Curve.__class__, Geom_Curve):
        twoD = False

    if twoD:
        projection = Geom2dAPI_ProjectPointOnCurve(Pnt, Curve)
    else:
        projection = GeomAPI_ProjectPointOnCurve(Pnt, Curve)
    Trigger = False
    for i in range(1, projection.NbPoints() + 1):
        if projection.Distance(i) <= tolerance:
            Trigger = True
        else:
            None
    return Trigger


def isPnt_on_BSplineLst(Pnt, BSplineLst, tolerance=1e-6):
    """
    checks if a point (either gp_Pnt or gp_Pnt2d) is on a BSplineLst (list of 
    either OCC.Geom_BSplineCurve or OCC.Geom2d_BSplineCurve) by projecting the 
    Point orthogonal on to the curve and checing if the distance is below a 
    tolerance tol. Using the fucntion 
    'SONATA.cbm.topo.BSplineLst_utils.isPnt_on_curve'.
    
    Parameter
    ---------
    Pnt : OCC.gp_Pnt or OCC.gp_Pnt2d
    BSplineLst : list
         list of either OCC.Geom_BSplineCurve or OCC.Geom2d_BSplineCurve
    tolerance : float, optional
        default tolerance 1e-6
    
    Return
    --------
    Bool : bool    
    """
    Bool = False
    for item in BSplineLst:
        if isPnt_on_curve(Pnt, item, tolerance):
            Bool = True
            break
    return Bool


# OBSOLETE
# def findPnt_on_2dcurve(Pnt2d,Curve2d):
#    projection = Geom2dAPI_ProjectPointOnCurve(Pnt2d,Curve2d)
#    u = projection.LowerDistanceParameter()
#    return u


def findPnt_on_curve(Pnt, Curve):
    if isinstance(Pnt, gp_Pnt2d) and issubclass(Curve.__class__, Geom2d_Curve):
        twoD = True
    elif isinstance(Pnt, gp_Pnt) and issubclass(Curve.__class__, Geom_Curve):
        twoD = False

    if twoD:
        projection = Geom2dAPI_ProjectPointOnCurve(Pnt, Curve)
    else:
        projection = GeomAPI_ProjectPointOnCurve(Pnt, Curve)
    u = projection.LowerDistanceParameter()
    return u


def findPnt_on_BSplineLst(Pnt, BSplineLst):
    for i, item in enumerate(BSplineLst):
        if isPnt_on_curve(Pnt, item):
            u = findPnt_on_curve(Pnt, item)
            coordinates = [i, u]
            break
        else:
            coordinates = None
    return coordinates


def get_BSpline_length(BSpline):
    """
    Returns the Length of the BSpline
    """
    if isinstance(BSpline, Geom_BSplineCurve):
        twoD = False
    elif isinstance(BSpline, Geom2d_BSplineCurve):
        twoD = True

    tolerance = 1e-10
    first = BSpline.FirstParameter()
    last = BSpline.LastParameter()

    if twoD:
        Adaptor = Geom2dAdaptor_Curve(BSpline)
    else:
        Adaptor = GeomAdaptor_Curve(BSpline)

    Length = GCPnts_AbscissaPoint().Length(Adaptor, first, last, tolerance)
    return Length


def copy_BSpline(BSpline):
    BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(BSpline.Copy())
    return BSplineCopy


def copy_BSplineLst(BSplineLst):
    copied_BSplineLst = []
    for i, item in enumerate(BSplineLst):
        copied_BSplineLst.append(copy_BSpline(item))
    return copied_BSplineLst


def get_BSplineLst_length(BSplineLst):
    """
    Returns the cummulated length of the BSplineLst
    """

    CummLength = 0
    for i, item in enumerate(BSplineLst):
        CummLength += get_BSpline_length(item)
    return CummLength


def equidistant_Points_on_BSplineLst(BSplineLst, minLen):
    """ minLen = the minimum distance between points. Should be determined by the segment0.BSplineLstLenght/Resolution
    """
    Pnt2dLst = []
    for item in BSplineLst:
        Adaptor = Geom2dAdaptor_Curve(item)
        length = get_BSpline_length(item)
        NbPoints = int(length // minLen) + 2
        discretization = GCPnts_QuasiUniformAbscissa(Adaptor, NbPoints)

        for j in range(1, NbPoints + 1):
            para = discretization.Parameter(j)
            Pnt = gp_Pnt2d()
            item.D0(para, Pnt)
            Pnt2dLst.append(Pnt)

    return Pnt2dLst


def equidistant_D1_on_BSplineLst(BSplineLst, NbPoints):
    """
    distributes NbPoints number of Points and NormalVectors onto equidistantly on 
    Geom2d_BSplineLst. 
    
    """
    PntLst = []
    VecLst = []
    for item in BSplineLst:

        if isinstance(item, Geom_BSplineCurve):
            print("ERROR: Argument must be a Geom2d_BSplineLst")

        Adaptor = Geom2dAdaptor_Curve(item)
        discretization = GCPnts_QuasiUniformAbscissa(Adaptor, NbPoints)

        for j in range(1, NbPoints + 1):
            para = discretization.Parameter(j)
            # print(para)
            Pnt2d = gp_Pnt2d()
            Vec2d = gp_Vec2d()
            item.D1(para, Pnt2d, Vec2d)
            PntLst.append(Pnt2d)
            VecLst.append(Vec2d.GetNormal())

    return (PntLst, VecLst)


def find_BSpline_coordinate(BSpline, s):
    # Be careful, s stands for the lenght coordinate of a single BSpline, while S represents the Global Coordinate!
    BSpline_length = get_BSpline_length(BSpline)
    tolerance = 1e-10
    Adaptor = Geom2dAdaptor_Curve(BSpline)
    if s <= BSpline_length:
        tmp = GCPnts_AbscissaPoint(tolerance, Adaptor, s, 0)
        u = tmp.Parameter()
    return u


def find_BSplineLst_coordinate(BSplineLst, S, start, end):
    """The function find_BSplineLst_coordinate returns the list index of 
    the bspline and its parameter where the coordinate is located.
    
    Parameters
    ----------
    S : float
        Coordinate to be found 
    start : float
        start of the interval where the BSplineLst is defined
    end : float
        end of the interval where the BSplineLst is defined
    
    Returns
    ----------
    [i,U]: list
        Return [index of the bspline, parameter U on the bspline]
        
    """

    if all(isinstance(s, Geom_BSplineCurve) for s in BSplineLst):
        twoD = False
    elif all(isinstance(s, Geom2d_BSplineCurve) for s in BSplineLst):
        twoD = True

    BSplineLstLength = get_BSplineLst_length(BSplineLst)
    x = BSplineLstLength * (S - start) / (end - start)
    CummLength = 0
    tolerance = 1e-8
    # print 'S: '+ str(S)
    for i, item in enumerate(BSplineLst):
        first = item.FirstParameter()
        last = item.LastParameter()

        if twoD:
            Adaptor = Geom2dAdaptor_Curve(item)
        else:
            Adaptor = GeomAdaptor_Curve(item)

        CummLength += get_BSpline_length(item)
        if x < CummLength or isclose(x, CummLength):
            dist = x - (CummLength - GCPnts_AbscissaPoint().Length(Adaptor, first, last, tolerance))
            tmp = GCPnts_AbscissaPoint(tolerance, Adaptor, dist, first)
            U = tmp.Parameter()
            break
    # print 'idx: '+ str(i) + '     U: '+ str(U)
    return [i, U]


def find_BSplineLst_pos(BSplineLst, para):
    idx = para[0]
    U = para[1]
    tol = 1e-7
    BSplineLst_length = get_BSplineLst_length(BSplineLst)
    CummLength = 0
    for x in range(0, idx):
        # print 'x',x
        CummLength += get_BSpline_length(BSplineLst[x])

    Adaptor = Geom2dAdaptor_Curve(BSplineLst[idx])
    First = BSplineLst[idx].FirstParameter()
    partial_length = GCPnts_AbscissaPoint().Length(Adaptor, First, U, tol)
    S = (CummLength + partial_length) / float(BSplineLst_length)
    return S


def reverse_BSplineLst(BSplineLst):
    BSplineLst = list(reversed(BSplineLst))
    for i, item in enumerate(BSplineLst):
        item.Reverse()
    return BSplineLst


def get_BSplineLst_Pnt2d(BSplineLst, S, start, end):
    if start > end:  # transfrom the interval if start>end to a basic interval [0..]
        D = 1 - start
        S = (S + D) % 1
        end = (end + D) % 1
        start = (start + D) % 1

    P = gp_Pnt2d()
    [idx, U] = find_BSplineLst_coordinate(BSplineLst, S, start, end)
    BSplineLst[idx].D0(U, P)
    return P


def get_BSplineLst_D2(BSplineLst, S, start, end, **kwargs):
    if start > end:  # transfrom the interval if start>end to a basic interval [0..]
        D = 1 - start
        S = (S + D) % 1
        end = (end + D) % 1
        start = (start + D) % 1

    if all(isinstance(s, Geom_BSplineCurve) for s in BSplineLst):
        twoD = False
    elif all(isinstance(s, Geom2d_BSplineCurve) for s in BSplineLst):
        twoD = True

    if twoD:
        P = gp_Pnt2d()
        V1 = gp_Vec2d()
        V2 = gp_Vec2d()
    else:
        P = gp_Pnt()
        V1 = gp_Vec()
        V2 = gp_Vec()

    [idx, U] = find_BSplineLst_coordinate(BSplineLst, S, start, end, **kwargs)
    BSplineLst[idx].D2(U, P, V1, V2)
    return P, V1, V2


def intersect_BSplineLst_with_BSpline(BSplineLst, BSpline, tolerance=1e-10):
    IntPnts = []
    IntPnts_Pnt2d = []
    for i, item in enumerate(BSplineLst):
        intersection = Geom2dAPI_InterCurveCurve(BSpline, item, tolerance)
        for j in range(1, intersection.NbPoints() + 1):
            IntPnt = intersection.Point(j)
            IntPnts_Pnt2d.append(IntPnt)
            u = findPnt_on_curve(IntPnt, item)
            IntPnts.append([i, u])
    return IntPnts, IntPnts_Pnt2d


def intersect_BSplineLst_with_plane(BSplineLst, plane):
    """
    intersects a 3D-BSplineLst (OCC.Geom_BSplineCurve) with a plane (Handle<Geom_Surface>) and returns the
    a list of intersection coordinates and the corresponding list of 
    intersection points
    
    Paramters
    ---------
    BSplineLst : list
         list of OCC.Geom_BSplineCurve
    plane : <Geom_Surface>
         Plane or surface - Geom_Surface, Geom_Plane is a subclass of Geom_Surface
    
    Return
    --------
    IntCoords : list
        list of intersection coordinates
        [index of the bspline, parameter U on the bspline]
    IntPnts : list 
        list of intersection points (OCC.gp_Pnt)
    """
    # plane = Geom_Plane(gp_Pnt(float(x),0,0), gp_Dir(1,0,0))
    IntCoords = []
    IntPnts = []
    for i, item in enumerate(BSplineLst):
        intersection = GeomAPI_IntCS(item, plane)
        for j in range(1, intersection.NbPoints() + 1):
            IntPnt = intersection.Point(j)
            IntPnts.append(IntPnt)
            u = findPnt_on_curve(IntPnt, item)
            IntCoords.append([i, u])

    return IntCoords, IntPnts


def BSplineLst_Orientation(BSplineLst, NbPoints=31):
    
    if isinstance(BSplineLst[0], Geom_BSplineCurve):
        twoD = False
    elif isinstance(BSplineLst[0], Geom2d_BSplineCurve):
        twoD = True    

    if twoD:
         P = gp_Pnt2d()
         V1 = gp_Vec2d()
         V2 = gp_Vec2d()
    else:
        P = gp_Pnt()
        V1 = gp_Vec()
        V2 = gp_Vec()

    # FASTER BUT NOT EQUIDISTANT APPROACH!
    Pnt2dLst = []
    for item in BSplineLst:
        First = item.FirstParameter()
        Pnt2dLst.append(item.StartPoint())
        Last = item.LastParameter()
        inc = (Last - First) / NbPoints
        para = First + inc
        while para <= Last:
            if twoD:
                p = gp_Pnt2d()
            else:
                p = gp_Pnt()
            
            item.D0(para, p)
            Pnt2dLst.append(p)
            para += inc
            
            
    b = Pnt2dLst_to_npArray(Pnt2dLst)

    #   Place Equidistant Points on BSplineLst:                             #THIS APPROACH WAS SUPER SLOW!
    #    a = np.linspace(0.0, 1.0, num=NbPoints, endpoint=True)
    #    Pnt2dLst = []
    #    for s in a:
    #        Pnt2dLst.append(get_BSplineLst_Pnt2d(BSplineLst,s, 0, 1))
    #    b = Pnt2dLst_to_npArray(Pnt2dLst)

    # Calculate of Polygon.
    area = PolygonArea(b)
    # print area
    if area > 0:  # counter-clockwise=True
        return True
    elif area < 0:  # clockwise=False
        return False
    else:
        return None


def trim_BSplineLst(BSplineLst, S1, S2, start, end):
    """ the trim_BSplineLst function trims the BSplineLst that is defined on 
    the interval between start and end to the interval defined by S1 and S2.
    
    Parameters
    ----------
    BSplineLst: list of geom2d_bsplines
        is the object to be modified.
    S1: float
        start of the interval to which the object is trimmed.
    S2: float
        end of the interval to which the object is trimmed.
    start: float
        start of the interval where the BSplineLst is defined
    end: float
        end of the interval where the BSplineLst is defined
            
    Returns
    ----------
    trimmed_BSplineLst: list of geom2d_bsplines
        is the modified object.
            
    Example
    ----------
    a BSplineLst is defined between 0.0 and 1.0 and shall be trimmed 
    to a interval S1=0.3 and S2=0.7
    """

    if start > end:  # transfrom the interval if start>end to a basic interval [0..]
        D = 1 - start
        S1 = (S1 + D) % 1
        S2 = (S2 + D) % 1
        end = (end + D) % 1
        start = (start + D) % 1

    if S1 > S2:
        trimmed_BSplineLst = []
        front_BSplineLst = []
        rear_BSplineLst = []
        para1 = find_BSplineLst_coordinate(BSplineLst, S1, start, end)
        para2 = find_BSplineLst_coordinate(BSplineLst, S2, start, end)

        # check if para1 and para2 is out of bounds!
        if isclose(para1[1], BSplineLst[para1[0]].FirstParameter()):
            para1[1] = BSplineLst[para1[0]].FirstParameter()
        if isclose(para2[1], BSplineLst[para2[0]].LastParameter()):
            para2[1] = BSplineLst[para2[0]].LastParameter()

        for i, item in enumerate(BSplineLst):
            First = item.FirstParameter()
            Last = item.LastParameter()
            if para1[0] == i and para2[0] != i:  # Okay
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                if isclose(para1[1], Last):
                    pass
                elif isclose(para1[1], First):
                    BSplineCopy.Segment(First, Last)
                    trimmed_BSplineLst.append(BSplineCopy)
                else:
                    BSplineCopy.Segment(para1[1], Last)
                    trimmed_BSplineLst.append(BSplineCopy)

            elif (para1[0] != i and para2[0] != i) and (para1[0] > i and para2[0] > i):
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                BSplineCopy.Segment(First, Last)
                front_BSplineLst.append(BSplineCopy)

            elif (para1[0] != i and para2[0] != i) and (para1[0] < i and para2[0] < i):
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                BSplineCopy.Segment(First, Last)
                rear_BSplineLst.append(BSplineCopy)

            elif para1[0] == i and para2[0] == i:  # Okay
                BSplineCopy1 = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                BSplineCopy1.Segment(para1[1], Last)
                trimmed_BSplineLst.append(BSplineCopy1)
                BSplineCopy2 = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                BSplineCopy2.Segment(First, para2[1])
                trimmed_BSplineLst.append(BSplineCopy2)
                break

            elif para1[0] != i and para2[0] == i:  # Okay
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                BSplineCopy.Segment(First, para2[1])
                front_BSplineLst.append(BSplineCopy)
                # break

        if front_BSplineLst and rear_BSplineLst:
            trimmed_BSplineLst = rear_BSplineLst + front_BSplineLst

        elif front_BSplineLst and not rear_BSplineLst:
            trimmed_BSplineLst = front_BSplineLst

        elif not front_BSplineLst and rear_BSplineLst:
            trimmed_BSplineLst = rear_BSplineLst

    elif S1 == start and S2 == end:
        trimmed_BSplineLst = copy_BSplineLst(BSplineLst)

    elif S2 == start or S1 == S2:
        trimmed_BSplineLst = []

    elif S2 > S1:
        trimmed_BSplineLst = []

        para1 = find_BSplineLst_coordinate(BSplineLst, S1, start, end)
        para2 = find_BSplineLst_coordinate(BSplineLst, S2, start, end)

        # check if para1 and para2 is out of bounds!
        if isclose(para1[1], BSplineLst[para1[0]].FirstParameter()):
            para1[1] = BSplineLst[para1[0]].FirstParameter()
        if isclose(para2[1], BSplineLst[para2[0]].LastParameter()):
            para2[1] = BSplineLst[para2[0]].LastParameter()

        for i, item in enumerate(BSplineLst):
            First = item.FirstParameter()
            Last = item.LastParameter()
            if para1[0] == i and para2[0] != i:
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                if isclose(para1[1], Last):
                    pass
                elif isclose(para1[1], First):
                    BSplineCopy.Segment(First, Last)
                    trimmed_BSplineLst.append(BSplineCopy)
                else:
                    BSplineCopy.Segment(para1[1], Last)
                    trimmed_BSplineLst.append(BSplineCopy)

            elif (para1[0] != i and para2[0] != i) and (para1[0] < i and para2[0] > i):
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                BSplineCopy.Segment(First, Last)
                trimmed_BSplineLst.append(BSplineCopy)

            elif para1[0] == i and para2[0] == i:
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                BSplineCopy.Segment(para1[1], para2[1])
                trimmed_BSplineLst.append(BSplineCopy)
                break

            elif para1[0] != i and para2[0] == i:
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                if isclose(para2[1], Last):
                    trimmed_BSplineLst.append(BSplineCopy)
                else:
                    BSplineCopy.Segment(First, para2[1])
                    trimmed_BSplineLst.append(BSplineCopy)
                break

    return trimmed_BSplineLst


def trim_BSplineLst_by_coordinates(BSplineLst, para1, para2):
    trimmed_BSplineLst = []
    #    para1 =  find_BSplineLst_coordinate(BSplineLst, S1, start, end)
    #    para2 =  find_BSplineLst_coordinate(BSplineLst, S2, start, end)
    for i, item in enumerate(BSplineLst):
        First = item.FirstParameter()
        Last = item.LastParameter()
        if para1[0] == i and para2[0] != i:
            BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
            BSplineCopy.Segment(para1[1], Last)
            trimmed_BSplineLst.append(BSplineCopy)

        elif (para1[0] != i and para2[0] != i) and (para1[0] < i and para2[0] > i):
            BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
            BSplineCopy.Segment(First, Last)
            trimmed_BSplineLst.append(BSplineCopy)

        elif para1[0] == i and para2[0] == i:
            BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
            BSplineCopy.Segment(para1[1], para2[1])
            trimmed_BSplineLst.append(BSplineCopy)
            break

        # NOTE: Somehow the execption if para2 is close to Last doesn't work properly!!!!!!!
        elif para1[0] != i and para2[0] == i:
            BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
            if isclose(para2[1], Last):
                trimmed_BSplineLst.append(BSplineCopy)
            else:
                BSplineCopy.Segment(First, para2[1])
                trimmed_BSplineLst.append(BSplineCopy)
            break

    #             elif para1[0] != i and para2[0] == i:
    #                 BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
    #                 BSplineCopy.Segment(First,para2[1])
    #                 trimmed_BSplineLst.append(BSplineCopy)
    #                 break
    #

    return trimmed_BSplineLst


def trim_BSplineLst_by_Pnt2d(BSplineLst, Pos1_Pnt2d, Pos2_Pnt2d):
    trimmed_BSplineLst = []
    front_BSplineLst = []
    rear_BSplineLst = []
    para1 = findPnt_on_BSplineLst(Pos1_Pnt2d, BSplineLst)
    para2 = findPnt_on_BSplineLst(Pos2_Pnt2d, BSplineLst)

    # if para2>para1
    if (para2[0] > para1[0]) or ((para2[0] == para1[0]) and (para2[1] > para1[1])):
        for i, item in enumerate(BSplineLst):
            First = item.FirstParameter()
            Last = item.LastParameter()
            if para1[0] == i and para2[0] != i:
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                BSplineCopy.Segment(para1[1], Last)
                trimmed_BSplineLst.append(BSplineCopy)

            elif (para1[0] != i and para2[0] != i) and (para1[0] < i and para2[0] > i):
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                BSplineCopy.Segment(First, Last)
                trimmed_BSplineLst.append(BSplineCopy)

            elif para1[0] == i and para2[0] == i:
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                BSplineCopy.Segment(para1[1], para2[1])
                trimmed_BSplineLst.append(BSplineCopy)
                break

            elif para1[0] != i and para2[0] == i:
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                BSplineCopy.Segment(First, para2[1])
                trimmed_BSplineLst.append(BSplineCopy)
                break

    # if para1>para2
    else:
        # para1, para2 = para2, para1
        for i, item in enumerate(BSplineLst):
            First = item.FirstParameter()
            Last = item.LastParameter()
            if para1[0] == i and para2[0] != i:  # Okay
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                BSplineCopy.Segment(para1[1], Last)
                rear_BSplineLst.append(BSplineCopy)

            elif (para1[0] != i and para2[0] != i) and (para1[0] > i and para2[0] > i):
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                BSplineCopy.Segment(First, Last)
                front_BSplineLst.append(BSplineCopy)

            elif (para1[0] != i and para2[0] != i) and (para1[0] < i and para2[0] < i):
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                BSplineCopy.Segment(First, Last)
                rear_BSplineLst.append(BSplineCopy)

            elif para1[0] == i and para2[0] == i:  # Okay
                BSplineCopy1 = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                BSplineCopy1.Segment(para1[1], Last)
                trimmed_BSplineLst.append(BSplineCopy1)
                BSplineCopy2 = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                BSplineCopy2.Segment(First, para2[1])
                trimmed_BSplineLst.append(BSplineCopy2)
                break

            elif para1[0] != i and para2[0] == i:  # Okay
                BSplineCopy = Handle_Geom2d_BSplineCurve_DownCast(item.Copy())
                BSplineCopy.Segment(First, para2[1])
                front_BSplineLst.append(BSplineCopy)
                # break

        if front_BSplineLst and rear_BSplineLst:
            trimmed_BSplineLst = rear_BSplineLst + front_BSplineLst

        elif front_BSplineLst and not rear_BSplineLst:
            trimmed_BSplineLst = front_BSplineLst

        elif not front_BSplineLst and rear_BSplineLst:
            trimmed_BSplineLst = rear_BSplineLst

    return trimmed_BSplineLst


def seg_boundary_from_dct(DCT_data, angular_deflection=30):
    # OUTDATED!!!
    # Check if DCT_Definition is closed, if not: close it
    if not np.allclose(DCT_data[0], DCT_data[-1]):
        print("INFO:\t Closing open discrete boundary definition")
        DCT_data = np.concatenate((DCT_data, DCT_data[0:1, :]), axis=0)

    # Find corners and edges of data
    DCT_angles = calc_DCT_angles(DCT_data)
    print(DCT_data.shape, DCT_angles.shape)
    corners = []
    for i in range(0, DCT_angles.shape[0]):
        if DCT_angles[i] < (180 - angular_deflection) or DCT_angles[i] > (180 + angular_deflection):
            corners.append(i)
    NbCorners = np.size(corners)

    # Segment by Corners
    print("NbCorners = ", NbCorners)
    DCT_Segments = []
    if NbCorners == 0:
        DCT_Segments.append(DCT_data)
    if NbCorners > 0:
        DCT_Segments.append(DCT_data[0 : corners[0] + 1])
        for i in range(0, NbCorners - 1):
            print("i, corners = ", i, corners[i])
            DCT_Segments.append(DCT_data[corners[i] : corners[i + 1] + 1])
        DCT_Segments.append(DCT_data[corners[-1] : -1])

    list_of_bsplines = []
    for i, item in enumerate(DCT_Segments):
        data = item.T

        tmp_harray = TColgp_HArray1OfPnt2d_from_nparray(data)
        if NbCorners == 0:
            # tmp_interpolation = Geom2dAPI_Interpolate(tmp_harray, True, 1e-6)             #Interpolate datapoints to bspline
            tmp_interpolation = Geom2dAPI_Interpolate(tmp_harray, False, 1e-6)  # Interpolate datapoints to bspline
        else:
            tmp_interpolation = Geom2dAPI_Interpolate(tmp_harray, False, 1e-4)  # Interpolate datapoints to bspline
        tmp_interpolation.Perform()
        tmp_bspline = tmp_interpolation.Curve()
        list_of_bsplines.append(tmp_bspline)

    return list_of_bsplines


def discretize_BSplineLst(BSplineLst, Deflection=2e-4, AngularDeflection=0.02, CurvatureDeflection=0.06, MinimumOfPoints=50):
    Pnt2dLst = []
    for i, item in enumerate(BSplineLst):
        Adaptor = Geom2dAdaptor_Curve(item)
        discretization = GCPnts_QuasiUniformDeflection(Adaptor, Deflection, 1)  # GeomAbs_Shape Continuity: 1=C0, 2=G1, 3=C1, 3=G2,...
        # discretization = GCPnts_TangentialDeflection(Adaptor,AngularDeflection,CurvatureDeflection,MinimumOfPoints)
        NbPoints = discretization.NbPoints()
        for j in range(1, NbPoints + 1):
            Pnt = discretization.Value(j)
            Pnt2dLst.append(Pnt)

    a = Pnt2dLst_to_npArray(Pnt2dLst)
    b = np.around(a, 10)  # Evenly round to the given number of decimals.

    # check if it is closed:
    if np.array_equal(b[0], b[-1]):
        closed = True
    else:
        closed = False

    npArray = unique_rows(b)  # Remove possible doubles!

    if closed:
        npArray = np.vstack((npArray, b[0]))
    else:
        None
    # print len(npArray)

    #    #new Interpolate Large spaces!
    #    refine = True
    #    Resolution = 1000
    #    seg_P2Plength = []
    #    cumm_length = 0
    #    data = npArray
    #
    #    for j in range(0,len(data)-1):
    #        seg_P2Plength.append(P2Pdistance(data[j],data[j+1]))
    #        cumm_length += P2Pdistance(data[j],data[j+1])
    #
    #    #Check if Refinement is necessary:
    #    if len(seg_P2Plength) > 0 and max(seg_P2Plength) > cumm_length/Resolution :
    #        Refinement = True
    #    else:
    #        Refinement = False
    #
    #
    #    while Refinement == True:
    #        print(Refinement)
    #
    #        temp_data = []
    #        for i in range(0,len(data)-1):
    #            if P2Pdistance(data[i],data[i+1]) > (cumm_length/Resolution):
    #                p0 = data[i]
    #                p1 = data[i+1]
    #                v1 = p1-p0
    #                p05 = p0+v1/2
    #                temp_data.append(p0)
    #                temp_data.append(p05)
    #            else:
    #                temp_data.append(data[i])
    #
    #        temp_data.append(data[-1])
    #        data = np.vstack(temp_data)
    #
    #        #Check if further Refinement is necessary
    #        seg_P2Plength = []
    #        cumm_length = 0
    #        for j in range(0,len(data)-1):
    #            seg_P2Plength.append(P2Pdistance(data[j],data[j+1]))
    #            cumm_length += P2Pdistance(data[j],data[j+1])
    #
    #        if max(seg_P2Plength) > cumm_length/Resolution:
    #            Refinement = True
    #        else:
    #            Refinement = False
    return npArray



# return data


def BSplineLst_from_dct(DCT_data, angular_deflection=15, closed=False, tol_interp=1e-6, twoD=True):
    # Close DCT_date if closed == True
    if closed and not np.allclose(DCT_data[0], DCT_data[-1]):
        print("INFO:\t Closing open discrete definition")
        DCT_data = np.concatenate((DCT_data, DCT_data[0:1, :]), axis=0)

    DCT_data = fuse_rows(DCT_data, 1e-6)  # check all datapoints and merge if allclose is true

    # Find corners and edges of data
    # print(DCT_data)
    DCT_angles = calc_DCT_angles(DCT_data)
    # print(DCT_angles)
    corners = []
    for i in range(0, DCT_angles.shape[0]):
        if DCT_angles[i] < (180 - angular_deflection) or DCT_angles[i] > (180 + angular_deflection):
            corners.append(i)
    NbCorners = np.size(corners)

    # Segmenting data according to corners
    # ====================================
    DCT_Segments = []

    # Closed:
    if np.allclose(DCT_data[0], DCT_data[-1]):
        closed = True
        if NbCorners == 0:
            DCT_Segments.append(DCT_data)

        elif NbCorners > 0:

            if corners[0] == 0 and corners[-1] == len(DCT_data) - 1:
                for i in range(0, NbCorners - 2):
                    DCT_Segments.append(DCT_data[corners[i] : corners[i + 1] + 1])
                DCT_Segments.append(DCT_data[corners[-2] : corners[-1] + 1])  # last segment

            # regular case:
            else:
                DCT_Segments.append(DCT_data[: corners[0] + 1])  # first
                if NbCorners > 1:  # intermediate segments:
                    for i in range(0, NbCorners - 1):
                        DCT_Segments.append(DCT_data[corners[i] : corners[i + 1] + 1])
                DCT_Segments.append(DCT_data[corners[-1] :])  # last
                DCT_Segments[0] = np.concatenate((DCT_Segments[-1][:-1], DCT_Segments[0]))
                del DCT_Segments[-1]

    else:  # Open:
        closed = False
        if NbCorners == 0:
            DCT_Segments.append(DCT_data)

        if NbCorners > 0:
            DCT_Segments.append(DCT_data[: corners[0] + 1])  # first
            if NbCorners > 1:  # intermediate segments:
                for i in range(0, NbCorners - 1):
                    DCT_Segments.append(DCT_data[corners[i] : corners[i + 1] + 1])
            DCT_Segments.append(DCT_data[corners[-1] :])  # last segment

    list_of_bsplines = []
    for i, item in enumerate(DCT_Segments):
        data = item.T

        if twoD:
            tmp_harray = TColgp_HArray1OfPnt2d_from_nparray(data)
        else:
            tmp_harray = TColgp_HArray1OfPnt_from_nparray(data)

        try:
            if twoD:
                tmp_interpolation = Geom2dAPI_Interpolate(tmp_harray, False, tol_interp)  # Interpolate datapoints to bspline
            else:
                tmp_interpolation = GeomAPI_Interpolate(tmp_harray, False, tol_interp)
        except:
            print("ERROR: BSplineLst_from_dct did not perform the Interpolation with tol_interp", tol_interp)  # RAISED ERROR
            plt.figure()
            plt.clf()
            plt.axis("equal")
            plt.plot(*item.T, marker=".")
            for i, it in enumerate(item):
                plt.annotate(i, (it[0], it[1]), color="black")
            plt.show()

        tmp_interpolation.Perform()
        tmp_bspline = tmp_interpolation.Curve()
        # print(tmp_bspline.Degree(), tmp_bspline.NbKnots())
        list_of_bsplines.append(tmp_bspline)

    return list_of_bsplines


def set_BSplineLst_to_Origin(BSplineLst, Theta=0):
    """
    The Origin is determined by the most right Intersection Point of the X-Axis
    with the segment boundary, the axis is rotated backwards by Theta (in degree). 
    """

    Theta = Theta / float(180) * np.pi
    x = math.cos(Theta)
    y = -math.sin(Theta)
    OPnt = gp_Pnt2d(0, 0)
    dir2d = gp_Dir2d(x, y)
    axis = Geom2d_Line(OPnt, dir2d)  # x-axis
    v0 = gp_Vec2d(dir2d)
    tolerance = 1e-10
    IntPnts = []

    for i, item in enumerate(BSplineLst):
        First = item.FirstParameter()
        Last = item.LastParameter()
        intersection = Geom2dAPI_InterCurveCurve(axis, item, tolerance)
        for j in range(1, intersection.NbPoints() + 1):
            IntPnt = intersection.Point(j)
            v1 = gp_Vec2d(OPnt, IntPnt)
            vres = v0.Dot(v1)
            # distance = IntPnt.Distance(OPnt)
            u = findPnt_on_curve(IntPnt, item)
            IntPnts.append([i, u, vres])

    # Determine Origin as point
    IntPntsarray = unique_rows(np.asarray(IntPnts))  # idx,W,XValue
    OriEdgePnt = IntPntsarray[np.argmax(IntPntsarray[:, 2]), :]
    # print(IntPntsarray,OriEdgePnt)

    # Reorder Sequence of BSplines of BSplinesLst
    OBSplineLst = []
    CorrectOrigin = False

    for i, item in enumerate(BSplineLst):
        if i == OriEdgePnt[0]:
            First = item.FirstParameter()
            Last = item.LastParameter()

            if isclose(OriEdgePnt[1], First) == True:
                OBSplineLst.append(item)
                CorrectOrigin = True

            elif isclose(OriEdgePnt[1], Last) == True:
                CorrectOrigin = False
                BSplineCurve2 = item

            else:
                CorrectOrigin = False
                BSplineCurve1 = copy_BSpline(item)
                BSplineCurve1.Segment(OriEdgePnt[1], Last)
                BSplineCurve2 = copy_BSpline(item)
                BSplineCurve2.Segment(First, OriEdgePnt[1])
                OBSplineLst.append(BSplineCurve1)

        elif i > OriEdgePnt[0]:
            OBSplineLst.append(item)
        else:
            None

    for i, item in enumerate(BSplineLst):
        if i < OriEdgePnt[0]:
            OBSplineLst.append(item)
        else:
            None

    if CorrectOrigin == False:
        OBSplineLst.append(BSplineCurve2)
    else:
        None

    return OBSplineLst


def set_BSplineLst_to_Origin2(BSplineLst, gp_Pnt2d, tol=1e-3):
    """
    this procedure reorders the self.BSplineLst to and origin if the layer 
    is closed. The Origin is detected by searching for an orthogonal 
    projection of the StartPoint of the self.Boundary_BSplineLst. If no 
    projection is found it takes the closest neighbor of the discrete 
    offlinepts (Offset Line Points).
            
    """
    # print(BSplineLst, gp_Pnt2d)
    # Determine Origin as point
    if BSplineLst[0].StartPoint().IsEqual(BSplineLst[-1].EndPoint(), 1e-5):
        proj = ProjectPointOnBSplineLst(BSplineLst, gp_Pnt2d, tol)

        if len(proj) > 0:
            OriPnt = proj[0]
        elif len(proj) == 0:
            print("WARNING : No Projection found to produce origin!")

    else:
        print("ERRROR: Only apply set_BSplineLst_to_Origin2 function when BSplineLst is Closed!")

    # Reorder Sequence of BSplines of BSplinesLst
    OriPara = findPnt_on_BSplineLst(OriPnt, BSplineLst)
    OBSplineLst = []
    CorrectOrigin = False

    for i, item in enumerate(BSplineLst):
        if i == OriPara[0]:
            First = item.FirstParameter()
            Last = item.LastParameter()

            if isclose(OriPara[1], First) == True:
                OBSplineLst.append(item)
                CorrectOrigin = True

            elif isclose(OriPara[1], Last) == True:
                CorrectOrigin = False
                BSplineCurve2 = item

            else:
                CorrectOrigin = False
                BSplineCurve1 = copy_BSpline(item)
                BSplineCurve1.Segment(OriPara[1], Last)
                BSplineCurve2 = copy_BSpline(item)
                BSplineCurve2.Segment(First, OriPara[1])
                OBSplineLst.append(BSplineCurve1)

        elif i > OriPara[0]:
            OBSplineLst.append(item)
        else:
            None

    for i, item in enumerate(BSplineLst):
        if i < OriPara[0]:
            OBSplineLst.append(item)
        else:
            None

    if CorrectOrigin == False:
        OBSplineLst.append(BSplineCurve2)
    else:
        None

    BSplineLst = OBSplineLst

    return BSplineLst


if __name__ == "__main__":
    pass
