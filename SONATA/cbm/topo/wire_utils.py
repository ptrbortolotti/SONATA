# Basic Libraries:
# Third party modules
import matplotlib.pyplot as plt
import numpy as np
from OCC.Core.BRepAdaptor import BRepAdaptor_CompCurve, BRepAdaptor_Curve
from OCC.Core.BRepBuilderAPI import (BRepBuilderAPI_MakeEdge,
                                     BRepBuilderAPI_MakeWire,
                                     BRepBuilderAPI_Transform,)
from OCC.Core.BRepLib import BRepLib_MakeFace
from OCC.Core.BRepTools import BRepTools_WireExplorer
from OCC.Core.GCPnts import (GCPnts_AbscissaPoint, GCPnts_QuasiUniformAbscissa,
                             GCPnts_QuasiUniformDeflection,
                             GCPnts_TangentialDeflection,
                             GCPnts_UniformDeflection,)
from OCC.Core.Geom import Geom_Plane
# PythonOCC Libraries
from OCC.Core.gp import (gp_Ax1, gp_Dir, gp_Pln, gp_Pnt,
                         gp_Pnt2d, gp_Trsf, gp_Vec,)
from OCC.Core.IntCurvesFace import IntCurvesFace_Intersector
from OCC.Core.TopoDS import topods
from OCC.Core.TopTools import (TopTools_ListIteratorOfListOfShape,
                               TopTools_ListOfShape,)

# First party modules
# Own Libraries:
from SONATA.cbm.topo.utils import (P2Pdistance, Pnt2dLst_to_npArray,
                                   PolygonArea, unique_rows,)

###############################################################################
# Wire Utilities
###############################################################################


def Unique_EdgeLst(EdgeLst):
    Preci = 1e-3
    BSplineLst3D = []
    NonUniques = []
    for edg in EdgeLst:
        Adaptor = BRepAdaptor_Curve(edg)
        BSplineLst3D.append(Adaptor.BSpline())

    for i, bspline in enumerate(BSplineLst3D):
        for j, item in enumerate(BSplineLst3D):
            if bspline.IsEqual(item, Preci) and i != j:
                NonUniques.append(j)

    return NonUniques


def Wire_Orientation(wire, NbPoints=31):
    # Place Equidistant Points on BSplineLst:
    # VERY SLOW
    a = np.linspace(0.0, 1.0, num=NbPoints, endpoint=True)
    Pnt2dLst = []
    for s in a:
        Pnt2d = get_wire_Pnt2d(wire, s)
        Pnt2dLst.append(Pnt2d)
    b = Pnt2dLst_to_npArray(Pnt2dLst)

    # plt.figure(3)
    # plt.clf()
    # plt.plot(*b.T, color='black', marker='.')
    # plt.axis('equal')
    # plt.show()

    # Calculate of Polygon.
    area = PolygonArea(b)
    # print area
    if area > 0:  # counter-clockwise=True
        return True
    elif area < 0:  # clockwise=False
        return False
    else:
        return None


def NbEdges_in_wire(Wire):
    ex = BRepTools_WireExplorer(Wire)
    counter = 0
    while ex.More():
        counter += 1
        ex.Next()
    return counter


def get_wire_length(TopoDS_wire):  # std_real L = get_wire_length(TopoDS_Wire)
    AdaptorComp = BRepAdaptor_CompCurve(TopoDS_wire, True)
    length = AdaptorComp.LastParameter() - AdaptorComp.FirstParameter()
    return length


def get_wire_Pnt2d(TopoDS_wire, S):  # gp_Pnt2d P= get_wire_point(TopoDS_Wire, std_real s)
    length = get_wire_length(TopoDS_wire)
    AdaptorComp = BRepAdaptor_CompCurve(TopoDS_wire, True)
    P = AdaptorComp.Value(S * length)
    return gp_Pnt2d(P.X(), P.Y())


def get_wire_Pnt(TopoDS_wire, S):  # gp_Pnt2d P= get_wire_point(TopoDS_Wire, std_real s)
    length = get_wire_length(TopoDS_wire)
    AdaptorComp = BRepAdaptor_CompCurve(TopoDS_wire, True)
    P = AdaptorComp.Value(S * length)
    return P


def equidistant_Points_on_wire(wire, NbPoints):
    AdaptorComp = BRepAdaptor_CompCurve(wire, True)
    discretization = GCPnts_QuasiUniformAbscissa(AdaptorComp, NbPoints)

    PntLst = []
    for j in range(1, NbPoints + 1):
        para = discretization.Parameter(j)
        Pnt = gp_Pnt()
        AdaptorComp.D0(para, Pnt)
        PntLst.append(Pnt)
    return PntLst


def rotate_wire(wire, ax1, angle, copy=False):  # for top
    aTrsf = gp_Trsf()
    aTrsf.SetRotation(ax1, angle)
    aBRespTrsf = BRepBuilderAPI_Transform(wire, aTrsf)
    RotatedWire = aBRespTrsf.Shape()
    rotWire = topods.Wire(RotatedWire)
    return rotWire


def translate_wire(wire, gp_Pnt1, gp_Pnt2, copy=False):
    aTrsf = gp_Trsf()
    aTrsf.SetTranslation(gp_Pnt1, gp_Pnt2)
    aBRespTrsf = BRepBuilderAPI_Transform(wire, aTrsf, copy)
    atraslatedShape = aBRespTrsf.Shape()
    translatedWire = topods.Wire(atraslatedShape)
    return translatedWire


def scale_wire(wire, gp_Pnt1, factor, copy=False):
    aTrsf = gp_Trsf()
    aTrsf.SetScale(gp_Pnt1, factor)
    aBRespTrsf = BRepBuilderAPI_Transform(wire, aTrsf, copy)
    aBRespTrsf.Build()
    aScaledShape = aBRespTrsf.Shape()
    scaledWire = topods.Wire(aScaledShape)
    return scaledWire


def trsf_wire(wire, gp_Trsf, copy=False):
    aBRespTrsf = BRepBuilderAPI_Transform(wire, gp_Trsf, copy)
    aBRespTrsf.Build()
    atrsfedShape = aBRespTrsf.Shape()
    trsfedWire = topods.Wire(atrsfedShape)
    return trsfedWire


def mirror_wire_pnt_dir(brep, pnt, direction, copy=False):
    """
    from pythonocc-utils
    @param brep:
    @param line:
    """
    trns = gp_Trsf()
    trns.SetMirror(gp_Ax1(pnt, direction))
    brep_trns = BRepBuilderAPI_Transform(brep, trns, copy)
    brep_trns.Build()
    return topods.Wire(brep_trns.Shape())


def mirror_wire_axe2(brep, axe2, copy=False):
    """
    from pythonocc-utils
    @param brep:
    @param line:
    """
    trns = gp_Trsf()
    trns.SetMirror(axe2)
    brep_trns = BRepBuilderAPI_Transform(brep, trns, copy)
    brep_trns.Build()
    return topods.Wire(brep_trns.Shape())


def discretize_wire(wire, Resolution=500, Deflection=5e-3, refine=True, display=None):
    Pnt2dLst = []
    AdaptorComp = BRepAdaptor_CompCurve(wire, True)
    # discretization = GCPnts_TangentialDeflection(AdaptorComp,AngularDeflection,CurvatureDeflection,MinimumOfPoints, UTol)
    discretization = GCPnts_QuasiUniformDeflection(AdaptorComp, Deflection, 1)
    # print 'NbPoints3 = ', discretization.NbPoints()

    NbPoints = discretization.NbPoints()
    # print(NbPoints)
    for j in range(1, NbPoints + 1):
        Pnt = discretization.Value(j)
        Pnt2dLst.append(Pnt)
        # display.DisplayShape(Pnt)

    a = Pnt2dLst_to_npArray(Pnt2dLst)
    b = np.around(a, 10)  # Evenly round to the given number of decimals.

    # check if it is closed:
    if np.array_equal(b[0], b[-1]):
        closed = True
    else:
        closed = False
    # print 'Closed: ',closed

    npArray = unique_rows(b)  # Remove possible doubles!
    if closed:
        npArray = np.vstack((npArray, b[0]))
    else:
        None

    # Interpolate Large spaces!
    seg_P2Plength = []
    cumm_length = 0
    data = npArray

    for j in range(0, len(data) - 1):
        seg_P2Plength.append(P2Pdistance(data[j], data[j + 1]))
        cumm_length += P2Pdistance(data[j], data[j + 1])

    # Check if Refinement is necessary:
    if len(seg_P2Plength) > 0 and max(seg_P2Plength) > cumm_length / Resolution:
        Refinement = True
    else:
        Refinement = False

    Refinement = refine

    while Refinement == True:
        temp_data = []
        for i in range(0, len(data) - 1):
            if P2Pdistance(data[i], data[i + 1]) > (cumm_length / Resolution):
                p0 = data[i]
                p1 = data[i + 1]
                v1 = p1 - p0
                p05 = p0 + v1 / 2
                temp_data.append(p0)
                temp_data.append(p05)
            else:
                temp_data.append(data[i])

        temp_data.append(data[-1])
        data = np.vstack(temp_data)

        # Check if further Refinement is necessary
        seg_P2Plength = []
        cumm_length = 0
        for j in range(0, len(data) - 1):
            seg_P2Plength.append(P2Pdistance(data[j], data[j + 1]))
            cumm_length += P2Pdistance(data[j], data[j + 1])

        if max(seg_P2Plength) > cumm_length / Resolution:
            Refinement = True
        else:
            Refinement = False

    return data


def build_wire_from_BSplineLst2(BSplineLst):  # Builds TopoDS_Wire from connecting BSplineSegments and returns it
    tmp_wire = BRepBuilderAPI_MakeWire()
    for i, item in enumerate(BSplineLst):
        P = gp_Pnt(0, 0, 0)
        V = gp_Dir(gp_Vec(0, 0, 1))
        Plane = Geom_Plane(P, V)
        tmp_edge = BRepBuilderAPI_MakeEdge(item, Plane)
        # print(tmp_edge.IsDone())
        tmp_wire.Add(tmp_edge.Edge())
    return tmp_wire.Wire()


def build_wire_from_BSplineLst(BSplineLst, twoD=True):  # Builds TopoDS_Wire from connecting BSplineSegments and returns it
    P = gp_Pnt(0, 0, 0)
    V = gp_Dir(gp_Vec(0, 0, 1))
    Plane = Geom_Plane(P, V)
    TopTools_EdgeLst = TopTools_ListOfShape()
    for i, item in enumerate(BSplineLst):
        if twoD:
            tmp_edge = BRepBuilderAPI_MakeEdge(item, Plane)
        else:
            tmp_edge = BRepBuilderAPI_MakeEdge(item)
        # print(tmp_edge.IsDone())
        TopTools_EdgeLst.Append(tmp_edge.Edge())

    wire = BRepBuilderAPI_MakeWire()
    wire.Add(TopTools_EdgeLst)
    wire.Build()
    return wire.Wire()
