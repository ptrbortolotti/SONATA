## 2016 Tobias Pflumm (tobias.pflumm@tum.de)
## Note that this is adapted from pythonocc-Utils.
""" This Module describes all displying functionalietes of the code"""

# Core Library modules
import math
import os
import sys

# Third party modules
import matplotlib as plt
import numpy as np
from OCC.Core.AIS import AIS_Shape
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.Core.gp import (gp_Ax1, gp_Ax3, gp_Dir, gp_Pln,
                         gp_Pnt, gp_Pnt2d, gp_Trsf, gp_Vec,)
from OCC.Core.Quantity import Quantity_Color
from OCC.Display.backend import get_qt_modules
from OCC.Display.SimpleGui import init_display

# First party modules
from SONATA.cbm.topo.wire_utils import rotate_wire, translate_wire, trsf_wire
from SONATA.utl.trsf import trsf_cbm_to_blfr


def display_config(DeviationAngle=1e-5, DeviationCoefficient=1e-5, bg_c=((20, 6, 111), (200, 200, 200)), cs_size=25):
    """
    CBM method that initializes and configures the pythonOcc 3D Viewer 
    and adds Menues to the toolbar. 
    
    Parameters
    ----------
    DeviationAngle : float, optional 
        default = 1e-5
    DeviationCoefficient : float, optional
        default = 1e-5 
    bg_c : tuple, optional
        Background Gradient Color ((RBG Tuple),(RBG Tuple)) the default 
        values are a CATIA style gradient for better 3D visualization. 
        for a white background use: ((255,255,255,255,255,255))
    cs_size : float 
        coordinate system size in [mm]
        
    Returns
    ----------
    tuple :
        (self.display: the display handler for the pythonOcc 3D Viewer
        self.start_display: function handle
        self.add_menu: function handle
        self.add_function_to_menu: function handle)
    
    
    See Also
    ----------
    OCC.Display.SimpleGui : PyhtonOcc wrapper provides more details on this 
        method
    """

    def export_png():
        return export_to_PNG(display)

    def export_jpg():
        return export_to_JPEG(display)

    #    def export_pdf(): return export_to_PDF(display)
    #    def export_svg(): return export_to_SVG(display)
    #    def export_ps(): return export_to_PS(display)

    # ===========DISPLAY CONFIG:===============
    display, start_display, add_menu, add_function_to_menu = init_display()
    display.Context.SetDeviationAngle(DeviationAngle)  # 0.001 default. Be careful to scale it to the problem.
    display.Context.SetDeviationCoefficient(DeviationCoefficient)  # 0.001 default. Be careful to scale it to the problem.
    display.set_bg_gradient_color([bg_c[0][0], bg_c[0][1], bg_c[0][2]], [bg_c[1][0], bg_c[1][1], bg_c[1][2]])
    show_coordinate_system(display, cs_size)

    add_menu("File")
    add_function_to_menu("File", close)

    add_menu("View")
    add_function_to_menu("View", display.FitAll)
    add_function_to_menu("View", display.View_Bottom)
    add_function_to_menu("View", display.View_Top)
    add_function_to_menu("View", display.View_Left)
    add_function_to_menu("View", display.View_Right)
    add_function_to_menu("View", display.View_Front)
    add_function_to_menu("View", display.View_Rear)
    add_function_to_menu("View", display.View_Iso)

    add_menu("Screencapture")
    add_function_to_menu("Screencapture", export_png)
    add_function_to_menu("Screencapture", export_jpg)
    #    add_function_to_menu('Screencapture', export_pdf)
    #    add_function_to_menu('Screencapture', export_svg)
    #    add_function_to_menu('Screencapture', export_ps)

    add_menu("Export")

    return (display, start_display, add_menu, add_function_to_menu)


def export_to_PDF(display, event=None):
    f = display.View.View()
    # display.set_bg_gradient_color([255,255,255],[255,255,255])
    i = 0
    while os.path.exists("img/capture_pdf%s.pdf" % i):
        i += 1
    f.Export("img/capture_pdf%s.pdf" % i, Graphic3d_EF_PDF)
    print("EXPORT: \t Screencapture exported to img/capture_pdf%s.pdf" % i)
    # display.set_bg_gradient_color([20,6,111],[200,200,200])


def export_to_SVG(display, event=None):
    f = display.View.View()
    # display.set_bg_gradient_color([255,255,255],[255,255,255])
    i = 0
    while os.path.exists("img/capture_svg%s.svg" % i):
        i += 1
    f.Export("img/capture_svg_%s.svg" % i, Graphic3d_EF_SVG)
    print("EXPORT: \t Screencapture exported to img/capture_svg%s.svg" % i)
    # display.set_bg_gradient_color([20,6,111],[200,200,200])


def export_to_PS(display, event=None):
    f = display.View.View()
    # display.set_bg_gradient_color([255,255,255],[255,255,255])
    i = 0
    while os.path.exists("img/capture_ps%s.ps" % i):
        i += 1
    f.Export("img/capture_ps%s.ps" % i, Graphic3d_EF_PostScript)
    print("EXPORT: \t Screencapture exported to img/capture_ps%s.ps" % i)
    # display.set_bg_gradient_color([20,6,111],[200,200,200])


def export_to_EnhPS(display, event=None):
    f = display.View.View()
    # display.set_bg_gradient_color([255,255,255],[255,255,255])
    i = 0
    while os.path.exists("img/capture_Enh_ps%s.ps" % i):
        i += 1
    f.Export("img/capture_Enh_ps%s.ps" % i, Graphic3d_EF_EnhPostScript)
    print("EXPORT: \t Screencapture exported to img/capture_Enh_ps%s.ps" % i)
    # display.set_bg_gradient_color([20,6,111],[200,200,200])


def export_to_TEX(display, event=None):
    f = display.View.View()
    display.set_bg_gradient_color([255, 255, 255], [255, 255, 255])
    i = 0
    while os.path.exists("img/capture_tex%s.tex" % i):
        i += 1
    f.Export("img/capture_tex%s.tex" % i, Graphic3d_EF_TEX)
    print("EXPORT: \t Screencapture exported to img/capture_tex%s.tex" % i)
    display.set_bg_gradient_color([20, 6, 111], [200, 200, 200])


def export_to_BMP(display, event=None):
    display.set_bg_gradient_color([255, 255, 255], [255, 255, 255])
    i = 0
    while os.path.exists("img/capture_bmp%s.bmp" % i):
        i += 1
    display.View.Dump("img/capture_bmp%s.bmp" % i)
    print("EXPORT: \t Screencapture exported to img/capture_bmp%s.bmp" % i)
    display.set_bg_gradient_color([20, 6, 111], [200, 200, 200])


def export_to_PNG(display, event=None):
    display.set_bg_gradient_color([255, 255, 255], [255, 255, 255])
    i = 0
    while os.path.exists("img/capture_png%s.png" % i):
        i += 1
    display.View.Dump("img/capture_png%s.png" % i)
    print("EXPORT: \t Screencapture exported to img/capture_png%s.bmp" % i)
    display.set_bg_gradient_color([20, 6, 111], [200, 200, 200])


def export_to_JPEG(display, event=None):
    display.set_bg_gradient_color([255, 255, 255], [255, 255, 255])
    i = 0
    while os.path.exists("img/capture_jpeg%s.jpeg" % i):
        i += 1
    display.View.Dump("img/capture_jpeg%s.jpeg" % i)
    print("EXPORT: \t Screencapture exported to img/capture_jpeg%s.jpeg" % i)
    display.set_bg_gradient_color([20, 6, 111], [200, 200, 200])


def export_to_TIFF(display, event=None):
    display.set_bg_gradient_color([255, 255, 255], [255, 255, 255])
    i = 0
    while os.path.exists("img/capture_tiff%s.tiff" % i):
        i += 1
    display.View.Dump("img/capture_tiff%s.tiff" % i)
    print("EXPORT: \t Screencapture exported to img/capture_tiff%s.tiff" % i)
    display.set_bg_gradient_color([20, 6, 111], [200, 200, 200])


def print_xy_click(SHP, *kwargs):
    for shape in SHP:
        print(("Shape selected: ", shape))
    print(kwargs)


def show_coordinate_system(display, length=1, event=None):
    """CREATE AXIS SYSTEM for Visualization"""
    O = gp_Pnt(0.0, 0.0, 0.0)
    p1 = gp_Pnt(length, 0.0, 0.0)
    p2 = gp_Pnt(0.0, length, 0.0)
    p3 = gp_Pnt(0.0, 0.0, length)

    h1 = BRepBuilderAPI_MakeEdge(O, p1).Shape()
    h2 = BRepBuilderAPI_MakeEdge(O, p2).Shape()
    h3 = BRepBuilderAPI_MakeEdge(O, p3).Shape()

    display.DisplayShape(O, color="BLACK")
    display.DisplayShape(h1, color="RED")
    display.DisplayShape(h2, color="GREEN")
    display.DisplayShape(h3, color="BLUE")
    display.DisplayMessage(p1, "x", message_color=(0, 0, 0))
    display.DisplayMessage(p2, "y", message_color=(0, 0, 0))
    display.DisplayMessage(p3, "z", message_color=(0, 0, 0))


def display_points_of_array(array, display):
    for j in range(array.Lower(), array.Upper() + 1):
        p = array.Value(j)
        display.DisplayShape(p, update=False)


def display_custome_shape(display, shape, linewidth=1.0, transparency=0.0, RGB=[0, 0, 0]):
    s = shape
    ais_shp = AIS_Shape(s)
    ais_shp.SetWidth(linewidth)
    ais_shp.SetTransparency(transparency)
    ais_shp.SetColor(Quantity_Color(RGB[0], RGB[1], RGB[2], 0))
    ais_context = display.GetContext()
    ais_context.Display(ais_shp)
    return None


def transform_wire_2to3d(display, wire, coord=(0, 0, 0), alpha=0, beta=0, color="BLACK", show=True):
    wire = rotate_wire(wire, gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), alpha)
    wire = rotate_wire(wire, gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 1, 0)), beta)
    wire = translate_wire(wire, gp_Pnt(0, 0, 0), gp_Pnt(float(coord[0]), float(coord[1]), float(coord[2])))
    if show:
        display.DisplayShape(wire, color=color)
    return wire


# =======================SONATA DISPLAY FUCTIONS===================================
def display_SONATA_SegmentLst(display, SegmentLst, coord=(0, 0, 0), alpha=0, beta=0):
    # transfer shapes and display them in the viewer
    if SegmentLst:
        # transform_wire_2to3d(display,SegmentLst[0].wire,coord,alpha,beta,color='BLACK')

        for i, seg in enumerate(SegmentLst):
            wire = transform_wire_2to3d(display, seg.wire, coord, alpha, beta,)

            k = 0
            for j, layer in enumerate(seg.LayerLst):
                [R, G, B, T] = plt.cm.jet(k * 50)

                wire = transform_wire_2to3d(display, layer.wire, coord, alpha, beta, show=False)
                display.DisplayColoredShape(wire, Quantity_Color(R, G, B, 0), update=True)
                # display Start Point
                #                string = 'Layer:'+str(layer.ID)+'(S1='+str(layer.S1)+')'
                #                P = gp_Pnt(layer.StartPoint.X(),layer.StartPoint.Y(),0)
                #                display.DisplayShape(P,color="BLUE")
                #                display.DisplayMessage(P,string,message_color=(0.0,0.0,0.0))
                #
                #                #display End Point
                #                string = 'Layer:'+str(layer.ID)+'(S1='+str(layer.S2)+')'
                #                P = gp_Pnt(layer.EndPoint.X(),layer.EndPoint.Y(),0)
                #                display.DisplayShape(P,color="RED")
                #                display.DisplayMessage(P,string,message_color=(0.0,0.0,0.0))
                k = k + 1
                if k > 5:
                    k = 0
    return None


def display_cbm_SegmentLst(display, SegmentLst, Ax2_blfr, Ax2_befr):
    """
    replaces the display SONATA_SegmentLst in the future!
    
    Parameters
    ---------
    display : OCC.Display.OCCViewer.Viewer3d 
        OCC 3d Viewer instance
    
    SegmentLst : list
        list of Segments of the cbm 
        
    fromAx2 : gp_Ax2
        OCC gp_Ax2 coordinate system
        
    toAx2: : gp_Ax2
        OCC gp_Ax2 coordinate system
        
    """
    Trsf = trsf_cbm_to_blfr(Ax2_blfr, Ax2_befr)

    # transfer shapes and display them in the viewer
    if SegmentLst:
        for i, seg in enumerate(SegmentLst):
            wire = trsf_wire(seg.wire, Trsf)
            display.DisplayColoredShape(wire, Quantity_Color(0, 0, 0, 0), update=True)
            k = 0
            for j, layer in enumerate(seg.LayerLst):
                [R, G, B, T] = plt.cm.jet(k * 50)
                wire = trsf_wire(layer.wire, Trsf)
                display.DisplayColoredShape(wire, Quantity_Color(R, G, B, 0), update=True)
                k = k + 1
                if k > 5:
                    k = 0
    return None


def display_Ax2(display, Ax2, length=1):
    """
    
    
    """
    p0 = Ax2.Location()
    px = p0.Translated(gp_Vec(Ax2.XDirection()).Normalized().Multiplied(length))
    py = p0.Translated(gp_Vec(Ax2.YDirection()).Normalized().Multiplied(length))
    pz = p0.Translated(gp_Vec(Ax2.Direction()).Normalized().Multiplied(length))

    e1 = BRepBuilderAPI_MakeEdge(p0, px).Shape()
    e2 = BRepBuilderAPI_MakeEdge(p0, py).Shape()
    e3 = BRepBuilderAPI_MakeEdge(p0, pz).Shape()

    display.DisplayShape(p0, color="BLACK")
    display.DisplayShape(e1, color="RED")
    display.DisplayShape(e2, color="GREEN")
    display.DisplayShape(e3, color="BLUE")

    return None


def findMainWindow():
    # Global function to find the (open) QMainWindow in application

    QtCore, QtGui, QtWidgets, QtOpenGL = get_qt_modules()
    app = QtWidgets.QApplication.instance()
    for widget in app.topLevelWidgets():
        if isinstance(widget, QtWidgets.QMainWindow):
            return widget
    return None


def close():
    win = findMainWindow()
    win.close()


if __name__ == "__main__":
    (display, start_display, add_menu, add_function_to_menu) = display_config()
