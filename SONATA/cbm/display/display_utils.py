## 2016 Tobias Pflumm (tobias.pflumm@tum.de)
## Note that this is adapted from pythonocc-Utils.
''' This Module describes all displying functionalietes of the code'''

import sys
import os
import matplotlib as plt
import math

from OCC.Display.SimpleGui import init_display
from OCC.gp import gp_Pnt2d, gp_Pnt, gp_Pln, gp_Dir, gp_Vec, gp_Trsf, gp_Ax3,gp_Ax1
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.Graphic3d import (Graphic3d_EF_PDF,
                           Graphic3d_EF_SVG,
                           Graphic3d_EF_TEX,
                           Graphic3d_EF_PostScript,
                           Graphic3d_EF_EnhPostScript)

from OCC.Quantity import Quantity_Color
from OCC.AIS import AIS_Shape

from SONATA.cbm.topo.wire_utils import rotate_wire, translate_wire

#===========================================================================
# MENU FUNCTIONALITIES
#===========================================================================


#===========================================================================
# DISPLAY
#===========================================================================
def init_viewer():
    display, start_display, add_menu, add_function_to_menu = init_display('wx')
    display.Context.SetDeviationAngle(0.000001)       # 0.001 default. Be careful to scale it to the problem.
    display.Context.SetDeviationCoefficient(0.000001) # 0.001 default. Be careful to scale it to the problem. 

    add_menu('screencapture')
    add_function_to_menu('screencapture', export_to_BMP)
    add_function_to_menu('screencapture', export_to_PNG)
    add_function_to_menu('screencapture', export_to_JPEG)
    add_function_to_menu('screencapture', export_to_TIFF)
    add_function_to_menu('screencapture', exit)              
    return display    

    
def export_to_PDF(display,event=None):
    f = display.View.View().GetObject()
    #display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('img/capture_pdf%s.pdf' % i):
        i += 1
    f.Export('img/capture_pdf%s.pdf' % i, Graphic3d_EF_PDF)
    print("EXPORT: \t Screencapture exported to img/capture_pdf%s.pdf" % i)
    #display.set_bg_gradient_color(20,6,111,200,200,200)
    
def export_to_SVG(display,event=None):
    f = display.View.View().GetObject()
    #display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('img/capture_svg%s.svg' % i):
        i += 1
    f.Export('img/capture_svg_%s.svg' % i, Graphic3d_EF_SVG)
    print("EXPORT: \t Screencapture exported to img/capture_svg%s.svg" % i)
    #display.set_bg_gradient_color(20,6,111,200,200,200)
    
def export_to_PS(display,event=None):
    f = display.View.View().GetObject()
    #display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('img/capture_ps%s.ps' % i):
        i += 1
    f.Export('img/capture_ps%s.ps' % i, Graphic3d_EF_PostScript)
    print("EXPORT: \t Screencapture exported to img/capture_ps%s.ps" % i)
    #display.set_bg_gradient_color(20,6,111,200,200,200)

def export_to_EnhPS(display,event=None):
    f = display.View.View().GetObject()
    #display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('img/capture_Enh_ps%s.ps' % i):
        i += 1
    f.Export('img/capture_Enh_ps%s.ps' % i, Graphic3d_EF_EnhPostScript)
    print("EXPORT: \t Screencapture exported to img/capture_Enh_ps%s.ps" % i)
    #display.set_bg_gradient_color(20,6,111,200,200,200)
    
def export_to_TEX(display,event=None):
    f = display.View.View().GetObject()
    display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('img/capture_tex%s.tex' % i):
        i += 1
    f.Export('img/capture_tex%s.tex' % i, Graphic3d_EF_TEX)
    print("EXPORT: \t Screencapture exported to img/capture_tex%s.tex" % i)
    display.set_bg_gradient_color(20,6,111,200,200,200)
    
def export_to_BMP(display,event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('img/capture_bmp%s.bmp' % i):
        i += 1
    display.View.Dump('img/capture_bmp%s.bmp' % i)
    print("EXPORT: \t Screencapture exported to img/capture_bmp%s.bmp" % i)
    display.set_bg_gradient_color(20,6,111,200,200,200)

def export_to_PNG(display,event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('img/capture_png%s.png' % i):
        i += 1
    display.View.Dump('img/capture_png%s.png' % i)
    print("EXPORT: \t Screencapture exported to img/capture_png%s.bmp" % i)
    display.set_bg_gradient_color(20,6,111,200,200,200)

def export_to_JPEG(display,event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('img/capture_jpeg%s.jpeg' % i):
        i += 1
    display.View.Dump('img/capture_jpeg%s.jpeg' % i)
    print("EXPORT: \t Screencapture exported to img/capture_jpeg%s.jpeg" % i)
    display.set_bg_gradient_color(20,6,111,200,200,200)

def export_to_TIFF(display,event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('img/capture_tiff%s.tiff' % i):
        i += 1
    display.View.Dump('img/capture_tiff%s.tiff' % i)
    print("EXPORT: \t Screencapture exported to img/capture_tiff%s.tiff" % i)
    display.set_bg_gradient_color(20,6,111,200,200,200)


def print_xy_click(SHP, *kwargs):
    for shape in SHP:
        print(("Shape selected: ", shape))
    print(kwargs)

def exit():
    sys.exit()
    
def show_coordinate_system(display,length,event=None):
    '''CREATE AXIS SYSTEM for Visualization'''
    O  = gp_Pnt(0., 0., 0.)
    p1 = gp_Pnt(length,0.,0.)
    p2 = gp_Pnt(0.,length,0.)
    p3 = gp_Pnt(0.,0.,length)
    
    h1 = BRepBuilderAPI_MakeEdge(O,p1).Shape()
    h2 = BRepBuilderAPI_MakeEdge(O,p2).Shape()
    h3 = BRepBuilderAPI_MakeEdge(O,p3).Shape()

    display.DisplayShape(O,color='BLACK')
    display.DisplayShape(h1,color='RED')
    display.DisplayShape(h2,color='GREEN')
    display.DisplayShape(h3,color='BLUE')
    display.DisplayMessage(p1,'x',message_color=(0,0,0))
    display.DisplayMessage(p2,'y',message_color=(0,0,0))
    display.DisplayMessage(p3,'z',message_color=(0,0,0))


def display_points_of_array(array,display):
    for j in range(array.Lower(),array.Upper()+1):
        p = array.Value(j)
        display.DisplayShape(p, update=False)         


def display_custome_shape(display,shape,linewidth=1.0,transparency=0.0,RGB=[0,0,0]):
    s = shape
    ais_shp = AIS_Shape(s)
    ais_shp.SetWidth(linewidth)
    ais_shp.SetTransparency(transparency)
    ais_shp.SetColor(Quantity_Color(RGB[0], RGB[1], RGB[2], 0))
    ais_context = display.GetContext().GetObject()
    ais_context.Display(ais_shp.GetHandle())
    return None


def transform_wire_2to3d(display,wire,coord=(0,0,0),alpha=0,beta=0,color='BLACK',show=True):
    wire = rotate_wire(wire,gp_Ax1(gp_Pnt(0,0,0),gp_Dir(0,0,1)),alpha)
    wire = rotate_wire(wire,gp_Ax1(gp_Pnt(0,0,0),gp_Dir(0,1,0)),beta)
    wire = translate_wire(wire,gp_Pnt(0,0,0),gp_Pnt(float(coord[0]),float(coord[1]),float(coord[2])))
    if show:
        display.DisplayShape(wire, color=color)
    return wire

#=======================SONATA DISPLAY FUCTIONS===================================
def display_SONATA_SegmentLst(display,SegmentLst,coord=(0,0,0),alpha=0,beta=0):
    # transfer shapes and display them in the viewer   
    transform_wire_2to3d(display,SegmentLst[0].wire,coord,alpha,beta,color='BLACK')
    
    for i,seg in enumerate(SegmentLst):
        wire = transform_wire_2to3d(display,seg.wire,coord,alpha,beta,)
    
        k = 0
        for j,layer in enumerate(seg.LayerLst):
            [R,G,B,T] =  plt.cm.jet(k*50)
            
            wire = transform_wire_2to3d(display,layer.wire, coord, alpha, beta, show=False)
            display.DisplayColoredShape(wire, Quantity_Color(R, G, B, 0),update=True)
#            #display Start Point
#            string = 'Layer:'+str(layer.ID)+'(S1='+str(layer.S1)+')'
#            P = gp_Pnt(layer.StartPoint.X(),layer.StartPoint.Y(),0)
#            display.DisplayShape(P,color="BLUE")
#            display.DisplayMessage(P,string,message_color=(0.0,0.0,0.0))
#            
#            #display End Point
#            string = 'Layer:'+str(layer.ID)+'(S1='+str(layer.S2)+')'
#            P = gp_Pnt(layer.EndPoint.X(),layer.EndPoint.Y(),0)
#            display.DisplayShape(P,color="RED")
#            display.DisplayMessage(P,string,message_color=(0.0,0.0,0.0))
            k = k+1;
            if k>5:
                k = 0
    return None




if __name__ == '__main__':   
    pass
