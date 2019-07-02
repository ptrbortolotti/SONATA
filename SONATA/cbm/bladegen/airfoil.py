# -*- coding: utf-8 -*-
"""
Created on Mon Sep 04 10:45:43 2017

@author: TPflumm
"""
import matplotlib.pyplot as plt

#PythonOCC Libraries
from OCC.Core.gp import gp_Pnt, gp_Vec
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire
from OCC.Core.GeomAPI import GeomAPI_Interpolate
from OCC.Core.Geom import Geom_BezierCurve

#SONATA modules:
from SONATA.cbm.fileIO.readinput  import UIUCAirfoil
from SONATA.cbm.topo.utils import TColgp_HArray1OfPnt_from_nparray,point_list_to_TColgp_Array1OfPnt

class Airfoil(object):
    class_counter= 0
    def __init__(self, name = None, trailing_edge_style=0):
        self.id = self.__class__.class_counter
        self.__class__.class_counter += 1
        
        if name:
            self.name = name  
            self.trailing_edge_style = trailing_edge_style
            self.array = UIUCAirfoil(self.name)
            self.wire = self.gen_wire()
    
    def restore_counter(self):
        self.__class__.class_counter = 0

    def gen_wire(self):
        harray = TColgp_HArray1OfPnt_from_nparray(self.array)
        anInterpolation = GeomAPI_Interpolate(harray.GetHandle(), False, 0.001)
        anInterpolation.Perform()
        bspline = anInterpolation.Curve().GetObject()
        edge = BRepBuilderAPI_MakeEdge(bspline.GetHandle())
        #check if bspline is closed:      
        if not bspline.IsClosed():
            #print self.name, 'not closed'
            First = bspline.FirstParameter()
            Last = bspline.LastParameter()
            
            StartPnt = gp_Pnt()
            EndPnt = gp_Pnt()
            V1 = gp_Vec()
            V2 = gp_Vec()
            bspline.D1(First, StartPnt, V1)
            bspline.D1(Last, EndPnt, V2)
            
            V1.Normalize()
            V2.Normalize()
            
            P1 = StartPnt.Translated(V1.Multiplied(-2e-3))
            P2 = EndPnt.Translated(V2.Multiplied(2e-3))
            
            if self.trailing_edge_style == 1:
                bezier = Geom_BezierCurve(point_list_to_TColgp_Array1OfPnt([StartPnt,P1,P2,EndPnt]))
                #display.DisplayShape(bezier, color='BLUE')
                trailing_edge = BRepBuilderAPI_MakeEdge(bezier.GetHandle())
            
            else:
                trailing_edge = BRepBuilderAPI_MakeEdge(StartPnt,EndPnt)
            
            wire = BRepBuilderAPI_MakeWire(edge.Edge(),trailing_edge.Edge()).Wire()
        else:
            wire = BRepBuilderAPI_MakeWire(edge.Edge()).Wire()
        return wire
    
        
    def plot(self):
        plt.figure()
        plt.plot(self.array[0,:], self.array[1,:], 'r-o')
        plt.axis('equal')
        plt.title(self.name, fontsize=10)
        plt.xlabel('x')
        plt.ylabel('y')
    
    