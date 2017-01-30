# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 12:49:54 2017

@author: TPflumm
"""
import numpy as np

from OCC.gp import gp_Circ2d, gp_Pnt2d, gp_Ax2d, gp_Dir2d, gp_Vec2d
from OCC.Geom2d import Geom2d_Circle, Geom2d_BezierCurve
from OCC.Display.SimpleGui import init_display

from utils import point2d_list_to_TColgp_Array1OfPnt2d


def WeightStyle1_BezierPoint(P1,P2,R,W):
            Vec12 =  	gp_Vec2d(P1,P2)
            n12 = Vec12.GetNormal()
            n12.Normalize()
            P_tmp = P1.Translated(Vec12.Multiplied(0.5))
            P12 = P_tmp.Translated(n12.Multiplied(W*R))
            return P12

class Weight(object):
   
    def __init__(self, ID, XPos, YPos, Diameter, MatID, weight_style=0):
        self.ID = ID
        self.X = XPos
        self.Y = YPos
        self.D = Diameter
        self.MatID = MatID
        self.style = weight_style
        
        Center = gp_Pnt2d(self.X,self.Y)
        Ax2d = gp_Ax2d(Center,gp_Dir2d(1,0))
        #Circ = gp_Circ2d(Ax2d,R)
        self.Curve = Geom2d_Circle(Ax2d,self.D/2)
        
        
#        elif style == 1:
#        W = 0.4
#        P1 = Center.Translated(gp_Vec2d(-R,0)) 
#        ANG = float(120)/180*np.pi
#        P2 = P1.Rotated(Center,ANG)
#        P3 = P2.Rotated(Center,ANG)
#        
#        P12 = WeightStyle1_BezierPoint(P1,P2,W)
#        P23 = WeightStyle1_BezierPoint(P2,P3,W)
#        P31 = WeightStyle1_BezierPoint(P3,P1,W)
#        
#        Bezier1 = Geom2d_BezierCurve(point2d_list_to_TColgp_Array1OfPnt2d([P1,P12,P2]))
#        Bezier2 = Geom2d_BezierCurve(point2d_list_to_TColgp_Array1OfPnt2d([P2,P23,P3]))
#        Bezier3 = Geom2d_BezierCurve(point2d_list_to_TColgp_Array1OfPnt2d([P3,P31,P1]))
        
        
       
#======================================================
#       MAIN
#======================================================
if __name__ == '__main__':
    #==========================
    #DISPLAY CONFIG:
    display, start_display, add_menu, add_function_to_menu = init_display('wx')
    display.Context.SetDeviationAngle(0.000001)       # 0.001 default. Be careful to scale it to the problem.
    display.Context.SetDeviationCoefficient(0.000001) # 0.001 default. Be careful to scale it to the problem. 
    #show_coordinate_system(display) #CREATE AXIS SYSTEM for Visualization

    BW = Weight(1,0.1,0.1,0.25,1,0)
    display.DisplayShape(BW.Curve, color="BLACK")
               
    

    display.set_bg_gradient_color(20,6,111,200,200,200)

    display.View_Top()
    display.FitAll()
    
      
    start_display()