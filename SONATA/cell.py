# -*- coding: utf-8 -*-
"""
Created on Fri Mar 03 13:39:55 2017

@author: TPflumm
"""
import numpy as np

from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeEdge
from OCC.gp import gp_Vec2d,gp_Pnt2d

from utils import PolygonArea, calc_angle_between



def calc_cell_angles(cell):
    corners = []
    for node in cell.nodes:
        corners.append(node.coordinates)         
    corners = np.asarray(corners)   
    temp = []
    for i in range(0,corners.shape[0]):
            if i == corners.shape[0]-1: #last point
                v1 = corners[i-1]-corners[i] 
                v2 = corners[0]-corners[i]
            else:
                v1 = corners[i-1]-corners[i]
                v2 = corners[i+1]-corners[i]
            temp.append(calc_angle_between(v1,v2))
    return np.array(temp)


class Cell(object):
    class_counter= 1
    def __init__(self,nodeLst):                  #int
        self.id = self.__class__.class_counter
        self.__class__.class_counter += 1
        
        self.nodes = nodeLst                #[node,node,node,nodes]      !!!counterclockwise direction!!!
        self.face  = []                     #[rear,top,front,bottom]        !!!counterclockwise direction!!!       
        #self.wire  = self.build_wire()      #TopoDs_wire
        self.wire  = None     #TopoDs_wire
        self.neighbours = []                #-[Cell_ID,CELL_ID... ]
        #self.theta_1 = self.calc_theta_1()  #Ply coordinate system is formed by rotating the global coordinate system in the right-hand sense about the amount 0<Theta_1<260.
                                            #Theta_1[0:9] is a list storing nine real numbers for the layer plane angles at the nodes of ths element. For simplification, if the 
                                            #ply orinetation can be considered as uniform this element. Theta_1[0] stores the layer plane angles and Theta_1[1] = 540, and all the 
                                            #rest can be zeroes or other real numbers because they do not enter the calculation. If the elements''' 
        self.theta_3 = None                 #The Ply coordiate system is rotated about y3 in the right hand sense by the amount -90<Theta_3<90 to for the material system. 
        self.MatID  = None                 #material id, int
        
        #Element quality critiria
        self.area = self.calc_area()
        self.minimum_angle = self.minimum_angle()
        self.maximum_angle = self.maximum_angle()
        self.minimum_edge = None
        self.maximum_edge = None
        self.shape_quality = None
        self.minimum_jacobinan = None
        #AREA RATIO to neighbors
        
    def __repr__(self): 
        #we can tell Python how to prepresent an object of our class (when using a print statement) for general purposes use  __repr__(self): 
        STR = ''
        STR += str('Cell %s w. nodes:\t' % (self.id))
        for n in self.nodes:
            STR += str('%i, ' % (n.id))
            
        return  STR
        
    def calc_theta_1(self):  
        theta_1 = [0] * 9
        v0 = gp_Vec2d(gp_Pnt2d(0,0),gp_Pnt2d(1,0))
        v1 = gp_Vec2d(self.nodes[1].Pnt2d,self.nodes[2].Pnt2d)
        theta_11 = (v0.Angle(v1))*180/np.pi
        if theta_11<0:
            theta_11 = 360+theta_11
        theta_1[0] = theta_11
        theta_1[1] = 540
        return theta_1

    def calc_area(self):  
        corners = []
        for node in self.nodes:
            corners.append(node.coordinates)      
        return PolygonArea(corners)     
    
    def minimum_angle(self):  
        #print np.amin(calc_cell_angles(self))
        return np.amin(calc_cell_angles(self))
    
    def maximum_angle(self):  
        #print np.amax(calc_cell_angles(self))
        return np.amax(calc_cell_angles(self))       

    def build_wire(self):
        WireBuilder = BRepBuilderAPI_MakeWire()
        for i in range(0,len(self.nodes)-1):
            me = BRepBuilderAPI_MakeEdge(self.nodes[i].Pnt, self.nodes[i+1].Pnt)
            if me.IsDone():
                WireBuilder.Add(me.Edge())
        
        me = BRepBuilderAPI_MakeEdge(self.nodes[-1].Pnt, self.nodes[0].Pnt)
        if me.IsDone():
            WireBuilder.Add(me.Edge())         
        
        return WireBuilder.Wire()