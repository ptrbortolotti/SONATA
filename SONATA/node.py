# -*- coding: utf-8 -*-
"""
Created on Fri Mar 03 13:37:55 2017

@author: TPflumm
"""
from OCC.gp import gp_Pnt

class Node(object):
    class_counter= 1
    def __init__(self, Pnt2d, parameters=['0',0,0]):
        self.id= self.__class__.class_counter
        self.__class__.class_counter += 1
        
        self.Pnt2d = Pnt2d  #gp_Pnt2d
        self.coordinates = [self.Pnt2d.X(),self.Pnt2d.Y()]  #[x,y]
        self.parameters = parameters    #[LayerID, idx, U]
        self.corner = False
        self.cornerstyle = None
        self.face_pointer = []
        
    def __repr__(self): 
        #we can tell Python how to prepresent an object of our class (when using a print statement) for general purposes use  __repr__(self): 
        return  str('Node: %s @ [%.3f,%.3f]' % (self.id, self.coordinates[0],self.coordinates[1]))
            
    @property
    def Pnt(self):
        return gp_Pnt(self.Pnt2d.X(),self.Pnt2d.Y(),0)  #gp_Pnt


def Pnt2dLst_to_NodeLst(Pnt2dLst):
    NodeLst = []
    for Pnt2d in Pnt2dLst:
        NodeLst.append(Node(Pnt2d))
    return NodeLst