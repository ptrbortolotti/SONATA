# -*- coding: utf-8 -*-
"""
Created on Fri Mar 03 13:37:55 2017

@author: TPflumm
"""
from OCC.gp import gp_Pnt,gp_Pnt2d

class Node(object):
    class_counter= 1
    def __init__(self, Pnt2d, parameters=['0',0,0]):
        self.id = self.__class__.class_counter
        self.__class__.class_counter += 1
        
        self.Pnt2d = Pnt2d  #gp_Pnt2d
        self.parameters = parameters    #[LayerID, idx, U]
        self.corner = False
        self.regular_corner = None
        self.cornerstyle = None
        self.displacement = [None,None,None]

    @property
    def Pnt(self):
        return gp_Pnt(self.Pnt2d.X(),self.Pnt2d.Y(),0)  #gp_Pnt
    
    @property
    def coordinates(self):
        return [self.Pnt2d.X(),self.Pnt2d.Y()]  #[x,y]
    
    def __repr__(self): 
        return  str('Node: %s @ [%.3f,%.3f]' % (self.id, self.coordinates[0],self.coordinates[1]))
        
    def __eq__(self,other):
        #return self.Pnt2d.IsEqual(other.Pnt2d,1e-6)    #slow but robust
        return self.id == other.id    #faster, but be careful not to assign the id's otherwise in the code
    
    def __getstate__(self):
        """Return state values to be pickled."""
        return (self.id, self.parameters, self.corner, self.cornerstyle, self.coordinates)   
    
    def __setstate__(self, state):
        """Restore state from the unpickled state values."""
        self.id, self.parameters, self.corner, self.cornerstyle, tmp_coords   = state
        self.Pnt2d = gp_Pnt2d(tmp_coords[0],tmp_coords[1])
    
    def Distance(self,other):
       return self.Pnt2d.Distance(other.Pnt2d)
   
    
        

def Pnt2dLst_to_NodeLst(Pnt2dLst):
    NodeLst = []
    for Pnt2d in Pnt2dLst:
        NodeLst.append(Node(Pnt2d))
    return NodeLst