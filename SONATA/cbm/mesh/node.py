# -*- coding: utf-8 -*-
"""
Created on Fri Mar 03 13:37:55 2017

@author: TPflumm
"""
# Core Library modules
from functools import total_ordering

# Third party modules
from OCC.Core.gp import gp_Pnt, gp_Pnt2d

# First party modules
from SONATA.cbm.topo.BSplineLst_utils import (find_BSplineLst_pos,
                                              findPnt_on_BSplineLst,)


@total_ordering
class Node(object):

    class_counter = 1  # class attribute
    __slots__ = ("id", "Pnt2d", "parameters", "corner", "cornerstyle", "regular_corner", "displacement")
    # tell Python, to manage class attributes memorypreserving.

    def __init__(self, Pnt2d, parameters=["0", 0, 0]):
        self.id = self.__class__.class_counter
        self.__class__.class_counter += 1

        self.Pnt2d = Pnt2d  # gp_Pnt2d
        self.parameters = parameters  # [LayerID, idx, U]  'a_0010','b_0010',
        self.corner = False
        self.cornerstyle = None
        self.regular_corner = None
        self.displacement = [None, None, None]

    @property
    def Pnt(self):
        return gp_Pnt(self.Pnt2d.X(), self.Pnt2d.Y(), 0)  # gp_Pnt

    @property
    def coordinates(self):
        return [self.Pnt2d.X(), self.Pnt2d.Y()]  # [x,y]

    def b_BSplineLst_pos(self, b_BSplineLst):
        para = findPnt_on_BSplineLst(self.Pnt2d, b_BSplineLst)  # slow but robust
        # para = [self.parameters[1],self.parameters[2]]          #fast but not so sure if it's robust
        try:
            return find_BSplineLst_pos(b_BSplineLst, para)

        except:
            print("Node", self.id, "not on b_BSplineLst")
            return None

    def __repr__(self):
        return str("Node: %s @ [%.3f,%.3f]" % (self.id, self.coordinates[0], self.coordinates[1]))

    def __eq__(self, other):
        # return self.Pnt2d.IsEqual(other.Pnt2d,1e-6)    #slow but robust
        return self.id == other.id  # faster, but be careful not to assign the id's otherwise in the code

    def __hash__(self):
        """if you define __eq__, the default __hash__ (namely, hashing the address 
        of the object in memory) goes away. This is important because hashing 
        needs to be consistent with equality: equal objects need to hash the same."""
        return id(self)  # faster, but be careful not to assign the id's otherwise in the code

    def __lt__(self, other):
        """with the definition of__lt__ magic methon an the decorator 
        @total_ordering it is possible to define a sorting methodology for the 
        instaces within a list eg.: sorted(b_nodes) rather than 
        b_nodes =  sorted(b_nodes, key=lambda n: (n.parameters[1], n.parameters[2]) )"""
        if self.parameters[1] > other.parameters[1]:
            return True
        elif self.parameters[1] < other.parameters[1]:
            return False
        elif self.parameters[1] == other.parameters[1]:
            # print (self.parameters[2],other.parameters[2])
            # print (self.id, self.parameters)
            if self.parameters[2] > other.parameters[2]:
                return True
            else:
                return False

    def __getstate__(self):
        """Return state values to be pickled."""
        return (self.id, self.parameters, self.corner, self.cornerstyle, self.regular_corner, self.displacement, self.coordinates)

    def __setstate__(self, state):
        """Restore state from the unpickled state values."""
        self.id, self.parameters, self.corner, self.cornerstyle, self.regular_corner, self.displacement, tmp_coords = state
        self.Pnt2d = gp_Pnt2d(tmp_coords[0], tmp_coords[1])

    def Distance(self, other):
        return self.Pnt2d.Distance(other.Pnt2d)


def Pnt2dLst_to_NodeLst(Pnt2dLst):
    NodeLst = []
    for Pnt2d in Pnt2dLst:
        NodeLst.append(Node(Pnt2d))
    return NodeLst
