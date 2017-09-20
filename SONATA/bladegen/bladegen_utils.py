# -*- coding: utf-8 -*-
"""
Created on Mon Sep 04 10:45:43 2017

@author: TPflumm
"""

#PythonOCC Libraries
from OCC.TopoDS import TopoDS_Wire,  TopoDS_Vertex
from OCC.BRepOffsetAPI import BRepOffsetAPI_ThruSections

def make_loft(elements, ruled=False, tolerance=1e-6, continuity=4, check_compatibility=True):
    'adapted from OCCUtils'
    sections = BRepOffsetAPI_ThruSections(False, ruled, tolerance)
    for i in elements:
        if isinstance(i, TopoDS_Wire):
            sections.AddWire(i)
        elif isinstance(i, TopoDS_Vertex):
            sections.AddVertex(i)
        else:
            raise TypeError('elements is a list of TopoDS_Wire or TopoDS_Vertex, found a %s fool' % i.__class__)

    sections.CheckCompatibility(check_compatibility)
    sections.SetContinuity(continuity)
    sections.Build()
    loft = sections.Shape()
    return loft
    