# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 13:37:00 2017

@author: TPflumm
"""
from OCC.STEPControl import STEPControl_Writer, STEPControl_AsIs, STEPControl_GeometricCurveSet
from OCC.Interface import Interface_Static_SetCVal
from OCC.IFSelect import IFSelect_RetDone

def export_to_step(itemLst, filename):
    #====================STEP-EXPORT=============================================
    # initialize the STEP exporter
    step_writer = STEPControl_Writer()
    #Interface_Static_SetCVal("write.step.schema", "AP203")
    #step_writer.SetTolerance(1e-6)
    
    # transfer shapes and write file
    for item in itemLst:
        step_writer.Transfer(item, STEPControl_AsIs)   
    status = step_writer.Write(filename)    
    assert(status == IFSelect_RetDone)
    
    return None