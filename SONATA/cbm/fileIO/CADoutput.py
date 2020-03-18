# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 13:37:00 2017

@author: TPflumm
"""
# Third party modules
from OCC.Core.IFSelect import IFSelect_RetDone
from OCC.Core.Interface import Interface_Static_SetCVal
from OCC.Core.STEPControl import (STEPControl_AsIs,
                                  STEPControl_GeometricCurveSet,
                                  STEPControl_Writer,)


def export_to_step(itemLst, filename):
    # ====================STEP-EXPORT=============================================
    # initialize the STEP exporter
    step_writer = STEPControl_Writer()
    # Interface_Static_SetCVal("write.step.schema", "AP203")
    # step_writer.SetTolerance(1e-6)

    # transfer shapes and write file
    for item in itemLst:
        step_writer.Transfer(item, STEPControl_AsIs)
    status = step_writer.Write(filename)
    assert status == IFSelect_RetDone

    return None
