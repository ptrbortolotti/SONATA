# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 13:37:00 2017

@author: TPflumm
"""
import os

# Third party modules
from OCC.Core.IFSelect import IFSelect_RetDone
from OCC.Core.Interface import Interface_Static_SetCVal
from OCC.Core.STEPControl import (STEPControl_AsIs,
                                  STEPControl_GeometricCurveSet,
                                  STEPControl_Writer,)
from OCC.Core.IGESControl import IGESControl_Reader, IGESControl_Writer


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


def export_shape(itemLst, filename="export.step"):
        """
        Adapted from Airconics.Base
        Writes the Components in this Airconics shape to filename using
        file format specified in extension of filename.

        Parameters
        ----------
        filename : string
            the BASE.ext name of the file e.g. 'airliner.stp'.
            Note the Component name will be prepended to the base name of each
            output file

        single_export : bool
            Writes a single output file if true, otherwise writes one file
            per component

        Returns
        -------
        status : list of int
            error status of the file output of EACH component

        Notes
        -----
        File format is extracted from filename.

        stl file write will prepend filename onto the Component name to be
        written to file (cannot write multiple files )
        """
        path, ext = os.path.splitext(filename)

        # Default to a step writer if no extension type was provided:
        if not ext:
            ext = '.stp'

        status = []
        print('Writing to extension {}'.format(ext))
        print('Writing to file {}'.format(filename))

        if ext == '.stl':
            stl_ascii_format = False
            if single_export:
                stl_writer = StlAPI_Writer()
                for component in itemLst:
                    shape = component
                    status.append(stl_writer.Write(shape, filename))
            else:
                for component in itemLst:
                    stl_writer = StlAPI_Writer()
                    f = path + '_' + ext
                    shape = component
                    status.append(stl_writer.Write(shape, f))

        elif ext in ['.stp', '.step']:
            step_writer = STEPControl_Writer()
            application_protocol="AP214IS"
            Interface_Static_SetCVal("write.step.schema", application_protocol)
            for component in itemLst:
                step_writer.Transfer(component, STEPControl_AsIs)
            status = step_writer.Write(filename)

        elif ext in ['igs', '.iges']:
            iges_writer = IGESControl_Writer()
            for component in itemLst:
                iges_writer.AddShape(component)
            status = iges_writer.Write(filename)
        else:
            raise ValueError('Unexpected file extension {}'.format(ext))

        return status