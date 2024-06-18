#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 09:38:49 2018

@author: Tobias Pflumm
"""
# Core Library modules
import itertools
import numbers
import os

# Third party modules
import matplotlib.pyplot as plt
import numpy as np
# PythonOCC Libraries
from OCC.Core.gp import gp_Pnt, gp_Vec

# First party modules
from SONATA.cbm.display.display_utils import (display_config,)
from SONATA.cbm.topo.BSplineLst_utils import (BSplineLst_from_dct,)
from SONATA.cbm.topo.utils import (PntLst_to_npArray,)
from SONATA.cbm.topo.wire_utils import (build_wire_from_BSplineLst,
                                        equidistant_Points_on_wire,
                                        trsf_wire,)
from SONATA.utl.trsf import trsf_af_to_blfr

# SONATA modules:
if __name__ == "__main__":
    os.chdir("..")



class Airfoil(object):
    """
    Airfoil Class object.
    
    Attributes
    ----------
    name : string
        The name of the airfoil
    
    coordinates : np.array
        The x and y coordinates of the airfoil as np.array of shape (...,2)
        x coordinates should be defined between 0 and 1.
        
    polars: list
        A list of Polar instances. That store c_l, c_d, and c_m together with 
        Reynolds Number re and Machnumber ma.
        
    relative thickness : float
        Relative Thickness of the airfoil
        
    Notes
    ----------
    - A very brief parameter check is only performed at __init__ and is not 
    implemented via a getter and setter functionality to save function calls 
    and to preserve simplicity
     
    """

    class_counter = 1  # class attribute
    __slots__ = ("name", "id", "coordinates", "polars", "relative_thickness", "wire", "BSplineLst", "display", "start_display", "add_menu", "add_function_to_menu")

    def __init__(self, yml=None, name="NONAME", coordinates=None, polars=None, relative_thickness=None):
        self.name = "NONAME"
        self.id = self.__class__.class_counter
        self.__class__.class_counter += 1

        self.coordinates = None
        self.polars = None
        self.relative_thickness = None
        self.wire = None
        self.BSplineLst = None

        if isinstance(yml, dict):
            self.read_yaml_airfoil(yml)

        if isinstance(name, str) and not "NONAME":
            self.name = name

        if isinstance(coordinates, np.ndarray) and coordinates.shape[1] == 2:
            self.coordinates = coordinates

        if isinstance(relative_thickness, numbers.Real) and relative_thickness >= 0:
            self.relative_thickness = relative_thickness

    def __repr__(self):
        """__repr__ is the built-in function used to compute the "official" 
        string reputation of an object, """
        return "Airfoil: " + self.name

    @property
    def te_coordinates(self):
        """
        returns the calculated trailing edge coordinates. Mean of the first and
        the last coordinate point
        """
        te = np.mean(np.vstack((self.coordinates[0], self.coordinates[-1])), axis=0)
        return te

    def read_yaml_airfoil(self, yml):
        """
        reads the Airfoil dictionary from the yaml dictionary and assigns them to
        the class attributes
        """
        self.name = yml["name"]
        self.relative_thickness = yml.get("relative_thickness")

        self.coordinates = np.asarray([yml["coordinates"]["x"], yml["coordinates"]["y"]], dtype=float).T
    def gen_OCCtopo(self, angular_deflection = 30 ):
        """
        generates a Opencascade TopoDS_Wire and BSplineLst from the airfoil coordinates.
        This can be used for interpolation and surface generation
        
        
        Returns
        ---------
        self.wire : TopoDS_Wire
            wire of the airfoil
        
        """
        data = np.hstack((self.coordinates, np.zeros((self.coordinates.shape[0], 1))))
        self.BSplineLst = BSplineLst_from_dct(data, angular_deflection=angular_deflection, closed=True, tol_interp=1e-5, twoD=False)
        # print("STATUS:\t CHECK Head2Tail: \t\t ", Check_BSplineLst_Head2Tail(self.BSplineLst))
        # print("STATUS:\t CHECK Counterclockwise: \t ", BSplineLst_Orientation(self.BSplineLst, 11))
        
        self.wire = build_wire_from_BSplineLst(self.BSplineLst, twoD=False)
        return self.wire

    def trsf_to_blfr(self, loc, pa_loc, chord, twist):
        """
        transforms the nondim. airfoil to the blade reference frame location
        and pitch-axis information, scales it with chord information and rotates 
        it with twist information
        
        Parameters
        ----------
        loc : array
            [x,y,z] position in blade reference coordinates
        pa_loc : float
            nondim. pitch axis location
        chord : float
            chordlength
        twist : float
            twist angle about x in radians
        
        Retruns
        ---------
        wire : TopoDS_Wire
        
        """
        Trsf = trsf_af_to_blfr(loc, pa_loc, chord, twist)
        if self.wire == None or self.BSplineLst == None:
            self.gen_OCCtopo()

        wire = trsf_wire(self.wire, Trsf)
        tmp_pnt = gp_Pnt(self.te_coordinates[0], self.te_coordinates[1], 0)
        te_pnt = tmp_pnt.Transformed(Trsf)
        return (wire, te_pnt)  # bspline, nodes, normals

    def transformed(self, airfoil2, k=0.5, n=200):
        """
        Performs and linear interpolation of the airfoil with another airfoil 
        by translating equidistant points in the direction of vector v. 
        The magnitude of the translation is the vector's magnitude multiplied 
        by factor k.
    
        Parameters
        ----------
        airfoil2 : Airfoil
            The airfoil the user whats the current airfoil to be transformed to
        k : float
            the vectors magnitude factor. k=0: the transformed airfoil remains 
            the airfoil. k=1: the transformed airfoil will become airfoil2
        n : int
            number of discretization points.
        
        Returns:
        ----------
        trf_af : Airfoil
            with the name = TRF_airfoil1_airfoil2_k
            
        """
        # check if wire exists, else create it
        if self.wire == None:
            self.gen_OCCtopo()
        if airfoil2.wire == None:
            airfoil2.gen_OCCtopo()

        p1_lst = equidistant_Points_on_wire(self.wire, n)
        p2_lst = equidistant_Points_on_wire(airfoil2.wire, n)
        v_lst = [gp_Vec(p1, p2) for p1, p2 in zip(p1_lst, p2_lst)]

        pres = []
        for i, p in enumerate(p1_lst):
            pres.append(p.Translated(v_lst[i].Multiplied(k)))

        trf_af = Airfoil()
        str_k = "%.3f" % k
        trf_af.name = "TRF_" + self.name + "_" + airfoil2.name + "_" + str_k.replace(".", "")
        trf_af.coordinates = PntLst_to_npArray(pres)[:, :2]
        # Shifting the largest x value coordinate to the first position in the airfoil coordinates so the TE origin can be defined
        max_x_index = max(range(len(trf_af.coordinates)), key = lambda i: trf_af.coordinates[i][0])+1
        trf_af.coordinates = np.vstack((trf_af.coordinates[max_x_index:], trf_af.coordinates[: max_x_index]))
        return trf_af

if __name__ == "__main__":
    plt.close("all")
    import yaml

    with open("jobs/VariSpeed/UH-60A_adv.yml", "r") as f:
        data = yaml.load(f.read())

    airfoils = [Airfoil(af) for af in data["airfoils"]]

    for af in airfoils:
        af.gen_OCCtopo()

    af1 = airfoils[0]
    af2 = airfoils[1]
    res = af1.transformed(af2, 1.0)
    res.gen_OCCtopo()

    with open("data1.yml", "w") as outfile:
        yaml.dump(af1.write_yaml_airfoil(), outfile)

    with open("data2.yml", "w") as outfile:
        yaml.dump(af2.write_yaml_airfoil(), outfile)
