#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 09:39:06 2018

@author: Tobias Pflumm
"""
# Core Library modules
import numbers

# Third party modules
import numpy as np


class Polar(object):
    """
    Airfoil Polar Class object. 
    
    Attributes
    ----------
    re : float
        The Reynolds-Number of the polar
        
    ma : float
        The Mach-Number of the polar
    
    c_l, c_d, c_m : np.ndarray
        Section aerodynamic lift,drag and moment coefficients
        np.array with shape ( ... , 2)
        The first collumn is the angle of attack, alpha between ~ [-pi,+pi]
        And the second collumn stores the lift, drag and moment coefficent values
        
    configuration : string
        Attribute to store additional information about the polar.    
        
    Notes
    ----------
    - The  __slots__ magic method  tells Python not to use a dict to store an 
    objects instance attributes, and only allocate space for a fixed set of 
    attributes.
    - A very brief parameter check is only performed at __init__ and is not 
    implemented via a getter and setter functionality to save function calls 
    and to preserve simplicity
     
    """

    __slots__ = ("re", "ma", "c_l", "c_d", "c_m", "configuration")

    def __init__(self, yml=None, c_l=None, c_d=None, c_m=None, re=None, ma=None, configuration=None):
        """initialization method that is called every time an instance is created
        
        """
        self.c_l = None
        self.c_d = None
        self.c_m = None
        self.re = None
        self.ma = None
        self.configuration = ""

        if isinstance(yml, dict):
            self.read_IEA37(yml)

        if self._check_coefficients(c_l):
            self.c_l = c_l

        if self._check_coefficients(c_d):
            self.c_d = c_d

        if self._check_coefficients(c_m):
            self.c_m = c_m

        if isinstance(re, numbers.Real) and re >= 0:
            self.re = re

        if isinstance(ma, numbers.Real) and ma >= 0:
            self.ma = ma

        if isinstance(configuration, str):
            self.configuration = configuration

    def read_IEA37(self, yml):
        """
        procedure that reads the IEA Wind Task 37 style Polar dictionary and assigns them to
        the class attributes
        """
        self.configuration = yml.get("configuration")
        self.re = yml.get("re")
        self.ma = yml.get("ma")
        self.c_l = np.asarray([yml["c_l"]["grid"], yml["c_l"]["values"]]).T
        self.c_d = np.asarray([yml["c_d"]["grid"], yml["c_d"]["values"]]).T
        self.c_m = np.asarray([yml["c_m"]["grid"], yml["c_m"]["values"]]).T

    def write_IEA37(self):
        tmp = {}
        tmp["configuration"] = self.configuration
        tmp["re"] = self.re
        tmp["ma"] = self.ma
        tmp["c_l"] = dict(zip(("grid", "values"), self.c_l.T.tolist()))
        tmp["c_d"] = dict(zip(("grid", "values"), self.c_d.T.tolist()))
        tmp["c_m"] = dict(zip(("grid", "values"), self.c_m.T.tolist()))

        return tmp

    def _check_coefficients(self, x):
        """checks the aerodynamic coefficients for correctness"""
        if isinstance(x, np.ndarray) and x.shape[1] == 2:
            if all((x[:, 0] >= -np.pi - 1e-4) & (x[:, 0] <= np.pi + 1e-4)):
                return True
            else:
                return False
                print("WARNING: Attribute not changed; Polar.coefficients" + " AoA values have to be between ~ [-pi,+pi] ")
        else:
            return False
            print("WARNING: Attribute not changed; Polar.coefficients has" + " to be a np.array with shape ( ... , 2)")

    def interpolate(self, grid):
        """procedure that interpolates data to a specific grid and returns the  the 
        c_l, c_d, and c_m values."""
        pass


if __name__ == "__main__":
    import os

    os.chdir("..")
    from jsonschema import validate
    import yaml

    with open("jobs/PBortolotti/IEAonshoreWT.yaml", "r") as myfile:
        inputs = myfile.read()
    with open("jobs/PBortolotti/IEAturbine_schema.yaml", "r") as myfile:
        schema = myfile.read()
    validate(yaml.load(inputs), yaml.load(schema))
    wt_data = yaml.load(inputs)

    af = wt_data["airfoils"][0]
    polars = [Polar(p) for p in af["polars"]]
