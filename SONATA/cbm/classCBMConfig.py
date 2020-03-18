# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 14:53:09 2018

@author: TPflumm
"""
# Core Library modules
import os
from collections import OrderedDict

# Third party modules
import numpy as np
import yaml

# First party modules
from SONATA.cbm.fileIO.read_yaml_input import clean_filestring
from SONATA.classMaterial import find_material, read_IEA37_materials
from SONATA.vabs.classVABSConfig import VABSConfig

if __name__ == "__main__":
    os.chdir("/media/gu32kij/work/TPflumm/SONATA")


class CBMConfig(object):
    """
    Second Generation Configuration Class for the SONATA_CBM Disciplin  
    
    Attributes
    ----------
    filename : str
        filename, when config filename is given.
        
    setup: dict
        contains the following fields: input_type, datasource, material_db, 
        radial_station, Theta, scale_factor, BalanceWeight and mesh_resolution
    
    webs: OrderedDict
        Ordered Dictionary that contains subdictionaries of with Pos1 and Pos2 
        keys
    
    segments: OrderedDict
         Ordered Dictionary that contains subdictionaries with Corematerial, 
         Layup and Layup_names as keys
    
    bw: dict
    
    flags: dict

    vabs_cfg: VABSConfig

    
    """

    __slots__ = ("filename", "setup", "webs", "segments", "bw", "flags", "vabs_cfg")

    def __init__(self, inputdata=None, materials=None, iea37=False):
        self.setup, self.webs, self.segments, self.bw = {}, {}, {}, {}
        self.filename = ""

        if isinstance(inputdata, str) and iea37 == False:
            self.filename = inputdata
            self._read_yaml_config(self.filename)

        elif isinstance(inputdata, dict) and iea37 == True:
            yml = inputdata
            self._read_IEA37(yml, materials)

        self.vabs_cfg = VABSConfig()
        self.flags = {"mesh_core": True}

    def _read_IEA37(self, yml, materials):
        """ read the IEA37 style yml dictionary of a section and assign class 
        attributes to this configuration object
        
        Parameters:
        ----------
        yml : dict
           dictionary of the yaml style IEA37 section input    
        """

        # Setup:
        self.setup = {}
        self.setup["input_type"] = 5
        self.setup["datasource"] = None
        self.setup["material_db"] = None
        self.setup["radial_station"] = None
        self.setup["Theta"] = None  # GET FROM BLADE DEFINITION
        self.setup["scale_factor"] = 1
        self.setup["BalanceWeight"] = False
        self.setup["mesh_resolution"] = 280  # GET FROM function

        if yml.get("mesh_resolution"):
            self.setup["mesh_resolution"] = yml.get("mesh_resolution")

        # Webs:
        foo = {}
        if yml.get("webs"):
            for item in yml.get("webs"):
                w = {}
                w["id"] = item.get("id")
                w["Pos1"] = item["position"][0]
                w["Pos2"] = item["position"][1]
                w["curvature"] = item.get("curvature")
                foo[item.get("id")] = w
            self.webs = OrderedDict(sorted(foo.items(), key=lambda x: x[1]["id"]))

        # Segments
        self.segments = OrderedDict()
        for i, s in enumerate(yml.get("segments")):
            d = {}
            key = s.get("id")
            if s.get("filler") == None:
                d["CoreMaterial"] = 0

            elif isinstance(s.get("filler"), int):
                d["CoreMaterial"] = materials[s.get("filler")].id

            else:
                d["CoreMaterial"] = find_material(materials, "name", s.get("filler")).id

            layerlst = s.get("layup")
            if layerlst and all(isinstance(l, list) for l in layerlst):
                layerlst = s.get("layup")
                d["Layup_names"] = np.asarray(layerlst)[:, 5].tolist()
                d["Layup"] = np.asarray(layerlst)[:, :5].astype(np.float)

            elif layerlst and all(isinstance(l, dict) for l in layerlst):
                d["Layup"] = np.asarray([[l.get("start"), l.get("end"), l.get("thickness"), l.get("orientation"), find_material(materials, "name", l.get("material_name")).id] for l in layerlst])
                d["Layup_names"] = [l.get("name") for l in layerlst]

            else:
                d["Layup"] = np.empty((0, 0))
                d["Layup_names"] = np.empty((0, 0))

            self.segments[key] = d

        # BalanceWeight
        if yml.get("trim_mass"):
            self.bw = yml.get("trim_mass")  #
            self.setup["BalanceWeight"] = True

    def _read_yaml_config(self, fname):
        """ read the CBM config file and assign class attributes
        
        Parameters:
        ----------
        fname : str
            filename of the yaml style cbm input configuration    

        """
        b_string = clean_filestring(fname, comments="#")
        yDict = yaml.load(b_string)
        self.setup = yDict["Setup"]

        if "Webs" in yDict.keys():
            D = {int(k.split()[-1]): v for (k, v) in yDict["Webs"].items()}
            self.webs = OrderedDict(sorted(D.items()))

        # print(self.setup['BalanceWeight'])
        if self.setup["BalanceWeight"] == True:
            self.bw = yDict["BalanceWeight"]

        # read segments:
        D = {int(k.split()[-1]): v for (k, v) in yDict["Segments"].items()}

        for k in D:
            if D[k]["Layup"] != None:
                D[k]["Layup_names"] = np.asarray(D[k]["Layup"])[:, 5].tolist()
                D[k]["Layup"] = np.asarray(D[k]["Layup"])[:, :5].astype(np.float)
            else:
                D[k]["Layup"] = np.empty((0, 0))
                D[k]["Layup_names"] = np.empty((0, 0))

        self.segments = OrderedDict(sorted(D.items()))


if __name__ == "__main__":

    # classic configuration file:
    os.chdir("/media/gu32kij/work/TPflumm/SONATA")
    # fname = 'jobs/VariSpeed/advanced/sec_config.yml'
    fname = "examples/sec_config.yml"
    cfg = CBMConfig(fname)

    # IEA37 Style configuration:
    with open("jobs/PBortolotti/IEAonshoreWT.yaml", "r") as myfile:
        inputs = myfile.read()

    with open("jobs/VariSpeed/UH-60A_adv.yml", "r") as myfile:
        inputs = myfile.read()

    yml = yaml.load(inputs)
    materials = read_IEA37_materials(yml.get("materials"))
    yml = yml.get("components").get("blade").get("2d_fem").get("sections")[0]
    wt_cfg = CBMConfig(yml, materials, iea37=True)
