#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 09:38:49 2018

@author: Tobias Pflumm
"""
import numpy as np
import numbers
import matplotlib.pyplot as plt
import itertools
import os

#PythonOCC Libraries
from OCC.gp import gp_Pnt, gp_Vec,gp_Dir
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire
from OCC.GeomAPI import GeomAPI_Interpolate
from OCC.Geom2dAPI import Geom2dAPI_Interpolate
from OCC.Geom import Geom_BezierCurve, Geom_Plane

#SONATA modules:
if __name__ == '__main__':
    os.chdir('..')
from SONATA.classPolar import Polar
from SONATA.cbm.topo.utils import TColgp_HArray1OfPnt_from_nparray, point_list_to_TColgp_Array1OfPnt, PntLst_to_npArray, TColgp_HArray1OfPnt2d_from_nparray
from SONATA.cbm.topo.BSplineLst_utils import BSplineLst_from_dct
from SONATA.cbm.topo.wire_utils import build_wire_from_BSplineLst, build_wire_from_BSplineLst2, get_wire_length, equidistant_Points_on_wire
    

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
    - The  __slots__ magic method  tells Python not to use a dict to store an 
    objects instance attributes, and only allocate space for a fixed set of 
    attributes.
    - A very brief parameter check is only performed at __init__ and is not 
    implemented via a getter and setter functionality to save function calls 
    and to preserve simplicity
     
    """
    class_counter= 1    #class attribute
    __slots__ = ( 'name', 'id', 'coordinates', 'polars', 'relative_thickness', 'wire', 'BSplineLst')
    
    def __init__(self, yml=None, name = 'NONAME', coordinates = None, polars = None, relative_thickness = None):
        self.name = 'NONAME'
        self.id = self.__class__.class_counter
        self.__class__.class_counter += 1
        
        self.coordinates = None
        self.polars = None
        self.relative_thickness = None
        self.wire = None
        self.BSplineLst = None
        
        if isinstance(yml, dict): 
            self.read_IAE37(yml)
            
        if isinstance(name, str) and not 'NONAME': 
            self.name = name

        if isinstance(coordinates, np.ndarray) and coordinates.shape[1] == 2 : 
            self.coordinates = coordinates

        if isinstance(polars, list) and all(isinstance(p, Polar) for p in polars):
            self.polars = polars

        if isinstance(relative_thickness, numbers.Real) and relative_thickness>=0:
            self.relative_thickness = relative_thickness


    def __repr__(self):
        """__repr__ is the built-in function used to compute the "official" 
        string reputation of an object, """
        return 'Airfoil: '+ self.name
    
    
    def read_IAE37(self, yml):
        """
        reads the IAE Wind Task 37 style Airfoil dictionary and assigsn them to
        the class attributes
        """
        self.name = yml['name']
        self.relative_thickness = yml['relative_thickness']
        self.coordinates = np.asarray([yml['coordinates']['x'],yml['coordinates']['y']]).T
        self.polars = [Polar(p) for p in yml['polars']]
        

    def gen_OCCtopo(self):
        """
        generates a Opencascade TopoDS_Wire and BSplineLst from the airfoil coordinates.
        This can be used for interpolation and surface generation
        """
        data = np.hstack((self.coordinates,np.zeros((self.coordinates.shape[0],1))))
        self.BSplineLst = BSplineLst_from_dct(data, angular_deflection = 30, closed=True, tol_interp=1e-6, twoD = False)
        #print(self.BSplineLst)
        self.wire = build_wire_from_BSplineLst(self.BSplineLst, twoD=False)
        return self.wire
        
    
    def transformed(self, airfoil2, k=0.5, n=500):
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
        #check if wire exists, else create it
        if self.wire == None:
            self.gen_OCCtopo()
        if airfoil2.wire == None:
            airfoil2.gen_OCCtopo()
        
        p1_lst = equidistant_Points_on_wire(self.wire, n)       
        p2_lst = equidistant_Points_on_wire(airfoil2.wire, n)
        v_lst = [gp_Vec(p1,p2) for p1,p2 in zip(p1_lst, p2_lst)]
        
        pres = []
        for i,p in enumerate(p1_lst):
            pres.append(p.Translated(v_lst[i].Multiplied(k)))
        
        trf_af = Airfoil()
        str_k = '%.3f' % k
        trf_af.name = 'TRF_'+self.name+'_'+airfoil2.name+'_'+str_k.replace('.','')
        trf_af.coordinates = PntLst_to_npArray(pres)[:,:2]      
        return trf_af
    
    
    def plot_polars(self, xlim = (-24,32), markercycle='.>^+*',):
        """
        plots the airfoil coordinates and the stored polars. 
        
        Parameters
        ----------
        xlim : tuple
            The user can define a lower and upper x-axis limit in deg.
        markercycle : str
            A string of marker symbols to be used when plotting multiple polars
            the markers would be repeated, if more polars that markers are set.

        """
        fig, ax = plt.subplots(2,2)
        fig.suptitle(self.name, fontsize=16)
        fig.subplots_adjust(wspace=0.25, hspace=0.25)
        
        #plot airfoil coordinates
        ax[0][0].plot(self.coordinates[:,0], self.coordinates[:,1], '.k--')
        ax[0][0].axis('equal')
        ax[0][0].set_xlabel('x/c')
        ax[0][0].set_ylabel('y/c')
        
        if self.polars:
            for po, ms in zip(self.polars, itertools.cycle(markercycle)):
                label = po.configuration + ', Re: ' + '{:.2e}'.format(po.re) + ', Ma: ' + str(po.ma)
                
                ax[0][1].plot(np.rad2deg(po.c_l[:,0]), po.c_l[:,1], label=label, marker = ms)
                ax[0][1].set_ylabel(r'Section lift coefficient, $c_l$')
    
                ax[1][0].plot(np.rad2deg(po.c_d[:,0]), po.c_d[:,1], label=label, marker = ms)
                ax[1][0].set_ylabel(r'Section drag coefficient, $c_d$')
    
                ax[1][1].plot(np.rad2deg(po.c_m[:,0]), po.c_m[:,1], label=label, marker = ms)
                ax[1][1].set_ylabel(r'Section drag coefficient, $c_d$')
    
                #Formatting
                for i,j in ((0,1),(1,0),(1,1)):
                    ax[i][j].set_xlabel(r'Section angle of attack, $\alpha$, deg')
                    ax[i][j].set_xlim(xlim[0],xlim[1])
                    ax[i][j].grid(b=True, which='major', color='k', linestyle='-')
                    ax[i][j].minorticks_on()
                    ax[i][j].grid(b=True, which='minor', color='k', linestyle=':', alpha=0.4)
                    ax[i][j].axhline(xmin=xlim[0],xmax=xlim[1], color='k', linestyle='-', linewidth=1.5)
                    ax[i][j].axvline(ymin=0,ymax=1, color='k', linestyle='-', linewidth=1.5)
                    ax[i][j].legend()
       
    
    def run_mses(self,re,ma):
        """
        run mses to calculate the polars.
        
        Notes
        ----------
        A possiblity in the future.
        """
        pass
    
    def run_xfoil(self,re,ma):
        """
        run xfoil to calculate the polars.
        
        Notes
        ----------
        A possiblity in the future.
        """
        pass
    
    
if __name__ == '__main__':
    plt.close('all')
    from jsonschema import validate
    import yaml


    
    with open('jobs/PBortolotti/IEAonshoreWT.yaml', 'r') as myfile:
        inputs  = myfile.read()
    with open('jobs/PBortolotti/IEAturbine_schema.yaml', 'r') as myfile:
        schema  = myfile.read()
    validate(yaml.load(inputs), yaml.load(schema))
    wt_data     = yaml.load(inputs)    
    
    airfoils = [Airfoil(af) for af in wt_data['airfoils']]
                        
    for af in airfoils:
        af.gen_OCCtopo()

    af1 = airfoils[0]
    af2 = airfoils[6]
    res = af1.transformed(af2, 0.5)
    res.gen_OCCtopo()

