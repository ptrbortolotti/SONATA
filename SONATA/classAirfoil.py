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
from urllib.request import urlopen
from collections import OrderedDict

#PythonOCC Libraries
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir, gp_Ax1, gp_Trsf, gp_Pln
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire
from OCC.Core.GeomAPI import GeomAPI_Interpolate
from OCC.Core.Geom2dAPI import Geom2dAPI_Interpolate
from OCC.Core.Geom import Geom_BezierCurve, Geom_Plane
from OCC.Display.SimpleGui import init_display


#SONATA modules:
if __name__ == '__main__':
    os.chdir('..')
from SONATA.utl.trsf import trsf_af_to_blfr
from SONATA.classPolar import Polar
from SONATA.cbm.topo.utils import TColgp_HArray1OfPnt_from_nparray, point_list_to_TColgp_Array1OfPnt, PntLst_to_npArray, TColgp_HArray1OfPnt2d_from_nparray
from SONATA.cbm.topo.wire_utils import rotate_wire, translate_wire, scale_wire, trsf_wire
from SONATA.cbm.topo.BSplineLst_utils import BSplineLst_from_dct, copy_BSplineLst, \
                                            equidistant_D1_on_BSplineLst
from SONATA.cbm.topo.wire_utils import build_wire_from_BSplineLst, build_wire_from_BSplineLst2, get_wire_length, equidistant_Points_on_wire

from SONATA.cbm.display.display_utils import export_to_JPEG, export_to_PNG, export_to_PDF, \
                                        export_to_SVG, export_to_PS, export_to_EnhPS, \
                                        export_to_TEX, export_to_BMP,export_to_TIFF, \
                                        show_coordinate_system, display_SONATA_SegmentLst,\
                                        display_custome_shape, transform_wire_2to3d, display_config

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
    class_counter= 1    #class attribute
    __slots__ = ( 'name', 'id', 'coordinates', 'polars', 'relative_thickness', 'wire', 'BSplineLst', 
                 'display', 'start_display', 'add_menu', 'add_function_to_menu')
    
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
            self.read_IEA37(yml)
            
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
    

    @property
    def te_coordinates(self):
        """
        returns the calculated trailing edge coordinates. Mean of the first and
        the last coordinate point
        """
        te = np.mean(np.vstack((self.coordinates[0],self.coordinates[-1])),axis=0)
        return te
    
    
    def read_IEA37(self, yml):
        """
        reads the IEA Wind Task 37 style Airfoil dictionary and assigns them to
        the class attributes
        """
        self.name = yml['name']
        self.relative_thickness = yml.get('relative_thickness')
        
        if yml['coordinates']:
            self.coordinates = np.asarray([yml['coordinates']['x'],yml['coordinates']['y']]).T
        else:
            self.get_UIUCCoordinates()        
        
        if self.polars:
            self.polars = [Polar(p) for p in yml['polars']]
        else:
            self.polars = None
        
    def get_UIUCCoordinates(self):   
        url = 'http://m-selig.ae.illinois.edu/ads/coord_seligFmt/%s.dat' % self.name
        try:
            with urlopen(url) as f:
                self.coordinates = np.loadtxt(f, skiprows=1)
        except:
            print('HTTPError: Not Found')
    
    def write_IEA37(self):
        """
        writes the class attributes to a dictionary conform with the IEA Wind 
        Task 37 style
        """ 
        tmp = {}
        tmp['name'] = self.name
        tmp['coordinates'] = dict(zip(('x','y'), self.coordinates.T.tolist()))
        tmp['relative_thickness'] = self.relative_thickness
        if self.polars != None:
            tmp['polars'] = [p.write_IEA37() for p in self.polars]
        else:
            tmp ['polars'] = None
        return tmp
        
    
    def gen_OCCtopo(self):
        """
        generates a Opencascade TopoDS_Wire and BSplineLst from the airfoil coordinates.
        This can be used for interpolation and surface generation
        
        
        Returns
        ---------
        self.wire : TopoDS_Wire
            wire of the airfoil
        
        """
        data = np.hstack((self.coordinates,np.zeros((self.coordinates.shape[0],1))))
        self.BSplineLst = BSplineLst_from_dct(data, angular_deflection = 30, closed=True, tol_interp=1e-5, twoD = False)
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
        
        wire = trsf_wire(self.wire,Trsf)
        tmp_pnt = gp_Pnt(self.te_coordinates[0], self.te_coordinates[1], 0)
        te_pnt = tmp_pnt.Transformed(Trsf)
#        print(self.BSplineLst)
#        BSplineLst_tmp = copy_BSplineLst(self.BSplineLst)
#        print(BSplineLst_tmp)
#        [s.Transform(Trsf) for s in BSplineLst_tmp]
        #print(BSplineLst)
        return (wire, te_pnt) #bspline, nodes, normals 
            
    
    def gen_wopwop_dist(self, NbPoints=50, divide_surf=True):
        """
        distributes points and normal vectors on the upper and lower part of 
        the airfoil. Subsequently those points are transformed to the blade 
        reference frame.
        
        
        Parameters
        ----------
        
        
        """
        
        data = self.coordinates
        
        le_idx = np.argmin(np.linalg.norm(data, axis=1))
        up_data = data[0:le_idx+1]
        lo_data = data[le_idx:-1]
        
        up_BSplineLst =  BSplineLst_from_dct(up_data, angular_deflection = 45, closed=False, tol_interp=1e-5, twoD = True)
        lo_BSplineLst = BSplineLst_from_dct(lo_data, angular_deflection = 45, closed=False, tol_interp=1e-5, twoD = True)
        
        up_pnts2d, up_vecs2d = equidistant_D1_on_BSplineLst(up_BSplineLst, NbPoints)
        lo_pnts2d, lo_vecs2d = equidistant_D1_on_BSplineLst(lo_BSplineLst, NbPoints)
        
        BSplineLst2d = up_BSplineLst+lo_BSplineLst
        pnts2d = up_pnts2d+lo_pnts2d
        vecs2d = up_vecs2d+lo_vecs2d
        
        return (BSplineLst2d, pnts2d, vecs2d)
        
    
    
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
     

    def post_3dviewer(self):
        (self.display, self.start_display, self.add_menu, self.add_function_to_menu) = display_config(cs_size = 0.3, DeviationAngle = 1e-7,  DeviationCoefficient = 1e-7)
        if self.wire == None:
            self.gen_OCCtopo()
            
        for s in self.BSplineLst:
            self.display.DisplayShape(s)
            
        #self.display.DisplayShape(self.wire)
        self.display.View_Top()
        self.display.FitAll()
        self.start_display()  
    
    
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
    import yaml

    #with open('jobs/PBortolotti/IEAonshoreWT.yaml', 'r') as f:
    with open('jobs/VariSpeed/UH-60A_adv.yml', 'r') as f:
        data = yaml.load(f.read())    
    
    airfoils = [Airfoil(af) for af in data['airfoils']]
                        
    for af in airfoils:
        af.gen_OCCtopo()

    af1 = airfoils[0]
    af2 = airfoils[1]
    res = af1.transformed(af2, 1.0)
    res.gen_OCCtopo()
    
    with open('data1.yml', 'w') as outfile:
        yaml.dump(af1.write_IEA37(), outfile)
        
    with open('data2.yml', 'w') as outfile:
        yaml.dump(af2.write_IEA37(), outfile)