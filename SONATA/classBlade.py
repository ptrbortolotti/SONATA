#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 09:38:19 2018

@author: Tobias Pflumm
"""
import os

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

from OCC.gp import gp_Ax2, gp_Pnt, gp_Dir, gp_Ax1
from OCC.Display.SimpleGui import init_display

if __name__ == '__main__':
    os.chdir('..')
    
from SONATA.classComponent import Component
from SONATA.classAirfoil import Airfoil
from SONATA.classMaterial import read_IAE37_materials

from SONATA.cbm.classCBM import CBM
from SONATA.cbm.classCBMConfig import CBMConfig

from SONATA.blade_utl import interp_airfoil_position, make_loft
from SONATA.cbm.topo.wire_utils import rotate_wire, translate_wire, scale_wire

from SONATA.cbm.display.display_utils import export_to_JPEG, export_to_PNG, export_to_PDF, \
                                        export_to_SVG, export_to_PS, export_to_EnhPS, \
                                        export_to_TEX, export_to_BMP,export_to_TIFF, \
                                        show_coordinate_system, display_SONATA_SegmentLst,\
                                        display_custome_shape, transform_wire_2to3d  




class Blade(Component):

    __slots__ = ('coordinates', 'chord', 'twist', 'blade_matrix', 'airfoilLst', 'sections')
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args,**kwargs)
    
    
    def __repr__(self):
        """__repr__ is the built-in function used to compute the "official" 
        string reputation of an object, """
        return 'Blade: '+ self.name
    
    
    def read_IAE37(self, yml, airfoils, materials):
        """
        reads the IAE Wind Task 37 style Blade dictionary 
        generates the blade matrix and airfoilLst to represent all given 
        information at every grid point by interpolating the input data 
        and assigsn them to the class attribute twist, choord, coordinates 
        and airfoil_positions with the first column representing the 
        non-dimensional x-location  

        Parameters
        ----------
        airfoils : list
            Is the database of airfoils
        
        Returns
        ----------
        blade_matrix : np.ndarray
            The blade_matrix with the follwowing colloums 
            (1: non-dimensionalized x-station, 2:chord, 3:twist, x-coordinate, 
            y-coordinate, z-coordinate)
        
        airfoilLst : list
            list of airfoils at every grid point (row of blade_matrix) with 
            lin. interpolated airfoils. 
        
        """
        #Read information from DataDictionary
        tmp_coords = {}
        tmp = byml.get('bem_aero')
        tmp_coords['x'] = np.asarray((tmp.get('coordinates').get('x').get('grid'),tmp.get('coordinates').get('x').get('values'))).T
        tmp_coords['y'] = np.asarray((tmp.get('coordinates').get('y').get('grid'),tmp.get('coordinates').get('y').get('values'))).T
        tmp_coords['z'] = np.asarray((tmp.get('coordinates').get('z').get('grid'),tmp.get('coordinates').get('z').get('values'))).T
        tmp_tw = np.asarray((tmp.get('twist').get('grid'),tmp.get('twist').get('values'))).T
        tmp_chord = np.asarray((tmp.get('chord').get('grid'),tmp.get('chord').get('values'))).T
        airfoil_position = (tmp.get('airfoil_position').get('grid'),tmp.get('airfoil_position').get('labels'))
        
        #Generate Blade Matrix 
        tmp = []
        for an in airfoil_position[1]: 
            tmp.append(next((x for x in airfoils if x.name == an), None).id)
        arr = np.asarray([airfoil_position[0],tmp]).T
            
        f_chord = interp1d(tmp_chord[:,0], tmp_chord[:,1], bounds_error=False, fill_value='extrapolate')
        f_twist = interp1d(tmp_tw[:,0], tmp_tw[:,1], bounds_error=False, fill_value='extrapolate')
        f_coordinates_x = interp1d(tmp_coords['x'][:,0], tmp_coords['x'][:,1], bounds_error=False, fill_value='extrapolate')
        f_coordinates_y = interp1d(tmp_coords['y'][:,0], tmp_coords['y'][:,1], bounds_error=False, fill_value='extrapolate')
        f_coordinates_z = interp1d(tmp_coords['z'][:,0], tmp_coords['z'][:,1], bounds_error=False, fill_value='extrapolate')
        
    
        cs_pos = np.asarray([cs.get('position') for cs in byml.get('2d_fem').get('sections')])
        x = np.unique(np.sort(np.hstack((tmp_chord[:,0], tmp_tw[:,0], tmp_coords['x'][:,0], tmp_coords['y'][:,0], tmp_coords['z'][:,0], arr[:,0], cs_pos))))
        
        self.blade_matrix = np.transpose(np.unique(np.array([x, f_coordinates_x(x), f_coordinates_y(x), f_coordinates_z(x), f_chord(x), f_twist(x)]), axis=1))
        self.airfoilLst = [interp_airfoil_position(airfoil_position, airfoils, x) for x in x]
        self.coordinates = self.blade_matrix[:,0:4]
        self.chord = self.blade_matrix[:,[0,4]]
        self.twist = self.blade_matrix[:,[0,5]]
                
        #get sections information and init the CBM instances.
        tmp = byml.get('2d_fem').get('sections')
        self.sections = {}
        for cs in tmp:
            cbm_config = CBMConfig(cs, materials, iae37=True)
            x = cs.get('position')
            bm = np.array([x, f_coordinates_x(x), f_coordinates_y(x), f_coordinates_z(x), f_chord(x), f_twist(x)])
            af = interp_airfoil_position(airfoil_position, airfoils, x)
            self.sections[x] = CBM(cbm_config, materials = materials, blade_matrix = bm, airfoil=af)
            #get blade_matrix and airfoil at position and pass to CBM init instance!
            
        return (self.blade_matrix, self.airfoilLst)


    def plot_blade_matrix(self):
        """
        illustrates the coordinates, chord, twist of the blade
        """
        
        fig, ax = plt.subplots(3,2)
        fig.suptitle(self.name, fontsize=16)
        fig.subplots_adjust(wspace=0.25, hspace=0.25)
        
        ax[0][0].plot(self.blade_matrix[:,0], self.blade_matrix[:,1], 'k.-')
        ax[0][0].set_ylabel('x-coordinate [m]')
        
        ax[1][0].plot(self.blade_matrix[:,0], self.blade_matrix[:,2], 'k.-')
        ax[1][0].set_ylabel('y-coordinate [m]')
        
        ax[2][0].plot(self.blade_matrix[:,0], self.blade_matrix[:,3], 'k.-')
        ax[2][0].set_ylabel('z-coordinate [m]')
        
        ax[0][1].plot(self.blade_matrix[:,0], self.blade_matrix[:,4], 'k.-')
        ax[0][1].set_ylabel('chord [m]')
        
        ax[1][1].plot(self.blade_matrix[:,0], self.blade_matrix[:,5], 'k.-')
        ax[1][1].set_ylabel('twist [rad]')
        
#        ax3d = fig.add_subplot(326, projection='3d')
#        for bm, af in zip(self.blade_matrix, self.airfoilLst):
#            tmp_shape = af.coordinates[:,0].shape
#            arr = af.coordinates*bm[4]
#            ax3d.plot(np.ones(tmp_shape)*bm[1],arr[:,0],arr[:,1])
        plt.show()
    
    
    def gen_surface(self):
        """
        generates the wireframe and the loft surface of the blade

        Returns
        ----------
        loft : OCC.TopoDS_surface
            the 3D surface of the blade
        
        wireframe : list
            list of every airfoil_wire scaled and rotated at every grid point
        """
        
        wireframe = []
        for bm, afl in zip(self.blade_matrix, self.airfoilLst):
            if afl.wire == None:
                afl.gen_OCCtopo()
            
            wire = afl.wire
            wire = rotate_wire(wire, gp_Ax1(gp_Pnt(0,0,0), gp_Dir(1,0,0)), np.pi/2)
            wire = rotate_wire(wire, gp_Ax1(gp_Pnt(0,0,0), gp_Dir(0,0,1)), -np.pi/2)
            wire = translate_wire(wire, gp_Pnt(0,0,0), gp_Pnt(0,0.25,0))
            wire = scale_wire(wire, gp_Pnt(0,0,0), bm[4])
            wire = rotate_wire(wire, gp_Ax1(gp_Pnt(0,0,0), gp_Dir(1,0,0)), bm[5])  
            wire = translate_wire(wire, gp_Pnt(0,0,0), gp_Pnt(bm[1],bm[2],bm[3]))
            wireframe.append(wire)
            display.DisplayShape(wire, color='BLACK')
            #loft = make_loft(wireframe[:12], ruled=False, tolerance=1e-3, continuity=1, check_compatibility=True)
        return wireframe         
        

if __name__ == '__main__':
    import yaml  
        
    with open('jobs/PBortolotti/IEAonshoreWT.yaml', 'r') as myfile:
        inputs  = myfile.read()
    with open('jobs/PBortolotti/IEAturbine_schema.yaml', 'r') as myfile:
        schema  = myfile.read()
    
    yml = yaml.load(inputs)
    yml_schema = yaml.load(schema)
    
    airfoils = [Airfoil(af) for af in yml.get('airfoils')]
    materials = read_IAE37_materials(yml.get('materials'))
    
    byml = yml.get('components').get('blade')
    B = Blade(name='TestBlade')
    B.read_IAE37(byml, airfoils, materials)

    for cs in B.sections.values():
        cs.cbm_gen_topo()
        cs.cbm_gen_mesh()
        cs.cbm_post_2dmesh()

    # ====== DISPLAY ==============
    display, start_display, add_menu, add_function_to_menu = init_display()
    p = gp_Pnt(0,0,0)
    display.DisplayShape(p)
    show_coordinate_system(display,3)
    B.gen_surface()
    display.View_Iso()
    display.FitAll()
    start_display()   