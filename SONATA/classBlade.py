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

    __slots__ = ('coordinates', 'chord', 'twist', 'pitch_axis', 'blade_matrix', 'airfoilLst',  \
                 'sections', 'f_chord', 'f_twist', 'f_coordinates_x',  \
                 'f_coordinates_y', 'f_coordinates_z', 'f_pa')
    
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
        print('STATUS:\t Reading IAE37 Definition for Blade: %s' % (self.name))
        #Read information from DataDictionary
        tmp_coords = {}
        tmp = byml.get('bem_aero')
        tmp_coords['x'] = np.asarray((tmp.get('coordinates').get('x').get('grid'),tmp.get('coordinates').get('x').get('values'))).T
        tmp_coords['y'] = np.asarray((tmp.get('coordinates').get('y').get('grid'),tmp.get('coordinates').get('y').get('values'))).T
        tmp_coords['z'] = np.asarray((tmp.get('coordinates').get('z').get('grid'),tmp.get('coordinates').get('z').get('values'))).T
        tmp_tw = np.asarray((tmp.get('twist').get('grid'),tmp.get('twist').get('values'))).T
        tmp_chord = np.asarray((tmp.get('chord').get('grid'),tmp.get('chord').get('values'))).T
        tmp_pa = np.asarray((tmp.get('pitch_axis').get('grid'),tmp.get('pitch_axis').get('values'))).T
        airfoil_position = (tmp.get('airfoil_position').get('grid'),tmp.get('airfoil_position').get('labels'))
        
        #Generate Blade Matrix 
        tmp = []
        for an in airfoil_position[1]: 
            tmp.append(next((x for x in airfoils if x.name == an), None).id)
        arr = np.asarray([airfoil_position[0],tmp]).T
            
        self.f_chord = interp1d(tmp_chord[:,0], tmp_chord[:,1], bounds_error=False, fill_value='extrapolate')
        self.f_twist = interp1d(tmp_tw[:,0], tmp_tw[:,1], bounds_error=False, fill_value='extrapolate')
        self.f_coordinates_x = interp1d(tmp_coords['x'][:,0], tmp_coords['x'][:,1], bounds_error=False, fill_value='extrapolate')
        self.f_coordinates_y = interp1d(tmp_coords['y'][:,0], tmp_coords['y'][:,1], bounds_error=False, fill_value='extrapolate')
        self.f_coordinates_z = interp1d(tmp_coords['z'][:,0], tmp_coords['z'][:,1], bounds_error=False, fill_value='extrapolate')
        self.f_pa = interp1d(tmp_pa[:,0], tmp_pa[:,1], bounds_error=False, fill_value='extrapolate')
        
        cs_pos = np.asarray([cs.get('position') for cs in byml.get('2d_fem').get('sections')])
        x = np.unique(np.sort(np.hstack((tmp_chord[:,0], tmp_tw[:,0], tmp_coords['x'][:,0], tmp_coords['y'][:,0], tmp_coords['z'][:,0], tmp_pa[:,0], arr[:,0], cs_pos))))
        
        self.blade_matrix = np.transpose(np.unique(np.array([x, self.f_coordinates_x(x), self.f_coordinates_y(x), self.f_coordinates_z(x), self.f_chord(x), self.f_twist(x), self.f_pa(x)]), axis=1))
        self.airfoilLst = [interp_airfoil_position(airfoil_position, airfoils, x) for x in x]
        self.coordinates = self.blade_matrix[:,0:4]
        self.chord = self.blade_matrix[:,[0,4]]
        self.twist = self.blade_matrix[:,[0,5]]
        self.pitch_axis = self.blade_matrix[:,[0,6]]         
        
        #get sections information and init the CBM instances.
        tmp = byml.get('2d_fem').get('sections')
        self.sections = {}
        for cs in tmp:            
            cbm_config = CBMConfig(cs, materials, iae37=True)
            x = cs.get('position')
            bm = np.array([x, self.f_coordinates_x(x), self.f_coordinates_y(x), self.f_coordinates_z(x), self.f_chord(x), self.f_twist(x), self.f_pa(x)])
            af = interp_airfoil_position(airfoil_position, airfoils, x)
            self.sections[x] = CBM(cbm_config, materials = materials, blade_matrix = bm, airfoil=af)
            #get blade_matrix and airfoil at position and pass to CBM init instance!
            
        return (self.blade_matrix, self.airfoilLst)


    def plot_blade_matrix(self):
        """
        plot the the coordinates, chord, twist of the blade
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
    
    
    def post_3dtopo(self, flag_wf = True, flag_lft = False, flag_topo = True, flag_mesh = False, **kwargs):
        """
        generates the wireframe and the loft surface of the blade

        Returns
        ----------
        loft : OCC.TopoDS_surface
            the 3D surface of the blade
        
        wireframe : list
            list of every airfoil_wire scaled and rotated at every grid point
        """
        display, start_display, add_menu, add_function_to_menu = init_display()
        display.Context.SetDeviationAngle(1e-5) # 0.001 default. Be careful to scale it to the problem.
        display.Context.SetDeviationCoefficient(1e-5) # 0.001 default. Be careful to scale it to the problem. 
        bg_c = ((20,6,111),(200,200,200))
        display.set_bg_gradient_color(bg_c[0][0],bg_c[0][1],bg_c[0][2],bg_c[1][0],bg_c[1][1],bg_c[1][2])
        show_coordinate_system(display,4)
        
        #display WIREFRAME
        if flag_wf:
            wireframe = []
            for bm, afl in zip(self.blade_matrix, self.airfoilLst):
                if afl.wire == None:
                    afl.gen_OCCtopo()
                wire = afl.wire
                
                wire = rotate_wire(wire, gp_Ax1(gp_Pnt(bm[6],0,0), gp_Dir(0,0,1)), -bm[5], copy=True)
                wire = translate_wire(wire, gp_Pnt(bm[6],0,0), gp_Pnt(0,0,0))
                wire = scale_wire(wire, gp_Pnt(0,0,0), bm[4])
                
                wire = rotate_wire(wire, gp_Ax1(gp_Pnt(0,0,0), gp_Dir(0,0,1)), -np.pi/2)
                wire = rotate_wire(wire, gp_Ax1(gp_Pnt(0,0,0), gp_Dir(0,1,0)), -np.pi/2)
                wire = translate_wire(wire, gp_Pnt(0,0,0), gp_Pnt(bm[1],bm[2],bm[3]))
                
                wireframe.append(wire)
                display.DisplayShape(wire, color='WHITE', transparency=0.3)
            
        if flag_lft:
            loft = make_loft(wireframe, ruled=False, tolerance=1e-3, continuity=1, check_compatibility=True)
            #display.DisplayShape(loft)
        
        if flag_topo:
            for k,cs in B.sections.items():
                coord = self.f_coordinates_x(k), self.f_coordinates_y(k), self.f_coordinates_z(k)
                display_SONATA_SegmentLst(display, cs.SegmentLst, coord, -np.pi/2, -np.pi/2)
        
        display.View_Iso()
        display.FitAll()
        start_display()   
        
        return None         
        

if __name__ == '__main__':
    import yaml  
        
    with open('jobs/PBortolotti/IEAonshoreWT.yaml', 'r') as myfile:
        inputs  = myfile.read()
        
    yml = yaml.load(inputs)
    airfoils = [Airfoil(af) for af in yml.get('airfoils')]
    materials = read_IAE37_materials(yml.get('materials'))
    
    byml = yml.get('components').get('blade')
    B = Blade(name='TestBlade')
    B.read_IAE37(byml, airfoils, materials)

#    for key, cs in B.sections.items():
#        print('STATUS:\t Building Section at grid location %s' % (key))
#        cs.cbm_gen_topo()
#        cs.cbm_gen_mesh()
#        cs.cbm_run_vabs()
#        cs.cbm_post_2dmesh()

    #%%====== DISPLAY ==============
    #B.sections[0].cbm_gen_topo()
    B.post_3dtopo()
