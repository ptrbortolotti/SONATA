#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 09:38:19 2018

@author: Tobias Pflumm
"""
import os

import numpy as np
import yaml  
from jsonschema import validate

from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

from OCC.gp import gp_Ax2, gp_Pnt, gp_Dir, gp_Ax1
from OCC.Display.SimpleGui import init_display

if __name__ == '__main__':
    os.chdir('..')
    
from SONATA.classComponent import Component
from SONATA.classAirfoil import Airfoil
from SONATA.classMaterial import read_IEA37_materials

from SONATA.cbm.classCBM import CBM
from SONATA.cbm.classCBMConfig import CBMConfig

from SONATA.vabs.classVABSConfig import VABSConfig

from SONATA.utl.blade_utl import interp_airfoil_position, make_loft, interp_loads
from SONATA.utl.plot import plot_beam_properties
from SONATA.utl.converter import iea37_converter
from SONATA.cbm.topo.wire_utils import rotate_wire, translate_wire, scale_wire

from SONATA.cbm.display.display_utils import export_to_JPEG, export_to_PNG, export_to_PDF, \
                                        export_to_SVG, export_to_PS, export_to_EnhPS, \
                                        export_to_TEX, export_to_BMP,export_to_TIFF, \
                                        show_coordinate_system, display_SONATA_SegmentLst,\
                                        display_custome_shape, transform_wire_2to3d, display_config

class Blade(Component):
    """
    SONATA Blade component object.
    
    Attributes
    ----------                      
    coordinates :  ndarray
        Describes the axis LE coordinates in meters along the span.
        nparray([[grid, x, y, z]]).
        The grid represents the nondimensional x position along the Blade from 
        0 to 1
    
    chord : ndarray
        Describes the blades chord lenght in meters in spanwise direction. 
        nparray([[grid, chord]]) 

    twist : ndarray
        Describes the blades twist angles in !radians! in spanwise direction. 
        nparray([[grid, twist]]) 

    pitch_axis : ndarray
        Describes the blades pitch-axis location in 1/chord lengths from the 
        leading edge. nparray([[grid, pitch_axis]]) 

    airfoils : ndarray
        array of grid location and airfoil instance 
        nparray([[grid, airfoil instance]],dtype = object)
        
    sections : ndarray
        array of CBM cross-sections 
        nparray([[grid, CBM instance]],dtype = object)
        
    beam_properties : ndarray
        array of grid location and VABSSectionalProp instance
        nparray([[grid, beam_properties]],dtype = object)
              
        
    Methods
    -------
    blade_matrix : ndarray
        Summons all the blades global properties in one array
        nparray([[grid, x, y, z, chord, twist, pitch_axis,....]])
    
    
    Notes
    --------
    Units: meter (m), Newton (N), kilogramm (kg), degree (deg), Kelvin (K),



    See Also
    --------
    Component,
    

    ToDo
    -----
    - Include the possibity to rotate the beam_properties non-twisted frame. 
        Default is the twisted frame
    -
    

    Examples
    --------
    Initialize Blade Instance:
    
    >>> job = Blade(name='UH-60A_adv')
    
    >>> job.read_IEA37(yml.get('components').get('blade'), airfoils, materials)  
    
    >>> job.blade_gen_section()
    >>> job.blade_run_vabs()
    >>> job.blade_plot_sections()
    >>> job.blade_post_3dtopo(flag_lft = True, flag_topo = True)

    """ 
    
    __slots__ = ('coordinates', 'chord', 'twist', 'pitch_axis', 'airfoils',  \
                 'sections', 'beam_properties', 'f_chord', 'f_twist', 'f_coordinates_x',  \
                 'f_coordinates_y', 'f_coordinates_z', 'f_pa', \
                 'display', 'start_display', 'add_menu', 'add_function_to_menu')
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args,**kwargs)
        self.beam_properties = None
        
    def __repr__(self):
        """__repr__ is the built-in function used to compute the "official" 
        string reputation of an object, """
        return 'Blade: '+ str(self.name)
    
    
    def read_IEA37(self, yml, airfoils, materials, stations = [], npts = 11, wt_flag=False):
        """
        reads the IEA Wind Task 37 style Blade dictionary 
        generates the blade matrix and airfoil to represent all given 
        information at every grid point by interpolating the input data 
        and assigsn them to the class attribute twist, choord, coordinates 
        and airfoil_positions with the first column representing the 
        non-dimensional x-location  

        Parameters
        ----------
        airfoils : list
            Is the database of airfoils
        
        
        """
        self.name = yml.get('name')
        print('STATUS:\t Reading IAE37 Definition for Blade: %s' % (self.name))
        #Read information from DataDictionary

        tmp_coords = {}

        tmp_coords['x'] = np.asarray((yml.get('outer_shape_bem').get('reference_axis').get('x').get('grid'),yml.get('outer_shape_bem').get('reference_axis').get('x').get('values'))).T
        tmp_coords['y'] = np.asarray((yml.get('outer_shape_bem').get('reference_axis').get('y').get('grid'),yml.get('outer_shape_bem').get('reference_axis').get('y').get('values'))).T
        tmp_coords['z'] = np.asarray((yml.get('outer_shape_bem').get('reference_axis').get('z').get('grid'),yml.get('outer_shape_bem').get('reference_axis').get('z').get('values'))).T
        tmp_tw = np.asarray((yml.get('outer_shape_bem').get('twist').get('grid'),yml.get('outer_shape_bem').get('twist').get('values'))).T
        tmp_chord = np.asarray((yml.get('outer_shape_bem').get('chord').get('grid'),yml.get('outer_shape_bem').get('chord').get('values'))).T
        tmp_pa = np.asarray((yml.get('outer_shape_bem').get('pitch_axis').get('grid'),yml.get('outer_shape_bem').get('pitch_axis').get('values'))).T
        airfoil_position = (yml.get('outer_shape_bem').get('airfoil_position').get('grid'),yml.get('outer_shape_bem').get('airfoil_position').get('labels'))
        
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
        
        if wt_flag:
            if stations == []:
                cs_pos = np.linspace(0.0, 1.0, npts)
            else:
                cs_pos = stations
        else:
            cs_pos = np.asarray([cs.get('position') for cs in yml.get('internal_structure_2d_fem').get('sections')])
            
        x = np.unique(np.sort(np.hstack((tmp_chord[:,0], tmp_tw[:,0], tmp_coords['x'][:,0], tmp_coords['y'][:,0], tmp_coords['z'][:,0], tmp_pa[:,0], arr[:,0], cs_pos))))
        
        blade_matrix = np.transpose(np.unique(np.array([x, self.f_coordinates_x(x), self.f_coordinates_y(x), self.f_coordinates_z(x), self.f_chord(x), self.f_twist(x), self.f_pa(x)]), axis=1))
        self.airfoils = np.asarray([[x, interp_airfoil_position(airfoil_position, airfoils, x)] for x in x])
        self.coordinates = blade_matrix[:,0:4]
        self.chord = blade_matrix[:,[0,4]]
        self.twist = blade_matrix[:,[0,5]]
        self.pitch_axis = blade_matrix[:,[0,6]]
        
        #Generate CBMConfigs
        if wt_flag:
            cbmconfigs = iea37_converter(self, cs_pos, yml, materials)
            
        else:
            lst = [[cs.get('position'), CBMConfig(cs, materials, iea37=True)] for cs in yml.get('internal_structure_2d_fem').get('sections')]
            cbmconfigs = np.asarray(lst)
 
        #Generate CBMs
        tmp = []
        for x, cfg in cbmconfigs:
            bm = np.array([x, self.f_coordinates_x(x), self.f_coordinates_y(x), self.f_coordinates_z(x), self.f_chord(x), self.f_twist(x), self.f_pa(x)])
            af = interp_airfoil_position(airfoil_position, airfoils, x)
            tmp.append([x, CBM(cfg, materials = materials, blade_matrix = bm, airfoil=af)])
        self.sections = np.asarray(tmp)
        

        
        return None
    
    
    @property
    def blade_matrix(self):
        """ getter method for the property blade_matrix to retrive the full
        information set of the class in one reduced array"""
        return np.column_stack((self.coordinates, self.chord[:,1], self.twist[:,1], self.pitch_axis[:,1]))

    @property
    def x(self):
        """ getter method for the property grid to retrive only the 
        nondimensional grid values """
        return self.coordinates[:,0]

    def blade_gen_section(self, topo_flag=True, mesh_flag=True):
        """
        generates and meshes all sections of the blade
        """
        for (x, cs) in self.sections:
            if topo_flag:
                print('STATUS:\t Building Section at grid location %s' % (x))
                cs.cbm_gen_topo()
            if mesh_flag:
                print('STATUS:\t Meshing Section at grid location %s' % (x))
                cs.cbm_gen_mesh()
        return None
               

    def blade_run_vabs(self, loads = None, **kwargs):
        """
        runs vabs for every section
        
        Parameters
        ----------
        loads : dict, optional
            dictionary of the following keys and values, (default=None)
            for detailed information see the VABSConfig documentation or the 
            VABS user manual
            F : nparray([[grid, F1, F2, F3]]) 
            M : nparray([[grid, M1, M2, M3]]) 
            f : nparray([[grid, f1, f2, f2]])
            df : nparray([[grid, f1', f2', f3']])
            dm :  nparray([[grid, m1', m2', m3']])
            ddf : nparray([[grid, f1'', f2'', f3'']])
            ddm : nparray([[grid, m1'', m2'', m3'']])
            dddf : nparray([[grid, f1''', f2''', f3''']])
            dddm : nparray([[grid, m1''', m2''', m3''']])
            
        """
        vc = VABSConfig()
        lst = []
        for (x, cs) in self.sections:
            if loads:
                vc.recover_flag = 1
                load = interp_loads(loads, x)
                for k,v in load.items():
                    setattr(vc,k,v)
            cs.config.vabs_cfg = vc
            
            print('STATUS:\t Running VABS at grid location %s' % (x))
            cs.cbm_run_vabs(**kwargs)
            lst.append([x, cs.BeamProperties])
        self.beam_properties = np.asarray(lst)
        return None        


    def blade_plot_attributes(self):
        """
        plot the coordinates, chord, twist and pitch axis location of the blade
        """
        fig, ax = plt.subplots(3,2)
        fig.suptitle(self.name, fontsize=16)
        fig.subplots_adjust(wspace=0.25, hspace=0.25)
        
        ax[0][0].plot(self.coordinates[:,0], self.coordinates[:,1], 'k.-')
        ax[0][0].set_ylabel('x-coordinate [m]')
        
        ax[1][0].plot(self.coordinates[:,0], self.coordinates[:,2], 'k.-')
        ax[1][0].set_ylabel('y-coordinate [m]')
        
        ax[2][0].plot(self.coordinates[:,0], self.coordinates[:,3], 'k.-')
        ax[2][0].set_ylabel('z-coordinate [m]')
        
        ax[0][1].plot(self.chord[:,0], self.chord[:,1], 'k.-')
        ax[0][1].set_ylabel('chord [m]')
        
        ax[1][1].plot(self.twist[:,0], self.twist[:,1], 'k.-')
        ax[1][1].set_ylabel('twist [rad]')
        
        ax[2][1].plot(self.pitch_axis[:,0], self.pitch_axis[:,1], 'k.-')
        ax[2][1].set_ylabel('pitch axis location [1/chord]')
        
#        ax3d = fig.add_subplot(326, projection='3d')
#        for bm, af in zip(self.blade_matrix, self.airfoil):
#            tmp_shape = af.coordinates[:,0].shape
#            arr = af.coordinates*bm[4]
#            ax3d.plot(np.ones(tmp_shape)*bm[1],arr[:,0],arr[:,1])
        plt.show()
 


    
    
    def blade_exp_beam_props(self, cosy='global', style='DYMORE', eta_offset=0):
        """
        Exports the beam_properties in the 
        
        Parameters
        ----------
        cosy : str
            either 'global' for the global beam coordinate system or 
            'local' for a coordinate system that is always pointing with 
            the chord-line (in the twisted frame)
        
        style : str
            select the style you want the beam_properties to be exported
            'DYMORE' will return an array of the following form:
            [[Massterms(6) (m00, mEta2, mEta3, m33, m23, m22) 
            Stiffness(21) (k11, k12, k22, k13, k23, k33,... k16, k26, ...k66)
            Viscous Damping(1) mu, Curvilinear coordinate(1) eta]]
        
            ...
            
        eta_offset : float
            if the beam eta coordinates from start to end of the beam doesn't 
            coincide with the global coorinate system of the blade. The unit
            is in nondimensional r coordinates (x/Radius)
                
            
        Returns
        ----------
        """
        
        lst = []
        for cs in self.sections:
            #define rotation angle in rad
            if cosy=='global':
                rot = 0
            elif cosy=='local':
                rot = float(job.f_twist(cs[0])) 
                        
            #export data for each section
            if style=='DYMORE':
                lst.append(cs[1].cbm_exp_dymore_beamprops(eta=cs[0], rotate=rot))
            elif style == 'CAMRADII':
                pass
            elif style == 'CPLambda':
                pass        
                        
        arr = np.asarray(lst)
        return arr
        
    
    def blade_plot_beam_props(self):
        """
        plots the beam properties of the blade
        
        self.beam_properties()
        """
        plot_beam_properties(self.blade_exp_beam_props())
        
    
    def blade_plot_sections(self, **kwargs):
        """
        plots the different sections of the blade
        """      
        for (x,cs) in self.sections:
            string = 'Blade: '+ str(self.name) + '; Section : '+str(x)
            cs.cbm_post_2dmesh(title=string, **kwargs)
        return None    
    
    
    def blade_post_3dtopo(self, flag_wf = True, flag_lft = False, flag_topo = False, flag_mesh = False):
        """
        generates the wireframe and the loft surface of the blade

        Returns
        ----------
        loft : OCC.TopoDS_surface
            the 3D surface of the blade
        
        wireframe : list
            list of every airfoil_wire scaled and rotated at every grid point
        """
        (self.display, self.start_display, self.add_menu, self.add_function_to_menu) = display_config(cs_size = 0.3)
        
        if flag_wf:
            wireframe = []
            for bm, afl in zip(self.blade_matrix, self.airfoils[:,1]):
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
                self.display.DisplayShape(wire, color='BLACK')
            
        if flag_lft:
            loft = make_loft(wireframe, ruled=True, tolerance=1e-3, continuity=1, check_compatibility=True)
            self.display.DisplayShape(loft, transparency=0.5, update=True)
        
        if flag_topo:
            for (x,cs) in self.sections:
                coord = self.f_coordinates_x(x), self.f_coordinates_y(x), self.f_coordinates_z(x)
                display_SONATA_SegmentLst(self.display, cs.SegmentLst, coord, -np.pi/2, -np.pi/2)
                
        self.display.View_Iso()
        self.display.FitAll()
        self.start_display()   
        return None         


#%% ====== M A I N ==============
if __name__ == '__main__':
    plt.close('all')
    
    #%% ====== WindTurbine ============== 
    with open('./jobs/PBortolotti/IEAonshoreWT.yaml', 'r') as myfile:
        inputs  = myfile.read()
    with open('jobs/PBortolotti/IEAontology_schema.yaml', 'r') as myfile:
        schema  = myfile.read()
    validate(yaml.load(inputs), yaml.load(schema))    
    yml = yaml.load(inputs)
    
    airfoils = [Airfoil(af) for af in yml.get('airfoils')]
    materials = read_IEA37_materials(yml.get('materials'))
    
    job = Blade(name='IEAonshoreWT')
    job.read_IEA37(yml.get('components').get('blade'), airfoils, materials, wt_flag=True)

    job.blade_gen_section(mesh_flag = True)
    job.blade_run_vabs()
    job.blade_plot_sections()

    job.blade_post_3dtopo(flag_lft = False, flag_topo = True)


#   %% ====== Helicopter ============== 
#    #with open('jobs/VariSpeed/UH-60A_adv.yml', 'r') as f:
#    with open('jobs/VHeuschneider/blade_config.yml', 'r') as f:
#         yml = yaml.load(f.read())
#        
#    airfoils = [Airfoil(af) for af in yml.get('airfoils')]
#    materials = read_IEA37_materials(yml.get('materials'))
#    
#    job = Blade()
#    job.read_IEA37(yml.get('components').get('blade'), airfoils, materials, wt_flag=False)     
#    job.blade_gen_section()
#    
#    loads = {}
#    loads['F'] = np.array([[0.3, 88500, 0, 0],
#                          [1.0, 0, 0, 0]])
#    loads['M'] = np.array([[0.3, -12, 640, -150],
#                          [1.0, 0, 0, 0]])
#    
#    job.blade_run_vabs(loads)  
#    job.blade_plot_sections(attribute='stressM.sigma11')
#    job.blade_plot_sections(attribute='sf', vmin=0, vmax=3)
#    job.beam_properties[0,1].MpUS
#    #job.blade_post_3dtopo(flag_lft = True, flag_topo = True)
#    print((job.sections[0,1].BeamProperties.Xm2/0.130)*100, '%')
#    #job.blade_plot_sections(attribute='sf', vmin=0, vmax=3)