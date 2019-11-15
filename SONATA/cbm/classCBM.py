# -*- coding: utf-8 -*-
"""Defines the Crosssectional Beam Model (CBM) class
Created on Wed Jan 03 13:56:37 2018
@author: TPflumm

https://numpydoc.readthedocs.io/en/latest/format.html
 """

#Basic PYTHON Modules:
import pickle as pkl
import matplotlib.pyplot as plt
from datetime import datetime
import subprocess
import os 
import getpass 
import math
import numpy as np
import platform
import time
import copy


from SONATA.utl.converter import anbax_converter


#PythonOCC Modules
try:
    from OCC.Display.SimpleGui import init_display
except:
    pass
from OCC.Core.gp import gp_Ax2, gp_Pnt, gp_Dir, gp_Ax1, gp_Trsf, gp_Ax3 

if __name__ == '__main__':
    os.chdir('../..')

#SONATA modules:
from SONATA.cbm.fileIO.CADoutput import export_to_step
from SONATA.cbm.fileIO.CADinput import load_3D, import_2d_stp, import_3d_stp
from SONATA.classMaterial import read_yml_materials
from SONATA.cbm.classCBMConfig import CBMConfig

from SONATA.cbm.bladegen.blade import Blade

from SONATA.cbm.topo.segment import Segment
from SONATA.cbm.topo.web import Web
from SONATA.cbm.topo.utils import  getID
from SONATA.cbm.topo.weight import Weight
from SONATA.cbm.topo.BSplineLst_utils import get_BSplineLst_length, BSplineLst_from_dct, set_BSplineLst_to_Origin
from SONATA.cbm.topo.wire_utils import rotate_wire, translate_wire, scale_wire, discretize_wire, get_wire_length, trsf_wire

from SONATA.cbm.mesh.node import Node
from SONATA.cbm.mesh.cell import Cell
from SONATA.cbm.mesh.mesh_utils import grab_nodes_of_cells_on_BSplineLst, sort_and_reassignID, merge_nodes_if_too_close
from SONATA.cbm.mesh.consolidate_mesh import consolidate_mesh_on_web
from SONATA.cbm.mesh.mesh_intersect import map_mesh_by_intersect_curve2d
from SONATA.cbm.mesh.mesh_core import gen_core_cells

from SONATA.vabs.classVABSConfig import VABSConfig
from SONATA.vabs.vabs_utl import export_cells_for_VABS

from SONATA.vabs.classStrain import Strain
from SONATA.vabs.classStress import Stress
from SONATA.vabs.failure_criteria import von_Mises, tsaiwu_2D, maxstress_2D, maxstrain_2D, hashin_2D

from SONATA.cbm.classBeamSectionalProps import BeamSectionalProps

from SONATA.cbm.cbm_utl import trsf_sixbysix
from SONATA.cbm.display.display_mesh import plot_cells
from SONATA.cbm.display.display_utils import export_to_JPEG, export_to_PNG, export_to_PDF, \
                                        export_to_SVG, export_to_PS, export_to_EnhPS, \
                                        export_to_TEX, export_to_BMP,export_to_TIFF, \
                                        show_coordinate_system, display_SONATA_SegmentLst,\
                                        display_custome_shape, transform_wire_2to3d, display_config


try:
    from SONATA.anbax.anbax_utl import build_dolfin_mesh 
    from anba4 import anbax
    #from SONATA.anbax.anba_v4.anba4.anbax import anbax
except:
    pass

class CBM(object):
    """ 
    This Class includes the SONATA Dicipline Module for Structural 
    Composite Beam Modelling (CBM). 
    Design Variables are passed in form of the configuration Object or a 
    configuration file.          

    Attributes
    ----------
    config : configuration
        Pointer to the Configuration object.
        
    materials: list
        List of Materials object instances.
    
    SegmentLst: list
        list of Segment object instances
    
    refL : float, default: 1.0
        reference length to account for different dimensions and sizes in the 
        cross-section. This length is approximately the circumference of the 
        outer curve.  

    Methods
    -------
    cbm_save(output_filename=None)
         saves the complete cbm instance as pickle
    
    cbm_load(input_filename=None)
        loads the complete cbm instance from pickle
    
    cbm_save_topo(output_filename=None)
        saves the topology (SegmentLst, WebLst, and BW) as pickle
        
    cbm_load_topo(input_filename=None)
        loads the topology (SegmentLst, WebLst, and BW) from pickle
        
    cbm_save_mesh(output_filename=None)
        saves the mesh (self.mesh) as pickle
        
    cbm_load_mesh(input_filename=None)
        loads the mesh (self.mesh) as pickle
        
    cbm_save_res(output_filename=None)
        saves the configuration and the VABS BeamProperties results as pickle
    
    cbm_load_res(input_filename=None)
        saves the configuration and the VABS BeamProperties results from pickle
        
    cbm_stpexport_topo(export_filename=None)
        exports all Layer wires and the Segment0 Boundary Wire as .step
       
    __cbm_generate_SegmentLst(**kwargs)
        psydo private method of the cbm class to generate the list of 
        Segments in the instance. 
        
    cbm_gen_topo(**kwargs)
        generates the topology.
         
    cbm_gen_mesh(**kwargs)
        generates the dicretization of topology 
        
    cbm_review_mesh()
        prints a summary of the mesh properties to the screen 
        
    cbm_run_vabs()
        runs the solver VABS (Variational Asymptotic Beam Sectional Analysis)
    
    cbm_run_anbax()
        runs the solver anbax from macro morandini
        
    cbm_post_2dmesh(attribute='MatID', title='NOTITLE', **kw)
        displays the mesh with specified attributes with matplotlib
      
    cbm_post_3dtopo()
        displays the topology with the pythonocc 3D viewer
    
    cbm_post_3dmesh()
        displays the 2D mesh with the pythonocc 3D viewer
        
    cbm_set_DymoreMK(x_offset=0)
        Converts the Units of CBM to DYMORE/PYMORE/MARC units and returns the 
        array of the beamproperties with Massterms(6), Stiffness(21), 
        damping(1) and curvilinear coordinate(1)
    
    
    Notes
    ----------
    1 
        For computational efficiency it is sometimes not suitable to recalculate 
        the topology or the crosssection every iteration, maybe design flags 
        to account for that.
    2 
        make it possible to construct an instance by passing the topology 
        and/or mesh


    Examples
    --------
    >>> config = CBMConfig(fname)
    >>> job = CBM(config)
    >>> job.cbm_gen_topo()
    >>> job.cbm_gen_mesh()
    >>> job.cbm_review_mesh()
    >>> job.cbm_run_vabs()
    >>> job.cbm_post_2dmesh(title='Hello World!')

    """

        
    #__slots__ = ('config' , 'materials' , 'SegmentLst' , 'WebLst' , 'BW' , 'mesh', 'BeamProperties', 'display', )
    def __init__(self, Configuration, materials = None, **kwargs):
        """
        Initialize attributes.

        Parameters
        ----------
        Configuration : <Configuration>
            Pointer to the <Configuration> object.
        materials : dict(id, Materials)
        """
        
        self.config = Configuration
        if isinstance(materials, dict):
            self.materials = materials
        else:
            self.materials = read_yml_materials(self.config.setup['material_db'])
        
        self.name = 'cbm_noname'
        if kwargs.get('name'):
            self.name = kwargs.get('name')
        
        self.Ax2 = gp_Ax2()
        self.SegmentLst = []
        self.WebLst = []    
        self.BW = None
        
        self.mesh = []
        self.BeamProperties = None
        self.display = None
        
        self.refL = 1.0       
        
        self.startTime = datetime.now()
        self.exportLst = [] #list that contains all objects to be exported as step   
        self.surface3d = None #TODO: Remove definition and set it up in classBlade
        self.Blade = None #TODO: Remove definition and set it up in classBlade
        
        if self.config.setup['input_type'] == 3:
            self.surface3d = load_3D(self.config.setup['datasource'])
            
        elif self.config.setup['input_type'] == 4:
            self.blade =  Blade(self.config.setup['datasource'])
            self.surface3d = self.blade.surface
      
        elif self.config.setup['input_type'] == 5:
            #wire = kwargs.get('wire') #in the blade reference frame!
            self.Ax2 = kwargs.get('Ax2')
            self.BoundaryBSplineLst = kwargs.get('BSplineLst')
            self.Theta = 0
            self.refL = get_BSplineLst_length(self.BoundaryBSplineLst)
                        
    def __getstate__(self):
        """Return state values to be pickled."""
        return (self.config, self.materials, self.SegmentLst, self.WebLst, self.BW, self.mesh, self.BeamProperties)   
    
    def __setstate__(self, state):
        """Restore state from the unpickled state values."""
        (self.config, self.materials, self.SegmentLst, self.WebLst, self.BW, self.mesh, self.BeamProperties)  = state
        if self.config.setup['input_type'] == 3:
            self.surface3d = load_3D(self.config.setup['datasource'])
        elif self.config.setup['input_type']  == 4:
            self.blade =  Blade(self.config.setup['datasource'])
            self.surface3d = self.blade.surface
    
    
    def cbm_save(self, output_filename=None):
        ''' saves the complete CBM instance as pickle
        
        Parameters
        ----------
        output_filename : string, optional 
            path/filename to save to.
            The Default uses the config.filename and replaces .yml with .pkl
        '''
        if output_filename is None:
            output_filename = self.config.filename
            output_filename = output_filename.replace('.yml', '.pkl')
   
        with open(output_filename, 'wb') as output:
            pkl.dump(self, output, protocol=pkl.HIGHEST_PROTOCOL)        
        return None


    def cbm_load(self, input_filename=None):
        ''' loads the complete CBM instance from pickled file
        
        Parameters
        ----------
        input_filename : string, optional 
            path/filename of the .pkl file
            The Default uses the config.filename and replaces .yml with .pkl
        '''
        if input_filename is None:
            input_filename = self.config.filename
            input_filename = input_filename.replace('.yml', '.pkl') 
        
        with open(input_filename, 'rb') as handle:
            tmp_dict = pkl.load(handle, encoding='latin1').__dict__
            self.__dict__.update(tmp_dict)    
        return None


    def cbm_save_topo(self, output_filename=None):
        '''saves the topology (SegmentLst, WebLst, and BW) as pickle
        
        Parameters
        ----------
        output_filename : string, optional 
            path/filename, 
            The Default uses the config.filename and replaces .yml with .pkl
        '''
        if output_filename is None:
            output_filename = self.config.filename
            output_filename = output_filename.replace('.yml', '_topo.pkl')
            
        with open(output_filename, 'wb') as output:
            pkl.dump((self.SegmentLst, self.WebLst, self.BW), output, protocol=pkl.HIGHEST_PROTOCOL)
        return None


    def cbm_load_topo(self, input_filename=None):
        '''loads the topology (SegmentLst, WebLst, and BW) as pickle
       
        Parameters
        ----------
        input_filename : string, optional 
            path/filename of the .pkl file
            The Default uses the config.filename and replaces .yml with .pkl
        '''
        if input_filename is None:
            input_filename = self.config.filename
            input_filename = input_filename.replace('.yml', '_topo.pkl')
        
        with open(input_filename, 'rb') as handle:
            (self.SegmentLst, self.WebLst, self.BW) = pkl.load(handle)
    
        #Build wires for each layer and segment
        for seg in self.SegmentLst:
            seg.build_wire()
            for layer in seg.LayerLst:
                layer.build_wire()
        return None
          
    
    def cbm_save_mesh(self, output_filename=None):
        '''saves the mesh (self.mesh) as pickle 
        
        Parameters
        ----------
        output_filename : string, optional 
            path/filename, 
            The default uses the config.filename and replaces .yml with 
            _mesh.pkl
        '''
        if output_filename is None:
            output_filename = self.config.filename
            output_filename = output_filename.replace('.yml', '_mesh.pkl')
            print('STATUS:\t Saving Mesh to: ',output_filename)
            
        with open(output_filename, 'wb') as output:
            pkl.dump(self.mesh, output, protocol=pkl.HIGHEST_PROTOCOL)
        return None
    
    
    def cbm_load_mesh(self, input_filename=None):
        '''loads the mesh (self.mesh) as pickle
        
        Parameters
        ----------
        input_filename : string, optional 
            path/filename of the .pkl file
            The default uses the config.filename and replaces .yml with 
            _mesh.pkl
        '''

        if input_filename is None:
            input_filename = self.config.filename
            input_filename = input_filename.replace('.yml', '_mesh.pkl')
        
        with open(input_filename, 'rb') as handle:
            mesh = pkl.load(handle)   
        (self.mesh, nodes) = sort_and_reassignID(mesh)
        return None


    def cbm_save_res(self, output_filename=None):
        '''saves the configuration and the VABS BeamProperties results as 
        pickle
        
        Parameters
        ----------
        output_filename : string, optional 
            The default uses the config.filename and replaces .yml with 
            _res.pkl
        '''
        
        if output_filename is None: 
            output_filename = self.config.filename 
            output_filename = output_filename.replace('.yml', '_res.pkl') 
          
        with open(output_filename, 'wb') as output: 
            pkl.dump((self.config, self.BeamProperties), output, protocol=pkl.HIGHEST_PROTOCOL) 
     
        
     
    def cbm_load_res(self, input_filename=None): 
        '''saves the configuration and the VABS BeamProperties results from 
        pickle
        
        Parameters
        ----------
        input_filename : string, optional 
            The default uses the config.filename and replaces .yml with 
            _res.pkl
        '''
        
        if input_filename is None: 
            input_filename = self.config.filename 
            input_filename = input_filename.replace('.yml', '_res.pkl') 
         
        with open(input_filename, 'rb') as handle: 
            (self.config, self.BeamProperties) = pkl.load(handle)


    def cbm_stpexport_topo(self, export_filename=None):
        '''exports all Layer wires and the Segment0 Boundary Wire as .step
        
        Parameters
        ----------
        export_filename : string, optional 
            The default uses the config.filename and replaces .yml with 
            .stp
            
        Notes
        ----------
        If the results are imported into CATIA or a similar CAD software they 
        are not smooth. Future improvements are needed.
        '''
        
        if export_filename is None:
            export_filename = self.config.filename
            export_filename = export_filename.replace('.yml', '.stp')
        
        self.exportLst.append(self.SegmentLst[0].wire)
        for seg in self.SegmentLst:
            for layer in seg.LayerLst:
                self.exportLst.append(layer.wire)            
            
        print('STATUS:\t Exporting Topology to: ', export_filename)   
        export_to_step(self.exportLst, export_filename)
        return None


    def _cbm_generate_SegmentLst(self, **kwargs):
        '''
        psydo private method of the cbm class to generate the list of 
        Segments in the instance. 
        
        '''
        self.SegmentLst = []   #List of Segment Objects
        
        #TODO cleanup this mess!
        for k, seg in self.config.segments.items():
            if k == 0:        
                if self.config.setup['input_type'] == 0:   #0) Airfoil from UIUC Database  --- naca23012
                    self.SegmentLst.append(Segment(k, **seg, **self.config.setup, OCC=False, airfoil = self.config.setup['datasource']))
                    
                elif self.config.setup['input_type'] == 1: #1) Geometry from .dat file --- AREA_R250.dat
                    self.SegmentLst.append(Segment(k, **seg, **self.config.setup, OCC=False, filename = self.config.setup['datasource']))
                
                elif self.config.setup['input_type'] == 2: #2)2d .step or .iges  --- AREA_R230.stp
                    BSplineLst = import_2d_stp(self.config.setup['datasource'], self.config.setup['scale_factor'], self.config.setup['Theta'])
                    self.SegmentLst.append(Segment(k, **seg, **self.config.setup, OCC=True, Boundary = BSplineLst))
                
                elif self.config.setup['input_type'] == 3: #3)3D .step or .iges and radial station of crosssection --- AREA_Blade.stp, R=250
                    BSplineLst = import_3d_stp(self.config.setup['datasource'],self.config.setup['radial_station'], self.config.setup['scale_factor'], self.config.setup['Theta'])
                    self.SegmentLst.append(Segment(k, **seg, **self.config.setup, OCC=True, Boundary = BSplineLst))  
        
                elif self.config.setup['input_type'] == 4: #4)generate 3D-Shape from twist,taper,1/4-line and airfoils, --- examples/UH-60A, R=4089, theta is given from twist distribution
                    BSplineLst = self.blade.get_crosssection(self.config.setup['radial_station'], self.config.setup['scale_factor'])
                    self.SegmentLst.append(Segment(k, **seg, Theta = self.blade.get_Theta(self.config.setup['radial_station']), OCC=True, Boundary = BSplineLst))  
                
                elif self.config.setup['input_type'] == 5: #5) IEA37 Formulation, everything is passed internally!
                    self.SegmentLst.append(Segment(k, **seg, Theta = self.Theta, OCC=True, Boundary = self.BoundaryBSplineLst))
                        
                    #BSplineLst get BSplineLst from IEA37 definition of blade. By generating the crosssection in the blade class and passing the BSplineLst to the section!
                    #Get Theta from the IEA37 definition of the blade !
                    #self.SegmentLst.append()
                    
                else:
                    print('ERROR:\t WRONG input_type')
         
            else:
                if self.config.setup['input_type'] == 4:
                    self.SegmentLst.append(Segment(k, **seg, Theta = self.blade.get_Theta(self.config.setup['radial_station'])))
                
                elif self.config.setup['input_type'] == 5:
                    self.SegmentLst.append(Segment(k, **seg, Theta = self.Theta))
                
                else:
                    self.SegmentLst.append(Segment(k, **seg, **self.config.setup))
        
        sorted(self.SegmentLst, key=getID)
        self.refL = get_BSplineLst_length(self.SegmentLst[0].BSplineLst)
        return None


    def cbm_gen_topo(self, **kwargs):
        '''
        CBM Method that generates the topology. It starts by generating the 
        list of Segments. It continous to gen all layers for Segment 0. 
        Subsequently the webs are defined and afterwards the layers of the 
        remaining Segments are build. The Balance Weight is defined at the end
        of this method.        
        '''               
        #Generate SegmentLst from config:
        self.SegmentLst = []
        self._cbm_generate_SegmentLst(**kwargs)
        #Build Segment 0:
        self.SegmentLst[0].build_wire()
        self.SegmentLst[0].build_layers(l0 = self.refL, **kwargs)
        self.SegmentLst[0].determine_final_boundary()
        
        #Build Webs:    
        self.WebLst = []
        if len(self.config.webs) > 0:
            for k, w in self.config.webs.items(): 
                print('STATUS:\t Building Web %s' %(k+1))
                self.WebLst.append(Web(k, w['Pos1'], w['Pos2'], self.SegmentLst))
            sorted(self.SegmentLst, key=getID)  
            
        #Build remaining Segments:
        if len(self.config.webs) > 0:
            for i,seg in enumerate(self.SegmentLst[1:], start=1):
                seg.Segment0 = self.SegmentLst[0]
                seg.WebLst = self.WebLst
                seg.build_segment_boundary_from_WebLst(self.WebLst,self.SegmentLst[0])            
                seg.build_layers(self.WebLst,self.SegmentLst[0], l0 = self.refL)
                seg.determine_final_boundary(self.WebLst,self.SegmentLst[0])
                seg.build_wire()
              
        self.BW = None
        #Balance Weight:
        if self.config.setup['BalanceWeight'] == True:
            #print('STATUS:\t Building Balance Weight')   
            #self.BW = Weight(0, self.config.bw['XPos'], self.config.bw['YPos'], self.config.bw['Diameter'], self.config.bw['Material'])
            p = self.SegmentLst[0].det_weight_Pnt2d(self.config.bw['s'], self.config.bw['t'])
            self.BW = Weight(0,p,self.config.bw['Diameter'], self.config.bw['Material'])
            #print(p.Coord())
        return None


    def cbm_gen_mesh(self, **kwargs):
        '''
        CBM Method that generates the dicretization of topology and stores the 
        cells and nodes in both the <Layer> instances, the <Segment> instances 
        and the attribute self.mesh that is a list of <Cell> instances


        Parameters:
        ----------
        split_quads : bool, optional
            This option can be passed as keyword argument and splits the quads 
            (4 node cells) in mesh into two cells of 3 nodes      
        
        
        Notes:
        ----------  
        More option keyword arguments shall be possible in the future      
        
        Examples:
        ----------  
        >>> job.cbm_gen_mesh(splitquads=True)
        
        '''
        
        split_quads  = False
        if 'split_quads' in kwargs:
            if type(kwargs['split_quads']) == bool:
                split_quads = kwargs['split_quads']
            else:
                print ('split_quads must provide a boolean value')
        
                
        self.mesh = []
        Node.class_counter = 1
        Cell.class_counter = 1
        #meshing parameters:  
        Resolution = self.config.setup['mesh_resolution'] # Nb of Points on Segment0
        global_minLen = round(self.refL/Resolution,5)
            
        core_cell_area = 1.0*global_minLen**2
        bw_cell_area = 0.7*global_minLen**2
        web_consolidate_tol = 0.5*global_minLen

        #===================MESH SEGMENT
        for j,seg in enumerate(reversed(self.SegmentLst)):
            self.mesh.extend(seg.mesh_layers(self.SegmentLst, global_minLen, self.WebLst, display=self.display, l0=self.refL))
            #mesh,nodes = sort_and_reassignID(mesh)
            
        #===================MESH CORE  
        if self.config.flags['mesh_core']:
            for j,seg in enumerate(reversed(self.SegmentLst)):
                #if seg.ID == 1:
                    #core_cell_area = 1.6*global_minLen**2
                    #print(core_cell_area)
                self.mesh.extend(seg.mesh_core(self.SegmentLst, self.WebLst, core_cell_area, display=self.display))
                
        #===================consolidate mesh on web interface         
        for web in self.WebLst:
            #print web.ID,  'Left:', SegmentLst[web.ID].ID, 'Right:', SegmentLst[web.ID+1].ID,
            print('STATUS:\t Consolidate Mesh on Web Interface', web.ID)  
            (web.wl_nodes, web.wl_cells) = grab_nodes_of_cells_on_BSplineLst(self.SegmentLst[web.ID].cells, web.BSplineLst)            
            (web.wr_nodes, web.wr_cells) = grab_nodes_of_cells_on_BSplineLst(self.SegmentLst[web.ID+1].cells, web.BSplineLst)
                      
            newcells = consolidate_mesh_on_web(web, web_consolidate_tol, self.display)
            self.mesh.extend(newcells)
        
        #=====================split quad cells into trias:
        if split_quads == True:
            print('STATUS:\t Splitting Quads into Trias') 
            tmp = []
            for c in self.mesh:
                tmp.extend(c.split_quads())   
            self.mesh = tmp
            
            
        #============= BALANCE WEIGHT - CUTTING HOLE ALGORITHM
        if self.config.setup['BalanceWeight'] == True:
            print('STATUS:\t Meshing Balance Weight')   
            
            self.mesh, boundary_nodes = map_mesh_by_intersect_curve2d(self.mesh, self.BW.Curve, self.BW.wire, global_minLen)
            #boundary_nodes = merge_nodes_if_too_close(boundary_nodes,self.BW.Curve,global_minLen,tol=0.05))
            [bw_cells,bw_nodes] = gen_core_cells(boundary_nodes, bw_cell_area)
            
            for c in bw_cells:
                c.structured = False
                c.theta_3 = 0
                c.MatID = self.config.bw['Material']
                c.calc_theta_1()
            
            self.mesh.extend(bw_cells)
       
        #invert nodes list of all cell to make sure they are counterclockwise for vabs in the right coordinate system!
        for c in self.mesh:
            if c.orientation == False:
                c.invert_nodes()     
        (self.mesh, nodes) = sort_and_reassignID(self.mesh)
        return None
    

    def cbm_review_mesh(self):
        '''
        CBM method that prints a summary of the mesh properties to the 
        screen 
        '''
        
        print('STATUS:\t Review Mesh:')
        print('\t   - Total Number of Cells: %s' %(len(self.mesh)))
        print('\t   - Duration: %s' % (datetime.now() - self.startTime))
        #print '\t   - Saved as: %s' % filename 
        minarea = min([c.area for c in self.mesh])
        print('\t   - smallest cell area: %s' % minarea) 
        minimum_angle = min([c.minimum_angle for c in self.mesh])
        print('\t   - smallest angle [deg]: %s' % minimum_angle) 
        #orientation = all([c.orientation for c in mesh])
        #print '\t   - Orientation [CC]: %s' % orientation 
        
        return None

   
    def cbm_run_vabs(self, jobid=None, rm_vabfiles=True, ramdisk=False):
        '''CBM method to run the solver VABS (Variational Asymptotic Beam 
        Sectional Analysis). Note that this method is designed to work if 
        VABSIII is set in the PATH variable. For Users at the TUM-HT please load 
        the vabs module beforehand.
                
        Parameters
        ----------
        jobid : string, optional
                assign a unique ID for the job. If no jobid is assigned the 
                isoformat of datetime with microseconds is used
        rm_vabfiles : bool, optional
                removes VABS files after the calculation is completed and 
                the results are stored.
        ramdisk : bool, optional, 
            Instead of storing the writing and reading the vabs job directory, 
            the ramdisk "/tmpfs/username" is used. This options is currently 
            designed for linux users make sure to mount it beforehand with to 
            assign 200MB of Memory to the virtual drive.
            >>> sudo mount -t tmpfs -o size=200M none /tmpfs/username
            
        Returns
        ----------
        None : everything is stored within the CBM instance
        
        Examples
        ----------
        >>> job.cbm_run_vabs(rm_vabfiles=True, ramdisk=True)

        '''
        (self.mesh, nodes) = sort_and_reassignID(self.mesh)
        
        if jobid == None:
            s = datetime.now().isoformat(sep='_',timespec='microseconds')
            jobid =  s.replace(':','').replace('.','')
        fstring = '_'+jobid+'.vab'
        
        if ramdisk == True:
            rm_vabfiles = True
            if os.path.exists('/tmpfs'):
                user = getpass.getuser()
                path = '/tmpfs/'+user
                if not os.path.exists(path):
                    os.makedirs(path)
                #tmp = path+'/'+self.config.filename.split('/')[-1]
                vabs_filename = path+'/'+self.name+fstring
                    
            else:
                print('ERROR: ramdisk directory "/tmpfs" does not exist!')

        elif self.config.filename == '':
            vabs_filename = self.name+fstring
            
        else:
            print('config_filename')
            vabs_filename = self.config.filename.replace('.yml', fstring)
        
        print('STATUS:\t Running VABS for constitutive modeling:')
        
        if platform.system() == 'Linux':
            #executable = 'SONATA/vabs/bin/VABSIII'
            #check if module vabs is loaded, if not load it!
            #vabs_cmd = 'VABSIII '+vabs_filename
            #cmd = ['/bin/bash', '-i', '-c', vabs_cmd]
            cmd = ['VABSIII', vabs_filename]
            
        elif platform.system == 'Windows':
            cmd = ['VABS/bin/VABSIII.exe', vabs_filename]

        # Mac
        elif  platform.system() == 'Darwin':
            cmd = ['wine', '/Users/rfeil/work/8_VABS/vabs_WIN/AnalySwift/VABS/VABSIII.exe', vabs_filename]
            
        #print('vabs_fname',vabs_filename)
        result = None
        counter = 0
        stdout = ''
        while result is None and counter<1000:
            #EXECUTE VABS:
            try:
                if self.config.vabs_cfg.recover_flag == 1:
                    self.config.vabs_cfg.recover_flag=0
                    export_cells_for_VABS(self.mesh, nodes, vabs_filename, self.config.vabs_cfg, self.materials)
                    stdout = subprocess.run(cmd, stdout=subprocess.PIPE).stdout.decode('utf-8')
                    self.config.vabs_cfg.recover_flag=1
                    export_cells_for_VABS(self.mesh, nodes, vabs_filename, self.config.vabs_cfg, self.materials)
                    print('STATUS:\t Running VABS for 3D Recovery:')
                    stdout = subprocess.run(cmd, stdout=subprocess.PIPE).stdout.decode('utf-8')
                    
                else:
                    export_cells_for_VABS(self.mesh, nodes, vabs_filename, self.config.vabs_cfg, self.materials)
                    stdout = subprocess.run(cmd, stdout=subprocess.PIPE).stdout.decode('utf-8')
                    subprocess.call(cmd, shell=True)

                stdout = stdout.replace('\r\n\r\n','\n\t   -')
                stdout = stdout.replace('\r\n','\n\t   -')
                stdout = stdout.replace('\n\n','\n\t   -')
                stdout = stdout[:-2]
                                
                if ' VABS finished successfully' in stdout:
                    stdout = 'STATUS:\t VABS Calculations Completed: \n\t   -' + stdout
                else:
                    stdout = 'ERROR:\t VABS Calculations Incomplete: \n\t   -' + stdout
               
                #print('STATUS:\t Total Elapsed Time: %s' % (datetime.now() - self.startTime))
                print(stdout)
                #VABS Postprocessing:
                result = BeamSectionalProps(vabs_filename+'.K')
            except Exception as e:
                    if 'All "vabsiii" licenses in us' in stdout or 'Something wrong with the license file' in stdout:
                        time.sleep(2)
                        counter += 1
                    else:
                        print(e)
                        break
                    
        self.BeamProperties = result
        
        if self.config.vabs_cfg.recover_flag == 1:
            self.BeamProperties.read_all_VABS_Results()
            #ASSIGN Stress and strains to elements:
            for i,c in enumerate(self.mesh):
                c.strain = Strain(self.BeamProperties.ELE[i][1:7])
                c.stress = Stress(self.BeamProperties.ELE[i][7:13])
                c.strainM = Strain(self.BeamProperties.ELE[i][13:19])
                c.stressM = Stress(self.BeamProperties.ELE[i][19:25])
            
            #ASSIGN Displacement U to nodes:
            for i,n in enumerate(nodes):
                n.displacement = self.BeamProperties.U[i][3:6]
            
            #Calculate standart failure criterias
            self.cbm_calc_failurecriteria()
        
        #print(vabs_filename)
        #REMOVE VABS FILES:
        if rm_vabfiles:
            folder = '/'.join(vabs_filename.split('/')[:-1])
            fstring = vabs_filename.split('/')[-1]
            if not folder:
                folder = None
                
            for file in os.listdir(folder):
                if fstring in file:
                    #print('removing: '+folder+'/'+file)
                    if folder:
                        os.remove(folder+'/'+file)
                    else:
                        os.remove(file)
        
        return None
            
    
    def cbm_run_anbax(self):
        """interface method to run the solver anbax from marco.morandini 
        
        Notes
        ----------
        To be defined.

        """



        self.mesh, nodes = sort_and_reassignID(self.mesh)

        # nodes = anbax_converter(nodes_SONATA=nodes)
        # x_coord = np.zeros(len(nodes))
        # y_coord = np.zeros(len(nodes))
        # for i in range(len(nodes)):
        #     x_coord[i] = nodes[i].coordinates[0]
        #     y_coord[i] = nodes[i].coordinates[1]
        #
        # plt.plot(x_coord, y_coord)

        try:
            (mesh, matLibrary, materials, plane_orientations, fiber_orientations, maxE) = \
                build_dolfin_mesh(self.mesh, nodes, self.materials)
        except:
            print('\n')
            print('==========================================\n\n')
            print('Error, Anba4 wrapper called, but ')
            print('Anba4 _or_ Dolfin are not installed\n\n')
            print('==========================================\n\n')
        #TBD: pass it to anbax and run it!
        anba = anbax(mesh, 1, matLibrary, materials, plane_orientations, fiber_orientations, maxE)
   
        tmp_TS = anba.compute().getValues(range(6),range(6))
        tmp_MM = anba.inertia().getValues(range(6),range(6))
        
        #Define transformation T (from ANBA to SONATA/VABS coordinates)
        # B = np.array([[0,0,1],[1,0,0],[0,1,0]])
        B = np.array([[0,0,1],[-1,0,0],[0,1,0]])  # new
        T = np.dot(np.identity(3),np.linalg.inv(B))
        
        self.AnbaBeamProperties = BeamSectionalProps()
        self.AnbaBeamProperties.TS = trsf_sixbysix(tmp_TS,T)
        self.AnbaBeamProperties.MM = trsf_sixbysix(tmp_MM,T)
        return self.AnbaBeamProperties


    def cbm_calc_failurecriteria(self, criteria='tsaiwu_2D', iso_criteria = 'nocriteria'):
        """
        Applies the selected failure criteria for all cells. It is necessary 
        that the strains and stresses are calculated for all cells and the 
        strength characteristics are defined for the used materials
    
        
        Notes
        ----------
        'von_Mises' can only be applied for isotropic materials and the others
        can only be applied for orthotropic materials
        
        
        Parameters
        ----------
        criteria : string
            current options are: 'tsaiwu_2D', 'maxstress_2D', 'maxstrain_2D',
            'hashin_2D', 'von_Mises'
                
        """
        
        for c in self.mesh:
            sf = 99
            mode = 'nocriteria'
            mat = self.materials[c.MatID]
            
            if mat.orth == 0:
                if iso_criteria == 'von_Mises':
                    (sf, mode) = von_Mises(mat, c.stressM, c.strainM)
        
            if mat.orth == 1:
                if criteria == 'tsaiwu_2D':
                    (sf, mode) = tsaiwu_2D(mat, c.stressM, c.strainM)
                elif criteria == 'maxstress_2D':
                    (sf, mode) = maxstress_2D(mat, c.stressM, c.strainM)
                elif criteria == 'maxstrain_2D':
                    (sf, mode) = maxstrain_2D(mat, c.stressM, c.strainM)
                elif criteria == 'hashin_2D':
                    (sf, mode) = hashin_2D(mat, c.stressM, c.strainM)

            c.sf = sf
            c.failure_mode = mode
        
        return None 

    def cbm_post_2dmesh(self, attribute='MatID', title='NOTITLE', **kw):
        '''
        CBM Postprocessing method that displays the mesh with matplotlib.
        
        Parameters
        ----------
        attribute : string, optional
            Uses the string to look for the cell attributes. 
            The default attribute is MatID. Possible other attributes can be 
            fiber orientation (theta_3) or strains and stresses. 
            If BeamProperties are already calculated by VABS or something
            similar, elastic-axis, center-of-gravity... are displayed.
        title : string, optional
            Title to be placed over the plot.
        **kw : keyword arguments, optional
            are passed to the lower "plot_cells" function. Such options are: 
            VABSProperties=None, title='None', plotTheta11=False, 
            plotDisplacement=False, savepath
            
        Returns
        ----------
        (fig, ax) : tuple
            figure and axis handler are returned by the method
        
        Examples
        ----------
        >>> job.cbm_post_2dmesh(title='Hello World!', attribute='theta_3', plotTheta11=True)

        '''
        mesh,nodes = sort_and_reassignID(self.mesh)
        fig,ax = plot_cells(self.mesh, nodes, attribute, self.BeamProperties, title, **kw)
        return fig,ax
    
    def cbm_post_3dtopo(self):
        '''
        CBM Postprocessing method that displays the topology with the pythonocc
        3D viewer. If the input_type is 3 (3d *.stp + radial station) or 4 
        (generic blade definiton) the 3D Surface is displayed, 
        else only the 2D topology is shown.
                
        Notes
        ----------
        Be careful to set the DeviationAngle/Coeficient and scale it to the
        prolem or else it might crash. Remeber that is is in absolute 
        values (mm).
        
        '''
        (self.display, self.start_display, self.add_menu, self.add_function_to_menu) = display_config(DeviationAngle = 1e-3, DeviationCoefficient=1e-3, cs_size = self.refL/5)

        #display_custome_shape(self.display,self.SegmentLst[0].wire,2,0,[0,0,0])

        if self.config.setup['input_type'] == 3 or self.config.setup['input_type'] == 4:
            #self.display.Context.SetDeviationAngle(1e-6)       
            #self.display.Context.SetDeviationCoefficient(1e-6) 
            
            display_SONATA_SegmentLst(self.display,self.SegmentLst,(self.config.setup['radial_station'],0,0),-math.pi/2,-math.pi/2)
            self.display.DisplayShape(self.surface3d, color=None, transparency=0.7, update=True)
            
            if self.config.setup['BalanceWeight']:
                transform_wire_2to3d(self.display,self.BW.wire,(self.config.setup['radial_station'],0,0),-math.pi/2,-math.pi/2)

        else:
            display_SONATA_SegmentLst(self.display,self.SegmentLst)
            if self.config.setup['BalanceWeight']:
                self.display.DisplayShape(self.BW.Curve, color="BLACK")
        
        self.display.View_Iso()
        self.display.FitAll()
        self.start_display()   
        
        return None
    
    
    def cbm_post_3dmesh(self):
        '''
        CBM Postprocessing method that displays the 2D mesh with the pythonocc
        3D viewer. Similar functionality can be used when debuggin in the 
        meshing routines. 
        '''
                
        (self.display, self.start_display, self.add_menu, self.add_function_to_menu) = display_config(DeviationAngle = 1e-6, DeviationCoefficient=1e-6, cs_size = self.refL/5)
        for c in self.mesh:
            self.display.DisplayShape(c.wire, color='BLACK', transparency=0.7)    
        self.display.View_Top()
        self.display.FitAll()
        self.start_display()   
        return None
    
        
    def cbm_exp_dymore_beamprops(self, eta, Theta=0, solver='vabs', units={'mass':'kg', 'length':'m', 'force': 'N'}):
        '''
        Converts the Units of CBM to DYMORE/PYMORE/MARC units and returns the 
        array of the beamproperties with Massterms(6), Stiffness(21), 
        damping(1) and curvilinear coordinate(1)

        Parameters
        ----------
        
        eta : float, 
            is the beam curvilinear coordinate of the beam from 0 to 1. 
        
        Theta: float
            is the angle of rotation of the coordinate system in "radians"

        Returns
        ----------
        arr : ndarray
            [Massterms(6) (m00, mEta2, mEta3, m33, m23, m22) 
            Stiffness(21) (k11, k12, k22, k13, k23, k33,... k16, k26, ...k66)
            Viscous Damping(1) mu, Curvilinear coordinate(1) eta]
            
            
        Notes
        ----------
        - Unit Convertion takes sooo much time. Commented out for now!
        
        '''
        if solver == 'vabs':
            if Theta != 0:
                tmp_bp = self.BeamProperties.rotate(Theta)        
            else:
                tmp_bp = self.BeamProperties

        elif solver == 'anbax':
            if Theta != 0:
                tmp_bp = self.AnbaBeamProperties.rotate(Theta)        
            else:
                tmp_bp = self.AnbaBeamProperties

        MM = tmp_bp.MM
        MASS = np.array([MM[0,0], MM[2,3], MM[0,4], MM[5,5], MM[4,5], MM[4,4]])
        STIFF = tmp_bp.TS[np.tril_indices(6)[1],np.tril_indices(6)[0]]
        mu = 0.0
        return np.hstack((MASS,STIFF,mu,eta))   
    
    def cbm_exp_BeamDyn_beamprops(self, Theta=0, solver='vabs'):
        """ 
        Converts the Beam Properties of CBM to the correct coordinate System of
        BeamDyn and returns the 6x6 Stiffness matrix, the 6x6 MassMatrix.
        
        The geometry of the blade is defined by key-point coordinates and initial
        twist angles (in units of degree) in the blade local coordinate system
        (IEC standard blade system where Zr is along blade axis from root to
        tip, Xr directs normally toward the suction side, and Yr directs 
        normally toward the trailing edge).
        https://openfast.readthedocs.io/en/master/source/user/beamdyn/input_files.html
        
        Parameters
        ----------
        Theta: float, optional
            is the angle of rotation of the coordinate system in "radians"
        solver: str, optional
        
        Returns
        ----------
            tuple of arrays
            (6x6 StiffnessMatrix, 6x6MassMatrix)
            
            
        Notes:
        ----------
        - Following the station location parameter η, there are two 
        6×6 matrices providing the structural and inertial properties for this
        cross-section. First is the stiffness matrix and then the mass matrix. 
        We note that these matrices are defined in a local coordinate system 
        along the blade axis with Zl directing toward the unit tangent vector 
        of the blade reference axis.
        - Does this create an oblique cross-section!?
        
        
        """
        if solver == 'vabs':
            if Theta != 0:
                tmp_bp = self.BeamProperties.rotate(Theta)        
            else:
                tmp_bp = self.BeamProperties

        elif solver == 'anbax':
            if Theta != 0:
                tmp_bp = self.AnbaBeamProperties.rotate(Theta)        
            else:
                tmp_bp = self.AnbaBeamProperties
        
        tmp_bp = copy.deepcopy(tmp_bp)
        
        #transform to BeamDYN Coordinates
        B = np.array([[0,0,1],[0,-1,0],[1,0,0]])
        T = np.dot(np.identity(3),np.linalg.inv(B))
        
        tmp_TS = trsf_sixbysix(tmp_bp.TS,T)
        tmp_MM = trsf_sixbysix(tmp_bp.MM,T)
        return (tmp_TS, tmp_MM) 
        
        
#%%############################################################################
#                           M    A    I    N                                  #
###############################################################################    
if __name__ == '__main__':   
    plt.close('all')
    fname = 'jobs/debug/issue20/sec_config.yml'
    #fname = 'jobs/VariSpeed/uh60a_cbm_simple/sec_config.yml'
    fname = 'jobs/AREA/R250/sec_config.yml'
    #fname = 'jobs/PBortolotti/sec_config.yml'
    config = CBMConfig(fname)
    
    job = CBM(config)
    
    job.cbm_gen_topo()
    job.cbm_gen_mesh(split_quads=True)
    
    job.cbm_review_mesh()
    job.cbm_run_vabs(rm_vabfiles=False)
    #AnbaBeamProperties = job.cbm_run_anbax()
    
    
    #job.cbm_post_2dmesh(title='Hello World!')
    job.cbm_post_3dtopo()
#    job.config.vabs_cfg.recover_flag = 1
#    job.config.vabs_cfg.M = [0,2000e4,0]
