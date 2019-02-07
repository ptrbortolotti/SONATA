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

#PythonOCC Modules
from OCC.Display.SimpleGui import init_display
from OCC.gp import gp_Ax2, gp_Pnt, gp_Dir, gp_Ax1

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
from SONATA.cbm.topo.wire_utils import rotate_wire, translate_wire, scale_wire, discretize_wire

from SONATA.cbm.mesh.node import Node
from SONATA.cbm.mesh.cell import Cell
from SONATA.cbm.mesh.mesh_utils import grab_nodes_of_cells_on_BSplineLst, sort_and_reassignID, merge_nodes_if_too_close
from SONATA.cbm.mesh.consolidate_mesh import consolidate_mesh_on_web
from SONATA.cbm.mesh.mesh_intersect import map_mesh_by_intersect_curve2d
from SONATA.cbm.mesh.mesh_core import gen_core_cells

from SONATA.vabs.VABS_interface import VABS_config, export_cells_for_VABS, XSectionalProperties
from SONATA.vabs.strain import Strain
from SONATA.vabs.stress import Stress



from SONATA.cbm.display.display_mesh import plot_cells
from SONATA.cbm.display.display_utils import export_to_JPEG, export_to_PNG, export_to_PDF, \
                                        export_to_SVG, export_to_PS, export_to_EnhPS, \
                                        export_to_TEX, export_to_BMP,export_to_TIFF, \
                                        show_coordinate_system, display_SONATA_SegmentLst,\
                                        display_custome_shape, transform_wire_2to3d  

#from SONATA.anbax.anbax_utl import build_dolfin_mesh
#from SONATA.anbax.anba_v4 import anbax

class CBM(object):
    ''' 
    This Class includes the SONATA Dicipline Module for Structural 
    Composite Beam Modelling (CBM). 
    Design Variables are passed in form of the configuration Object or a 
    configuration file.          

    Attributes
    ----------
    config : configuration
        Pointer to the <Configuration> object.
        
    materials: list
        List of <Materials> object instances.
    
    SegmentLst: list
        list of <Segment> object instances
        
   
    Methods
    -------
    cbm_save : saves the complete cbm instance as pickle
    
    cbm_load : loads the complete cbm instance as pickle
    
    cbm_save_topo :
        
    cbm_load_topo : 
        
    cbm_save_mesh : 
        
    cbm_load_mesh :
        
    cbm_save_res :
    cbm_load_res :
    cbm_stpexport_topo :
    cbm_gen_topo :
    cbm_review_mesh :
    cbm_run_vabs :
    cbm_post_2dmesh :
    cbm_display_config :
    cbm_post_3dtopo : 
    cbm_post_3dmesh :
    cbm_set_DymoreMK :       
        
    
    Notes
    ----------
    1 : For computational efficiency it is sometimes not suitable to recalculate 
    the topology or the crosssection every iteration, maybe design flags to account for that.
    2 : make it possible to construct an instance by passing the 
     topology and/or mesh

    '''

        
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
        
        self.SegmentLst = []
        self.WebLst = []    
        self.BW = None
                
        self.mesh = []
        self.BeamProperties = None
        self.display = None
              
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
            bm = kwargs.get('blade_matrix')
            af = kwargs.get('airfoil')
            self.Theta = np.degrees(bm[5])
            wire = af.gen_OCCtopo()
            
            wire = rotate_wire(wire, gp_Ax1(gp_Pnt(bm[6],0,0), gp_Dir(0,0,1)), -bm[5], copy=True)
            wire = translate_wire(wire, gp_Pnt(bm[6],0,0), gp_Pnt(0,0,0))
            wire = scale_wire(wire, gp_Pnt(0,0,0), bm[4])
            #wire = translate_wire(wire, gp_Pnt(0,0,0), gp_Pnt(bm[1],bm[2],bm[3]))
            
            npArray = discretize_wire(wire, Deflection = 1e-4)
            BSplineLst = BSplineLst_from_dct(npArray, angular_deflection = 20, tol_interp=1e-4)
            self.BoundaryBSplineLst = set_BSplineLst_to_Origin(BSplineLst,self.Theta) 
            
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


    def __cbm_generate_SegmentLst(self, **kwargs):
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
                
                else:
                    self.SegmentLst.append(Segment(k, **seg, **self.config.setup))
        
        sorted(self.SegmentLst, key=getID)  
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
        self.__cbm_generate_SegmentLst(**kwargs)
        #Build Segment 0:
        self.SegmentLst[0].build_wire()
        l0 = get_BSplineLst_length(self.SegmentLst[0].BSplineLst)
        #l0 = 300
        self.SegmentLst[0].build_layers(l0 = l0)
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
                seg.build_layers(self.WebLst,self.SegmentLst[0], l0 = l0)
                seg.determine_final_boundary(self.WebLst,self.SegmentLst[0])
                seg.build_wire()
              
        self.BW = None
        #Balance Weight:
        if self.config.setup['BalanceWeight'] == True:
            print('STATUS:\t Building Balance Weight')   
            self.BW = Weight(0, self.config.bw['XPos'], self.config.bw['YPos'], self.config.bw['Diameter'], self.config.bw['Material'])
            
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
        length = get_BSplineLst_length(self.SegmentLst[0].BSplineLst)
        #length = 300
        global_minLen = round(length/Resolution,5)
            
        core_cell_area = 1.25*global_minLen**2
        bw_cell_area = 0.25*global_minLen**2
        web_consolidate_tol = 0.5*global_minLen

        #===================MESH SEGMENT
        for j,seg in enumerate(reversed(self.SegmentLst)):
            self.mesh.extend(seg.mesh_layers(self.SegmentLst, global_minLen, self.WebLst, display=self.display))
            #mesh,nodes = sort_and_reassignID(mesh)
            
        #===================MESH CORE  
        if self.config.flags['mesh_core']:
            for j,seg in enumerate(reversed(self.SegmentLst)):
                if seg.ID == 1:
                    core_cell_area = 0.25*global_minLen**2
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


    def cbm_run_anbax(self):
        """interface method to run the solver anbax from marco.morandini 
        
        Notes:
            To be defined.
        
        """
        self.mesh, nodes = sort_and_reassignID(self.mesh)
        (mesh, matLibrary, materials, plane_orientations, fiber_orientations, maxE) = \
            build_dolfin_mesh(self.mesh, nodes, self.materials)
        #TBD: pass it to anbax and run it!
        anba = anbax(mesh, 1, matLibrary, materials, plane_orientations, fiber_orientations, maxE)
        stiff = anba.compute()
        stiff.view()
        
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
                tmp = path+'/'+self.config.filename.split('/')[-1]
                vabs_filename = tmp.replace('.yml', fstring)
                    
            else:
                print('ERROR: ramdisk directory "/tmpfs" does not exist!')

        elif self.config.filename == '':
            vabs_filename = self.name+fstring
            
        else:
            print('config_filename')
            vabs_filename = self.config.filename.replace('.yml', fstring)
        
        print('STATUS:\t Running VABS for Constitutive modeling:')
        
        if platform.system() == 'Linux':
            #executable = 'SONATA/vabs/bin/VABSIII'
            #check if module vabs is loaded, if not load it!
            #vabs_cmd = 'VABSIII '+vabs_filename
            #cmd = ['/bin/bash', '-i', '-c', vabs_cmd]
            cmd = ['VABSIII', vabs_filename]
            
        elif platform.system == 'Windows':
            cmd = ['SONATA/vabs/bin/VABSIII.exe', vabs_filename]
            
        print('vabs_fname',vabs_filename)
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
                result = XSectionalProperties(vabs_filename+'.K')
            
            except Exception as e:
                    if 'All "vabsiii" licenses in us' in stdout:
                        time.sleep(0.01)
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
    
        if rm_vabfiles:
            folder = '/'.join(vabs_filename.split('/')[:-1])
            fstring = vabs_filename.split('/')[-1]
            #print(vabs_filename)
            if folder != '':
                for file in os.listdir(folder):
                    if fstring in file:
                        #print('removing: '+folder+'/'+file)
                        os.remove(folder+'/'+file)
            else:
                print(folder)
        
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
            If BeamProperties are allready calculated by VABS or something 
            similar, elastic-axis, center-of-gravity... are displayed.
        title : string, optional
            Title to be placed over the plot.
        **kw : keyword arguments, optional
            are passed to the lower "plot_cells" function. Such options are: 
            VABSProperties=None, title='None', plotTheta11=False, 
            plotDisplacement=False
            
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
                
    
    def cbm_display_config(self, DeviationAngle = 1e-5, DeviationCoefficient = 1e-5, bg_c = ((20,6,111),(200,200,200)), cs_size = 25):
        '''
        CBM method that initializes and configures the pythonOcc 3D Viewer 
        and adds Menues to the toolbar. 
        
        Parameters
        ----------
        DeviationAngle : float, optional 
            default = 1e-5
        DeviationCoefficient : float, optional
            default = 1e-5 
        bg_c : tuple, optional
            Background Gradient Color ((RBG Tuple),(RBG Tuple)) the default 
            values are a CATIA style gradient for better 3D visualization. 
            for a white background use: ((255,255,255,255,255,255))
        cs_size : float 
            coordinate system size in [mm]
            
        Returns
        ----------
        tuple :
            (self.display: the display handler for the pythonOcc 3D Viewer
            self.start_display: function handle
            self.add_menu: function handle
            self.add_function_to_menu: function handle)
        
        
        See Also
        ----------
        OCC.Display.SimpleGui : PyhtonOcc wrapper provides more details on this 
            method
        '''
        
        def export_png(): return export_to_PNG(self.display)
        def export_jpg():return export_to_JPEG(self.display)
        def export_pdf(): return export_to_PDF(self.display)
        def export_svg(): return export_to_SVG(self.display)
        def export_ps(): return export_to_PS(self.display)
        
        #===========DISPLAY CONFIG:===============
        self.display, self.start_display, self.add_menu, self.add_function_to_menu = init_display()
        self.display.Context.SetDeviationAngle(DeviationAngle) # 0.001 default. Be careful to scale it to the problem.
        self.display.Context.SetDeviationCoefficient(DeviationCoefficient) # 0.001 default. Be careful to scale it to the problem. 
        self.display.set_bg_gradient_color(bg_c[0][0],bg_c[0][1],bg_c[0][2],bg_c[1][0],bg_c[1][1],bg_c[1][2])
        show_coordinate_system(self.display,cs_size)
        
        self.add_menu('View')
        self.add_function_to_menu('View', self.display.FitAll)
        self.add_function_to_menu('View', self.display.View_Bottom)
        self.add_function_to_menu('View', self.display.View_Top)
        self.add_function_to_menu('View', self.display.View_Left)
        self.add_function_to_menu('View', self.display.View_Right)
        self.add_function_to_menu('View', self.display.View_Front)
        self.add_function_to_menu('View', self.display.View_Rear)
        self.add_function_to_menu('View', self.display.View_Iso)
        
        self.add_menu('Screencapture')
        self.add_function_to_menu('Screencapture', export_png)
        self.add_function_to_menu('Screencapture', export_jpg)
        self.add_function_to_menu('Screencapture', export_pdf)
        self.add_function_to_menu('Screencapture', export_svg)
        self.add_function_to_menu('Screencapture', export_ps)
       
        return (self.display, self.start_display, self.add_menu, self.add_function_to_menu)
        
    
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
        
        self.cbm_display_config()
        #display_custome_shape(display,seöf.SegmentLst[0].wire,2,0,[0,0,0])

        if self.config.setup['input_type'] == 3 or self.config.setup['input_type'] == 4:
            self.display.Context.SetDeviationAngle(1e-6)       
            self.display.Context.SetDeviationCoefficient(1e-6) 
            
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
                
        self.cbm_display_config()
        for c in self.mesh:
            self.display.DisplayColoredShape(c.wire, 'BLACK')    
        self.display.View_Top()
        self.display.FitAll()
        self.start_display()   
        return None
    
    
    def cbm_set_DymoreMK(self, x_offset = 0):
        '''
        Converts the Units of CBM to DYMORE/PYMORE/MARC units and returns the 
        array of the beamproperties with Massterms(6), Stiffness(21), 
        damping(1) and curvilinear coordinate(1)

        Parameters
        ----------
        x_offset : float, optional
            the x_offset parameter can be used if the DYMORE coordinate system 
            is different to the SONATA coordinate system of the rotor blade. 
            For example for the UH-60A Rotor model it is x_offset = 0.81786984, 
            because the curvilinear_coordinates "eta" are starting at the blade 
            root attachment, and not at the center of rotation.

        Returns
        ----------
        arr : ndarray
            [Massterms(6) (m00, mEta2, mEta3, m33, m23, m22) 
            Stiffness(21) (k11, k12, k22, k13, k23, k33,... k16, k26, ...k66)
            Viscous Damping(1) mu, Curvilinear coordinate(1) eta]
        
        '''
        MM = self.BeamProperties.MM_convert_units()
        MASS = np.array([MM[0,0], MM[2,3], MM[0,4], MM[5,5], MM[4,5], MM[4,4]])
        TS_u = np.triu(self.BeamProperties.TS_convert_units())
        TS_f = TS_u.flatten('F')
        STIFF = TS_f[TS_f != 0]
        mu = 0.0
        eta = self.config.setup['radial_station']/1000 - x_offset 
        return np.hstack((MASS,STIFF,mu,eta))
        
#%%############################################################################
#                           M    A    I    N                                  #
###############################################################################    
if __name__ == '__main__':   
    plt.close('all')
    fname = 'jobs/VariSpeed/uh60a_cbm_advanced/sec_config.yml'
    #fname = 'jobs/VariSpeed/uh60a_cbm_simple/sec_config.yml'
    #fname = 'jobs/AREA/R250/sec_config.yml'
    #fname = 'jobs/PBortolotti/sec_config.yml'
    config = CBMConfig(fname)
    
    job = CBM(config)
    
    job.cbm_gen_topo()
    job.cbm_gen_mesh()
    
    job.cbm_review_mesh()
    job.cbm_run_vabs()
    job.cbm_post_2dmesh(title='Hello World!')
    
    
#    #job.cbm_post_3dtopo()
#    job.config.vabs_cfg.recover_flag = 1
#    job.config.vabs_cfg.M = [0,2000e4,0]
