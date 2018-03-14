# -*- coding: utf-8 -*-
"""Define the Composite Beam Model (CBM) class
Created on Wed Jan 03 13:56:37 2018
@author: TPflumm
 """

#Basic PYTHON Modules:
import pickle as pkl
import matplotlib.pyplot as plt
from datetime import datetime
import subprocess
import sys,os,math
import numpy as np


#PythonOCC Modules
from OCC.Display.SimpleGui import init_display

#SONATA modules:
from SONATA.cbm.fileIO.CADoutput import export_to_step
from SONATA.cbm.fileIO.CADinput import load_3D, import_2d_stp, import_3d_stp
from SONATA.cbm.fileIO.readinput import read_material_input
from SONATA.cbm.fileIO.configuration import Configuration

from SONATA.cbm.bladegen.blade import Blade

from SONATA.cbm.topo.segment import Segment
from SONATA.cbm.topo.web import Web
from SONATA.cbm.topo.utils import  getID
from SONATA.cbm.topo.weight import Weight
from SONATA.cbm.topo.BSplineLst_utils import get_BSplineLst_length

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
                                        display_custome_shape   
                                        


class CBM(object):
    ''' This Class includes the SONATA Dicipline Module for Structural 
    Composite Beam Modelling (CBM). 
    Design Variables are passed in form of the configuration Object or a 
    configuration file.          

    Attributes
    ----------
    config : <Configuration>
        Pointer to the <Configuration> object.
    MaterialLst: list of <Materials> object instances.
    SegmentLst: list of <Segment> object instances
    
    Notes
    ----------
    -For computational efficiency it is sometimes not suitable to 
     recalculate the topology or the crosssection every iteration, 
     maybe design flags to account for that.
    -make it possible to construct an instance by passing the 
     topology and/or mesh

    '''
        
    #__slots__ = ('Configuration','MaterialLst','__tel','__email','__alter','__partner')
    def __init__(self, Configuration):
        """
        Initialize attributes.

        Parameters
        ----------
        Configuration : <Configuration>
            Pointer to the <Configuration> object.
        MaterialLst : <MaterialLst>
        Pointer to the  <MaterialLst> object.

        """
        
        self.config = Configuration
        self.MaterialLst = read_material_input(self.config.SETUP_mat_filename)
        
        self.SegmentLst = []
        self.WebLst = []    
        self.BW = None
                
        self.mesh = []
        self.BeamProperties = None
        self.display = None
              
        self.startTime = datetime.now()
        self.exportLst = [] #list that contains all objects to be exported as step   
        self.surface3d = None
        self.Blade = None
        
        if self.config.SETUP_input_type == 3:
            self.surface3d = load_3D(self.config.SETUP_datasource)
        elif self.config.SETUP_input_type == 4:
            self.blade =  Blade(self.config.SETUP_datasource, self.config.SETUP_datasource, False, False)
            self.surface3d = self.blade.surface
    

    def __getstate__(self):
        """Return state values to be pickled."""
        return (self.config, self.MaterialLst, self.SegmentLst, self.WebLst, self.BW, self.mesh, self.BeamProperties)   
    
    
    def __setstate__(self, state):
        """Restore state from the unpickled state values."""
        (self.config, self.MaterialLst, self.SegmentLst, self.WebLst, self.BW, self.mesh, self.BeamProperties)  = state
        if self.config.SETUP_input_type == 3:
            self.surface3d = load_3D(self.config.SETUP_datasource)
        elif self.config.SETUP_input_type == 4:
            self.blade =  Blade(self.config.SETUP_datasource, self.config.SETUP_datasource, False, False)
            self.surface3d = self.blade.surface
    

    def cbm_save(self, output_filename=None):
        '''saves the complete <CBM> object as pickle
        Args:
            filename: <string> of the configuration file with path.        
        '''
        if output_filename is None:
            output_filename = self.config.filename
            output_filename = output_filename.replace('.input', '.pkl')
   
        with open(output_filename, 'wb') as output:
            pkl.dump(self, output, protocol=pkl.HIGHEST_PROTOCOL)        
        return None


    def cbm_load(self, input_filename=None):
        '''loads the complete <CBM> object as pickle
        Args:
            filename: <string> of the configuration file with path.    
        '''
        if input_filename is None:
            input_filename = self.config.filename
            input_filename = input_filename.replace('.input', '.pkl') 
        
        with open(input_filename, 'rb') as handle:
            tmp_dict = pkl.load(handle, encoding='latin1').__dict__
            self.__dict__.update(tmp_dict)    
        return None


    def cbm_save_topo(self, output_filename=None):
        '''saves the topology (SegmentLst, WebLst, and BW) as pickle
        Args:
            filename: <string> of the configuration file with path.        
        '''
        if output_filename is None:
            output_filename = self.config.filename
            output_filename = output_filename.replace('.input', '_topo.pkl')
            
        with open(output_filename, 'wb') as output:
            pkl.dump((self.SegmentLst, self.WebLst, self.BW), output, protocol=pkl.HIGHEST_PROTOCOL)
        return None


    def cbm_load_topo(self, input_filename=None):
        '''loads the topology (SegmentLst, WebLst, and BW) as pickle
        Args:
            filename: <string> of the configuration file with path.        
        '''
        if input_filename is None:
            input_filename = self.config.filename
            input_filename = input_filename.replace('.input', '_topo.pkl')
        
        with open(input_filename, 'rb') as handle:
            (self.SegmentLst, self.WebLst, self.BW) = pkl.load(handle)
    
        #Build wires for each layer and segment
        for seg in self.SegmentLst:
            seg.build_wire()
            for layer in seg.LayerLst:
                layer.build_wire()
        return None
        
    
    def cbm_stpexport_topo(self, export_filename=None):
        if export_filename is None:
            export_filename = self.config.filename
            export_filename = export_filename.replace('.input', '.stp')
        
        self.exportLst.append(self.SegmentLst[0].wire)
        for seg in self.SegmentLst:
            for layer in seg.LayerLst:
                self.exportLst.append(layer.wire)            
            
        print('STATUS:\t Exporting Topology to: ', export_filename)   
        export_to_step(self.exportLst, export_filename)
        return None
    
    
    def cbm_save_mesh(self, output_filename=None):
        '''saves the mesh (self.mesh) as pickle 
        '''
        if output_filename is None:
            output_filename = self.config.filename
            output_filename = output_filename.replace('.input', '_mesh.pkl')
            print('STATUS:\t Saving Mesh to: ',output_filename)
            
        with open(output_filename, 'wb') as output:
            pkl.dump(self.mesh, output, protocol=pkl.HIGHEST_PROTOCOL)
        return None
    
    
    def cbm_load_mesh(self, input_filename=None):
        '''loads the mesh (self.mesh) as pickle
        '''
        if input_filename is None:
            input_filename = self.config.filename
            input_filename = input_filename.replace('.input', '_mesh.pkl')
        
        with open(input_filename, 'rb') as handle:
            mesh = pkl.load(handle)   
        (self.mesh, nodes) = sort_and_reassignID(mesh)
        return None

    def cbm_generate_SegmentLst(self):
        ''' generate Segment Lst'''
        self.SegmentLst = []   #List of Segment Objects
        for i,item in enumerate(self.config.SEG_ID):
            if item == 0:        
                if self.config.SETUP_input_type == 0:   #0) Airfoil from UIUC Database  --- naca23012
                    self.SegmentLst.append(Segment(item, Layup = self.config.SEG_Layup[i], CoreMaterial = self.config.SEG_CoreMaterial[i], scale_factor = self.config.SETUP_scale_factor, Theta = self.config.SETUP_Theta, OCC=False, airfoil = self.config.SETUP_datasource))
                
                elif self.config.SETUP_input_type == 1: #1) Geometry from .dat file --- AREA_R250.dat
                    self.SegmentLst.append(Segment(item, Layup = self.config.SEG_Layup[i], CoreMaterial = self.config.SEG_CoreMaterial[i], scale_factor = self.config.SETUP_scale_factor,  Theta = self.config.SETUP_Theta, OCC=False, filename = self.config.SETUP_datasource))
                
                elif self.config.SETUP_input_type == 2: #2)2d .step or .iges  --- AREA_R230.stp
                    BSplineLst = import_2d_stp(self.config.SETUP_datasource, self.config.SETUP_scale_factor,self.config.SETUP_Theta)
                    self.SegmentLst.append(Segment(item, Layup = self.config.SEG_Layup[i], CoreMaterial = self.config.SEG_CoreMaterial[i], Theta = self.config.SETUP_Theta, OCC=True, Boundary = BSplineLst))
                
                elif self.config.SETUP_input_type == 3: #3)3D .step or .iges and radial station of crosssection --- AREA_Blade.stp, R=250
                    BSplineLst = import_3d_stp(self.config.SETUP_datasource, self.config.SETUP_radial_station, self.config.SETUP_scale_factor, self.config.SETUP_Theta)
                    self.SegmentLst.append(Segment(item, Layup = self.config.SEG_Layup[i], CoreMaterial = self.config.SEG_CoreMaterial[i], Theta = self.config.SETUP_Theta, OCC=True, Boundary = BSplineLst))  
        
                elif self.config.SETUP_input_type == 4: #4)generate 3D-Shape from twist,taper,1/4-line and airfoils, --- examples/UH-60A, R=4089, theta is given from twist distribution
                    BSplineLst = self.blade.get_crosssection(self.config.SETUP_radial_station, self.config.SETUP_scale_factor)
                    self.SegmentLst.append(Segment(item, Layup = self.config.SEG_Layup[i], CoreMaterial = self.config.SEG_CoreMaterial[i], Theta = self.blade.get_Theta(self.config.SETUP_radial_station), OCC=True, Boundary = BSplineLst))  
                    
                else:
                    print('ERROR:\t WRONG input_type')
         
            else:
                if self.config.SETUP_input_type == 4:
                    self.SegmentLst.append(Segment(item, Layup = self.config.SEG_Layup[i], CoreMaterial = self.config.SEG_CoreMaterial[i], Theta = self.blade.get_Theta(self.config.SETUP_radial_station)))
                
                else:
                    self.SegmentLst.append(Segment(item, Layup = self.config.SEG_Layup[i], CoreMaterial = self.config.SEG_CoreMaterial[i], Theta = self.config.SETUP_Theta))
        
        sorted(self.SegmentLst, key=getID)  
        return None


    def cbm_gen_topo(self):
        '''generates the topology 
        '''               
        #Generate SegmentLst from config:
        self.SegmentLst = []
        self.cbm_generate_SegmentLst()
        #Build Segment 0:
        self.SegmentLst[0].build_wire()
        self.SegmentLst[0].build_layers()
        self.SegmentLst[0].determine_final_boundary()
        
        #Build Webs:    
        self.WebLst = []
        if self.config.SETUP_NbOfWebs > 0:
            for i in range(0,self.config.SETUP_NbOfWebs):
                print('STATUS:\t Building Web %s' %(i+1))
                self.WebLst.append(Web(self.config.WEB_ID[i],self.config.WEB_Pos1[i],self.config.WEB_Pos2[i],self.SegmentLst))
            sorted(self.SegmentLst, key=getID)  
            
        #Build remaining Segments:
        if self.config.SETUP_NbOfWebs > 0:
            for i,seg in enumerate(self.SegmentLst[1:],start=1):
                seg.Segment0 = self.SegmentLst[0]
                seg.WebLst = self.WebLst
                seg.build_segment_boundary_from_WebLst(self.WebLst,self.SegmentLst[0])            
                seg.build_layers(self.WebLst,self.SegmentLst[0])
                seg.determine_final_boundary(self.WebLst,self.SegmentLst[0])
                seg.build_wire()
              
        self.BW = None
        #Balance Weight:
        if self.config.SETUP_BalanceWeight == True:
            print('STATUS:\t Building Balance Weight')   
            self.BW = Weight(0,self.config.BW_XPos,self.config.BW_YPos,self.config.BW_Diameter,self.config.BW_MatID)
            
        return None



    def cbm_gen_mesh(self):
        """
        generates the dicretization of topology and stores the cells and nodes
        in both the <Layer> instances, the <Segment> instances and the attribute 
        self.mesh that is a list of <Cell> instances

        Returns
        -------
        None
        """
        self.mesh = []
        Node.class_counter = 1
        Cell.class_counter = 1
        #meshing parameters:  
        Resolution = self.config.SETUP_mesh_resolution # Nb of Points on Segment0
        length = get_BSplineLst_length(self.SegmentLst[0].BSplineLst)
        global_minLen = round(length/Resolution,5)
            
        core_cell_area = 1.25*global_minLen**2
        bw_cell_area = 0.6*global_minLen**2
        web_consolidate_tol = 0.5*global_minLen

        #===================MESH SEGMENT
        for j,seg in enumerate(reversed(self.SegmentLst)):
            self.mesh.extend(seg.mesh_layers(self.SegmentLst, global_minLen, self.WebLst, display=self.display))
            #mesh,nodes = sort_and_reassignID(mesh)
            
        #===================MESH CORE  
        if self.config.flag_mesh_core:
            for j,seg in enumerate(reversed(self.SegmentLst)):
                if seg.ID == 1:
                    core_cell_area = 0.5*global_minLen**2
                self.mesh.extend(seg.mesh_core(self.SegmentLst, self.WebLst, core_cell_area, display=self.display))
                
        #===================consolidate mesh on web interface         
        for web in self.WebLst:
            #print web.ID,  'Left:', SegmentLst[web.ID].ID, 'Right:', SegmentLst[web.ID+1].ID,
            print('STATUS:\t Consolidate Mesh on Web Interface ', web.ID)  
            web.wl_nodes = grab_nodes_of_cells_on_BSplineLst(self.SegmentLst[web.ID].cells, web.BSplineLst)
            web.wr_nodes = grab_nodes_of_cells_on_BSplineLst(self.SegmentLst[web.ID+1].cells, web.BSplineLst)
            self.mesh = consolidate_mesh_on_web(self.mesh,web.BSplineLst, web.wl_nodes, web.wr_nodes, web_consolidate_tol,self.display)
            
        #============= BALANCE WEIGHT - CUTTING HOLE ALGORITHM
        if self.config.SETUP_BalanceWeight == True:
            print('STATUS:\t Meshing Balance Weight')   
            
            self.mesh, boundary_nodes = map_mesh_by_intersect_curve2d(self.mesh,self.BW.Curve,self.BW.wire,global_minLen)
            #boundary_nodes = merge_nodes_if_too_close(boundary_nodes,self.BW.Curve,global_minLen,tol=0.05)
            
            
            [bw_cells,bw_nodes] = gen_core_cells(boundary_nodes, bw_cell_area)
            
            for c in bw_cells:
                c.structured = False
                c.theta_3 = 0
                c.MatID = self.config.BW_MatID
                c.calc_theta_1()
            
            self.mesh.extend(bw_cells)
       
        (self.mesh, nodes) = sort_and_reassignID(self.mesh)
        return None
    

    def cbm_review_mesh(self):
        '''prints a summary of the mesh properties to the screen '''
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

    
    def cbm_run_vabs(self):
        '''runs the solver
        
        Args:
            filename: <string> of the configuration file with path.        
        '''
        self.mesh,nodes = sort_and_reassignID(self.mesh)
        #TODO: BE CAREFUL TO USE THE RIGHT COORDINATE SYSTEM FOR THE CALCULATIONS!!!!  
        vabs_filename = self.config.filename.replace('.input', '.vab')
        print('STATUS:\t RUNNING VABS for Constitutive modeling:')
        #EXECUTE VABS:
        if self.config.VABS.recover_flag == 1:
            self.config.VABS.recover_flag=0
            export_cells_for_VABS(self.mesh,nodes,vabs_filename,self.config.VABS,self.MaterialLst)
            command = 'VABSIII.exe '+ vabs_filename
            stdout = subprocess.check_output(command, shell=True)
            self.config.VABS.recover_flag=1
            export_cells_for_VABS(self.mesh,nodes,vabs_filename,self.config.VABS,self.MaterialLst)
            print('STATUS:\t RUNNING VABS for 3D Recovery:')
            command = 'VABSIII.exe '+ vabs_filename
            stdout = stdout + subprocess.check_output(command, shell=True)
            
        else:
            export_cells_for_VABS(self.mesh,nodes,vabs_filename,self.config.VABS,self.MaterialLst)
            command = 'VABSIII.exe '+ vabs_filename
            stdout = subprocess.check_output(command, shell=True)
            
        stdout = stdout.decode('utf-8')
        stdout = stdout.replace('\r\n\r\n','\n\t   -')
        stdout = stdout.replace('\r\n','\n\t   -')
        stdout = 'STATUS:\t VABS CALCULATIONS COMPLETED: \n\t   -' + stdout
        print(stdout) 
        print('STATUS:\t Total Elapsed Time: %s' % (datetime.now() - self.startTime))
        
        #VABS Postprocessing:
        self.BeamProperties = XSectionalProperties(vabs_filename+'.K')
        if self.config.VABS.recover_flag == 1:
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
    
    
    def cbm_post_2dmesh(self, attribute='MatID',title='NOTITLE', **kw):
        mesh,nodes = sort_and_reassignID(self.mesh)
        fig,ax = plot_cells(self.mesh, nodes, attribute, self.BeamProperties, title, **kw)
        return fig,ax
        
        
    def cbm_display_config(self):
        #===========DISPLAY CONFIG:===============
        self.display, self.start_display, self.add_menu, self.add_function_to_menu = init_display()
        self.display.Context.SetDeviationAngle(1e-6)       # 0.001 default. Be careful to scale it to the problem.
        self.display.Context.SetDeviationCoefficient(1e-6) # 0.001 default. Be careful to scale it to the problem. 
        self.display.set_bg_gradient_color(20,6,111,200,200,200)
        show_coordinate_system(self.display,5)
        return (self.display, self.start_display, self.add_menu, self.add_function_to_menu)
        
    
    def cbm_post_3dtopo(self):
        self.cbm_display_config()
        
        if self.config.SETUP_BalanceWeight:
            self.display.DisplayShape(self.BW.Curve, color="BLACK")
        #display.DisplayShape(SegmentLst[0].BSplineLst[0].StartPoint())
        #display_custome_shape(display,SegmentLst[0].wire,2,0,[0,0,0])

        if self.config.SETUP_input_type == 3 or self.config.SETUP_input_type == 4:
            self.display.Context.SetDeviationAngle(1e-5)       # 0.001 default. Be careful to scale it to the problem, or else it will crash :) 
            self.display.Context.SetDeviationCoefficient(1e-5) 
            
            display_SONATA_SegmentLst(self.display,self.SegmentLst,self.config.SETUP_radial_station,-math.pi/2,-math.pi/2)
            self.display.DisplayShape(self.surface3d, color=None, transparency=0.7, update=True)
               
        else:
            display_SONATA_SegmentLst(self.display,self.SegmentLst)
        
        #self.display.View_ISO()
        self.display.FitAll()
        self.start_display()   
        return None
    
    
    def cbm_post_3dmesh(self):
        self.cbm_display_config()
        for c in self.mesh:
            self.display.DisplayColoredShape(c.wire, 'BLACK')    
            
        self.display.View_Top()
        self.display.FitAll()
        self.start_display()   
        return None
    
    def cbm_set_DymoreMK(self, x_offset = 0):
        '''
        Returns the tablerow with the Massterms(6), Stiffness(21),damping(1) 
        and coordinate(1) stacked horizontally for the use in DYMORE/PYMORE/MARC        
        '''
        MM = self.BeamProperties.MM_convert_units()
        MASS = np.array([MM[0,0], MM[2,3], MM[0,4], MM[5,5], MM[4,5], MM[4,4]])
        #print '@MASS_TERMS: {m00, mEta2, mEta3, m33, m23, m22}'
        #print MASS
        TS_u = np.triu(self.BeamProperties.TS_convert_units())
        TS_f = TS_u.flatten('F')
        STIFF = TS_f[TS_f != 0]
        #print '@STIFFNESS_MATRIX {k11, k12, k22, k13, k23, k33,... k16, k26, ...k66}: '
        #print STIFF
        mu = 0.0
        #print '@VISCOUS DAMPING:', mu
        eta = self.config.SETUP_radial_station/1000 - x_offset 
        #print'@CURVILINEAR_COORDINATE:', eta
        return np.hstack((MASS,STIFF,mu,eta))
        

        
#%%############################################################################
#                           M    A    I    N                                  #
###############################################################################    
if __name__ == '__main__':   
    plt.close('all')    
    filename = 'jobs/VHeuschneider/sec_config.input'
    
    config = Configuration(filename)
    #TODO: include optionflags and Vabs_setup in Configuration
    MaterialLst = read_material_input(config.SETUP_mat_filename)
    
    job1 = CBM(config,MaterialLst)
    job1.cbm_gen_topo()
    job1.cbm_gen_mesh()
    #job1.cbm_run_vabs(filename)
    job1.cbm_post_2dmesh()
    job1.cbm_post_3dtopo()
    job1.cbm_post_3dmesh()
