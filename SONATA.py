# -*- coding: utf-8 -*-
"""
This is currently the main excecution script for SONATA.
Because it is currently still under extensive development the code isn't as 
clean as it could be but allows easy debugging,printing and 3D displaying.

SONATA is a preprocessor for parametric analysis and design of composite beam 
cross-sections in a multidisciplinary rotor design environment. A helicopter 
rotor blade represents a classical aeroelastic problem, where the aerodynamic 
behavior, the structural elasticity and vibrational dynamics have to be studied 
simultaneously.  While a geometric definition of a rotorblade with CAD tools is 
simple, the transfer to a meshed cross-sectional representation may prohibit 
automated design optimization. Consequently, most researches have developed 
individual parametric mesh generators for the cross-sectional analysis, that 
reduces their structural model to few design variables in the process. 
SONATA represents such a preprocessor. SONATA is written in python and is using
for a lot of operations the Opencascade (CAD) kernel with its python wrapper 
(pythonocc).

Because it is currently still under extensive development the code isn't as 
clean as it could be but allows easy debugging, printing and 3D displaying. In 
the future, the SONATA execution script shall inlude an openmdao structure 
which can call the then unlying functionalities. 

#def SONATA_CBM(Configuration,Flags,)
    ''' This function includes the SONATA Dicipline Module for Structural 
    Composite Beam Modelling (CBM).
        Design Variables are passed in form of the Configuration file. 
        
        NOTE: For computational efficiency it is sometimes not suitable to 
        recalculate the topology or the crosssection every iteration,
              -maybe design flags to account for that.
        
        return: BeamProperties(allready inlcude Postprocessed parameters such 
        as Failure Critirion and Safety Margin...)
        '''

Created on Thu Nov 02 14:36:29 2017
@author: TPflumm
"""

#Basic PYTHON Modules:
import numpy as np       
import pickle
import sys
import math
import subprocess
import itertools
import matplotlib.pyplot as plt
from functools import partial
from datetime import datetime

#PythonOCC Modules
from OCC.Display.SimpleGui import init_display
from OCC.gp import gp_Vec2d

#SONATA modules:
from SONATA.fileIO.readinput import section_config, read_material_input,str2bool
from SONATA.fileIO.CADinput import import_2d_stp, import_3d_stp, load_3D
from SONATA.fileIO.CADoutput import export_to_step
 
from SONATA.topo.segment import Segment, generate_SegmentLst
from SONATA.topo.web import Web
from SONATA.topo.weight import Weight
from SONATA.topo.utils import  getID                 
from SONATA.bladegen.blade import Blade

from SONATA.topo.BSplineLst_utils import reverse_BSplineLst, BSplineLst_Orientation, \
                            get_BSplineLst_length,copy_BSplineLst, trim_BSplineLst
                            
from SONATA.mesh.mesh_byprojection import mesh_by_projecting_nodes_on_BSplineLst
from SONATA.mesh.mesh_core import gen_core_cells
from SONATA.mesh.mesh_utils import modify_cornerstyle_one, modify_sharp_corners,second_stage_improvements,grab_nodes_of_cells_on_BSplineLst,\
                                 determine_a_nodes, equidistant_nodes_on_BSplineLst, sort_and_reassignID, find_cells_that_contain_node                                        
from SONATA.mesh.mesh_intersect import map_mesh_by_intersect_curve2d
from SONATA.mesh.consolidate_mesh import consolidate_mesh_on_web

from SONATA.vabs.VABS_interface import VABS_config, export_cells_for_VABS, XSectionalProperties
from SONATA.vabs.strain import Strain
from SONATA.vabs.stress import Stress

from SONATA.display.display_mesh import plot_cells
from SONATA.display.display_utils import export_to_JPEG, export_to_PNG, export_to_PDF, export_to_SVG, export_to_PS, export_to_EnhPS, export_to_TEX, \
                                          export_to_BMP,export_to_TIFF, show_coordinate_system, display_SONATA_SegmentLst, display_custome_shape                        
#%%#############################################################################
#                           M    A    I    N                                  #
###############################################################################
plt.close('all')

#filename = str(sys.argv[1])#           #to run SONATA from command       
filename = 'sec_config.input'

FLAG_TOPO = True
FLAG_MESH = True
FLAG_VABS = True
FLAG_SHOW_3D_TOPO = False
FLAG_SHOW_2D_MESH = True
FLAG_SHOW_3D_MESH = False
FLAG_EXPORT_STEP = False

startTime = datetime.now()
#=========READ INPUT:===============
print "STATUS:\t Reading Crossection Configuration File"
Configuration = section_config(filename)
MaterialLst = read_material_input(Configuration.SETUP_mat_filename)
    
#===========DISPLAY CONFIG:===============
if FLAG_SHOW_3D_TOPO or FLAG_SHOW_3D_MESH:
    display, start_display, add_menu, add_function_to_menu = init_display()
    display.Context.SetDeviationAngle(1e-6)       # 0.001 default. Be careful to scale it to the problem.
    display.Context.SetDeviationCoefficient(1e-6) # 0.001 default. Be careful to scale it to the problem. 

#%%============================================================================ 
#               TOPOLOGY AND GEOMETRY:
#=============================================================================           
if FLAG_TOPO:
    #Generate SegmentLst from Configuration:
    SegmentLst = generate_SegmentLst(Configuration)
    #Build Segment 0:
    SegmentLst[0].build_wire()
    SegmentLst[0].build_layers()
    SegmentLst[0].determine_final_boundary()        
    #Build Webs:
        #TODO: CHECK IF WEB DEFINITION INTERSECT EACH OTHER
        #TODO: SORT WEBS BY POS1 VALUES:
    WebLst = []
    if Configuration.SETUP_NbOfWebs > 0:
        for i in range(0,Configuration.SETUP_NbOfWebs):
            print 'STATUS:\t Building Web %s' %(i+1)
            WebLst.append(Web(Configuration.WEB_ID[i],Configuration.WEB_Pos1[i],Configuration.WEB_Pos2[i],SegmentLst[0].BSplineLst, SegmentLst[0].final_Boundary_BSplineLst))
        sorted(SegmentLst, key=getID)  
        
    #Build remaining Segments:
    if Configuration.SETUP_NbOfWebs > 0:
        for i,seg in enumerate(SegmentLst[1:],start=1):
            seg.build_segment_boundary_from_WebLst(WebLst,SegmentLst[0].final_Boundary_BSplineLst)
            seg.build_layers()
    
    #Balance Weight:
    if Configuration.SETUP_BalanceWeight == True:
        print 'STATUS:\t Building Balance Weight'   
        BW = Weight(0,Configuration.BW_XPos,Configuration.BW_YPos,Configuration.BW_Diameter,Configuration.BW_MatID)

    #====================STEP-EXPORT=============================================
    if FLAG_EXPORT_STEP:
        exportLst = [SegmentLst[0].wire]
        export_filename = filename.replace('.input', '.stp')
        print 'STATUS:\t Exporting Topology to: ', filename   
        export_to_step(exportLst,"SONATA.stp")
    
    #=====================PICKLE TOPOLOGY:======================================
    output_filename = filename.replace('.input', '.pkl')
    with open(output_filename, 'wb') as output:
        pickle.dump(SegmentLst, output, protocol=pickle.HIGHEST_PROTOCOL)
 
else: 
    #LOAD PICKLED TOPO
    input_filename = filename.replace('.input', '.pkl')
    with open(input_filename, 'rb') as handle:
        SegmentLst = pickle.load(handle)

    #Build wires for each layer and segment
    for seg in SegmentLst:
        seg.build_wire()
        for layer in seg.LayerLst:
            layer.build_wire()

#%%============================================================================ 
#                           M E S H
#==============================================================================
if FLAG_MESH: 
    #meshing parameters:  
    Resolution = Configuration.SETUP_mesh_resolution # Nb of Points on Segment0
    length = get_BSplineLst_length(SegmentLst[0].BSplineLst)
    global_minLen = round(length/Resolution,5)
       
    proj_tol_1 = 5e-2
    proj_tol_2 = 5e-2
    non_dct_factor = 2.6
    crit_angle_1 = 115
    alpha_crit_2 = 60
    growing_factor = 1.8   #critical growing factor of cell before splitting 
    shrinking_factor = 0.10  #critical shrinking factor for cells before merging nodes
    
    core_cell_area = 0.9*global_minLen**2
    web_consolidate_tol = 0.5*global_minLen
    
    mesh = []
    #===================MESH SEGMENT===============================================
    disco_nodes = []
    for j,seg in enumerate(reversed(SegmentLst)):
        for i,layer in enumerate(reversed(seg.LayerLst)):
            print 'STATUS:\t Meshing Segment %s, Layer %s' %(seg.ID,len(seg.LayerLst)-i)
            
            a_BSplineLst = layer.BSplineLst       
            b_BSplineLst = trim_BSplineLst(layer.Boundary_BSplineLst, layer.S1, layer.S2, 0, 1)
            if BSplineLst_Orientation(b_BSplineLst,11) == False:
                b_BSplineLst = reverse_BSplineLst(b_BSplineLst)  
             
            if (i == 0 and seg.ID != 0) or i == 0 and len(SegmentLst) == 1:
                a_nodes = equidistant_nodes_on_BSplineLst(a_BSplineLst, True, True, True, minLen = global_minLen, LayerID = layer.ID[0])
            else: 
                a_nodes = determine_a_nodes(mesh,a_BSplineLst,global_minLen,layer.ID[0],non_dct_factor)
            
            #TODO: Scale tolerance to problem size!    
            if FLAG_SHOW_3D_MESH: 
                displaymesh=display 
            else: 
                displaymesh=None   
                
            a_nodes, b_nodes, cells = mesh_by_projecting_nodes_on_BSplineLst(a_BSplineLst,a_nodes,b_BSplineLst,layer.thickness, proj_tol_1,crit_angle_1, display=displaymesh) 
            #enhanced_cells = modify_cornerstyle_one(cells,b_BSplineLst)
            enhanced_cells = modify_sharp_corners(cells,b_BSplineLst,global_minLen,layer.thickness, proj_tol_2,alpha_crit_2,display=displaymesh)
            enhanced_cells = second_stage_improvements(enhanced_cells,b_BSplineLst,global_minLen,growing_factor,shrinking_factor)
            
            for c in enhanced_cells:
                c.calc_theta_1()
                c.theta_3 = layer.Orientation
                c.MatID = int(layer.MatID)
                c.structured = True

            layer.cells = enhanced_cells
            mesh.extend(enhanced_cells) 
            mesh,nodes = sort_and_reassignID(mesh)
        
        #===================MESH CORE================================================
        if seg.ID==0 and len(SegmentLst)>1:
            pass
        
        else:
            print 'STATUS:\t Meshing Segment %s, Core' %(seg.ID)
            core_Boundary_BSplineLst = []
            if seg.LayerLst[-1].S1<seg.LayerLst[-1].S2:
                core_Boundary_BSplineLst += trim_BSplineLst(seg.LayerLst[-1].Boundary_BSplineLst, 0, seg.LayerLst[-1].S1, 0, 1)  #start und ende der lage
                core_Boundary_BSplineLst += copy_BSplineLst(seg.LayerLst[-1].BSplineLst)
                core_Boundary_BSplineLst += trim_BSplineLst(seg.LayerLst[-1].Boundary_BSplineLst, seg.LayerLst[-1].S2, 1, 0, 1)  #start und ende der lage
            
            #FIX: Why doesn't this occure in Segments.buildLayer?????
            else:
                core_Boundary_BSplineLst += copy_BSplineLst(seg.LayerLst[-1].BSplineLst)
                core_Boundary_BSplineLst += trim_BSplineLst(seg.LayerLst[-1].Boundary_BSplineLst, seg.LayerLst[-1].S2, seg.LayerLst[-1].S1, 0, 1)  #start und ende der lage
                   
            a_nodes = determine_a_nodes(mesh,core_Boundary_BSplineLst,global_minLen,layer.ID[0])           
            [c_cells,c_nodes] = gen_core_cells(a_nodes,core_cell_area)
            
            for c in c_cells:
                c.structured = False
                c.theta_3 = 0
                c.MatID = int(seg.CoreMaterial)
                c.calc_theta_1()
               
            mesh.extend(c_cells)

    #===================consolidate mesh on web interface==================   
    for seg in SegmentLst:
        for layer in seg.LayerLst:
            for c in layer.cells:
                seg.cells.append(c)
    
    for web in WebLst:
        #print web.ID,  'Left:', SegmentLst[web.ID].ID, 'Right:', SegmentLst[web.ID+1].ID,
        web.wl_nodes = grab_nodes_of_cells_on_BSplineLst(SegmentLst[web.ID].cells,web.BSplineLst)
        web.wr_nodes = grab_nodes_of_cells_on_BSplineLst(SegmentLst[web.ID+1].cells,web.BSplineLst)
        mesh = consolidate_mesh_on_web(mesh,web.BSplineLst,web.wl_nodes,web.wr_nodes,web_consolidate_tol,display)
        
         
    #============= BALANCE WEIGHT - CUTTING HOLE ALGORITHM====================================
    if Configuration.SETUP_BalanceWeight == True:
        print 'STATUS:\t Meshing Balance Weight'   
        
        if FLAG_SHOW_3D_MESH: 
            mesh,boundary_nodes = map_mesh_by_intersect_curve2d(mesh,BW.Curve,BW.wire,display=display) 
        else:
            mesh,boundary_nodes = map_mesh_by_intersect_curve2d(mesh,BW.Curve,BW.wire)
        triangle_options = 'pa.3'
        [bw_cells,bw_nodes] = gen_core_cells(boundary_nodes,options=triangle_options)
        
        for c in bw_cells:
            c.structured = False
            c.theta_3 = 0
            c.MatID = Configuration.BW_MatID
            c.calc_theta_1()
        
        mesh.extend(bw_cells)
   
    mesh,nodes = sort_and_reassignID(mesh)
    
    #=====================PICKLE MESH ===========================================
    output_filename = filename.replace('.input', '_mesh.pkl')
    with open(output_filename, 'wb') as output:
        pickle.dump(mesh, output, protocol=pickle.HIGHEST_PROTOCOL)
    
    #====================REVIEW==================================================
    print 'STATUS:\t MESHING COMPLETED:'
    print '\t   - Total Number of Cells: %s' %(len(mesh))
    print '\t   - Duration: %s' % (datetime.now() - startTime)
    print '\t   - Saved as: %s' % filename 
    minarea = min([c.area for c in mesh])
    print '\t   - smallest cell area: %s' % minarea 
    minimum_angle = min([c.minimum_angle for c in mesh])
    print '\t   - smallest angle [deg]: %s' % minimum_angle 
    #orientation = all([c.orientation for c in mesh])
    #print '\t   - Orientation [CC]: %s' % orientation 


else:
    #LOAD PICKLED MESH 
    input_filename = filename.replace('.input', '_mesh.pkl')
    with open(input_filename, 'rb') as handle:
        mesh = pickle.load(handle)   
    mesh,nodes = sort_and_reassignID(mesh)
      
    #====================REVIEW==================================================
    print 'STATUS:\t MESH LOADED:'
    print '\t   - from file: %s' % input_filename 
    print '\t   - Total Number of Cells: %s' %(len(mesh))
    minarea = min([c.area for c in mesh])
    print '\t   - smallest cell area: %s' % minarea 
    minimum_angle = min([c.minimum_angle for c in mesh])
    print '\t   - smallest angle [deg]: %s' % minimum_angle 
    #orientation = all([c.orientation for c in mesh])
    #print '\t   - Orientation [CC]: %s' % orientation 


#%%============================================================================ 
#                           V A B S
#==============================================================================
if FLAG_VABS:
    #TODO: BE CAREFUL TO USE THE RIGHT COORDINATE SYSTEM FOR THE CALCULATIONS!!!!  
    vabs_filename = filename.replace('.input', '.vab')
    VABSsetup = VABS_config(recover_flag=1)
    VABSsetup.F = [0,0,0]    #in Newton
    VABSsetup.M = [0,220e3,0]     #in Newton/mm
    
    print 'STATUS:\t RUNNING VABS for Constitutive modeling:'
    #EXECUTE VABS:
    if VABSsetup.recover_flag == 1:
        VABSsetup.recover_flag=0
        export_cells_for_VABS(mesh,nodes,vabs_filename,VABSsetup,MaterialLst)
        command = 'VABSIII.exe '+ vabs_filename
        stdout = subprocess.check_output(command, shell=True)
        VABSsetup.recover_flag=1
        export_cells_for_VABS(mesh,nodes,vabs_filename,VABSsetup,MaterialLst)
        print 'STATUS:\t RUNNING VABS for 3D Recovery:'
        command = 'VABSIII.exe '+ vabs_filename
        stdout = stdout + subprocess.check_output(command, shell=True)
        
    else:
        export_cells_for_VABS(mesh,nodes,vabs_filename,VABSsetup,MaterialLst)
        command = 'VABSIII.exe '+ vabs_filename
        stdout = subprocess.check_output(command, shell=True)
    
    stdout = stdout.replace('\r\n\r\n','\n\t   -')
    stdout = stdout.replace('\r\n','\n\t   -')
    stdout = 'STATUS:\t VABS CALCULATIONS COMPLETED: \n\t   -' + stdout
    print stdout 
    
    
    #VABS Postprocessing:
    BeamProperties = XSectionalProperties(vabs_filename+'.K')
    if VABSsetup.recover_flag == 1:
        BeamProperties.read_all_VABS_Results()
        #ASSIGN Stress and strains to elements:
        for i,c in enumerate(mesh):
            c.strain = Strain(BeamProperties.ELE[i][1:7])
            c.stress = Stress(BeamProperties.ELE[i][7:13])
            c.strainM = Strain(BeamProperties.ELE[i][13:19])
            c.stressM = Stress(BeamProperties.ELE[i][19:25])
        
        #ASSIGN Displacement U to nodes:
        for i,n in enumerate(nodes):
            n.displacement = BeamProperties.U[i][3:6]


##%%============================================================================ 
##                 P O S T  -  P R O C E S S I N G
## =============================================================================
print 	BeamProperties.Xm2

#====================2D: MATPLOTLIB-DISPLAY======================
if FLAG_SHOW_2D_MESH:
    
    if FLAG_VABS:
        plot_cells(mesh, nodes, 'MatID', BeamProperties)
    
    
    #plt.savefig('SONATA_MESH.pdf', dpi=900, facecolor='w', edgecolor='w',
    #    orientation='landscape', papertype='a4', format='pdf')

#   from matplotlib2tikz import save as tikz_save
#   folder = 'publication/preprocessor_paper/img/'
#   tikz_filename = folder + filename.replace('.input', '.tex')
#   tikz_save(tikz_filename, figureheight = '\\figureheight', figurewidth = '\\figurewidth' )    
#    
#    import os  
#    os.system("img/lualatex minimal_latex.tex")
    
    if VABSsetup.recover_flag == 1:
        plot_cells(mesh, nodes, 'stress.sigma11', BeamProperties,'NACA0012, 150mm chord')


#====================3D: OCC-DISPLAY=====================
if FLAG_SHOW_3D_TOPO or FLAG_SHOW_3D_MESH:
    display.set_bg_gradient_color(20,6,111,200,200,200)
    show_coordinate_system(display,5)
#    add_menu('screencapture')
#    add_function_to_menu('screencapture','export to PDF', partial(export_to_PDF,display))
#    add_function_to_menu('screencapture','export to SVG', partial(export_to_SVG,display))
#    add_function_to_menu('screencapture','export to PS', partial(export_to_PS,display))
#    add_function_to_menu('screencapture','export to EnhPS', partial(export_to_EnhPS,display))
#    add_function_to_menu('screencapture','export to TEX', partial(export_to_TEX,display))
#    add_function_to_menu('screencapture','export to BMP', partial(export_to_BMP,display))
#    add_function_to_menu('export to PNG', partial(export_to_PNG,display))
#    add_function_to_menu('screencapture', 'export to JPEG', partial(export_to_JPEG,display))
#    add_function_to_menu('screencapture', 'export to TIFF', partial(export_to_TIFF,display))
    
    if FLAG_SHOW_3D_TOPO:
        display_SONATA_SegmentLst(display,SegmentLst)
        if Configuration.SETUP_BalanceWeight:
            display.DisplayShape(BW.Curve, color="BLACK")
        #display.DisplayShape(SegmentLst[0].BSplineLst[0].StartPoint())
        #display_custome_shape(display,SegmentLst[0].wire,2,0,[0,0,0])

        if Configuration.SETUP_input_type == 3:
            display.Context.SetDeviationAngle(1e-5)       # 0.001 default. Be careful to scale it to the problem, or else it will crash :) 
            display.Context.SetDeviationCoefficient(1e-5) 
            
            display_SONATA_SegmentLst(display,SegmentLst,Configuration.SETUP_radial_station,-math.pi/2,-math.pi/2)
            aResShape = load_3D(Configuration.SETUP_datasource)
            display.DisplayShape(aResShape,color=None, transparency=0.7, update=True)
        
        elif Configuration.SETUP_input_type == 4:
            display.Context.SetDeviationAngle(1e-5)       # 0.001 default. Be careful to scale it to the problem, or else it will crash :) 
            display.Context.SetDeviationCoefficient(1e-5) 
            
            genblade = Blade(Configuration.SETUP_datasource,Configuration.SETUP_datasource,False,False)
            display_SONATA_SegmentLst(display,SegmentLst,Configuration.SETUP_radial_station,-math.pi/2,-math.pi/2) 
            display.DisplayShape(genblade.surface,color=None, transparency=0.7, update=True)
       
        else:
            display_SONATA_SegmentLst(display,SegmentLst)
             
    if FLAG_SHOW_3D_MESH:        
        for c in mesh:
            display.DisplayColoredShape(c.wire, 'BLACK')    
  
    display.View_Top()
    display.FitAll()
    start_display()

