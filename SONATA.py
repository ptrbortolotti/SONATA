# -*- coding: utf-8 -*-
"""
THIS IS THE MAIN SONATA EXECUTION FILE!
@author: TPflumm
"""
#Basic PYTHON Modules:
import numpy as np       
import pickle
import sys
import subprocess
import matplotlib as plt
from functools import partial
from datetime import datetime

#PythonOCC Libraries
from OCC.Display.SimpleGui import init_display
from OCC.STEPControl import STEPControl_Writer, STEPControl_AsIs, STEPControl_GeometricCurveSet
from OCC.Interface import Interface_Static_SetCVal

#SONATA modules:
from SONATA.fileIO.readinput import section_config, read_material_input,str2bool
from SONATA.fileIO.CADinput import import_2d_stp, import_3d_stp, load_3D
 
from SONATA.topo.segment import Segment
from SONATA.topo.web import Web
from SONATA.topo.weight import Weight
from SONATA.topo.utils import  getID
                            
from SONATA.topo.BSplineLst_utils import reverse_BSplineLst, BSplineLst_Orientation, \
                            get_BSplineLst_length,copy_BSplineLst, trim_BSplineLst
                            
from SONATA.mesh.mesh_byprojection import mesh_by_projecting_nodes_on_BSplineLst
from SONATA.mesh.mesh_core import gen_core_cells
from SONATA.mesh.mesh_utils import first_stage_improvements,second_stage_improvements, determine_a_nodes, equidistant_nodes_on_BSplineLst, sort_and_reassignID                                        
from SONATA.mesh.mesh_intersect import map_mesh_by_intersect_curve2d

from SONATA.vabs.VABS_interface import VABS_config, export_cells_for_VABS, XSectionalProperties
from SONATA.vabs.strain import Strain
from SONATA.vabs.stress import Stress

from SONATA.display.display_mesh import plot_cells
from SONATA.display.display_utils import export_to_JPEG, export_to_PNG, export_to_PDF, export_to_SVG, export_to_PS, export_to_EnhPS, export_to_TEX, \
                                          export_to_BMP,export_to_TIFF, show_coordinate_system, display_SONATA_SegmentLst, display_custome_shape                        
#%%#############################################################################
#                           M    A    I    N                                  #
###############################################################################
#filename = str(sys.argv[1])#           #to run SONATA from command       
#SHOW = str2bool(sys.argv[2])
filename = 'sec_config.input'
SHOW = True
startTime = datetime.now()

#=========READ INPUT:===============
print "STATUS:\t Reading Crossection Configuration File"
Configuration = section_config(filename)
MaterialLst = read_material_input(Configuration.SETUP_mat_filename)

if SHOW:
    #===========DISPLAY CONFIG:===============
    display, start_display, add_menu, add_function_to_menu = init_display('wx')
    display.Context.SetDeviationAngle(0.000001)       # 0.001 default. Be careful to scale it to the problem.
    display.Context.SetDeviationCoefficient(0.000001) # 0.001 default. Be careful to scale it to the problem. 
    #show_coordinate_system(display) #CREATE AXIS SYSTEM for Visualization

#%%============================================================================ 
#               TOPOLOGY AND GEOMETRY:
# =============================================================================           
SegmentLst = []   #List of Segment Objects
for i,item in enumerate(Configuration.SEG_ID):
    if item == 0:        
        if Configuration.SETUP_input_type == 0:   #0) Airfoil from UIUC Database  --- naca23012
            SegmentLst.append(Segment(item, Layup = Configuration.SEG_Layup[i], CoreMaterial = Configuration.SEG_CoreMaterial[i], scale_factor = Configuration.SETUP_scale_factor, Theta = Configuration.SETUP_Theta, OCC=False, airfoil = Configuration.SETUP_datasource))
        
        elif Configuration.SETUP_input_type == 1: #1) Geometry from .dat file --- AREA_R250.dat
            SegmentLst.append(Segment(item, Layup = Configuration.SEG_Layup[i], CoreMaterial = Configuration.SEG_CoreMaterial[i], scale_factor = Configuration.SETUP_scale_factor,  Theta = Configuration.SETUP_Theta, OCC=False, filename = Configuration.SETUP_datasource))
        
        elif Configuration.SETUP_input_type == 2: #2)2d .step or .iges  --- AREA_R230.stp
            BSplineLst = import_2d_stp(Configuration.SETUP_datasource, Configuration.SETUP_scale_factor,Configuration.SETUP_Theta)
            SegmentLst.append(Segment(item, Layup = Configuration.SEG_Layup[i], CoreMaterial = Configuration.SEG_CoreMaterial[i], Theta = Configuration.SETUP_Theta, OCC=True, Boundary = BSplineLst))
        
        elif Configuration.SETUP_input_type == 3: #3)3D .step or .iges and radial station of crosssection --- AREA_Blade.stp, R=250
            BSplineLst = import_3d_stp(Configuration.SETUP_datasource,Configuration.SETUP_radial_station,Configuration.SETUP_scale_factor,Configuration.SETUP_Theta)
            SegmentLst.append(Segment(item, Layup = Configuration.SEG_Layup[i], CoreMaterial = Configuration.SEG_CoreMaterial[i], Theta = Configuration.SETUP_Theta, OCC=True, Boundary = BSplineLst))  

        else:
            print 'ERROR: \t WRONG input_type'
 
    else:
        SegmentLst.append(Segment(item, Layup = Configuration.SEG_Layup[i], CoreMaterial = Configuration.SEG_CoreMaterial[i],Theta = Configuration.SETUP_Theta))
sorted(SegmentLst, key=getID)  

#====================Build SEGMENT 0:========================================
SegmentLst[0].build_wire()
SegmentLst[0].build_layers()
SegmentLst[0].determine_final_boundary()    #Determine Boundary from Segment 0:
    
#==================== Build Webs:============================================
#TODO: CHECK IF WEB DEFINITION INTERSECT EACH OTHER
#TODO: SORT WEBS BY POS1 VALUES:
WebLst = []
if Configuration.SETUP_NbOfWebs > 0:
    for i in range(0,Configuration.SETUP_NbOfWebs):
        print 'STATUS: \t Building Web %s' %(i+1)
        WebLst.append(Web(Configuration.WEB_ID[i],Configuration.WEB_Pos1[i],Configuration.WEB_Pos2[i],SegmentLst[0].BSplineLst, SegmentLst[0].final_Boundary_BSplineLst))
    sorted(SegmentLst, key=getID)  
    
#======================Build remaining SEGMENTS =============================
if Configuration.SETUP_NbOfWebs > 0:
    for i,seg in enumerate(SegmentLst[1:],start=1):
        seg.build_segment_boundary_from_WebLst(WebLst,SegmentLst[0].final_Boundary_BSplineLst)
        seg.build_layers()


#======================Balance Weight========================================
if Configuration.SETUP_BalanceWeight == True:
    print 'STATUS: \t Building Balance Weight'   
    BW = Weight(0,Configuration.BW_XPos,Configuration.BW_YPos,Configuration.BW_Diameter,Configuration.BW_MatID)

    
#====================STEP-EXPORT=============================================
step_writer = STEPControl_Writer()  # initialize the STEP exporte
Interface_Static_SetCVal("write.step.schema", "AP203")
#step_writer.Transfer(SegmentLst[0].wire, STEPControl_AsIs)
#step_writer.Transfer(layer.wire, STEPControl_AsIs)
#step_writer.Transfer(BW.Curve, STEPControl_AsIs)      
#status = step_writer.Write("SONATA.stp")    
#assert(status == IFSelect_RetDone)

#=====================PICKLE TOPOLOGY:======================================
output_filename = filename.replace('.input', '.pkl')
with open(output_filename, 'wb') as output:
    pickle.dump(SegmentLst, output, protocol=pickle.HIGHEST_PROTOCOL)

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
# =============================================================================  
Resolution = Configuration.SETUP_mesh_resolution # Nb of Points on Segment0
length = get_BSplineLst_length(SegmentLst[0].BSplineLst)
global_minLen = round(length/Resolution,5)

#===================MESH SEGMENT===============================================
mesh = []
disco_nodes = []
for j,seg in enumerate(reversed(SegmentLst)):
    for i,layer in enumerate(reversed(seg.LayerLst)):
        print 'STATUS: \t Meshing Segment %s, Layer %s' %(seg.ID,len(seg.LayerLst)-i)
        
        a_BSplineLst = layer.BSplineLst       
        b_BSplineLst = trim_BSplineLst(layer.Boundary_BSplineLst, layer.S1, layer.S2, 0, 1)
        if BSplineLst_Orientation(b_BSplineLst,11) == False:
            b_BSplineLst = reverse_BSplineLst(b_BSplineLst)  
         
        if i==0:
            a_nodes = equidistant_nodes_on_BSplineLst(a_BSplineLst, True, True, True, minLen = global_minLen, LayerID = layer.ID[0])
        else: 
            a_nodes = determine_a_nodes(mesh,a_BSplineLst,global_minLen,layer.ID[0])
    
        #TODO: Scale tolerance to problem size!    
        a_nodes, b_nodes, cells = mesh_by_projecting_nodes_on_BSplineLst(a_BSplineLst,a_nodes,b_BSplineLst,layer.thickness,global_minLen,1e-1) 
        
#        if i==9:
#            for n in b_nodes:
#                display.DisplayShape(n.Pnt2d)
#            for c in cells:
#                display.DisplayShape(c.wire)
        

        enhanced_cells = first_stage_improvements(cells,b_BSplineLst,global_minLen)
        for c in enhanced_cells:
            c.calc_theta_1()
        enhanced_cells = second_stage_improvements(enhanced_cells,b_BSplineLst,global_minLen)

        for c in enhanced_cells:
            c.theta_3 = layer.Orientation
            c.MatID = int(layer.MatID)
            c.structured = True
            #print cell,'\t',cell.theta_3,cell.theta_1,cell.MatID,cell.area

        layer.cells = enhanced_cells
        mesh.extend(enhanced_cells)

    
    #===================MESH CORE================================================
    print 'STATUS: \t Meshing Segment %s, Core' %(seg.ID)
    core_Boundary_BSplineLst = []
    core_Boundary_BSplineLst += trim_BSplineLst(SegmentLst[-1].LayerLst[-1].Boundary_BSplineLst, 0, SegmentLst[-1].LayerLst[-1].S1, 0, 1)  #start und ende der lage
    core_Boundary_BSplineLst += copy_BSplineLst(SegmentLst[-1].LayerLst[-1].BSplineLst)
    core_Boundary_BSplineLst += trim_BSplineLst(SegmentLst[-1].LayerLst[-1].Boundary_BSplineLst, SegmentLst[-1].LayerLst[-1].S2, 1, 0, 1)  #start und ende der lage
        
    a_nodes = determine_a_nodes(mesh,core_Boundary_BSplineLst,global_minLen,layer.ID[0])
    [c_cells,c_nodes] = gen_core_cells(a_nodes,0.7)
    
    for c in c_cells:
        c.structured = False
        c.theta_3 = 0
        c.MatID = int(SegmentLst[0].CoreMaterial)
        c.calc_theta_1()
    
    mesh.extend(c_cells)
    mesh,nodes = sort_and_reassignID(mesh)


#%% BALANCE WEIGHT - CUTTING HOLE ALGORITHM====================================
#if Configuration.SETUP_BalanceWeight == True:
#    print 'STATUS: \t Meshing Balance Weight'   
#
#    mesh,boundary_nodes = map_mesh_by_intersect_curve2d(mesh,BW.Curve,BW.wire,display=display) 
#    triangle_options = 'pa.3'
#    [bw_cells,bw_nodes] = gen_core_cells(boundary_nodes,options=triangle_options)
#    
#    for c in bw_cells:
#        c.structured = False
#        c.theta_3 = 0
#        c.MatID = Configuration.BW_MatID
#        c.calc_theta_1()
#    
#    mesh.extend(bw_cells)
#
#    for c in mesh:
#         display.DisplayShape(c.wire,color='BLACK')         
#
#mesh,nodes = sort_and_reassignID(mesh)

#%%=====================PICKLE MESH ===========================================
output_filename = filename.replace('.input', '_mesh.pkl')
with open(output_filename, 'wb') as output:
    pickle.dump(mesh, output, protocol=pickle.HIGHEST_PROTOCOL)

input_filename = filename.replace('.input', '_mesh.pkl')
with open(input_filename, 'rb') as handle:
    mesh = pickle.load(handle)   
mesh,nodes = sort_and_reassignID(mesh)
#====================REVIEW==================================================
print 'STATUS: \t MESHING COMPLETED:'
print '\t   - Total Number of Cells: %s' %(len(mesh))
print '\t   - Duration: %s' % (datetime.now() - startTime)
print '\t   - Saved as: %s' % filename 
minarea = min([c.area for c in mesh])
print '\t   - smallest cell area: %s' % minarea 
minimum_angle = min([c.minimum_angle for c in mesh])
print '\t   - smallest angle [deg]: %s' % minimum_angle 
#orientation = all([c.orientation for c in mesh])
#print '\t   - Orientation [CC]: %s' % orientation 


#%%============================================================================ 
#                           V A B S
# =============================================================================
#TODO: BE CAREFUL TO USE THE RIGHT COORDINATE SYSTEM FOR THE CALCULATIONS!!!!  
filename = filename.replace('.input', '.vab')
VABSsetup = VABS_config(recover_flag=0)
VABSsetup.F = [0,0,0]    #in Newton
VABSsetup.M = [0,220e3,0]     #in Newton/mm

print 'STATUS: \t RUNNING VABS for Constitutive modeling:'
#EXECUTE VABS:
if VABSsetup.recover_flag == 1:
    VABSsetup.recover_flag=0
    export_cells_for_VABS(mesh,nodes,filename,VABSsetup,MaterialLst)
    command = 'VABSIII.exe '+ filename
    stdout = subprocess.check_output(command, shell=True)
    VABSsetup.recover_flag=1
    export_cells_for_VABS(mesh,nodes,filename,VABSsetup,MaterialLst)
    print 'STATUS: \t RUNNING VABS for 3D Recovery:'
    command = 'VABSIII.exe '+ filename
    stdout = stdout + subprocess.check_output(command, shell=True)
    
else:
    export_cells_for_VABS(mesh,nodes,filename,VABSsetup,MaterialLst)
    command = 'VABSIII.exe '+ filename
    stdout = subprocess.check_output(command, shell=True)

stdout = stdout.replace('\r\n\r\n','\n\t   -')
stdout = stdout.replace('\r\n','\n\t   -')
stdout = 'STATUS: \t VABS CALCULATIONS COMPLETED: \n\t   -' + stdout
print stdout 


#%%============================================================================ 
#                 P O S T  -  P R O C E S S I N G
# =============================================================================  
print 'STATUS: \t POST-PROCESSING:'
print 'INFO: \t Total duration: %s' % (datetime.now() - startTime)
filename_K = filename+'.K'
BeamProperties = XSectionalProperties(filename_K)
plot_cells(mesh, nodes, 'MatID',BeamProperties,'NACA0012, 150mm chord')


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
    
    plot_cells(mesh, nodes, 'stress.sigma11', BeamProperties,'NACA0012, 150mm chord')


#%%===================CROSS-CHECKING====================================================
#if VABSsetup.recover_flag == 1:
#    Integrator = 0
#    for c in mesh:
#        Integrator += c.stress.sigma11*c.area
#    print Integrator


#%%====================OCC-DISPLAY=============================================
if SHOW:
    
#    if Configuration.SETUP_input_type == 3:
#        aResShape = load_3D(Configuration.SETUP_datasource)
#        display.Context.SetDeviationAngle(0.0001)       # 0.001 default. Be careful to scale it to the problem.
#        display.Context.SetDeviationCoefficient(0.0001) # 0.001 default. Be careful to scale it to the problem. 
#        display.DisplayShape(aResShape,color=None, transparency=0.7, update=True)
        
    for c in mesh:
        display.DisplayColoredShape(c.wire, 'BLACK')
    
    display_custome_shape(display,SegmentLst[0].wire,2,0,[0,0,0])
    display.DisplayShape(SegmentLst[0].BSplineLst[0].StartPoint())
    display_SONATA_SegmentLst(display,SegmentLst)
    
    #display.DisplayShape(BW.Curve, color="BLACK")
    display.set_bg_gradient_color(20,6,111,200,200,200)
    show_coordinate_system(display,5)
    
    add_menu('screencapture')
    add_function_to_menu('screencapture','export to PDF', partial(export_to_PDF,display))
    add_function_to_menu('screencapture','export to SVG', partial(export_to_SVG,display))
    add_function_to_menu('screencapture','export to PS', partial(export_to_PS,display))
    add_function_to_menu('screencapture','export to EnhPS', partial(export_to_EnhPS,display))
    add_function_to_menu('screencapture','export to TEX', partial(export_to_TEX,display))
    add_function_to_menu('screencapture','export to BMP', partial(export_to_BMP,display))
    add_function_to_menu('screencapture', 'export to PNG', partial(export_to_PNG,display))
    add_function_to_menu('screencapture', 'export to JPEG', partial(export_to_JPEG,display))
    add_function_to_menu('screencapture', 'export to TIFF', partial(export_to_TIFF,display))
    
    display.View_Top()
    display.FitAll()
    start_display()

