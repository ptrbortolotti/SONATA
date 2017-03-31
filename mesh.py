#Basic PYTHON Modules:
import subprocess
import numpy as np
import math
import pickle
import matplotlib as plt
from functools import partial
from datetime import datetime

#PythonOCC Libraries
from OCC.Display.SimpleGui import init_display
from OCC.BRepExtrema import BRepExtrema_DistShapeShape 
from OCC.GCPnts import GCPnts_QuasiUniformAbscissa
from OCC.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.gp import gp_Pnt2d,gp_Lin2d, gp_Dir2d
from OCC.Geom2d import Geom2d_Line
from OCC.Geom2dAPI import Geom2dAPI_InterCurveCurve

#THE PYTHON SHAPELY MODULE:
import shapely.geometry as shp_geom

#SONATA modules:
from SONATA.topo.BSplineLst_utils import get_BSplineLst_length, get_BSpline_length, trim_BSplineLst, set_BSplineLst_to_Origin, \
                            BSplineLst_Orientation, reverse_BSplineLst, findPnt_on_BSplineLst, copy_BSplineLst,\
                            distance_on_BSplineLst, trim_BSplineLst_by_Pnt2d, trim_BSplineLst_by_coordinates, \
                            ProjectPointOnBSplineLst, findPnt_on_2dcurve
from SONATA.topo.weight import Weight

from SONATA.display.display_mesh import plot_cells

from SONATA.mesh.mesh_byprojection import mesh_by_projecting_nodes_on_BSplineLst
from SONATA.mesh.mesh_core import gen_core_cells
from SONATA.mesh.mesh_utils import first_stage_improvements,second_stage_improvements, determine_a_nodes, equidistant_nodes_on_BSplineLst, sort_and_reassignID
from SONATA.mesh.mesh_intersect import map_mesh_by_intersect_curve2d                        

from SONATA.display.display_utils import export_to_JPEG, export_to_PNG, export_to_PDF, export_to_SVG, export_to_PS, export_to_EnhPS, export_to_TEX, \
                          export_to_BMP,export_to_TIFF, show_coordinate_system, display_SONATA_SegmentLst
                          
from SONATA.fileIO.readinput import read_material_input,section_config
from SONATA.vabs.VABS_interface import VABS_config, export_cells_for_VABS, XSectionalProperties
from SONATA.vabs.strain import Strain
from SONATA.vabs.stress import Stress



#=================================================================================
#                               M A I N 
#=================================================================================
SHOW = True

if SHOW == True:
    #====================INIT DISPLAY:====================================
    display, start_display, add_menu, add_function_to_menu = init_display('wx')
    display.Context.SetDeviationAngle(0.000001)       # 0.001 default. Be careful to scale it to the problem.
    display.Context.SetDeviationCoefficient(0.000001) # 0.001 default. Be careful to scale it to the problem.   

#===================LOAD CROSSSECTION==========================================
print "STATUS:\t Reading Crossection Configuration File"
filename = 'sec_config.input'
Configuration = section_config(filename)

startTime = datetime.now()
#LOAD .pkl data with SegmentLst
filename = 'sec_config.pkl'
with open(filename, 'rb') as handle:
    SegmentLst = pickle.load(handle)
    

#Build wires for each layer and segment
for seg in SegmentLst:
    seg.build_wire()
    for layer in seg.LayerLst:
        layer.build_wire()

#
##%%===================MESH SEGMENT================================================
#Resolution = 300 # Nb of Points on Segment0
#length = get_BSplineLst_length(SegmentLst[0].BSplineLst)
#global_minLen = round(length/Resolution,5)
#
#
#mesh = []
#disco_nodes = []
#k = 0
#for i,layer in enumerate(reversed(SegmentLst[-1].LayerLst)):
#    print 'STATUS: \t Meshing Layer %s' %(i)
#    [R,G,B,T] =  plt.cm.jet(k*50)
#    
#    a_BSplineLst = layer.BSplineLst       
#    b_BSplineLst = trim_BSplineLst(layer.Boundary_BSplineLst, layer.S1, layer.S2, 0, 1)
#    if BSplineLst_Orientation(b_BSplineLst,11) == False:
#        b_BSplineLst = reverse_BSplineLst(b_BSplineLst)  
#     
#    if i==0:
#        a_nodes = equidistant_nodes_on_BSplineLst(a_BSplineLst, True, True, True, minLen = global_minLen, LayerID = layer.ID[0])
#    else: 
#        a_nodes = determine_a_nodes(mesh,a_BSplineLst,global_minLen,layer.ID[0])
#       
#    a_nodes, b_nodes, cells = mesh_by_projecting_nodes_on_BSplineLst(a_BSplineLst,a_nodes,b_BSplineLst,layer.thickness,global_minLen,1e-1)
#    enhanced_cells = first_stage_improvements(cells,b_BSplineLst,global_minLen)
#    for c in enhanced_cells:
#        c.calc_theta_1()
#    enhanced_cells = second_stage_improvements(enhanced_cells,b_BSplineLst,global_minLen)
#
#    for c in enhanced_cells:
#        c.theta_3 = layer.Orientation
#        c.MatID = int(layer.MatID)
#        c.structured = True
#        #print cell,'\t',cell.theta_3,cell.theta_1,cell.MatID,cell.area
#        
#    k = k+1;
#    if k>5:
#        k = 0
#    layer.cells = enhanced_cells
#    mesh.extend(enhanced_cells)
#
#    
##%%===================MESH CORE================================================
#print 'STATUS: \t Meshing Core %s' %(1)
#core_Boundary_BSplineLst = []
#core_Boundary_BSplineLst += trim_BSplineLst(SegmentLst[-1].LayerLst[-1].Boundary_BSplineLst, 0, SegmentLst[-1].LayerLst[-1].S1, 0, 1)  #start und ende der lage
#core_Boundary_BSplineLst += copy_BSplineLst(SegmentLst[-1].LayerLst[-1].BSplineLst)
#core_Boundary_BSplineLst += trim_BSplineLst(SegmentLst[-1].LayerLst[-1].Boundary_BSplineLst, SegmentLst[-1].LayerLst[-1].S2, 1, 0, 1)  #start und ende der lage
#    
#a_nodes = determine_a_nodes(mesh,core_Boundary_BSplineLst,global_minLen,layer.ID[0])
#[c_cells,c_nodes] = gen_core_cells(a_nodes,0.7)
#
#for c in c_cells:
#    c.structured = False
#    c.theta_3 = 0
#    c.MatID = int(SegmentLst[0].CoreMaterial)
#    c.calc_theta_1()
#
#mesh.extend(c_cells)
#mesh,nodes = sort_and_reassignID(mesh)
#
#
##%%=====================SAVE MESH as pickle:===================================
#output_filename = filename.replace('.pkl', '_mesh.pkl')
#with open(output_filename, 'wb') as output:
#    pickle.dump(mesh, output, protocol=pickle.HIGHEST_PROTOCOL)

filename = filename.replace('.pkl', '_mesh.pkl')
with open(filename, 'rb') as handle:
    mesh = pickle.load(handle)
    mesh,nodes = sort_and_reassignID(mesh)
    
# %%====================REVIEW=================================================
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


#%% BALANCE WEIGHT - CUTTING HOLE ALGORITHM====================================

if Configuration.SETUP_BalanceWeight == True:
    print 'STATUS: \t Building Balance Weight'   
    BW = Weight(0,Configuration.BW_XPos,Configuration.BW_YPos,Configuration.BW_Diameter,Configuration.BW_MatID)
    display.DisplayShape(BW.wire)
    
    for c in mesh:
        display.DisplayShape(c.wire,color='BLACK')         

    for x in nodes:
        if x.id == 112:
            display.DisplayShape(x.Pnt2d,color='GREEN')
        elif x.id == 113:
            display.DisplayShape(x.Pnt2d,color='YELLOW')

        
    mesh,boundary_nodes = map_mesh_by_intersect_curve2d(mesh,BW.Curve,BW.wire,display=display)
    
    for n in boundary_nodes:
        print n
        display.DisplayShape(n.Pnt2d,color='WHITE')         

    
    plot_cells(mesh, nodes, 'MatID')   
    triangle_options = 'pa.3'
    [bw_cells,bw_nodes] = gen_core_cells(boundary_nodes,options=triangle_options)
    
    for c in bw_cells:
        c.structured = False
        c.theta_3 = 0
        c.MatID = Configuration.BW_MatID
        c.calc_theta_1()
    
    mesh.extend(bw_cells)
    

mesh,nodes = sort_and_reassignID(mesh)
plot_cells(mesh, nodes, 'MatID')       



## %%===================VABS====================================================
#filename = filename.replace('.pkl', '.vab')
#MaterialLst = read_material_input('mat_database.input')
#VABSsetup = VABS_config(recover_flag=1)
#VABSsetup.F = [0,0,0]    #in Newton
#VABSsetup.M = [100e3,0,0]     #in Newton/mm
#
#print 'STATUS: \t RUNING VABS for Constitutive modeling:'
##EXECUTE VABS:
#if VABSsetup.recover_flag == 1:
#    VABSsetup.recover_flag=0
#    export_cells_for_VABS(mesh,nodes,filename,VABSsetup,MaterialLst)
#    command = 'VABSIII.exe '+ filename
#    stdout = subprocess.check_output(command, shell=True)
#    VABSsetup.recover_flag=1
#    export_cells_for_VABS(mesh,nodes,filename,VABSsetup,MaterialLst)
#    print 'STATUS: \t RUNING VABS for 3D Recovery:'
#    command = 'VABSIII.exe '+ filename
#    stdout = stdout + subprocess.check_output(command, shell=True)
#    
#else:
#    export_cells_for_VABS(mesh,nodes,filename,VABSsetup,MaterialLst)
#    command = 'VABSIII.exe '+ filename
#    stdout = subprocess.check_output(command, shell=True)
#
#stdout = stdout.replace('\r\n\r\n','\n\t   -')
#stdout = stdout.replace('\r\n','\n\t   -')
#stdout = 'STATUS: \t VABS CALCULATIONS COMPLETED: \n\t   -' + stdout
#print stdout 
#
#
##%%===================POST-PROCESSING====================================================
#print 'STATUS: \t POST-PROCESSING:'
#filename_K = filename+'.K'
#BeamProperties = XSectionalProperties(filename_K)
#plot_cells(mesh, nodes, 'MatID',BeamProperties,'NACA0012, 150mm chord',False)
#
#if VABSsetup.recover_flag == 1:
#    BeamProperties.read_all_VABS_Results()
#    #ASSIGN Stress and strains to elements:
#    for i,c in enumerate(mesh):
#        c.strain = Strain(BeamProperties.ELE[i][1:7])
#        c.stress = Stress(BeamProperties.ELE[i][7:13])
#        c.strainM = Strain(BeamProperties.ELE[i][13:19])
#        c.stressM = Stress(BeamProperties.ELE[i][19:25])
#    
#    #ASSIGN Displacement U to nodes:
#    for i,n in enumerate(nodes):
#        n.displacement = BeamProperties.U[i][3:6]
#    
#    plot_cells(mesh, nodes, 'stress.sigma11', BeamProperties,'NACA0012, 150mm chord')


#%%===================CROSS-CHECKING====================================================
#if VABSsetup.recover_flag == 1:
#    Integrator = 0
#    for c in mesh:
#        Integrator += c.stress.sigma11*c.area
#    print Integrator


#%%====================OCC-DISPLAY=================================================
if SHOW == True:
#    for c in mesh:
#        display.DisplayColoredShape(c.wire, 'BLACK')
    #    
    display_SONATA_SegmentLst(display, SegmentLst)      
    ##plot_cells(mesh)
    
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