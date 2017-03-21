import os
import numpy as np
import pickle
import matplotlib as plt
from functools import partial
from datetime import datetime

from OCC.gp import gp_Vec2d,gp_Pnt2d
from OCC.Geom2dAPI import Geom2dAPI_ProjectPointOnCurve
from OCC.Display.SimpleGui import init_display

from BSplineLst_utils import get_BSplineLst_length, get_BSpline_length, trim_BSplineLst, set_BSplineLst_to_Origin, \
                            BSplineLst_Orientation, reverse_BSplineLst, findPnt_on_BSplineLst, copy_BSplineLst,\
                            isPnt_on_BSplineLst, distance_on_BSplineLst, trim_BSplineLst_by_Pnt2d, trim_BSplineLst_by_coordinates, \
                            ProjectPointOnBSplineLst


from display_mesh import plot_mesh, plot_cells, plot_nodes
from mesh_utils import mesh_quality_enhancer, grab_nodes_on_BSplineLst, grab_nodes_of_cells_on_BSplineLst, determine_a_nodes, \
                        equidistant_nodes_on_BSplineLst, export_cells
                        
from mesh_core import gen_core_cells
from node import Node
from cell import Cell
from display_utils import export_to_JPEG, export_to_PNG, export_to_PDF, export_to_SVG, export_to_PS, export_to_EnhPS, export_to_TEX, \
                          export_to_BMP,export_to_TIFF, show_coordinate_system, display_SONATA_SegmentLst
from VABS_interface import VABS_config, export_cells_for_VABS
from readinput import read_material_input
#====================NOTES:==========================================================

def corners_of_BSplineLst(BSplineLst):
    corners = [] 
    for item in BSplineLst:
        corners.append(item.EndPoint())
        
    corners.pop(-1)
    return corners #gp_Pnt2d Lst
    

def mesh_by_projecting_nodes_on_BSplineLst(a_BSplineLst,a_nodes,b_BSplineLst,layer_thickness,minLen, tol=8e-3):
    LayerID = 'T_' + a_nodes[0].parameters[0]
    b_nodes = []
    cellLst = []
    distance = (1+tol)*layer_thickness
               
               
    #Is a_BSplineLst closed? 
    closed_a = False
    if a_BSplineLst[0].StartPoint().IsEqual(a_BSplineLst[-1].EndPoint(),1e-5):
        closed_a = True
        
        #print 'closed_a :', closed_a
        #print a_nodes[0],a_nodes[-2]

               
    #==================PROJECT POINTS ON LOWER BOUNDARY =======================            
    if closed_a == True:
        prj_nodes = a_nodes
    else:
        prj_nodes = a_nodes[1:-1]
    
    for i,node in enumerate(prj_nodes, start=1):
        Pnt2d = node.Pnt2d
        pPnts = []
        pPara = []
        pIdx = []
        
        for idx,item in enumerate(b_BSplineLst):
            projection = Geom2dAPI_ProjectPointOnCurve(Pnt2d,item.GetHandle())
            
            for j in range(1,projection.NbPoints()+1):
                if projection.Distance(j)<=distance:
                    pPnts.append(projection.Point(j))
                    pPara.append(projection.Parameter(j))
                    pIdx.append(idx)
                else: None   
            
        #==================DETECT CORNERS====================================== 
        if len(pPnts) == 1: 
            b_nodes.append(Node(pPnts[0],[LayerID,pIdx[0],pPara[0]]))    
            

        elif len(pPnts) == 2:
            
#            display.DisplayShape(node.Pnt2d,color='GREEN')
#            display.DisplayShape(pPnts[0],color='YELLOW')
#            display.DisplayShape(pPnts[1],color='RED')
            
            v1 = gp_Vec2d(Pnt2d,pPnts[0])
            v2 = gp_Vec2d(Pnt2d,pPnts[1])
            angle = (180-abs(v1.Angle(v2)*180/np.pi))
            crit_angle = 100
            
            node.cornerstyle = 1
            
            
            #print 'corner angle: ', angle 
            if angle < crit_angle:
                node.corner = True

#                display.DisplayShape(pPnts[0],color='YELLOW') 
#                display.DisplayShape(pPnts[1],color='RED')
                #=======================
                #print pIdx[0],pPara[0],'|',pIdx[1], pPara[1]
                #TODO: DETECT ALL EXTERIOR CORNERS WITHIN THAT INTERVAL pIdx[0],pPara[0],'|',pIdx[1], pPara[1]
                #TODO: IF NO EXTERIOR CORNER IS FOUND USE BISECTOR!
                #TODO: USE node.cornertype = 1,2,3,4 to detemine Shape of Elements
                #       cornertype = 0 no exterior corner, cornertype=1 one exterior corner, ....
                #TODO: THIS ALGORITHM DOESNT FIND ALL 
                exterior_corners = [] 
                exterior_corners_para = []
                
                #if b_nodes auf dem BSpline dazwischen existieren?
                try:
                    b_nodes[-1]
                    #print b_nodes[-1], b_nodes[-1].parameters
                    x =  b_nodes[-1].parameters[1]*1e5 + b_nodes[-1].parameters[2]
                    a = pIdx[0]*1e5 + pPara[0]
                    b = pIdx[1]*1e5 + pPara[1]
                    if a<x and x<b:
                        z_BSplineLst = b_BSplineLst[pIdx[1]:] + b_BSplineLst[:pIdx[0]+1]
                        regular_corner = False

                    else:
                        regular_corner = True

                except IndexError:
                    #print "IndexError", pIdx[0], pIdx[1], c, d
                    z_BSplineLst = b_BSplineLst[pIdx[1]:] + b_BSplineLst[:pIdx[0]+1]
                    regular_corner = False


                if regular_corner == True:
                    for j,item in enumerate(b_BSplineLst[pIdx[0]:pIdx[1]], start=pIdx[0]):
                        #print "j: ",j
                        spline1 = item
                        spline2 = b_BSplineLst[j+1]
                        u1,p1,v1 = spline1.LastParameter(),gp_Pnt2d(),gp_Vec2d()
                        u2,p2,v2  = spline2.FirstParameter(),gp_Pnt2d(),gp_Vec2d()
                        spline1.D1(u1,p1,v1)
                        spline2.D1(u2,p2,v2)
                        
                        Angle = abs(v1.Angle(v2))*180/np.pi
                        if Angle>0.05:
                            exterior_corners.append(item.EndPoint())
                            exterior_corners_para.append([LayerID,j,u1])  
                            #display.DisplayShape(item.EndPoint(),color='WHITE')

                else: 
                    for j,item in enumerate(z_BSplineLst[:-1]):
                        #print "j: ",j, 'len: ',len(z_BSplineLst) 
                        spline1 = item
                        spline2 = z_BSplineLst[j+1]
                        u1,p1,v1 = spline1.LastParameter(),gp_Pnt2d(),gp_Vec2d()
                        u2,p2,v2  = spline2.FirstParameter(),gp_Pnt2d(),gp_Vec2d()
                        spline1.D1(u1,p1,v1)
                        spline2.D1(u2,p2,v2)
                        
                        Angle = abs(v1.Angle(v2))*180/np.pi
                        if Angle>0.05:
                            exterior_corners.append(item.EndPoint())
                            if len(b_BSplineLst) > j+pIdx[1]:
                                idx = j+pIdx[1]
                            else: 
                                idx = j+pIdx[1]-len(b_BSplineLst)
                            exterior_corners_para.append([LayerID,idx,u1])  
                            #display.DisplayShape(item.EndPoint(),color='WHITE')
    
                
                #print 'Len:exterior_corners =',len(exterior_corners)
                
                if len(exterior_corners) == 0:
                    node.cornerstyle = 2
                    #print pIdx[0],pPara[0],pIdx[1],pPara[1]
                    b_nodes.append(Node(pPnts[0],[LayerID,pIdx[0],pPara[0]]))
                    newPnt = gp_Pnt2d()
                    newPara = (pPara[0]+pPara[1])/2                    
                    b_BSplineLst[pIdx[0]].D0(newPara,newPnt)
                    b_nodes.append(Node(newPnt,[LayerID,pIdx[0],newPara]))
                    b_nodes.append(Node(pPnts[1],[LayerID,pIdx[1],pPara[1]]))
                
                elif len(exterior_corners) == 1:
                    node.cornerstyle = 3
                    if regular_corner == True:
                        b_nodes.append(Node(pPnts[0],[LayerID,pIdx[0],pPara[0]]))
                        b_nodes.append(Node(exterior_corners[0],[exterior_corners_para[0][0],exterior_corners_para[0][1],exterior_corners_para[0][2]]))
                        b_nodes.append(Node(pPnts[1],[LayerID,pIdx[1],pPara[1]]))
                    
                    else:
                        b_nodes.append(Node(pPnts[1],[LayerID,pIdx[1],pPara[1]]))
                        b_nodes.append(Node(exterior_corners[0],[exterior_corners_para[0][0],exterior_corners_para[0][1],exterior_corners_para[0][2]]))
                        b_nodes.append(Node(pPnts[0],[LayerID,pIdx[0],pPara[0]]))
                        
                    #display.DisplayShape(exterior_corners[0],color='WHITE')
                
                elif len(exterior_corners) == 2:
                    node.cornerstyle = 4
                    if regular_corner == True:
                        #print 'R',[exterior_corners_para[0][0],exterior_corners_para[0][1],exterior_corners_para[0][2]],[exterior_corners_para[1][0],exterior_corners_para[1][1],exterior_corners_para[1][2]]
                        b_nodes.append(Node(pPnts[0],[LayerID,pIdx[0],pPara[0]]))
                        b_nodes.append(Node(exterior_corners[0],[exterior_corners_para[0][0],exterior_corners_para[0][1],exterior_corners_para[0][2]]))
                        
                        newPnt = gp_Pnt2d()
                        newPara = (exterior_corners_para[0][2]+exterior_corners_para[1][2])/2                    
                        b_BSplineLst[exterior_corners_para[0][1]].D0(newPara,newPnt)
                        b_nodes.append(Node(newPnt,[LayerID,exterior_corners_para[0][1],newPara]))
                        
                        b_nodes.append(Node(exterior_corners[1],[exterior_corners_para[1][0],exterior_corners_para[1][1],exterior_corners_para[1][2]]))
                        b_nodes.append(Node(pPnts[1],[LayerID,pIdx[1],pPara[1]]))
                        
                    else:
                        #print 'IR',[exterior_corners_para[0][0],exterior_corners_para[0][1],exterior_corners_para[0][2]],[exterior_corners_para[1][0],exterior_corners_para[1][1],exterior_corners_para[1][2]]
                        b_nodes.append(Node(pPnts[1],[LayerID,pIdx[1],pPara[1]]))
                        b_nodes.append(Node(exterior_corners[0],[exterior_corners_para[0][0],exterior_corners_para[0][1],exterior_corners_para[0][2]]))
                        
                        v = gp_Vec2d(exterior_corners[0],exterior_corners[1])
                        p = exterior_corners[0].Translated(v.Multiplied(0.5))
                        p2 = []
                        for idx_2,item_2 in enumerate(b_BSplineLst):
                            projection2 = Geom2dAPI_ProjectPointOnCurve(p,item_2.GetHandle())
                            p2.append([projection2.NearestPoint(),idx_2,projection2.LowerDistanceParameter(),projection2.LowerDistance()])
                        p2 = np.asarray(p2)
                        min_index = p2[:,3].argmin()                       
                        b_nodes.append(Node(p2[min_index,0],[LayerID,p2[min_index,1],p2[min_index,2]]))
                        
                        b_nodes.append(Node(exterior_corners[1],[exterior_corners_para[1][0],exterior_corners_para[1][1],exterior_corners_para[1][2]]))
                        b_nodes.append(Node(pPnts[0],[LayerID,pIdx[0],pPara[0]]))
                
                elif len(exterior_corners) > 2:
                    node.cornerstyle = 5
                    for p in exterior_corners:
                        display.DisplayShape(p,color='RED')
                
                    #b_nodes.append(Node(exterior_corners[0],[exterior_corners_para[0][0],exterior_corners_para[0][1],exterior_corners_para[0][2]]))  
                
                #=======================
                
                #b_nodes.append(Node(b_BSplineLst[pIdx[0]].EndPoint(),[LayerID,pIdx[0],b_BSplineLst[pIdx[0]].LastParameter()]))
                #b_nodes.append(Node(pPnts[1],[LayerID,pIdx[1],pPara[1]]))    

                
                
#                display.DisplayShape(pPnts[0],color='YELLOW')        
#                display.DisplayShape(b_BSplineLst[pIdx[0]].EndPoint(),color='ORANGE')
#                display.DisplayShape(pPnts[1],color='RED')
             
            else:
                if b_BSplineLst[pIdx[0]].EndPoint().IsEqual(b_nodes[-1].Pnt2d,1e-5):
                   b_nodes.append(Node(pPnts[1],[LayerID,pIdx[1],pPara[1]]))
                else:
                    b_nodes.append(Node(b_BSplineLst[pIdx[0]].EndPoint(),[LayerID,pIdx[0],b_BSplineLst[pIdx[0]].LastParameter()]))
                
                #display.DisplayShape(b_BSplineLst[pIdx[0]].EndPoint(),color='WHITE')
                #display.DisplayShape(Pnt2d,color='GREEN')
                
        elif len(pPnts) == 3:
            b_nodes.append(Node(pPnts[1],[LayerID,pIdx[1],pPara[1]]))
            
        else:
            print 'Projection Error, number of projection points: ', len(pPnts)
    
    
    #==============REVERSED PROJECTION=========================================
    leftover_exterior = [] 
    leftover_exterior_para = []
    for i,item in enumerate(b_BSplineLst[:-1]):
        spline1 = item
        spline2 = b_BSplineLst[i+1]
        u1,p1,v1 = spline1.LastParameter(),gp_Pnt2d(),gp_Vec2d()
        u2,p2,v2  = spline2.FirstParameter(),gp_Pnt2d(),gp_Vec2d()
        spline1.D1(u1,p1,v1)
        spline2.D1(u2,p2,v2)
        
        Angle = abs(v1.Angle(v2))*180/np.pi       
        if Angle>0.5:
            leftover_exterior.append(item.EndPoint())
            leftover_exterior_para.append([LayerID,i,u1])  
    
    #find exterior corner Points that are not part of b_nodes
    to_delete = []
    LinearTolerance = 1e-3
    for idx,corn in enumerate(leftover_exterior):
        for node in b_nodes:
            if node.Pnt2d.IsEqual(corn, LinearTolerance):
                to_delete.append(idx)
                break                
    
    for offset,idx in enumerate(to_delete):
        idx -= offset
        del leftover_exterior[idx]
        del leftover_exterior_para[idx]
    
    #print len(leftover_exterior)
    #do the reversed projection! -> the original Pnt2dLst must be modified and be returned as well!
    leftover_interior = []
    leftover_interior_para = []
    for Pnt2d in leftover_exterior:
        pPnts = []
        pIdx = []
        pPara = []
        for idx,item in enumerate(a_BSplineLst):
            projection = Geom2dAPI_ProjectPointOnCurve(Pnt2d,item.GetHandle())
            for i in range(1,projection.NbPoints()+1):
                if projection.Distance(i)<=distance:
                    pPnts.append(projection.Point(i))
                    pPara.append(projection.Parameter(i))
                    pIdx.append(idx)
              
        if len(pPnts) == 1:
            leftover_interior.append(pPnts[0])
            leftover_interior_para.append([a_nodes[0].parameters[0],pIdx[0],pPara[0]]) 


    
    leftover_exterior_nodes = []
    leftover_interior_nodes = []
    #print len(leftover_exterior), len(leftover_interior)
    for i,p in enumerate(leftover_interior):
        leftover_exterior_nodes.append(Node(leftover_exterior[i],leftover_exterior_para[i]))
        leftover_interior_nodes.append(Node(leftover_interior[i],leftover_interior_para[i]))
    

    #======INSERT LEFTOVER NODES ==============================================
#    newlist = a_nodes + leftover_interior_nodes       
#    a_nodes =  sorted(newlist, key=lambda Node: (Node.parameters[1],Node.parameters[2]))
#       
#    newlist = b_nodes + leftover_exterior_nodes       
#    b_nodes =  sorted(newlist, key=lambda Node: (Node.parameters[1],Node.parameters[2]))
    
    #Assosiate a_nodes[0] to b_nodes[0]
#    pNode = []
#    for i,node in enumerate(b_nodes):
#        if a_nodes[0].Pnt2d.Distance(node.Pnt2d)<=distance:
#            pNode.append(node)
#            break
#        
#    b_nodes_start = pNode[0]
    
#    for n in b_nodes:
#        print n

    #display.DisplayShape(b_nodes_start.Pnt2d,color='WHITE')  
    

    #==============CREATE CELLS PROJECTION=========================================
#    for n in b_nodes:
#        display.DisplayShape(n.Pnt2d,color='GREEN')
#    
#    for n in a_nodes:
#        display.DisplayShape(n.Pnt2d,color='ORANGE')  
#     
    #display.DisplayShape(a_nodes[0].Pnt2d,color='RED')  
#    for n in b_nodes:
#        display.DisplayShape(n.Pnt2d,color='YELLOW')  

    #Last Cell as Triangle:
    b = 0   #b_nodes idx
    if closed_a == True:
        start = 0
        end = len(a_nodes)
    else: 
        start = 1
        end = len(a_nodes)-1
    
    
    #for a,node in enumerate(a_nodes[1:-1], start=beginning):
    for a in range(start,end):
        #print 'Closed_a: ', closed_a, ', a: ', a, ', len(a_nodes): ', len(a_nodes),', b: ', b, ', len(b_nodes):', len(b_nodes), '\n',  
        if closed_a == False and a == 1: #Start Triangle
            cellLst.append(Cell([a_nodes[a],a_nodes[a-1],b_nodes[b]]))
        
        elif closed_a == False and a == len(a_nodes)-2: #End Triangle
            cellLst.append(Cell([a_nodes[a-1],b_nodes[b-1],b_nodes[b],a_nodes[a]]))
            cellLst.append(Cell([a_nodes[a],b_nodes[b],a_nodes[a+1]]))

        else: #Regular Cell Creation
            if a_nodes[a].cornerstyle == 2 or a_nodes[a].cornerstyle==3:
                #print a, a_nodes[a], a_nodes[a].cornerstyle
                cellLst.append(Cell([a_nodes[a-1],b_nodes[b-1],b_nodes[b],a_nodes[a]]))
                b += 2
                cellLst.append(Cell([a_nodes[a],b_nodes[b-2],b_nodes[b-1],b_nodes[b]]))
            
            elif a_nodes[a].cornerstyle == 4:
                cellLst.append(Cell([a_nodes[a-1],b_nodes[b-1],b_nodes[b],a_nodes[a]]))
                b += 2
                cellLst.append(Cell([a_nodes[a],b_nodes[b-2],b_nodes[b-1],b_nodes[b]]))
                b += 2
                cellLst.append(Cell([a_nodes[a],b_nodes[b-2],b_nodes[b-1],b_nodes[b]]))
            else:   
                cellLst.append(Cell([a_nodes[a-1],b_nodes[b-1],b_nodes[b],a_nodes[a]]))
        
        b += 1
        
    return a_nodes, b_nodes, cellLst


def theta_1_from_2nodes(node1,node2):
    #calc theta_1_angle for middle Triangle
    theta_1 = [0] * 9
    v0 = gp_Vec2d(gp_Pnt2d(0,0),gp_Pnt2d(1,0))
    v1 = gp_Vec2d(node1.Pnt2d,node2.Pnt2d)
    theta_11 = (v0.Angle(v1))*180/np.pi
    if theta_11<0:
        theta_11 = 360+theta_11
    theta_1[0] = theta_11
    theta_1[1] = 540
    return theta_1


def second_stage_improvements(cells,b_BSplineLst,global_minLen):
    enhanced_cells2 = []
    for i,c in enumerate(cells):
        if len(c.nodes)==4:
            v = gp_Vec2d(c.nodes[1].Pnt2d,c.nodes[2].Pnt2d)
            magnitude = v.Magnitude()
            cP = c.nodes[1].Pnt2d.Translated(v.Multiplied(0.5)) 
            #display.DisplayColoredShape(cP, 'GREEN')  
            p2 = ProjectPointOnBSplineLst(b_BSplineLst,cP,1)
            #display.DisplayColoredShape(p2[0], 'YELLOW')  
            
            #SPLIT CELLS INTO TRIANGLES AND ADD NODE!
            if magnitude>=1.8*global_minLen:
                #display.DisplayColoredShape(p2[0], 'ORANGE')
                nodeLst = c.nodes
                newNode = Node(p2[0],['test',p2[1],p2[2]])
                #MODIFY EXISTING CELL
                c.nodes = [nodeLst[0],nodeLst[1],newNode]
                enhanced_cells2.append(c)
                enhanced_cells2[-1].calc_theta_1()
                #ADD NEW CELLS
                enhanced_cells2.append(Cell([nodeLst[0],newNode,nodeLst[3]]))
                enhanced_cells2[-1].theta_1 = theta_1_from_2nodes(nodeLst[0],nodeLst[3])
                #Append last triangle
                enhanced_cells2.append(Cell([nodeLst[3],newNode,nodeLst[2]]))
                enhanced_cells2[-1].calc_theta_1()
                
            #MERGE NODES when to small
            elif magnitude<=0.15*global_minLen:
                #display.DisplayColoredShape(p2[0], 'RED')
                nodeLst = c.nodes
                #Modify Node 2
                nodeLst[2].Pnt2d = p2[0]
                nodeLst[2].parameters = ['modified',p2[1],p2[2]]
                #MODIFY EXISTING CELL
                c.nodes = [nodeLst[0],nodeLst[2],nodeLst[3]]
                c.theta_1 = theta_1_from_2nodes(nodeLst[0],nodeLst[3])
                enhanced_cells2.append(c)
                
                #MODIFY Last CELL
                cells[i-1].nodes[2] = nodeLst[2]
                
                
            else:
                enhanced_cells2.append(c)

        else:
            enhanced_cells2.append(c) 
    return enhanced_cells2

#=========================================================================
#                   M A I N 
#=========================================================================
SHOW = True

if SHOW == True:
    #====================INIT DISPLAY:==========================================================
    display, start_display, add_menu, add_function_to_menu = init_display('wx')
    display.Context.SetDeviationAngle(0.00001)       # 0.001 default. Be careful to scale it to the problem.
    display.Context.SetDeviationCoefficient(0.00001) # 0.001 default. Be careful to scale it to the problem.   



#===================LOAD CROSSSECTION==========================================
startTime = datetime.now()
#LOAD .pkl data with SegmentLst
filename = 'naca0012_cspar_topo.pkl'
with open(filename, 'rb') as handle:
    SegmentLst = pickle.load(handle)
    

#Build wires for each layer and segment
for seg in SegmentLst:
    seg.build_wire()
    for layer in seg.LayerLst:
        layer.build_wire()



#%%===================MESH SEGMENT================================================
#Resolution = 500 # Nb of Points on Segment0
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
#    enhanced_cells = mesh_quality_enhancer(cells,b_BSplineLst,global_minLen)
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

##%%=====================SAVE MESH as pickle:===================================
#output_filename = filename.replace('_topo.pkl', '_mesh.pkl')
#with open(output_filename, 'wb') as output:
#    pickle.dump(mesh, output, protocol=pickle.HIGHEST_PROTOCOL)

filename = 'naca0012_cspar_mesh.pkl'
with open(filename, 'rb') as handle:
    mesh = pickle.load(handle)
    
print 'STATUS: \t MESHING COMPLETED:'
print '\t   - Total Number of Cells: %s' %(len(mesh))
print '\t   - Duration: %s' % (datetime.now() - startTime)
print '\t   - Saved as: %s' % filename 

# %%====================REVIEW=================================================
minarea = min([c.area for c in mesh])
print '\t   - smallest cell area: %s' % minarea 
minimum_angle = min([c.minimum_angle for c in mesh])
print '\t   - smallest angle [deg]: %s' % minimum_angle 
orientation = all([c.orientation for c in mesh])
print '\t   - Orientation [CC]: %s' % orientation 

plot_cells(mesh,'theta_3')

# %%===================VABS====================================================
filename = filename.replace('.pkl', '_vabs.dat')
MaterialLst = read_material_input('mat_database.input')
VABSsetup = VABS_config(Timoshenko_flag=0)
export_cells_for_VABS(mesh,filename,VABSsetup,MaterialLst)
command = 'VABSIII.exe '+ filename
os.system(command)



#%%====================DISPLAY=================================================
if SHOW == True:
    for c in mesh:
        display.DisplayColoredShape(c.wire, 'BLACK')
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