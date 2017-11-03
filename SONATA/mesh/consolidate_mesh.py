# -*- coding: utf-8 -*-
"""
Created on Thu Nov 02 10:46:29 2017
@author: TPflumm
"""

import numpy as np

from OCC.gp import gp_Vec2d

from SONATA.mesh.mesh_utils import find_cells_that_contain_node
from SONATA.mesh.cell import Cell
from SONATA.display.display_utils import display_custome_shape    

def consolidate_mesh_on_web(mesh,w_BSplineLst,w1_nodes,w2_nodes,w_tol,display):
    '''Consolidates mesh on the web interface.
    
    After the mesh has been generated for every segment. This function makes 
    sure that no hanging nodes remain in the mesh. 
    
    Args:
        mesh: The overall mesh (list of cells) to consolidate and also to return
        w_BSplineLst: is the definition of the web as a BSplineLst, usually
                    it is a List with one spline that is a straight line.
        w1_nodes: are the nodes of the left segment that lie on the web
        w2_nodes: are the nodes of the right segment that lie on the web 
        w_tol (float) : tolerance to decide whether to match a left node to a right node 
               or if the cell should rather be devided into triangles.
        display: 
            
    Returns: 
        mesh: the updated mesh formulation (list of cells)
                
    '''
    
    #Generate NodeMatching Matrix [NM]
    tmp = []
    for w1 in w1_nodes:
        for w2 in w2_nodes:
            tmp.append([w1.Pnt2d.Distance(w2.Pnt2d),w1.id,w2.id])   
    tmp=np.asarray(tmp)
    NM = tmp[tmp[:,0]<w_tol]
    NM = NM[NM[:,0].argsort()]
    
    #remove possible double nodes in NodeMatching Matrix
    tmp, tmp_idx = np.unique(NM[:,1],True,axis=0) 
    tmp_idx.sort(axis=0)
    NM = NM[tmp_idx]
    tmp, tmp_idx = np.unique(NM[:,2],True,axis=0)
    tmp_idx.sort(axis=0)
    NM = NM[tmp_idx]
    
    #MERGE Nodes according to NodeMatching(NM) Matrix:
    for match in NM:
        #print int(match[1]),int(match[2])
        n1 = filter(lambda x: x.id == match[1], w1_nodes)[0]
        n2 = filter(lambda x: x.id == match[2], w2_nodes)[0]
        #move the first node to the middle!
        n1.Pnt2d
        v = gp_Vec2d(n1.Pnt2d, n2.Pnt2d)
        v.Multiply(0.5)
        n1.Pnt2d.Translate(v)
        
        #find all cells that contain n2 and replace it with n1
        disco_cells = find_cells_that_contain_node(mesh,n2)
        for c in disco_cells:
            c.nodes = [n1 if x==n2 else x for x in c.nodes]
    
    #determine ramaining w1_nodes. If w1_node.id is not in NM[:,1]
    rem_w1_nodes = []
    for w1 in w1_nodes:
        if w1.id not in NM[:,1]:
            rem_w1_nodes.append(w1)
            #display.DisplayShape(w1.Pnt2d, color="RED")   
            
    rem_w2_nodes = []
    for w2 in w2_nodes:
        if w2.id not in NM[:,2]:
            rem_w2_nodes.append(w2)
            #display.DisplayShape(w2.Pnt2d, color="GREEN")   
    
    mesh = split_cells_to_consolidate(mesh,rem_w1_nodes,display)
    mesh = split_cells_to_consolidate(mesh,rem_w2_nodes,display)
    return mesh                    
   

                 
def split_cells_to_consolidate(mesh,rem_nodes,display):
    '''Subfunction to split cells with hanging nodes.
    
    Subfunction of consolidate_mesh_on_web to split cells with hanging nodes.
    It first identifies cells that are beeing intersected by the hanging nodes
    and then splits quad cells into 1 triangle and one quad and a triangle cell
    into two triangles.
    
    Args:
        mesh: The overall mesh (list of cells) to consolidate and also to return
        rem_nodes: remaining hanging nodes, where no match was found
        display: passes OCC display environment
            
    Returns: 
        mesh: the updated mesh formulation (list of cells)
        
    TODO: Determine edge of cell on which the hanging node is and split cells 
          accordingly.
    TODO: Split cells according to edg_idx. edg_idx is the index of the edge of
             the cell, where the new node rn is on
                
    '''
    for rn in rem_nodes:
        #mark all cells that are very close to the node but do not contain it! 
        for c in mesh:
            #identify cells that are beeing intersected by the hanging node!  
            if c.cell_node_distance(rn) < 1e-4 and (rn not in c.nodes):
                #calculate the minimun distance between the wire and the node.Pnt2d
                display_custome_shape(display,c.wire,2.0, 0.0, [1,1,0]) #Yellow
                display.DisplayShape(rn.Pnt2d, color="YELLOW")   
                edg_idx = c.closest_cell_edge(rn)     
                #print rn, c, edg_idx 
                
                #if the cell is a quad split cell into 1 triangles and one quad 
                if len(c.nodes) == 4:
                    
                    #display.DisplayColoredShape(p2[0], 'ORANGE')
                    nodeLst = c.nodes
                    #MODIFY EXISTING CELL
                    c.nodes = [nodeLst[edg_idx-1],nodeLst[edg_idx],rn]
                    #ADD NEW CELL
                    newcell_1 = Cell([nodeLst[edg_idx-1],rn,nodeLst[edg_idx-2]])
                    newcell_1.theta_1 = c.theta_1
                    newcell_1.theta_3 = c.theta_3
                    newcell_1.MatID = c.MatID
                    mesh.append(newcell_1)
                    #ADD NEW CELL
                    newcell_2 = Cell([nodeLst[edg_idx-2],rn,nodeLst[edg_idx-3]])
                    newcell_2.theta_1 = c.theta_1
                    newcell_2.theta_3 = c.theta_3
                    newcell_2.MatID = c.MatID
                    mesh.append(newcell_2)         
                    
                #if cell is triangle split into 2 further triangles
                elif len(c.nodes) == 3: 
                    c.nodes = [nodeLst[edg_idx-1],nodeLst[edg_idx],rn]
                    newcell = Cell([nodeLst[edg_idx-1],rn,nodeLst[edg_idx]])
                    newcell.theta_1 = c.theta_1
                    newcell.theta_3 = c.theta_3
                    newcell.MatID = c.MatID
                    mesh.append(newcell)                 
            
    return mesh             

    