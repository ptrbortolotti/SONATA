# -*- coding: utf-8 -*-
"""
Created on Thu Nov 02 10:46:29 2017

@author: TPflumm
"""

import numpy as np

from OCC.gp import gp_Vec2d

from SONATA.mesh.mesh_utils import find_cells_that_contain_node
from SONATA.mesh.cell import Cell

def consolidate_mesh_on_web(mesh,w_BSplineLst,w1_nodes,w2_nodes,w_tol,display):
    '''Consolidates mesh on the web interface'''
       
    
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
            display.DisplayShape(w1.Pnt2d, color="RED")   
            
    rem_w2_nodes = []
    for w2 in w2_nodes:
        if w2.id not in NM[:,2]:
            rem_w2_nodes.append(w2)
            display.DisplayShape(w2.Pnt2d, color="GREEN")   
        
    #identify cells that are beeing intersected by the hanging node!
   
    for rw1 in rem_w1_nodes:
        #mark all cells that are very close to the node but do not contain it! 
        for c in mesh:
            if c.cell_node_distance(rw1) < 1e-4 and (rw1 not in c.nodes):
                #calculate the minimun distance between the wire and the node.Pnt2d
                #display_custome_shape(display,c.wire,2.0, 0.0, [1,0.5,0])
                #print rw1, c
                        
                #if the cell is a quad split cell into 2 triangles
                if len(c.nodes) == 4:
                    #display.DisplayColoredShape(p2[0], 'ORANGE')
                    nodeLst = c.nodes
                    newNode = rw1
                    #MODIFY EXISTING CELL
                    c.nodes = [nodeLst[0],nodeLst[1],rw1]
                    #ADD NEW CELL
                    newcell = Cell([nodeLst[0],newNode,nodeLst[2],nodeLst[3]])
                    newcell.theta_1 = c.theta_1
                    newcell.theta_3 = c.theta_3
                    newcell.MatID = c.MatID
                    mesh.append(newcell)
                #if cell is triangle split into 2 further triangles
                elif len(c.nodes) == 3: 
                    c.nodes = [nodeLst[0],nodeLst[1],rw1]
                    newcell = Cell([nodeLst[0],newNode,nodeLst[2]])
                    newcell.theta_1 = c.theta_1
                    newcell.theta_3 = c.theta_3
                    newcell.MatID = c.MatID
                    mesh.append(newcell)
                    
    for rw2 in rem_w2_nodes:
        #mark all cells that are very close to the node but do not contain it! 
        for c in mesh:
            if c.cell_node_distance(rw2) < 1e-4 and (rw2 not in c.nodes):
                #calculate the minimun distance between the wire and the node.Pnt2d
                #display_custome_shape(display,c.wire,2.0, 0.0, [1,1,0])
                #print rw2, c     
        
                #if the cell is a quad split cell into 2 triangles
                if len(c.nodes) == 4:
                    #display.DisplayColoredShape(p2[0], 'ORANGE')
                    nodeLst = c.nodes
                    newNode = rw2
                    #MODIFY EXISTING CELL
                    c.nodes = [nodeLst[0],nodeLst[1],rw2]
                    #ADD NEW CELL
                    newcell = Cell([nodeLst[0],newNode,nodeLst[2],nodeLst[3]])
                    newcell.theta_1 = c.theta_1
                    newcell.theta_3 = c.theta_3
                    newcell.MatID = c.MatID
                    mesh.append(newcell)
                #if cell is triangle split into 2 further triangles
                elif len(c.nodes) == 3: 
                    c.nodes = [nodeLst[0],nodeLst[1],rw2]
                    newcell = Cell([nodeLst[0],newNode,nodeLst[2]])
                    newcell.theta_1 = c.theta_1
                    newcell.theta_3 = c.theta_3
                    newcell.MatID = c.MatID
                    mesh.append(newcell)
    
    return mesh