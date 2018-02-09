# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 10:30:01 2017

@author: TPflumm
"""

import numpy as np

from OCC.gp import gp_Vec2d, gp_Pnt2d
from OCC.Geom2dAPI import Geom2dAPI_ProjectPointOnCurve
from OCC.Geom2dAPI import Geom2dAPI_PointsToBSpline

from SONATA.cbm.mesh.node import Node
from SONATA.cbm.mesh.cell import Cell
from SONATA.cbm.mesh.mesh_utils import move_node_on_BSplineLst, theta_1_from_2nodes
from SONATA.cbm.topo.BSplineLst_utils import ProjectPointOnBSplineLst

def modify_sharp_corners(cells, b_BSplineLst, global_minLen, layer_thickness, LayerID = 0, tol=1e-2,alpha_crit = 50,**kwargs):
    '''
    
    '''
    
     
    #KWARGS:
    if kwargs.get('display') !=  None:
        display = kwargs.get('display')
        
    enhanced_cells = []
    for i,c in enumerate(cells):
        if len(c.nodes) == 4:  
            
            #cs4_counter = 0            
            if c.nodes[0].cornerstyle == 2 or c.nodes[0].cornerstyle == 3:
                #display.DisplayShape(c.nodes[0].Pnt2d,color='RED')
                
                #v1 = gp_Vec2d(c.nodes[0].Pnt2d,c.nodes[1].Pnt2d)
                #v2 = gp_Vec2d(c.nodes[0].Pnt2d,c.nodes[3].Pnt2d)
                #angle = (180-abs(v1.Angle(v2)*180/np.pi))
                #if v2.Magnitude() == 0:
                    #print c.nodes[0].coordinates, c.nodes[3].coordinates
                
                v21 = gp_Vec2d(c.nodes[2].Pnt2d,c.nodes[1].Pnt2d)
                v23 = gp_Vec2d(c.nodes[2].Pnt2d,c.nodes[3].Pnt2d)
                angle = abs(v21.Angle(v23)*180/np.pi)
                
                #print c,angle
                if angle < alpha_crit:
#                    display.DisplayShape(c.nodes[0].Pnt2d,color='RED')
#                    display.DisplayShape(c.nodes[2].Pnt2d,color='GREEN')
                    L = c.nodes[0].Pnt2d.Distance(c.nodes[2].Pnt2d)*1.5
                    BS_Vec2d = gp_Vec2d(c.nodes[0].Pnt2d,c.nodes[2].Pnt2d)
                    MiddleNodes = []
                    for i in range(0,int(L//global_minLen)-1):
                        P = c.nodes[0].Pnt2d.Translated(BS_Vec2d.Multiplied((1+i)/float(int(L//global_minLen))))
                        MiddleNodes.append(Node(P))
                        #display.DisplayShape(P)
                    

                    FrontNodes = []
                    BackNodes= [] 
                    distance = (1+tol)*layer_thickness
                    for n in MiddleNodes:
                        pPnts = []
                        pPara = []
                        pIdx = []
                        for idx,item in enumerate(b_BSplineLst):
                            projection = Geom2dAPI_ProjectPointOnCurve(n.Pnt2d,item.GetHandle())
                            for j in range(1,projection.NbPoints()+1):
                                if projection.Distance(j)<=distance:
                                    if not any(item.IsEqual(projection.Point(j),1e-6) for item in pPnts):
                                        pPnts.append(projection.Point(j))
                                        pPara.append(projection.Parameter(j))
                                        pIdx.append(idx)
                                else: None
                        

                        
                        #print 'Nuber of Projected Middle Nodes pPnts:', len(pPnts)                        
                        trigger_f = True
                        trigger_b = True
                        for i,P in enumerate(pPnts):
                                v01 = gp_Vec2d(c.nodes[0].Pnt2d,c.nodes[1].Pnt2d)
                                #v03 = gp_Vec2d(c.nodes[0].Pnt2d,c.nodes[3].Pnt2d)
                                vnP = gp_Vec2d(n.Pnt2d,P)

                                if  len(pPnts)>2:
                                    print 'vnP.Dot(v01): ', vnP.Dot(v01)
                                    
                                    
                                if vnP.Dot(v01)>0 and trigger_f :
                                    trigger_f = False
                                    FrontNodes.append(Node(P,[LayerID, pIdx[i],pPara[i]]))
                                
                                elif vnP.Dot(v01)<0 and  trigger_b:
                                    trigger_b = False
                                    BackNodes.append(Node(P,[LayerID, pIdx[i],pPara[i]]))
                                
                                else:
                                    print 'ERROR: cannot determine FRONT and BACK nodes @ ',c.nodes[0], 'because vnp and v01 are orthogonal'
                                    print pPara
#                            if c.nodes[1].Pnt2d.Distance(P)<c.nodes[3].Pnt2d.Distance(P):
#                                FrontNodes.append(Node(P,['',pIdx[i],pPara[i]]))
#                            else: 
#                                BackNodes.append(Node(P,['',pIdx[i],pPara[i]]))
                         
#                        for mn in MiddleNodes:
#                              display.DisplayShape(mn.Pnt2d,color='BLUE')
#                              display.DisplayMessage(mn.Pnt,str(mn.id),message_color=(1,1,1))    
#                          
#                        for fn in FrontNodes:
#                              display.DisplayShape(fn.Pnt2d,color='ORANGE')
#                              display.DisplayMessage(fn.Pnt,str(fn.id),message_color=(1,1,1))
#                        
#                        for bn in BackNodes:
#                              display.DisplayShape(bn.Pnt2d,color='RED')
#                              display.DisplayMessage(bn.Pnt,str(bn.id),message_color=(1,1,1))
                              
                          
                    #=====================CREATE FRONT CELLS
                    FrontCellLst = []
                    #print '@', c.nodes[0],'  len(Middle):',len(MiddleNodes),'len(Front):',len(FrontNodes),'len(Back):',len(BackNodes)
                    
                    for i in range(0,len(MiddleNodes)):
        
                        if i == 0:  #FIRST
                            nodeLst = [c.nodes[0],c.nodes[1],FrontNodes[i],MiddleNodes[i]]                    
                        else:
                            nodeLst = [MiddleNodes[i-1],FrontNodes[i-1],FrontNodes[i],MiddleNodes[i]]    
                        FrontCellLst.append(Cell(nodeLst))
                    
                    if len(MiddleNodes)>0:   #LAST
                        nodeLst = [MiddleNodes[-1],FrontNodes[-1],c.nodes[2]]   
                        FrontCellLst.append(Cell(nodeLst))
                                               
                    #=====================CREATE BACK CELLS
                    BackCellLst = []
                    for i in range(0,len(MiddleNodes)):
        
                        if i == 0:  #FIRST
                            nodeLst = [MiddleNodes[i],BackNodes[i],c.nodes[3],c.nodes[0]]                                  
                        else:
                            nodeLst = [MiddleNodes[i],BackNodes[i],BackNodes[i-1],MiddleNodes[i-1]]    
                        BackCellLst.append(Cell(nodeLst))
                    
                    if len(MiddleNodes)>0: #LAST
                        nodeLst = [MiddleNodes[-1],c.nodes[2],BackNodes[-1]]  
                        BackCellLst.append(Cell(nodeLst))
                        
                    enhanced_cells.extend(FrontCellLst)
                    enhanced_cells.extend(reversed(BackCellLst))

                    if len(MiddleNodes) == 0:
                        enhanced_cells.append(c)
                else:
                    enhanced_cells.append(c)
            
            
           #elif c.nodes[0].cornerstyle == 4 and c.nodes[0].id != trigger_id_cs4:
               #print c.nodes[0].id
#              trigger_id_cs4 = c.nodes[0].id
#                
#                display.DisplayShape(c.nodes[0].Pnt2d,color='RED')
#                display.DisplayShape(c.nodes[2].Pnt2d,color='ORANGE')
#                display.DisplayShape(c.nodes[3].Pnt2d,color='YELLOW')
#                L = c.nodes[0].Pnt2d.Distance(c.nodes[2].Pnt2d)*1.5
#                BS_Vec2d = gp_Vec2d(c.nodes[0].Pnt2d,c.nodes[2].Pnt2d)
#                MiddleNodes = []
#                for i in range(0,int(L//global_minLen)-1):
#                    P = c.nodes[0].Pnt2d.Translated(BS_Vec2d.Multiplied((1+i)/float(int(L//global_minLen))))
#                    MiddleNodes.append(Node(P))
#        #                for P in PntLst:
#        #                    display.DisplayShape(P)
                
                
                
                #enhanced_cells.append(c)
 
            else:
                enhanced_cells.append(c)
        else:
            enhanced_cells.append(c)
    try:
        new_b_nodes = FrontNodes+BackNodes
    except:
        new_b_nodes = []
        
    return enhanced_cells,new_b_nodes





def modify_cornerstyle_one(cells,b_BSplineLst,**kwargs):
    
    #KWARGS:
    if kwargs.get('display') !=  None:
        display = kwargs.get('display')
    
    enhanced_cells = []
    for i,c in enumerate(cells):
        if len(c.nodes) == 4:
            if c.nodes[0].cornerstyle == 1:           
                #TODO: Wrap into single function that has the variable of affecting neighboring range
                #=====Neighbor to the Right==============
                #determine distance between node 0 and 3 and move node 2 by 2/3(y-x) closer to node 1
                x = c.nodes[0].Pnt2d.Distance(c.nodes[3].Pnt2d)
                y = c.nodes[1].Pnt2d.Distance(c.nodes[2].Pnt2d)
                delta = 2/float(3)*(y-x)
                move_node_on_BSplineLst(b_BSplineLst,c.nodes[2],-delta)             
        
                #=====Neighbor to the Right +1 ==============
                x = cells[i+1].nodes[0].Pnt2d.Distance(cells[i+1].nodes[3].Pnt2d)
                y = cells[i+1].nodes[1].Pnt2d.Distance(cells[i+1].nodes[2].Pnt2d)
                delta = 1/float(3)*(y-x)
                move_node_on_BSplineLst(b_BSplineLst,cells[i+1].nodes[2],-delta)   
                
            if c.nodes[3].cornerstyle == 1:
                #=====Neighbor to the LEFT ==============
                x = c.nodes[0].Pnt2d.Distance(c.nodes[3].Pnt2d)
                y = c.nodes[1].Pnt2d.Distance(c.nodes[2].Pnt2d)
                delta = 2/float(3)*(y-x)
                move_node_on_BSplineLst(b_BSplineLst,c.nodes[1],delta)   
        
                #=====Neighbor to the LEFT +1 ==============
                x = cells[i-1].nodes[0].Pnt2d.Distance(cells[i-1].nodes[3].Pnt2d)
                y = cells[i-1].nodes[1].Pnt2d.Distance(cells[i-1].nodes[2].Pnt2d)
                delta = 1/float(3)*(y-x)         
                move_node_on_BSplineLst(b_BSplineLst,cells[i-1].nodes[1],delta)   
                
            else:
                enhanced_cells.append(c)
        else:
            enhanced_cells.append(c)
    
    return enhanced_cells


def second_stage_improvements(cells,b_BSplineLst,global_minLen, LayerID = 0,factor1=1.8,factor2=0.15, **kw):
    
    if kw.get('display') !=  None:
        display = kw.get('display')

    enhanced_cells2 = []
    new_b_nodes = []
    for i,c in enumerate(cells):
        if len(c.nodes)==4:
            v = gp_Vec2d(c.nodes[1].Pnt2d,c.nodes[2].Pnt2d)
            magnitude = v.Magnitude()

            #SPLIT CELLS INTO TRIANGLES AND ADD NODE!
            if magnitude>=factor1*global_minLen:
                cP = c.nodes[1].Pnt2d.Translated(v.Multiplied(0.5)) 
                #display.DisplayColoredShape(cP, 'GREEN')  
                p2 = ProjectPointOnBSplineLst(b_BSplineLst,cP,1)
                #display.DisplayColoredShape(p2[0], 'YELLOW')  
                nodeLst = c.nodes
                newNode = Node(p2[0],[LayerID,p2[1],p2[2]])
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
            elif magnitude<=factor2*global_minLen and c.nodes[1].corner == False and c.nodes[2].corner == False:
                cP = c.nodes[1].Pnt2d.Translated(v.Multiplied(0.5)) 
                #display.DisplayColoredShape(cP, 'GREEN')  
                p2 = ProjectPointOnBSplineLst(b_BSplineLst,cP,1)
                #display.DisplayColoredShape(p2[0], 'RED')
                nodeLst = c.nodes
                #Modify Node 2
                nodeLst[2].Pnt2d = p2[0]
                nodeLst[2].parameters = [LayerID,p2[1],p2[2]]
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
        
        try:
            new_b_nodes.append(newNode)
        except: 
            None
            
    return enhanced_cells2, new_b_nodes




def integrate_leftover_interior_nodes(a_nodes, a_BSplineLst, b_nodes, b_BSplineLst, cells, thickness, LayerID, **kw):
    
    if kw.get('display') !=  None:
        display = kw.get('display')
    
    aglTol = 0.5
    linTol = 1e-6
    prjTol = 1e-2
    celTol = 0.1
    leftover_exterior_corners = []
    for i,item in enumerate(b_BSplineLst[:-1]):
        spline1 = item
        spline2 = b_BSplineLst[i+1]
        u1,p1,v1 = spline1.LastParameter(),gp_Pnt2d(),gp_Vec2d()
        u2,p2,v2  = spline2.FirstParameter(),gp_Pnt2d(),gp_Vec2d()
        spline1.D1(u1,p1,v1)
        spline2.D1(u2,p2,v2)
        
        Angle = abs(v1.Angle(v2))*180/np.pi       
        if Angle>aglTol and not any(n.Pnt2d.IsEqual(item.EndPoint(), linTol) for n in b_nodes):
            leftover_exterior_corners.append((item.EndPoint(),[LayerID,i,u1]))
                
    #reversed projection:
    for p1 in leftover_exterior_corners:
        new_b_node = Node(p1[0], p1[1])
        p2 = ProjectPointOnBSplineLst(a_BSplineLst, p1[0], (1+prjTol)*thickness)
        
        if len(p2) > 0:
            for c in cells:
                if c.cell_node_distance(new_b_node) < (1+celTol)*thickness:
                    print len(leftover_exterior_corners)
                    display.DisplayShape(p2[0], color='YELLOW')
                    display.DisplayShape(c.wire, color='ORANGE')
                    

    if kw.get('display') !=  None:
        for p in leftover_exterior_corners:
            display.DisplayShape(p[0], color='RED')



    return None