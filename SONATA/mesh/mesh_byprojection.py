# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 09:58:14 2017

@author: TPflumm
"""
import numpy as np
import math
import itertools

from OCC.gp import gp_Vec2d,gp_Pnt2d,gp_Pnt,gp_Vec
from OCC.Geom2dAPI import Geom2dAPI_ProjectPointOnCurve
from OCC.Display.SimpleGui import init_display


from SONATA.mesh.node import Node
from SONATA.mesh.cell import Cell
from SONATA.display.display_utils import display_custome_shape
from SONATA.topo.BSplineLst_utils import get_BSplineLst_length, find_BSplineLst_coordinate



def mesh_by_projecting_nodes_on_BSplineLst(a_BSplineLst,a_nodes,b_BSplineLst,layer_thickness, tol=1e-2, crit_angle = 95, **kwargs):
    """
    *function to mesh the SONATA topologies by projecting nodes onto the 
    generated BSplineLists.   
    
    :rtype: [list,list,list]
    
    Variables and arguments:
    ----------
    a_BSplineLst & a_nodes: the projection source
    b_BSplineLst & b_nodes: the projection destination curve and its resulting nodes
    layer_thickness & tol: is exactly that and will be used with the tolerance(tol) to determine
                     a distance, in which the resulting projection point has to be.
    crit_angle: is the critical angle to determine a corner if 2 projection points are found.
    display: the kwargs display object can be passed to plot/display within the main OCC3DViewer
                    
    
    Workflow:
    ----------        
      
    
    
    Notes & Comments:
    ----------  
    TODO: * scale distance not only to layerthickenss but also to min_len.
            or adapt the distance individually for each node. 


    """
    
    #KWARGS:
    if kwargs.get('display') !=  None:
        display = kwargs.get('display')
    
    LayerID = 'T_' + a_nodes[0].parameters[0]
    b_nodes = []
    cellLst = []
    distance = (1+tol)*layer_thickness
               
               
    #Is a_BSplineLst closed? 
    closed_a = False
    if a_BSplineLst[0].StartPoint().IsEqual(a_BSplineLst[-1].EndPoint(),1e-5):
        closed_a = True
                       
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
            first =  item.FirstParameter()
            last = item.LastParameter()
            tol = 1e-4
            Umin = first-(last-first)*tol
            Umax = last+(last-first)*tol
            projection = Geom2dAPI_ProjectPointOnCurve(Pnt2d,item.GetHandle(),Umin,Umax)
            #projection = Geom2dAPI_ProjectPointOnCurve(Pnt2d,item.GetHandle())
            
            for j in range(1,projection.NbPoints()+1):
                if projection.Distance(j)<=distance:
                    pPnts.append(projection.Point(j))
                    pPara.append(projection.Parameter(j))
                    pIdx.append(idx)
                else: None   
        
        
        #==================making sure the pPnts are unique:
        '''It happend that somehow the same points were found multiple times
        TODO: wrap into own function'''
        rm_idx=[]
        for a, b in itertools.combinations(enumerate(pPnts), 2):
            if a[1].IsEqual(b[1],1e-6):
                rm_idx.append(a[0])
        
        pPnts = [j for k, j in enumerate(pPnts) if k not in rm_idx]
        pPara = [j for k, j in enumerate(pPara) if k not in rm_idx]
        pIdx = [j for k, j in enumerate(pIdx) if k not in rm_idx]  
    
        #=========if 3 Points are found. Select the 2 points that create the larger angle.
        if len(pPnts) == 3:
            v0 = gp_Vec2d(Pnt2d,pPnts[0])
            v1 = gp_Vec2d(Pnt2d,pPnts[1])
            v2 = gp_Vec2d(Pnt2d,pPnts[2])
            angle01 = abs(v0.Angle(v1)*180/np.pi)
            angle02 = abs(v0.Angle(v2)*180/np.pi)
            angle12 = abs(v1.Angle(v2)*180/np.pi)
            aLst = [angle01,angle02,angle12]
            pPnts.pop(aLst.index(max(aLst)))
         
            
        #==================DETECT CORNERS====================================== 
        if len(pPnts) == 0: 
            print 'ERROR:\t No Projection found for node,', node.id
        elif len(pPnts) == 1: 
            b_nodes.append(Node(pPnts[0],[LayerID,pIdx[0],pPara[0]]))      

        elif len(pPnts) == 2:
            #=======================determine the angle of the potential corner 
            v1 = gp_Vec2d(Pnt2d,pPnts[0])
            v2 = gp_Vec2d(Pnt2d,pPnts[1])
            angle = (180-abs(v1.Angle(v2)*180/np.pi))
            vres = v1.Added(v2)

            if angle < crit_angle:
                node.corner = True
                
            if angle >= crit_angle:
                node.corner = False
            #print 'Node:', node.id, 'corner angle: ', angle 
            
            #=======================determine the corner character - regular_corner type:(Bool)
            node.regular_corner = True
            vp0 = gp_Vec2d() 
            vp1 = gp_Vec2d()
            p_tmp = gp_Pnt2d()
            b_BSplineLst[pIdx[0]].D1(pPara[0],p_tmp,vp0)
            b_BSplineLst[pIdx[1]].D1(pPara[1],p_tmp,vp1)
            dot0 = vres.Dot(vp0)
            dot1 = vres.Dot(vp1)
            
            if dot0<0:
                node.regular_corner = False
                z_BSplineLst = b_BSplineLst[pIdx[1]:] + b_BSplineLst[:pIdx[0]+1]
            
            elif dot0>0:
                node.regular_corner = True
            
            else:
                print 'ERROR: cannot determine regular_corner because vres and v0 are orthogonal'
            
            #=======================determine the exterior corners on z_BsplineLst
            #TODO: DETECT ALL EXTERIOR CORNERS WITHIN THAT INTERVAL pIdx[0],pPara[0],'|',pIdx[1], pPara[1]
            #TODO: IF NO EXTERIOR CORNER IS FOUND USE BISECTOR!
            #TODO: USE node.cornertype = 1,2,3,4 to detemine Shape of Elements
            #       cornertype = 0 no exterior corner, cornertype=1 one exterior corner, ....
            #TODO: THIS ALGORITHM DOESNT FIND ALL 
            exterior_corners = [] 
            exterior_corners_para = []
            
            if node.regular_corner == True:
                for j,item in enumerate(b_BSplineLst[pIdx[0]:pIdx[1]], start=pIdx[0]):
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
    
            #=======================generate b_nodes 
            #print node,'corner: ', node.corner, ', regular_corner = ',node.regular_corner, ',  Len:exterior_corners =',len(exterior_corners)
            if len(exterior_corners) == 0 and node.corner==True:
                node.cornerstyle = 2
                if node.regular_corner == True:
                    b_nodes.append(Node(pPnts[0],[LayerID,pIdx[0],pPara[0]]))
                    newPnt = gp_Pnt2d()
                    newPara = (pPara[0]+pPara[1])/2                    
                    b_BSplineLst[pIdx[0]].D0(newPara,newPnt)
                    b_nodes.append(Node(newPnt,[LayerID,pIdx[0],newPara]))
                    b_nodes.append(Node(pPnts[1],[LayerID,pIdx[1],pPara[1]]))
                
                else:
                    b_nodes.append(Node(pPnts[1],[LayerID,pIdx[1],pPara[1]]))
                    newPara = b_BSplineLst[pIdx[0]].FirstParameter()
                    newPnt = gp_Pnt2d()                  
                    b_BSplineLst[pIdx[0]].D0(newPara,newPnt)
                    b_nodes.append(Node(newPnt,[LayerID,pIdx[0],newPara]))
                    b_nodes.append(Node(pPnts[0],[LayerID,pIdx[0],pPara[0]]))
                    

            elif len(exterior_corners) == 1 and node.corner==True:
                node.cornerstyle = 3
                if node.regular_corner == True:
                    b_nodes.append(Node(pPnts[0],[LayerID,pIdx[0],pPara[0]]))
                    b_nodes.append(Node(exterior_corners[0],[exterior_corners_para[0][0],exterior_corners_para[0][1],exterior_corners_para[0][2]]))
                    b_nodes.append(Node(pPnts[1],[LayerID,pIdx[1],pPara[1]]))
                
                else:
                    b_nodes.append(Node(pPnts[1],[LayerID,pIdx[1],pPara[1]]))
                    b_nodes.append(Node(exterior_corners[0],[exterior_corners_para[0][0],exterior_corners_para[0][1],exterior_corners_para[0][2]]))
                    b_nodes.append(Node(pPnts[0],[LayerID,pIdx[0],pPara[0]]))
                
            #===4======
            elif len(exterior_corners) == 2 and node.corner==True:
                node.cornerstyle = 4
                if node.regular_corner == True:
                    print 'R',[exterior_corners_para[0][0],exterior_corners_para[0][1],exterior_corners_para[0][2]],[exterior_corners_para[1][0],exterior_corners_para[1][1],exterior_corners_para[1][2]]
                    b_nodes.append(Node(pPnts[0],[LayerID,pIdx[0],pPara[0]]))
                    b_nodes.append(Node(exterior_corners[0],[exterior_corners_para[0][0],exterior_corners_para[0][1],exterior_corners_para[0][2]]))
                    
                    newPnt = gp_Pnt2d()
                    newPara = (exterior_corners_para[0][2]+exterior_corners_para[1][2])/2   
                    
                    #TODO: Find Middle between the two exterior corners on b_BsplineLst
                    c_BSplineLst = b_BSplineLst[exterior_corners_para[0][1]+1:exterior_corners_para[1][1]+1]
                    [tmp_idx,tmp_u] =  find_BSplineLst_coordinate(c_BSplineLst,0.5,0,1)
                    newIdx = exterior_corners_para[0][1]+1+tmp_idx
                    newPara = tmp_u
                    b_BSplineLst[newIdx].D0(newPara,newPnt)
                    display.DisplayShape(newPnt,color='RED')
                    
                    
                    b_nodes.append(Node(newPnt,[LayerID,newIdx,newPara]))
                    
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
            
            elif len(exterior_corners) > 2 and node.corner==True:
                node.cornerstyle = 5
                #TODO: 
                  #for p in exterior_corners:
                        #display.DisplayShape(p,color='RED')
                #b_nodes.append(Node(exterior_corners[0],[exterior_corners_para[0][0],exterior_corners_para[0][1],exterior_corners_para[0][2]]))  
                #b_nodes.append(Node(b_BSplineLst[pIdx[0]].EndPoint(),[LayerID,pIdx[0],b_BSplineLst[pIdx[0]].LastParameter()]))
                #b_nodes.append(Node(pPnts[1],[LayerID,pIdx[1],pPara[1]]))    

            if len(exterior_corners) == 0 and node.corner==False:
                node.cornerstyle = 0
                #TODO: a more robust possibilit is to use a bisector vres
                if pIdx[0] == pIdx[1]:
                    newPnt = gp_Pnt2d()
                    newPara = (pPara[0]+pPara[1])/2                    
                    b_BSplineLst[pIdx[0]].D0(newPara,newPnt)
                    b_nodes.append(Node(newPnt,[LayerID,pIdx[0],newPara]))
                elif node.regular_corner == False:
                    print 'ERROR: this possibility has not been implemented yet.'
                
                else:
                    print 'ERROR: this possibility has not been implemented yet. pIdx[0] != pIdx[0]. Create Bisector and intersect with bsplinelist!'                     
                    
                
            if len(exterior_corners) == 1 and node.corner==False:
                node.cornerstyle = 1
                print 'Node ID: ',node.id,', Len(exterior_corners):', len(exterior_corners)
    
                try: 
                    if b_BSplineLst[pIdx[0]].EndPoint().IsEqual(b_nodes[-1].Pnt2d,1e-5):
                       b_nodes.append(Node(pPnts[1],[LayerID,pIdx[1],pPara[1]]))
                      
                    else:
                        b_nodes.append(Node(b_BSplineLst[pIdx[0]].EndPoint(),[LayerID,pIdx[0],b_BSplineLst[pIdx[0]].LastParameter()]))
                
                except:
                    b_nodes.append(Node(b_BSplineLst[pIdx[0]].StartPoint(),[LayerID,pIdx[0],b_BSplineLst[pIdx[0]].FirstParameter()]))
                #display.DisplayShape(b_BSplineLst[pIdx[0]].EndPoint(),color='WHITE')
                #display.DisplayShape(pPnts[0],color='GREEN')
                #display.DisplayShape(pPnts[1],color='ORANGE')    
                

            
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
        
        

    #==============OCC3DVIEWER========================================
    if kwargs.get('display') !=  None:
        for i,a in enumerate(a_nodes):
                if a.corner == True:
                    display.DisplayShape(a.Pnt,color='WHITE')  
                    string = str(a.id)+' (cs='+str(a.cornerstyle)+', rg='+str(a.regular_corner)+')'
                    display.DisplayMessage(a.Pnt,string,message_color=(1.0,0.0,0.0))
                    
                elif a.cornerstyle == 1 or a.cornerstyle == 10 :
                    display.DisplayShape(a.Pnt,color='WHITE')  
                    string = str(a.id)+' (cs='+str(a.cornerstyle)+', rg='+str(a.regular_corner)+')'
                    display.DisplayMessage(a.Pnt,string,message_color=(1.0,0.5,0.0))
                    
                    
                else: 
                    display.DisplayShape(a.Pnt,color='WHITE')  
                    display.DisplayMessage(a.Pnt,str(a.id))
                    
        for i,b in enumerate(b_nodes):
                display.DisplayShape(b.Pnt,color='GREEN')  
                #display.DisplayMessage(b.Pnt,str(b.id),message_color=(1.0,0.5,0.0))
    
    
        for i,a_spline in enumerate(a_BSplineLst):
            #display_custome_shape(display,a_spline,1.0,0.0,[0.2,0.9,0.8])
            display.DisplayShape(a_spline,color='CYAN')
            p = gp_Pnt2d()
            v = gp_Vec2d()
            u = (a_spline.LastParameter()-a_spline.FirstParameter())/2+a_spline.FirstParameter()
            a_spline.D1(u,p,v)
            display.DisplayMessage(p,str(i),height=30,message_color=(0,1,1))
            #display.DisplayVector(gp_Vec(v.X(),v.Y(),0), gp_Pnt(p.X(),p.Y(),0))
            
        for i,b_spline in enumerate(b_BSplineLst):
            #display_custome_shape(display,b_spline,1.0,0.0,[0.1,0.5,1.0 ])
            display.DisplayShape(b_spline,color='BLUE')
            p = gp_Pnt2d()
            v = gp_Vec2d()
            u = (b_spline.LastParameter()-b_spline.FirstParameter())/2+b_spline.FirstParameter()
            b_spline.D1(u,p,v)
            display.DisplayMessage(p,str(i),height=30,message_color=(0,0,1))    
            #display.DisplayVector(gp_Vec(v.X(),v.Y(),0), gp_Pnt(p.X(),p.Y(),0))
            
    return a_nodes, b_nodes, cellLst