# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 09:58:14 2017

@author: TPflumm
"""
import numpy as np

from OCC.gp import gp_Vec2d,gp_Pnt2d
from OCC.Geom2dAPI import Geom2dAPI_ProjectPointOnCurve


from SONATA.mesh.node import Node
from SONATA.mesh.cell import Cell



def mesh_by_projecting_nodes_on_BSplineLst(a_BSplineLst,a_nodes,b_BSplineLst,layer_thickness,minLen, tol=8e-3,**kwargs):
    
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
                    print pIdx[0],pPara[0],pIdx[1],pPara[1]
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