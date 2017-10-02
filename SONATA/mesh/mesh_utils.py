import numpy as np
import math
import itertools
import toolz

from OCC.AIS import AIS_Shape
from OCC.BRep import BRep_Builder, BRep_Tool
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeFace
from OCC.BRepAdaptor import BRepAdaptor_CompCurve
from OCC.gp import gp_Pnt2d, gp_Pnt, gp_Vec2d,gp_Lin2d, gp_Dir2d
from OCC.GCPnts import GCPnts_QuasiUniformAbscissa, GCPnts_AbscissaPoint
from OCC.Geom2d import Geom2d_Line
from OCC.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.Geom2dAPI import Geom2dAPI_ProjectPointOnCurve
from OCC.Display.SimpleGui import init_display
from OCC.Quantity import Quantity_Color
from OCC.Graphic3d import Graphic3d_EF_PDF, Graphic3d_EF_SVG, Graphic3d_EF_TEX, Graphic3d_EF_PostScript, Graphic3d_EF_EnhPostScript
from OCC.BRepLib import breplib_BuildCurves3d
from OCC.TopoDS import TopoDS_Compound, topods_Face, topods_Edge


from SONATA.topo.BSplineLst_utils import get_BSplineLst_length, get_BSpline_length, trim_BSplineLst, set_BSplineLst_to_Origin, \
                            BSplineLst_Orientation, reverse_BSplineLst, findPnt_on_BSplineLst, copy_BSplineLst, \
                            isPnt_on_BSplineLst, distance_on_BSplineLst, trim_BSplineLst_by_Pnt2d, trim_BSplineLst_by_coordinates, \
                            ProjectPointOnBSplineLst
from SONATA.fileIO.CADinput import order_BSplineLst_Head2Tail, Check_BSplineLst_Head2Tail
from SONATA.topo.wire_utils import build_wire_from_BSplineLst,get_wire_length
from SONATA.topo.utils import Pnt2dLst_to_npArray, unique_rows, PolygonArea, calc_DCT_angles,calc_angle_between
from SONATA.mesh.node import Node
from SONATA.mesh.cell import Cell




def sort_and_reassignID(mesh):
    #Get all nodes in cells
    temp = []
    for cell in mesh:
        temp.extend(cell.nodes)
        
    nodes = sorted(set(temp), key=lambda Node: (Node.id))
    for i,n in enumerate(nodes):
        n.id = i+1
        
    mesh = sorted(mesh, key=lambda Cell: (Cell.id))    
    for i,c in enumerate(mesh):
        c.id = i+1
    
    return mesh, nodes



def equidistant_nodes_on_BSplineLst(BSplineLst, IC=False, IncStart=True, IncEnd=True, **kwargs): #include corners
    ''' minLen = the minimum distance between points. Should be determined by the segment0.BSplineLstLenght/Resolution
        IC (Include Corners): True or False, Note: IC=True the the Points are not equidistant placed on BSplineLst!!!
        IncStart(Include StartPoint)
        IncEnd(Include EndPoint)
        (BSplineLst, IC=False, either NbPoints or MinLen)
        (BSplineLst, IC=True,  MinLen)
        
    ''' 
    if IC==True:
        closed = False 
        if BSplineLst[0].StartPoint().IsEqual(BSplineLst[-1].EndPoint(),1e-5):
            closed = True
        
        #KWARGS:
        if kwargs.get('minLen') !=  None:
            minLen = kwargs.get('minLen')

        if kwargs.get('LayerID') !=  None:
            LayerID = kwargs.get('LayerID')
        
        nodes = []
        for idx,item in enumerate(BSplineLst):
            Adaptor = Geom2dAdaptor_Curve(item.GetHandle())
            length = get_BSpline_length(item)
            if not math.isnan(length):
                NbPoints = int(length//minLen)+2  
                discretization = GCPnts_QuasiUniformAbscissa(Adaptor,NbPoints)
                
                for j in range(1, NbPoints):
                        para = discretization.Parameter(j)
                        Pnt = gp_Pnt2d()
                        item.D0(para,Pnt)
                        if j==1 and IncStart==False and idx==0:
                            pass
                            
                        else:
                            node = Node(Pnt,[LayerID,idx,para])
                            nodes.append(node)
            
        if closed == False and IncEnd==True: #add last point 
            para = discretization.Parameter(j+1)
            Pnt = gp_Pnt2d()
            item.D0(para,Pnt)
            node = Node(Pnt, [LayerID,idx,para])
            nodes.append(node)       
                    
            

    else:
        wire = build_wire_from_BSplineLst(BSplineLst)
        wire_length = get_wire_length(wire)
        
        #KWARGS:
        if  kwargs.get('NbPoints') !=  None:
            NbPoints = kwargs.get('NbPoints')
            
        elif kwargs.get('minLen') !=  None and kwargs.get('NbPoints') == None:
            NbPoints = int(wire_length//kwargs.get('minLen'))+2  
                          
        AdaptorComp = BRepAdaptor_CompCurve(wire, True)
        discretization = GCPnts_QuasiUniformAbscissa(AdaptorComp,NbPoints)
        
        nodes = []
        for j in range(1, NbPoints+1):
            para = discretization.Parameter(j)
            P = AdaptorComp.Value(para)
            node = Node(gp_Pnt2d(P.X(),P.Y()))
            nodes.append(nodes)
       
    return nodes

def move_node_on_BSplineLst(BSplineLst,node,dist,tol=1e-6):
    CRL = 0    #Cummulative Remaining Length
    direction = True
    idx = node.parameters[1]
    U = node.parameters[2]
    BSplineLst = copy_BSplineLst(BSplineLst)
    RO_BSplineLst = BSplineLst[idx:] + BSplineLst[:idx]

    #Reverse Reordered_BSplineLst if dist<0
    if dist < 0:
        direction = False
        RO_BSplineLst = reverse_BSplineLst(BSplineLst[:idx+1]) + reverse_BSplineLst(BSplineLst[idx+1:])
        dist = abs(dist)

    #Infinitely Loop BSplineLst until P is found
    for j,BSpline in enumerate(itertools.cycle(RO_BSplineLst)):  #Infinitely Loop BSplineLst!
        first = BSpline.FirstParameter()
        last = BSpline.LastParameter()
        Adaptor = Geom2dAdaptor_Curve(BSpline.GetHandle())
        L = GCPnts_AbscissaPoint().Length(Adaptor, first, last, tol) 
        #Add new Spline to CRL
        if j==0:
            if direction == False:
                U = last-(U-first)
            CRL += GCPnts_AbscissaPoint().Length(Adaptor, U, last, tol)
        else:
            CRL += L
        
        #Determin P on BSpline:
        if CRL >= dist and j==0:
            X = GCPnts_AbscissaPoint(Adaptor, dist, U)
            X = X.Parameter()
            P = gp_Pnt2d()
            BSpline.D0(X,P)
            break
            
        elif CRL >= dist and j!=0:
            RD = L - (CRL - dist)
            X = GCPnts_AbscissaPoint(Adaptor, RD, first)
            X = X.Parameter()
            P = gp_Pnt2d()
            BSpline.D0(X,P)
            break
    
    #Replace NodeValues!        
    node.Pnt2d = P
    node.parameters[1] = idx+j
    node.parameters[2] = X

    return None


def grab_nodes_of_cells_on_BSplineLst(cells,BSplineLst):
    disco_nodes = []
    for c in cells:
        disco_nodes.extend(grab_nodes_on_BSplineLst(c.nodes,BSplineLst))
        
    disco_nodes = list(set(disco_nodes))               
    disco_nodes = sorted(disco_nodes, key=lambda Node: (Node.parameters[1],Node.parameters[2]))
    return disco_nodes


def grab_nodes_on_BSplineLst(nodes,BSplineLst):
    disco_nodes = []
    tolerance = 1e-6
    
    for n in nodes:
        for idx,item in enumerate(BSplineLst):
            projection = Geom2dAPI_ProjectPointOnCurve(n.Pnt2d,item.GetHandle())
            for j in range(1,projection.NbPoints()+1):
                if projection.Distance(j) <= tolerance:
                    n.parameters[1] = idx
                    n.parameters[2] = projection.LowerDistanceParameter()
                    disco_nodes.append(n)
                else:
                    None           
                    
    disco_nodes = list(set(disco_nodes))               
    disco_nodes = sorted(disco_nodes, key=lambda Node: (Node.parameters[1],Node.parameters[2]))
    
#    splitter = None
#    for j in range(1,len(disco_nodes)):
#        if disco_nodes[j].parameters[1]-disco_nodes[j-1].parameters[1]>1:
#            splitter = j
#            break 
#        
#    #print 'splitter: ',splitter
#    disco_nodes = disco_nodes[splitter:] + disco_nodes[:splitter] 
    
    return disco_nodes

def determine_a_nodes(mesh,a_BSplineLst,global_minLen,LayerID,factor=5):
    disco_nodes = grab_nodes_of_cells_on_BSplineLst(mesh,a_BSplineLst)
    #determine distance between neighboring nodes and discover the remaining segments to discretize mit equidistant points!
    non_dct_segments = []
    disco_nodes = sorted(disco_nodes, key=lambda Node: (Node.parameters[1],Node.parameters[2]))
    para_start = [0,0]
    para_end = [len(a_BSplineLst)-1, a_BSplineLst[-1].LastParameter()]
    
    #print 'global_minLen:' , global_minLen 
    for j in range(0,len(disco_nodes)+1):
        if j==0:
            d = distance_on_BSplineLst(a_BSplineLst,para_start,disco_nodes[j].parameters[1:])
            if d>(factor*global_minLen):
                non_dct_segments.append([para_start,disco_nodes[j].parameters[1:]])
                #print 'distance on BSplineLst d:', d

        elif j==len(disco_nodes):
#            if closed:
#                d = distance_on_BSplineLst(a_BSplineLst,disco_nodes[j-1].Pnt2d, EndPoint, True) 
#                if d>(factor*global_minLen):
#                    non_dct_segments.append([disco_nodes[j-1].Pnt2d, EndPoint])

            d = distance_on_BSplineLst(a_BSplineLst,disco_nodes[j-1].parameters[1:], para_end)
            if d>(factor*global_minLen):
                non_dct_segments.append([disco_nodes[j-1].parameters[1:], para_end])
                #print 'distance on BSplineLst d:', d
        else:                             
            d = distance_on_BSplineLst(a_BSplineLst,disco_nodes[j-1].parameters[1:], disco_nodes[j].parameters[1:])
            if d>(factor*global_minLen):
                non_dct_segments.append([disco_nodes[j-1].parameters[1:], disco_nodes[j].parameters[1:]])
                #print 'distance on BSplineLst d:', d
                
    #print len(non_dct_segments), " non_dct_segments found"
    tmp_nodes = []
    for seg in non_dct_segments:
        tmp_BSplineLst = trim_BSplineLst_by_coordinates(a_BSplineLst,seg[0],seg[1])
        tmp_nodes.extend(equidistant_nodes_on_BSplineLst(tmp_BSplineLst, True, True, True, minLen = global_minLen, LayerID = LayerID))
    
    #print len(tmp_nodes), "tmp_nodes added"
    #disco_nodes.extend(grab_nodes_on_BSplineLst(tmp_nodes,a_BSplineLst))
    #print 'compare disco_nodes with tmp_nodes and return matches', nf
        
    disco_nodes.extend(grab_nodes_on_BSplineLst(tmp_nodes,a_BSplineLst))    
    disco_nodes = remove_dublicate_nodes(disco_nodes)
    disco_nodes = sorted(disco_nodes, key=lambda Node: (Node.parameters[1],Node.parameters[2]))
    
    return disco_nodes
    


def remove_dublicate_nodes(nodes,tol=1e-6):
    nodes = list(set(nodes))
    doublicated_nodes = [] 
    for i,a in enumerate(nodes):
        for b in nodes[i:]:
            if a.id != b.id and a.Pnt2d.IsEqual(b.Pnt2d,tol):
                doublicated_nodes.append(a)
    
    #print 'doublicated nodes:', doublicated_nodes
    doublicated_nodes = list(set(doublicated_nodes))
    for dn in doublicated_nodes:
        nodes.remove(dn)
    return nodes

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




def modify_sharp_corners(cells,b_BSplineLst,global_minLen,layer_thickness, tol=1e-2,alpha_crit = 50,**kwargs):
     
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
                    #display.DisplayShape(c.nodes[0].Pnt2d,color='RED')
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
                                    pPnts.append(projection.Point(j))
                                    pPara.append(projection.Parameter(j))
                                    pIdx.append(idx)
                                else: None
                        
                        print 'Nuber of Projected Middle Nodes pPnts:', len(pPnts)
                        trigger_f = True
                        trigger_b = True
                        for i,P in enumerate(pPnts):
                                v01 = gp_Vec2d(c.nodes[0].Pnt2d,c.nodes[1].Pnt2d)
                                v03 = gp_Vec2d(c.nodes[0].Pnt2d,c.nodes[3].Pnt2d)
                                vnP = gp_Vec2d(n.Pnt2d,P)

                                if  len(pPnts)>2:
                                    print vnP.Dot(v01)
                                
                                if vnP.Dot(v01)>0 and trigger_f :
                                    trigger_f = False
                                    FrontNodes.append(Node(P,['',pIdx[i],pPara[i]]))
                                
                                elif vnP.Dot(v01)<0 and  trigger_b:
                                    trigger_b = False
                                    BackNodes.append(Node(P,['',pIdx[i],pPara[i]]))
                                
                                else:
                                    print 'ERROR: cannot determine FRONT and BACK nodes because vnp and v01 are orthogonal'
            
                                
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
                    #print 'len(Middle):',len(MiddleNodes),'len(Front):',len(FrontNodes),'len(Back):',len(BackNodes)
                    
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
            
            
#            elif c.nodes[0].cornerstyle == 4 and c.nodes[0].id != trigger_id_cs4:
#                print c.nodes[0].id
#                trigger_id_cs4 = c.nodes[0].id
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
    
    return enhanced_cells



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


def second_stage_improvements(cells,b_BSplineLst,global_minLen,factor1=1.8,factor2=0.15):
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
            if magnitude>=factor1*global_minLen:
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
            elif magnitude<=factor2*global_minLen:
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







def export_cells(cells, filename):
    #Get all nodes in cells
    nodes = [] 
    for cell in cells:
        for node in cell.nodes:
            if node not in nodes:
                nodes.append(node)
                
    nodes = sorted(nodes, key=lambda Node: (Node.id))
    for i,n in enumerate(nodes):
        n.id = i+1
    cells = sorted(cells, key=lambda Cell: (Cell.id))    
    
    
    f = open(filename,'w+')
    f.write('! Number of Nodes, Number of Elements, Number of Materials \n')
    f.write('%i\t%i\t%i\n' % (len(nodes),len(cells),3))
    f.write('\n! Node number, coordinates x_2, coordinatex x_3 \n')
    
    for n in nodes:
        f.write('%i\t\t%f\t%f\n' % (n.id,n.coordinates[0],n.coordinates[1]))
    f.write('\n! Element number, connectivity \n') 

    for c in cells:
        f.write('%i\t\t' % (c.id))
        for i in range(0,9):
            if i<len(c.nodes):
                f.write('%i\t' % (c.nodes[i].id))  
            else:
                f.write('%i\t' % (0))
        f.write('\n')   
    
    f.write('\n! Element number, Layup orientation \n')    
    for c in cells:
        f.write('%i\t\t%i\t%.1f\t' % (c.id,c.MatID,c.theta_3))
        for t in c.theta_1:
            f.write('%.1f\t' % (t))
        f.write('\n')   


