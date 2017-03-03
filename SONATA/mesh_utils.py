import numpy as np
import itertools

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


from BSplineLst_utils import get_BSplineLst_length, get_BSpline_length, trim_BSplineLst, set_BSplineLst_to_Origin, \
                            BSplineLst_Orientation, reverse_BSplineLst, findPnt_on_BSplineLst, copy_BSplineLst
from CADinput import order_BSplineLst_Head2Tail, Check_BSplineLst_Head2Tail
from wire_utils import build_wire_from_BSplineLst,get_wire_length
from utils import Pnt2dLst_to_npArray, unique_rows, PolygonArea, calc_DCT_angles,calc_angle_between
from node import Node
from cell import Cell

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




def mesh_quality_enhancer(cells,b_BSplineLst,global_minLen):
    enhanced_cells = []
    for i,c in enumerate(cells):
        
        
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
            
            
        if c.nodes[0].cornerstyle == 2:
            #display.DisplayShape(c.nodes[0].Pnt2d,color='RED')
            
            v1 = gp_Vec2d(c.nodes[0].Pnt2d,c.nodes[1].Pnt2d)
            v2 = gp_Vec2d(c.nodes[0].Pnt2d,c.nodes[3].Pnt2d)
            angle = (180-abs(v1.Angle(v2)*180/np.pi))
    
            if angle < 60:
                #display.DisplayShape(c.nodes[0].Pnt2d,color='RED')
                L = c.nodes[0].Pnt2d.Distance(c.nodes[2].Pnt2d)
                BS_Vec2d = gp_Vec2d(c.nodes[0].Pnt2d,c.nodes[2].Pnt2d)
                PntLst = []
                for i in range(0,int(L//global_minLen)-1):
                    PntLst.append(c.nodes[0].Pnt2d.Translated(BS_Vec2d.Multiplied((1+i)/float(int(L//global_minLen)))))
                
    #                for P in PntLst:
    #                    display.DisplayShape(P)
                
                pPnts = []
                pPara = []
                pIdx = []
                distance = global_minLen
                for Pnt2d in PntLst: 
                    for idx,item in enumerate(b_BSplineLst):
                        projection = Geom2dAPI_ProjectPointOnCurve(Pnt2d,item.GetHandle())
                        for j in range(1,projection.NbPoints()+1):
                            if projection.Distance(j)<=distance:
                                pPnts.append(projection.Point(j))
                                pPara.append(projection.Parameter(j))
                                pIdx.append(idx)
                            else: None
                Front = []
                Front_parameters = []
                Middle = PntLst
                Back = [] 
                Back_parameters = []
                for i,P in enumerate(pPnts):     
                    if c.nodes[1].Pnt2d.Distance(P)<c.nodes[3].Pnt2d.Distance(P):
                        Front.append(P)
                        Front_parameters.append(([pIdx[i],pPara[i]]))
                    else: 
                        Back.append(P)
                        Back_parameters.append(([pIdx[i],pPara[i]]))
                        
    #                for P in Front:
    #                      display.DisplayShape(P,color='ORANGE')
    #                
    #                for P in Back:
    #                      display.DisplayShape(P,color='GREEN')
                      
                #CREATE FRONT NODES and CELLS
                FrontCellLst = []
                for i in range(0,len(Middle)):
    
                    if i == 0:  #FIRST
                        node0 = c.nodes[0]
                        node1 = c.nodes[1]
                        node2 = Node(Front[i],['',Front_parameters[i][0],Front_parameters[i][1]])
                        node3 = Node(Middle[i])
                        nodeLst = [node0,node1,node2,node3]                    
    
                    else:
                        node0 = Node(Middle[i-1])
                        node1 = Node(Front[i-1],['',Front_parameters[i-1][0],Front_parameters[i-1][1]])
                        node2 = Node(Front[i],['',Front_parameters[i][0],Front_parameters[i][1]])
                        node3 = Node(Middle[i])
                        nodeLst = [node0,node1,node2,node3]    
                
                    FrontCellLst.append(Cell(nodeLst))
                
                #LAST
                if len(Middle)>0:
                    node0 = Node(Middle[-1])
                    node1 = Node(Front[-1],['',Front_parameters[-1][0],Front_parameters[-1][1]])
                    node2 =  c.nodes[2]
                    nodeLst = [node0,node1,node2]   
                    FrontCellLst.append(Cell(nodeLst))
                    
                for fc in FrontCellLst:
                    fc.wire = fc.build_wire()
    #                    display.DisplayShape(fc.wire,color='ORANGE')
                    
                 #CREATE FRONT NODES and CELLS
                BackCellLst = []
                for i in range(0,len(Middle)):
    
                    if i == 0:  #FIRST
                        node0 = c.nodes[0]
                        node1 = c.nodes[3]
                        node2 = Node(Back[i],['',Back_parameters[i][0],Back_parameters[i][1]])
                        node3 = Node(Middle[i])
                        nodeLst = [node0,node1,node2,node3]                    
    
                    else:
                        node0 = Node(Middle[i-1])
                        node1 = Node(Back[i-1],['',Back_parameters[i-1][0],Back_parameters[i-1][1]])
                        node2 = Node(Back[i],['',Back_parameters[i][0],Back_parameters[i][1]])
                        node3 = Node(Middle[i])
                        nodeLst = [node0,node1,node2,node3]    
                
                    BackCellLst.append(Cell(nodeLst))
                
                #LAST
                if len(Middle)>0:
                    node0 = Node(Middle[-1])
                    node1 = Node(Back[-1],['',Back_parameters[-1][0],Back_parameters[-1][1]])
                    node2 =  c.nodes[2]
                    nodeLst = [node0,node1,node2]   
                    BackCellLst.append(Cell(nodeLst))
                    
    #                for bc in BackCellLst:
    #                    bc.wire = bc.build_wire()
    #                    display.DisplayShape(bc.wire,color='GREEN')     
                    
                enhanced_cells.extend(FrontCellLst)
                enhanced_cells.extend(reversed(BackCellLst))
            
            else:
                enhanced_cells.append(c)
    
        else:
            enhanced_cells.append(c)
    
    return enhanced_cells