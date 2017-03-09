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
                            BSplineLst_Orientation, reverse_BSplineLst, findPnt_on_BSplineLst, copy_BSplineLst, \
                            isPnt_on_BSplineLst, distance_on_BSplineLst, trim_BSplineLst_by_Pnt2d, trim_BSplineLst_by_coordinates
from CADinput import order_BSplineLst_Head2Tail, Check_BSplineLst_Head2Tail
from wire_utils import build_wire_from_BSplineLst,get_wire_length
from utils import Pnt2dLst_to_npArray, unique_rows, PolygonArea, calc_DCT_angles,calc_angle_between
from node import Node
from cell import Cell

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
            NbPoints = int(length//minLen)+2  
            discretization = GCPnts_QuasiUniformAbscissa(Adaptor,NbPoints)
            
            for j in range(1, NbPoints):
                    para = discretization.Parameter(j)
                    Pnt = gp_Pnt2d()
                    item.D0(para,Pnt)
                    if j==1:
                       if IncStart==False and idx==0:
                           pass
                       else: 
                           node = Node(Pnt,[LayerID,idx,para])
                           nodes.append(node)
                        
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

def determine_a_nodes(mesh,a_BSplineLst,global_minLen,LayerID):
    print "grab nodes of mesh"
    disco_nodes = grab_nodes_of_cells_on_BSplineLst(mesh,a_BSplineLst)
    print "dicover remaining segments"
    #determine distance between neighboring nodes and discover the remaining segments to discretize mit equidistant points!
    non_dct_segments = []
    disco_nodes = sorted(disco_nodes, key=lambda Node: (Node.parameters[1],Node.parameters[2]))
    para_start = [0,0]
    para_end = [len(a_BSplineLst)-1, a_BSplineLst[-1].LastParameter()]
    
    factor = 30
    for j in range(0,len(disco_nodes)+1):
        if j==0:
            d = distance_on_BSplineLst(a_BSplineLst,para_start,disco_nodes[j].parameters[1:])
            if d>(factor*global_minLen):
                non_dct_segments.append([para_start,disco_nodes[j].parameters[1:]])

        elif j==len(disco_nodes):
#            if closed:
#                d = distance_on_BSplineLst(a_BSplineLst,disco_nodes[j-1].Pnt2d, EndPoint, True) 
#                if d>(factor*global_minLen):
#                    non_dct_segments.append([disco_nodes[j-1].Pnt2d, EndPoint])

            d = distance_on_BSplineLst(a_BSplineLst,disco_nodes[j-1].parameters[1:], para_end)
            if d>(factor*global_minLen):
                non_dct_segments.append([disco_nodes[j-1].parameters[1:], para_end])
        else:                             
            d = distance_on_BSplineLst(a_BSplineLst,disco_nodes[j-1].parameters[1:], disco_nodes[j].parameters[1:])
            if d>(factor*global_minLen):
                non_dct_segments.append([disco_nodes[j-1].parameters[1:], disco_nodes[j].parameters[1:]])
        
    
    #print len(non_dct_segments), " non_dct_segments found"
    tmp_nodes = []
    for seg in non_dct_segments:
        tmp_BSplineLst = trim_BSplineLst_by_coordinates(a_BSplineLst,seg[0],seg[1])
        tmp_nodes.extend(equidistant_nodes_on_BSplineLst(tmp_BSplineLst, True, False, False, minLen = global_minLen, LayerID = LayerID))
    
    #print len(tmp_nodes), "tmp_nodes added"
    disco_nodes.extend(grab_nodes_on_BSplineLst(tmp_nodes,a_BSplineLst))
    disco_nodes = sorted(disco_nodes, key=lambda Node: (Node.parameters[1],Node.parameters[2]))
    return disco_nodes


def mesh_quality_enhancer(cells,b_BSplineLst,global_minLen):
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
                
                
            if c.nodes[0].cornerstyle == 2:
                #display.DisplayShape(c.nodes[0].Pnt2d,color='RED')
                
                v1 = gp_Vec2d(c.nodes[0].Pnt2d,c.nodes[1].Pnt2d)
                v2 = gp_Vec2d(c.nodes[0].Pnt2d,c.nodes[3].Pnt2d)
                angle = (180-abs(v1.Angle(v2)*180/np.pi))
        
                if angle < 60:
                    #display.DisplayShape(c.nodes[0].Pnt2d,color='RED')
                    L = c.nodes[0].Pnt2d.Distance(c.nodes[2].Pnt2d)*1.2
                    BS_Vec2d = gp_Vec2d(c.nodes[0].Pnt2d,c.nodes[2].Pnt2d)
                    MiddleNodes = []
                    for i in range(0,int(L//global_minLen)-1):
                        P = c.nodes[0].Pnt2d.Translated(BS_Vec2d.Multiplied((1+i)/float(int(L//global_minLen))))
                        MiddleNodes.append(Node(P))
        #                for P in PntLst:
        #                    display.DisplayShape(P)
                    
                    pPnts = []
                    pPara = []
                    pIdx = []
                    distance = global_minLen
                    for n in MiddleNodes: 
                        for idx,item in enumerate(b_BSplineLst):
                            projection = Geom2dAPI_ProjectPointOnCurve(n.Pnt2d,item.GetHandle())
                            for j in range(1,projection.NbPoints()+1):
                                if projection.Distance(j)<=distance:
                                    pPnts.append(projection.Point(j))
                                    pPara.append(projection.Parameter(j))
                                    pIdx.append(idx)
                                else: None
                    
                    
                    FrontNodes = []
                    BackNodes= [] 
                    for i,P in enumerate(pPnts):     
                        if c.nodes[1].Pnt2d.Distance(P)<c.nodes[3].Pnt2d.Distance(P):
                            FrontNodes.append(Node(P,['',pIdx[i],pPara[i]]))
                        else: 
                            BackNodes.append(Node(P,['',pIdx[i],pPara[i]]))
                            
        #                for P in Front:
        #                      display.DisplayShape(P,color='ORANGE')
        #                
        #                for P in Back:
        #                      display.DisplayShape(P,color='GREEN')
                          
                    #=====================CREATE FRONT NODES and CELLS
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
                        
#                    for fc in FrontCellLst:
#                        fc.wire = fc.build_wire()
#                        display.DisplayShape(fc.wire,color='ORANGE')
                        
                    #=====================CREATE BACK NODES and CELLS
                    BackCellLst = []
                    for i in range(0,len(MiddleNodes)):
        
                        if i == 0:  #FIRST
                            nodeLst = [c.nodes[0],MiddleNodes[i],BackNodes[i],c.nodes[3]]                           
                        else:
                            nodeLst = [MiddleNodes[i],BackNodes[i],BackNodes[i-1],MiddleNodes[i-1]]    
                        BackCellLst.append(Cell(nodeLst))
                    
                    if len(MiddleNodes)>0: #LAST
                        nodeLst = [c.nodes[2],BackNodes[-1],MiddleNodes[-1]]   
                        BackCellLst.append(Cell(nodeLst))
                        
        #                for bc in BackCellLst:
        #                    bc.wire = bc.build_wire()
        #                    display.DisplayShape(bc.wire,color='GREEN')     
                        
                    enhanced_cells.extend(FrontCellLst)
                    enhanced_cells.extend(reversed(BackCellLst))

                    if len(MiddleNodes) == 0:
                        enhanced_cells.append(c)
                else:
                    enhanced_cells.append(c)
        
            else:
                enhanced_cells.append(c)
        else:
            enhanced_cells.append(c)
    
    return enhanced_cells




def export_cells_to_patran(filename, cells):
    #Get all nodes in cells
    nodes = [] 
    for cell in cells:
        for node in cell.nodes:
            if node not in nodes:
                nodes.append(node)
    nodes = sorted(nodes, key=lambda Node: (Node.id))
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
