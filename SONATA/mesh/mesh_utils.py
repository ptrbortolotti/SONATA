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


def find_cells_that_contain_node(cells,n2find):
    #search cells that contain the node n2find
    disco_cells = []
    for c in cells: 
        if n2find in c.nodes:
            disco_cells.append(c)
    return disco_cells


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
    '''the grab_nodes_of_cells_on_BSplineLst fuction determines the nodes of 
    the cells that are located on the BSplineLst (list of geom2d_BSpline 
    objects) with a given tolerance and uses the subfunction 
    grab_nodes_on_BSplineLst
                    
    Args:
        cells: (list of cells)
        BSplineLst: (list of geom2d_BSpline objects) to be searched
            
   Returns: 
        disco_nodes: (list of nodes) the discovered nodes that are located on 
        the BSplineLst 
    '''
        
    disco_nodes = []
    tmp_nodes = []
    for c in cells:
        tmp_nodes.extend(c.nodes)
    
    tmp_nodes = list(set(tmp_nodes))    
    disco_nodes.extend(grab_nodes_on_BSplineLst(tmp_nodes,BSplineLst))
    disco_nodes = list(set(disco_nodes))               
    disco_nodes = sorted(disco_nodes, key=lambda Node: (Node.parameters[1],Node.parameters[2]))
    return disco_nodes


def grab_nodes_on_BSplineLst(nodes,BSplineLst, tolerance=1e-5):
    '''the grab_nodes_on_BSplineLst fuction determines the nodes of the input 
    argument list 'nodes' that are located on the BSplineLst. It loops though 
    all nodes and for each nodes it generates a projection on all bsplins in 
    the list of Bsplines (BSplineLst). If the projection distance is below a 
    certain tolerance, the node is to be returned in the disco_nodes list.
                
    Args:
        nodes: (list of nodes)
        BSplineLst: (list of geom2d_BSpline objects) to be searched
        tolerance: (float) to decide whether a node is on the BSpline, 
                    the default value is 1e-5.
            
   Returns: 
        disco_nodes: (list of nodes) the discovered nodes that are located on 
        the BSplineLst 
        
    Notes: Projection is slow. It would be nice to not use this function very 
        often.
    '''
    
    disco_nodes = []    
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
    return disco_nodes


def determine_a_nodes(mesh,a_BSplineLst,global_minLen,LayerID,factor=5):
    #is not needed anymore!
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
    #nodes = list(set(nodes))
    nodes = remove_duplicates_from_list_preserving_order(nodes)
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


def remove_duplicates_from_list_preserving_order(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


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


def merge_nodes_if_too_close(nodes,BSplineLst,global_minLen,tol=0.1):
    rm_idx=[]
    for i,n1 in enumerate(nodes[0:], start=0):
        n2 = nodes[i-1]
        v = gp_Vec2d(n1.Pnt2d,n2.Pnt2d)
        magnitude = v.Magnitude()

        if magnitude<=tol*global_minLen:
            cP = n1.Pnt2d.Translated(v.Multiplied(0.5))
            p2 = ProjectPointOnBSplineLst(BSplineLst,cP,1)
            n1.Pnt2d = p2[0]
            n1.parameters = ['modified',p2[1],p2[2]]
            rm_idx.append(i-1)         
   
    for index in sorted(rm_idx, reverse=True):
        del nodes[index]            
    return nodes


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