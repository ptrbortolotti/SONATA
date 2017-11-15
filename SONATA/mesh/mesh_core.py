
import numpy as np
import pickle
import matplotlib.pyplot as plt

from OCC.gp import gp_Pnt2d
from triangle import triangulate, plot as tplot 

from SONATA.topo.BSplineLst_utils import get_BSplineLst_length
from SONATA.mesh.mesh_utils import equidistant_nodes_on_BSplineLst                   
from SONATA.mesh.node import Node
from SONATA.mesh.cell import Cell

#===================FUNCTIONALITIES==========================================

def triangle_mesh(nodes,options):
    points = []
    for n in nodes:
        points.append([n.Pnt2d.X(),n.Pnt2d.Y()])
      
    old_vertices = np.asarray(points)
    tmp = []
    for i,v in enumerate(old_vertices):  
        if i == len(old_vertices)-1:
            tmp.append([i,0])
        else: tmp.append([i,i+1])
    segments_core = np.asarray(tmp)
        
    poly = {'vertices': None, 'holes': None, 'segments': None}
    poly['vertices'] = old_vertices
    poly['segments'] = segments_core

    mesh =  triangulate(poly, options)
    
    new_vertices = []
    for v in mesh['vertices']:
        if v not in old_vertices:
            new_vertices.append(v)
    return mesh


def find_node(nodeLst,ID):
    for n in nodeLst:
        if n.id == ID:
            tmp = n
            break
        else:
            tmp = None
            
    return tmp


def gen_core_cells(a_nodes,area=1.0,**kwargs):
    '''The gen_core_cells function generates the triagular mesh within the 
    a_nodes polygon.        
    
    Args:
        a_nodes: (list of nodes), the list contains all nodes that are on the
            innermost boundary of the generated topology.
        area: (float), is a resolution parameter / area constraint for the 
            triangle_mesh function.
    
    Kwargs: 
        options: the options can be passed to the triagle_mesh function,
            default is the optionstring 'pa%s' % (area)
            
        
    Returns: 
        [c_cells,c_nodes]: c_cells is a list of the newly generated cell 
            objects, while c_nodes is a list of the newly generated node 
            objects.              
    '''
    
    #KWARGS:
    if kwargs.get('options') !=  None:
        options = kwargs.get('options')
    else:
        options = 'pa%s' % (area)
    
    mesh = triangle_mesh(a_nodes,options)  
    
    #TODO: substract v from old_nodes, build norm, use argmin and check if argmin is below tol.!!!!! Schickies geniale idee!
    existence = False
    c_nodes = []
    connector = []      #connector stores the information [tri_vertex_id, old_node_id]
    for tri_id,v in enumerate(mesh['vertices']):
        for an in a_nodes:    
            if np.all(np.isclose(np.asarray(v),np.asarray([an.Pnt2d.X(),an.Pnt2d.Y()]))):
                #if the vertice allready exists as a_node, do not create a new node, and write id to connector
                existence=True
                connector.append([tri_id,an.id])
                break
            else:
                existence=False
         
        if existence == False:
            #create new node for every new vertice that is not in a_nodes:
            tmp_node = Node(gp_Pnt2d(v[0],v[1]))
            connector.append([tri_id,tmp_node.id])
            c_nodes.append(tmp_node)

                
    connector = np.asarray(connector)
    nodes = a_nodes+c_nodes
    
    c_cells =[]
    for ele in mesh['triangles']:
        nodeLst = []
        for tri_id in ele: 
            NodeID = connector[tri_id,1]
            nodeLst.append(find_node(nodes,NodeID))
    
        c_cells.append(Cell(nodeLst))
        
#    plt.figure(figsize=(20, 20))
#    ax = plt.subplot(111)
#    tplot.plot(ax, **mesh)
#    plt.show()
    
    return [c_cells,nodes]

#===================MAIN==========================================
if __name__ == '__main__':   
    #LOAD .pkl data with SegmentLst
    filename = 'C:\TPflumm_local\work\SONATA\sec_config.pkl'
    with open(filename, 'rb') as handle:
        SegmentLst = pickle.load(handle)
        
    #Build wires for each layer and segment
    for seg in SegmentLst:
        seg.build_wire()
        for layer in seg.LayerLst:
            layer.build_wire()
    
    
    layer = SegmentLst[-1].LayerLst[-1]
    Resolution = 100 # Nb of Points on Segment0
    length = get_BSplineLst_length(SegmentLst[0].BSplineLst)
    global_minLen = round(length/Resolution,5)
    
    a_nodes = equidistant_nodes_on_BSplineLst(layer.Boundary_BSplineLst, True, True, True, minLen = global_minLen, LayerID = layer.ID[0])
    [c_cells,c_nodes] = gen_core_cells(a_nodes,2)
    
    for c in c_cells:
        c.structured = False
        c.theta_3 = 0
        c.MatID = int(seg.CoreMaterial)
        c.calc_theta_1()
    
    from SONATA.display.display_mesh import plot_cells
    plot_cells(c_cells, c_nodes, 'MatID')
    
#    from OCC.Display.SimpleGui import init_display
#    display, start_display, add_menu, add_function_to_menu = init_display()
#    display.Context.SetDeviationAngle(1e-6)       # 0.001 default. Be careful to scale it to the problem.
#    display.Context.SetDeviationCoefficient(1e-6) # 0.001 default. Be careful to scale it to the problem. 
#    for n in c_nodes:
#        display.DisplayColoredShape(n.Pnt2d, 'BLACK')    