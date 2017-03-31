
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
    
    #KWARGS:
    if kwargs.get('options') !=  None:
        options = kwargs.get('options')
    else:
        options = 'pa%s' % (area)
    
    mesh = triangle_mesh(a_nodes,options)  
    
    tmp = []
    for n in a_nodes:
        tmp.append([n.Pnt2d.X(),n.Pnt2d.Y(),n.id])
    old_vertices = np.asarray(tmp)
    
    new_vertices = []
    c_nodes = []
    connector = []
    for tri_id,v in enumerate(mesh['vertices']):
        if v in old_vertices[:,:2]:
            new_vertices.append({'coords': v, 'tri_id':tri_id, 'new': False})
            idx = int(np.where(np.all(old_vertices[:,:2]==v,axis=1))[0][0])
            connector.append([tri_id,int(old_vertices[idx,2])])
        elif v not in old_vertices[:,:2]:
            new_vertices.append({'coords': v, 'tri_id':tri_id, 'new': True})
            tmp_node = Node(gp_Pnt2d(v[0],v[1]))
            NodeID = tmp_node.id
            connector.append([tri_id,NodeID])
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
    
    return [c_cells,c_nodes]

#===================MAIN==========================================
if __name__ == '__main__':   
    #LOAD .pkl data with SegmentLst
    filename = 'naca0012_cspar.pkl'
    with open(filename, 'rb') as handle:
        SegmentLst = pickle.load(handle)
        
    #Build wires for each layer and segment
    for seg in SegmentLst:
        seg.build_wire()
        for layer in seg.LayerLst:
            layer.build_wire()
    
    
    layer = SegmentLst[-1].LayerLst[-1]
    Resolution = 500 # Nb of Points on Segment0
    length = get_BSplineLst_length(SegmentLst[0].BSplineLst)
    global_minLen = round(length/Resolution,5)
    
    a_nodes = equidistant_nodes_on_BSplineLst(layer.Boundary_BSplineLst, True, True, True, minLen = global_minLen, LayerID = layer.ID[0])
    [c_cells,c_nodes] = gen_core_cells(a_nodes,22,1.4)
    
#    for c in c_cells:
#        print c
