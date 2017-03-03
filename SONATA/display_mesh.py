# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 11:01:06 2017

@author: TPflumm
"""
import scipy.io
import numpy as np

from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt



def centroid(points):
    x = [p[0] for p in points]
    y = [p[1] for p in points]
    centroid = (sum(x) / len(points), sum(y) / len(points))
    return centroid

def plot_mesh(nodes,elements,data,show_element_number=False,show_node_number=False):
    fig, ax = plt.subplots()
    patches = []
    centroids = []
    for i,ele in enumerate(elements):
        
        if int(0) in ele:
            array = np.vstack((nodes[ele[0]-1],nodes[ele[1]-1],nodes[ele[2]-1]))
        else:
            array = np.vstack((nodes[ele[0]-1],nodes[ele[1]-1],nodes[ele[2]-1],nodes[ele[3]-1]))
            centroids.append(centroid(array))
        polygon = Polygon(array, True)
        patches.append(polygon)
    
    p = PatchCollection(patches, alpha=0.5)
    p.set_array(data)
    ax.add_collection(p)
    cbar = fig.colorbar(p, ax=ax)
    cbar = cbar.ax.set_ylabel('data')
    
    ax.scatter(nodes[:,0],nodes[:,1],c='k',)
    plt.axis('equal')
    ax.set_title('example')
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')
    
    ##display element number
    if show_element_number == True:
        for i, item in enumerate(centroids):
           ax.annotate(i+1, (item[0],item[1]))
    
    ##display node number
    if show_node_number == True:
        for i, item in enumerate(nodes):
            ax.annotate(i+1, (item[0],item[1]), color='red')
    plt.show()
    

def plot_cells(cells):
    #Get all nodes in cells
    nodes = [] 
    for cell in cells:
        for node in cell.nodes:
            if node not in nodes:
                nodes.append(node)
    nodes = sorted(nodes, key=lambda Node: (Node.id))
    cells = sorted(cells, key=lambda Cell: (Cell.id))   
    
    nodes_array = []
    for n in nodes:
        nodes_array.append([n.coordinates[0],n.coordinates[1]])
    nodes_array = np.asarray(nodes_array)   
     
    element_array = []
    for c in cells:
        tmp = []
        for i in range(0,4):
            if i<len(c.nodes):
                tmp.append(c.nodes[i].id)
            else:
                tmp.append(0)
        element_array.append(tmp)
    element_array = np.asarray(element_array)  
    
    data = []
    for c in cells:
        data.append(c.minimum_angle)  
    data = np.asarray(data)  
    
    plot_mesh(nodes_array,element_array,data,False,False)    
    
    return None
    
    
    
    
#=============================================================
#       IF MAIN:
#============================================================
if __name__ == '__main__':

    nodes = scipy.io.loadmat('nodes.mat')
    nodes = nodes['nodes']
    elements = scipy.io.loadmat('elements.mat')
    elements = elements['elements']
    layup_angle = scipy.io.loadmat('layupAngle.mat')
    layup_angle = layup_angle['layup_angle']
    layup_angle = layup_angle[:,0]

    plot_mesh(nodes,elements,layup_angle,False,False)
    
    


