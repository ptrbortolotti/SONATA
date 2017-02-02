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
    p.set_array(data[:,0])
    ax.add_collection(p)
    fig.colorbar(p, ax=ax)
    
    ax.scatter(nodes[:,0],nodes[:,1],c='k',)
    plt.axis('equal')
    
    ##display element number
    if show_element_number == True:
        for i, item in enumerate(centroids):
           ax.annotate(i+1, (item[0],item[1]))
    
    ##display node number
    if show_node_number == True:
        for i, item in enumerate(nodes):
            ax.annotate(i+1, (item[0],item[1]), color='red')
    plt.show()
    

    
    
    
    
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
    

    plot_mesh(nodes,elements,layup_angle,False,False)
    
    


