# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 11:01:06 2017

@author: TPflumm
"""
import scipy.io
import numpy as np
import math
import datetime

from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt



def centroid(points):
    x = [p[0] for p in points]
    y = [p[1] for p in points]
    centroid = (sum(x) / len(points), sum(y) / len(points))
    return centroid

def plot_nodes(nodes):
    fig, ax = plt.subplots()
    nodes_array = []
    for n in nodes:
        nodes_array.append([n.coordinates[0],n.coordinates[1]])
    nodes_array = np.asarray(nodes_array)  
    plt.plot(nodes_array[:,0],nodes_array[:,1],'.-')
    plt.axis('equal')
    plt.show()


def plot_mesh(nodes,elements,theta_11,data,data_name,title=None,VABSProperties=None,show_element_number=False,show_node_number=False,):
    fig, ax = plt.subplots()
    patches = []
    centroids = []
    for i,ele in enumerate(elements):
        #print ele
        if int(0) in ele:
            array = np.vstack((nodes[ele[0]-1],nodes[ele[1]-1],nodes[ele[2]-1]))
            centroids.append(centroid(array))
        else:
            array = np.vstack((nodes[ele[0]-1],nodes[ele[1]-1],nodes[ele[2]-1],nodes[ele[3]-1]))
            centroids.append(centroid(array))
        polygon = Polygon(array, True, edgecolor='k')
        patches.append(polygon)
    
    p = PatchCollection(patches, alpha=0.5, edgecolors = 'k')
    p.set_array(data)
    ax.add_collection(p)
    cbar = fig.colorbar(p, ax=ax)
    cbar = cbar.ax.set_ylabel(data_name)
    
    if len(theta_11)==len(elements):
        for i,cent in enumerate(centroids):
            lfactor=0.5
            dx = lfactor*math.cos(math.radians(theta_11[i]))
            dy = lfactor*math.sin(math.radians(theta_11[i]))
            ax.arrow(cent[0], cent[1], dx, dy, head_width=0.05, head_length=0.1, fc='k', ec='k') 
        
    #ax.scatter(nodes[:,0],nodes[:,1],c='k',)
    plt.axis('equal')
    if title!=None:
        ax.set_title(title)
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
            
    if VABSProperties != None:
        pass
        CG, = plt.plot(VABSProperties.Xm2,VABSProperties.Xm3,'ro', label='CG: Mass Center')
        ax.annotate('CG', (VABSProperties.Xm2,VABSProperties.Xm3),fontsize=20)
        GC, = plt.plot(VABSProperties.Xg2,VABSProperties.Xg3,'bo', label='GC: Geometric Center')
        ax.annotate('GC', (VABSProperties.Xg2,VABSProperties.Xg3),fontsize=20)
        NA, = plt.plot(VABSProperties.Xt2,VABSProperties.Xt3,'go',  label='NA: Neutral Axes')
        ax.annotate('NA', (VABSProperties.Xt2,VABSProperties.Xt3),fontsize=20)
        plt.legend(handles=[CG,GC,NA])
        
        
        #place a text box in upper left in axes coords
        textstr = 'mass per unit span \t\t\t\t\t = $%.2f$ [g/mm]\nmass moments of intertia about x1 axis \t= $%.2f$ [g/mm2]'%(VABSProperties.MpUS,VABSProperties.I1)
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        ax.text(0.02, 0.97, textstr, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=props)
                
        
        #place a text box in upper left in axes coords
        textstr = 'Classical Stiffness Matrix:\n'+np.array2string(VABSProperties.CS, precision=2, separator=',  ')
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        ax.text(0.02, 0.12, textstr, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=props)
            
        

        if VABSProperties.Xs2 != None:
            SC, = plt.plot(VABSProperties.Xs2,VABSProperties.Xs3,'ko',label='SC: Generalized Shear Center')
            ax.annotate('SC', (VABSProperties.Xs2,VABSProperties.Xs3),fontsize=20)
            plt.legend(handles=[CG,GC,NA,SC])

    
            
    plt.show()
    return (fig,ax)
    

def plot_cells(cells,nodes,attr1,VABSProperties=None,title='None', savepath=None, plotTheta11=False, plotDisplacement=False,):
    nodes_array = []
    for n in nodes:
        if plotDisplacement:
            nodes_array.append([[n.coordinates[0]+n.displacement[1],n.coordinates[1]+n.displacement[2]]])
        else:
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
        if '.' in attr1:
            attrLst = attr1.split('.')
            tmp = getattr(c, attrLst[0])
            data.append(getattr(tmp,attrLst[1]))
        else:
            data.append(getattr(c, attr1))  
    data = np.asarray(data)
    data_name = attr1
    
    theta_11 = []
    if plotTheta11==True:
        for c in cells:
            theta_11.append(getattr(c, 'theta_11'))  
        theta_11 = np.asarray(theta_11)
    
    
    fig,ax = plot_mesh(nodes_array,element_array,theta_11,data,data_name,title,VABSProperties,False,False)    
    
    if savepath!=None:
        #savepath = 'jobs/VHeuschneider/figures/R90_config.svg'
        #datestr = datetime.datetime.now().strftime("%Y%m%d_%H%M")
        tmp_fig = plt.gcf()
        tmp_fig.set_size_inches(11.69, 8.27)    #a4 landscape
        tmp_fig.savefig(savepath, dpi=300, orientation='landscape', papertype='a4')
   
    
    return (fig, ax)


    
    
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
    
    


