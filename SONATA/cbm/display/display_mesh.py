# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 11:01:06 2017

@author: TPflumm
"""
import scipy.io
import numpy as np
import math
import datetime
import os

from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap



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


def plot_mesh(nodes, elements, theta_11, data, data_name, title=None, VABSProperties=None, 
              show_element_number=False, show_node_number=False, invert_xaxis = True, lfactor=0.5e-2, **kw):
    
    """
    To be continued...
    
    Parameters
    ----------
    nodes : list
    elements : list
    data : 
    data_name : string
    
    """
    #print(data_name, max(data))
    alpha = 0.7
    if 'cmap' in kw:
        cmap = plt.cm.get_cmap(kw['cmap'])
    
    elif data_name == 'sf':
        colors = [(0.6, 0, 0), (1, 1, 0), (0, 0.5, 0)]  # R -> G -> B
        cmap_name = 'my_list'
        cmap = LinearSegmentedColormap.from_list(
        cmap_name, colors, N=6)
        #cmap.set_over(color='white')
        #cmap.set_under(color='k')
        
    elif data_name == 'MatID':
        cmap = plt.cm.get_cmap()
        # extract all colors from the .jet map
        cmaplist = [cmap(i) for i in range(cmap.N)]
        # force the first color entry to be grey
        #cmaplist[0] = (.5, .5, .5, 1.0)
        # create the new map
        cmap = LinearSegmentedColormap.from_list('Custom cmap', cmaplist, max(data))
        
    else:
        cmap = plt.cm.get_cmap()
        
    if 'vmin' in kw:
        vmin = kw['vmin']
    else:
        vmin = None
        
    if 'vmax' in kw:
        vmax = kw['vmax']
    else:
        vmax = None
        
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
    
    p = PatchCollection(patches, alpha=alpha, cmap=cmap, edgecolors = 'k')
    p.set_array(data)
    p.set_clim(vmin, vmax)
    _ = ax.add_collection(p)
    
    
    cbar = fig.colorbar(p, ax=ax, drawedges=True)
    cbar.ax.set_ylabel(data_name)
    
    if data_name == 'MatID':
        cbar.set_ticks(np.linspace(1, max(data), max(data)))
        cbar.set_ticklabels(np.linspace(1, max(data), max(data)))
        p.set_clim(0.5, max(data)+0.5)
    
    if len(theta_11)==len(elements):
        for i,cent in enumerate(centroids):
            
            dx = lfactor*math.cos(math.radians(theta_11[i]))
            dy = lfactor*math.sin(math.radians(theta_11[i]))
            ax.arrow(cent[0], cent[1], dx, dy, width = 0.01e-2, head_width=0.1e-2, head_length=0.1e-2, fc='k', ec='k')
        
    #ax.scatter(nodes[:,0],nodes[:,1],c='k',)
    plt.axis('equal')
    if title!=None:
        ax.set_title(title)
    ax.set_xlabel(r'$x_2$ in m')
    ax.set_ylabel(r'$x_3$ in m')
    
    #plot coordinate system.
    cslength=0.015
    ax.arrow(0, 0, cslength, 0, color='lime')
    ax.arrow(0, 0, 0, cslength, color='deepskyblue')

    ax.annotate(r'$x_2$', (cslength,0),fontsize=10, color = 'lime')
    ax.annotate(r'$x_3$', (0,cslength),fontsize=10, color='deepskyblue')
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
        CG, = plt.plot(VABSProperties.Xm[0],VABSProperties.Xm[1],'ro', label='CG: Mass Center')
        #ax.annotate('CG', (VABSProperties.Xm2,VABSProperties.Xm3),fontsize=20)
        GC, = plt.plot(VABSProperties.Xg[0],VABSProperties.Xg[1],'b^', label='GC: Geometric Center')
        #ax.annotate('GC', (VABSProperties.Xg2,VABSProperties.Xg3),fontsize=20)
        NA, = plt.plot(VABSProperties.Xt[0],VABSProperties.Xt[1],'gs',  label='NA: Neutral Axes')
        #ax.annotate('NA', (VABSProperties.Xt2,VABSProperties.Xt3),fontsize=20)
        plt.legend(handles=[CG,GC,NA])
        
        
#        #place a text box in upper left in axes coords
#        textstr = 'mass per unit span \t\t\t\t\t = $%.2f$ [g/mm]\nmass moments of intertia about x1 axis \t= $%.2f$ [g/mm2]'%(VABSProperties.MpUS,VABSProperties.I1)
#        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
#        ax.text(0.02, 0.97, textstr, transform=ax.transAxes, fontsize=10,
#                verticalalignment='top', bbox=props)
#                
#        
#        #place a text box in upper left in axes coords
#        textstr = 'Classical Stiffness Matrix:\n'+np.array2string(VABSProperties.CS, precision=2, separator=',  ')
#        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
#        ax.text(0.02, 0.12, textstr, transform=ax.transAxes, fontsize=10,
#                verticalalignment='top', bbox=props)
            
        

        if isinstance(VABSProperties.Xs,np.ndarray):
            SC, = plt.plot(VABSProperties.Xs[0],VABSProperties.Xs[1],'kD',label='SC: Generalized Shear Center')
            #ax.annotate('SC', (VABSProperties.Xs2,VABSProperties.Xs3),fontsize=20)
            plt.legend(handles=[CG,GC,NA,SC])

    if invert_xaxis:
        ax.invert_xaxis()
        
    plt.show()
    return (fig,ax)
    

def plot_cells(cells,nodes,attr1, VABSProperties=None, title='None', plotTheta11=False, plotDisplacement=False, **kw):
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
    
    
    fig,ax = plot_mesh(nodes_array,element_array,theta_11,data,data_name,title, VABSProperties, False, False, **kw)    
   
    if 'savepath' in kw:

        if not os.path.exists(kw['savepath']+'figures'):  # create 'figures' Folder if not already existing
            os.mkdir(kw['savepath']+'figures')

        # datestr = datetime.datetime.now().strftime("%Y%m%d_%H%M")
        # fname = kw['savepath'].split('.')[0]+'_'+datestr+'.'+kw['savepath'].split('.')[1]
        fname = kw['savepath']+'/figures/blade_section_'+kw['section']+'.png'
        # print(fname)
        # tmp_fig = plt.gcf()
        tmp_fig = fig
        # tmp_fig.set_size_inches(11.69, 8.27)    #a4 landscape
        tmp_fig.set_size_inches(10, 5)    #a4 landscape
        tmp_fig.savefig(fname, dpi=300, orientation='landscape', papertype='a4')

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
    
    


