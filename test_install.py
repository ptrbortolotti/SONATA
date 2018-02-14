# -*- coding: utf-8 -*-
"""
Created on Mon May 29 09:30:03 2017

@author: TPflumm
"""

import OCC
import intervaltree 
import shapely
import openmdao

#==============================================================================
#TEST SHAPELY Parallel OFFSET:
#==============================================================================
from matplotlib import pyplot
from shapely.geometry import LineString

def plot_coords(ax, x, y, color='#999999', zorder=1):
    ax.plot(x, y, 'o', color=color, zorder=zorder)

def plot_line(ax, ob):
    parts = hasattr(ob, 'geoms') and ob or [ob]
    for part in parts:
        x, y = part.xy
        ax.plot(x, y, linewidth=3, solid_capstyle='round', zorder=1)

def set_limits(ax, x_range, y_range):
    ax.set_xlim(*x_range)
    ax.set_xticks(list(range(*x_range)) + [x_range[-1]])
    ax.set_ylim(*y_range)
    ax.set_yticks(list(range(*y_range)) + [y_range[-1]])
    ax.set_aspect(1)

line = LineString([(0, 0), (1, 1), (0, 2), (2, 2), (3, 1), (1, 0)])
line_bounds = line.bounds
ax_range = [int(line_bounds[0] - 1.0), int(line_bounds[2] + 1.0)]
ay_range = [int(line_bounds[1] - 1.0), int(line_bounds[3] + 1.0)]

fig = pyplot.figure(1)

# 1
ax = fig.add_subplot(221)

plot_line(ax, line)
x, y = list(line.coords)[0]
plot_coords(ax, x, y)
offset = line.parallel_offset(0.5, 'left', join_style=1)
plot_line(ax, offset)

ax.set_title('a) left, round')
set_limits(ax, ax_range, ay_range)

#2
ax = fig.add_subplot(222)

plot_line(ax, line)
x, y = list(line.coords)[0]
plot_coords(ax, x, y)

offset = line.parallel_offset(0.5, 'left', join_style=2)
plot_line(ax, offset)

ax.set_title('b) left, mitred')
set_limits(ax, ax_range, ay_range)

#3
ax = fig.add_subplot(223)

plot_line(ax, line)
x, y = list(line.coords)[0]
plot_coords(ax, x, y)
offset = line.parallel_offset(0.5, 'left', join_style=3)
plot_line(ax, offset)

ax.set_title('c) left, beveled')
set_limits(ax, ax_range, ay_range)

#4
ax = fig.add_subplot(224)

plot_line(ax, line)
x, y = list(line.coords)[0]
plot_coords(ax, x, y)
offset = line.parallel_offset(0.5, 'right', join_style=1)
plot_line(ax, offset)

ax.set_title('d) right, round')
set_limits(ax, ax_range, ay_range)

pyplot.show()



#==============================================================================
#TEST TRIANGLE Package
#==============================================================================
import triangle
import triangle.plot
import numpy as np
import matplotlib.pyplot as plt

theta = np.linspace(0, 2*np.pi, 33)[:-1]
pts = np.vstack((np.cos(theta), np.sin(theta))).T
A = dict(vertices=pts)
B = triangle.triangulate(A, 'q')
triangle.plot.compare(plt, A, B)
plt.show()


#==============================================================================
#TEST INTERVALTREE
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
import intervaltree 


#=====================================================
def timelines(y, xstart, xstop, c='b'):
    """Plot timelines at y from xstart to xstop with given color."""   
    if xstart > xstop:
        plt.hlines(y,0,xstop, colors=c, lw=4)
        plt.hlines(y,xstart,1,colors=c,lw=4)
    else:
        plt.hlines(y, xstart, xstop,colors=c, lw=4)

        
def datafunc1(interval, islower=None):
    if islower == True:
        return False
    elif islower == False:
        return True
    else: 
        return False
    
        
def datafunc2(interval, islower=None):
    if islower == True:
        return True
    elif islower == False:
        return False
    else: 
        return False
         

       
#=====================================================
def projection_of_layers2(layup,begin,end,idx):
    if layup.ndim == 1:
         layup = np.array([layup])     
    #print layup
        
    #Define gloab interval between 0 and 1
    t0 = intervaltree.IntervalTree()
    t0.addi(0,1)
    
    #Put layup format into Intervaltree datastructure
    t1 = intervaltree.IntervalTree()
    for i,item in enumerate(layup):
        if item[1]<item[0]:
            t1.addi(0,float(item[1]),i+1)
            t1.addi(float(item[0]),1,i+1)
        else:
            t1.addi(float(item[0]),float(item[1]),i+1)
    
    #Initialize to be projected interval between begin and end
    t2 = intervaltree.IntervalTree()
    t2.addi(begin,end)
    
    #Generate Negativ of Intervals
    for iv in t1:
        t0.chop(iv.begin,iv.end)
    
    #Intersect t0 with t2 and add to t1
    for iv in t0:
        t2.slice(iv.begin,datafunc1)
    
    for iv in t0:   
        t2.slice(iv.end,datafunc2)
        
    for iv in t1:
        t2.remove_overlap(iv.begin,iv.end)
        
    #Define Layer Number   
#    datalist = []
#    for iv in t1:
#        datalist.append(iv.data)
#    nextlayer = max(datalist)+1
    
    #Bring everything together
    for iv in t2: 
        if iv.data == True:
            t1.addi(iv.begin,iv.end,idx)
        elif iv.data == None:
            t1.addi(iv.begin,iv.end,idx)
    
    #transform back to python np.array
    tmp = []
    for iv in t1:
        tmp.append([iv.begin,iv.end,iv.data])
    Projection = np.asarray(tmp)       
    
    return Projection

def projection_of_layers(layup,begin,end):
    ''' projection_of_layers fills the current set of layers between begin and end with the new layer
    begin: float[0:1] 
    end: float[0:1] 
    '''
    idx = max(layup[:,2])+1
    t1 = intervaltree.IntervalTree()
    for i,item in enumerate(layup):
        t1.addi(float(item[0]),float(item[1]),int(item[2]))
    
    #Initialize to be projected interval between begin and end
    t2 = intervaltree.IntervalTree()
    if begin>end:
        t2.addi(0,end)
        t2.addi(begin,1)
    else:
        t2.addi(begin,end)
        
    #chop intervaltree t1
    for iv in t2:
        t1.chop(iv.begin,iv.end)
    
    #insert t2  
    for iv in t2:
        t1.addi(iv.begin,iv.end,idx)
    
    #transform back to python np.array
    tmp_pj = []
    for iv in t1:
        tmp_pj.append([iv.begin,iv.end,iv.data])
    Projection = np.asarray(tmp_pj)   
    return Projection

        

def layup_to_projection(Layup):
    nlayers = int(np.size(Layup,0))
    #clean up Layup array
    a = np.delete(Layup,[2,3,4],1)
    b = np.transpose(np.linspace(1,nlayers,nlayers))
    c = np.c_[a,b]  
    clean_layup = np.insert(c,0, np.array((0,1,0)),0)
    tmp = clean_layup[[0]]
    #print clean_layup
    
    #Iterate over Layup
    projectionlist =  []
    for i in range(1,nlayers,1):
        begin   = float(clean_layup[i][0])
        end     = float(clean_layup[i][1])  
        tmp = projection_of_layers(tmp,begin,end)
        tmp = tmp[np.lexsort(np.fliplr(tmp).T)]
        projectionlist.append(tmp)

    return projectionlist        
        

 
#==============================================================================       
if __name__ == '__main__':   
    
    Layup = np.array([[  0.  ,   1.  ,   0.23 , 45.  ,   1.  ],
                  [  0.   ,  1.  ,   0.23 , 45.  ,   1.  ],
                  [  0.3  ,  0.7  ,  0.6  , 45.  ,   2.  ],
                  [  0.05 ,  0.6 ,  0.6  , 45.  ,   2.  ],
                  [  0.85  , 0.2  ,  0.6  , 45.  ,   2.  ],
                  [  0.9 ,  0.4 ,  0.6  , 45.  ,   2.  ],
                  [  0.20 ,  0.532,  0.6  , 45.  ,   2.  ],
                  [  0.25 ,  0.3 ,  0.6  , 45.  ,   2.  ]],)
    
    nlayers = int(np.size(Layup,0))
    
    #clean up Layup array
    a = np.delete(Layup,[2,3,4],1)
    b = np.transpose(np.linspace(1,nlayers,nlayers))
    c = np.c_[a,b]
    clean_layup = np.insert(c,0, np.array((0,1,0)),0)
    tmp = clean_layup[[0]]

    
    projectionlist=layup_to_projection(Layup)
    #Projection = projection_of_layers2(Layup,0,1,7)

    #=====================================================
    #               PLOT!!!!
    #=====================================================
    
    plt.figure(1)
    #============================= 
    plt.subplot(121)
    plt.xlim(0,1)
    plt.ylim(0,10)
    plt.ylabel('Layup')    
    plt.xlabel('Coordinate')
    color=plt.cm.Paired(np.linspace(0,1,nlayers+4))
    
    for i,item in enumerate(clean_layup):
        timelines(i,item[0],item[1],color[int(item[2])])
        
        
    plt.subplot(122)
    plt.xlim(0,1)
    plt.ylim(0,10)
    plt.ylabel('Cummulated Boundaries')    
    plt.xlabel('Coordinate')

    
    for i,idx in enumerate(projectionlist):
        for j,item in enumerate(idx):
            timelines(i+1,item[0],item[1],color[int(item[2])])
    


