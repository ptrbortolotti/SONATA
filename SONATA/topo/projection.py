# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 16:32:15 2016

@author: TPflumm
"""
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
    
    #timelines(1,begin,end,color[4])

#    for i,item in enumerate(projection):
#        timelines(0,item[0],item[1],color[int(item[2])])