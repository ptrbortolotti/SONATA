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
def projection_of_layers(layup,begin,end):
    
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
    datalist = []
    for iv in t1:
        datalist.append(iv.data)
    nextlayer = max(datalist)+1
    
    #Bring everything together
    for iv in t2: 
        if iv.data == True:
            t1.addi(iv.begin,iv.end,nextlayer)
        elif iv.data == None:
            t1.addi(iv.begin,iv.end,nextlayer)
    
    #transform back to python np.array
    tmp = []
    for iv in t1:
        tmp.append([iv.begin,iv.end,iv.data])
    Projection = np.asarray(tmp)       
    
    return Projection

    
    
if __name__ == '__main__':   
     
    layup = np.array([[  0.15  ,   0.3  ,   0.23 , 45.  ,   1.  ],
                  [  0.5   ,  0.9  ,   0.23 , 45.  ,   1.  ],
                  [  0.0  ,   0.15  ,   0.23 , 45.  ,   1.  ],])    
    
    begin = 0.9
    end   = 0.92
    projection = projection_of_layers(layup,begin,end)
    
                
    #=====================================================
    #               PLOT!!!!
    #=====================================================
    
    plt.figure(1)
    #============================= 
    plt.subplot(211)
    plt.xlim(0,1)
    plt.ylim(0,4)
    plt.ylabel('Layup')    
    plt.xlabel('Coordinate')
    color=plt.cm.Paired(np.linspace(0,1,6))
    
    for i,item in enumerate(layup):
        timelines(2,item[0],item[1],color[i+1])
    
    timelines(1,begin,end,color[4])

    for i,item in enumerate(projection):
        timelines(0,item[0],item[1],color[int(item[2])])