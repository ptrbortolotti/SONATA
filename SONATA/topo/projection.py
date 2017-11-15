# -*- coding: utf-8 -*-
"""
The SONATA.topo.projection module uses the python invervaltree package 
(Editable interval tree data structure for Python) to determine the cummulated 
layup boundaries of the SONATA Layup definition.

 consits of the following functions:
    
    - timelines(y, xstart, xstop, c='b'): Plot timelines at y from xstart to 
            xstop with given color.
    - datafunc1(interval,islower=None): boolean function that can be used in 
            the intervaltree slice function 
    - datafunc2(interval,islower=None): boolean function that can be used in 
            the intervaltree slice function 
    - projection_of_layers2(layup,begin,end,idx): Truthfully I currently don't
            know the exact inded of this funtion.
    - insert_interval_in_layup(layup,begin,end): insert_interval_in_layup 
            generates a intervaltree structure t1 from layup 
            and inserts a new intervaltree set t2 [begin, end] or 
            [0,end][begin,1] into the existing the set t1
    - chop_interval_from_layup(layup,begin,end): chops an intervaltree 
            structure t1 from layup to a given interval [begin, end] or 
            [0,end][begin,1].
    - cummulated_layup_boundaries(Layup): fuctions generates a interval 
            structure from the Layup definitions so that for each layer the 
            information of the underlying layers that create and interval#
            between zero and one is given.

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

        
def datafunc1(interval,islower=None):
    '''boolean function that can be used in the intervaltree slice function '''
    if islower == True:
        return False
    elif islower == False:
        return True
    else: 
        return False
    
        
def datafunc2(interval,islower=None):
    '''boolean function that can be used in the intervaltree slice function '''
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

def insert_interval_in_layup(layup,begin,end):
    ''' insert_interval_in_layup generates a intervaltree structure t1 from layup 
    and inserts a new intervaltree set t2 [begin, end] or [0,end][begin,1] into
    the existing the set t1. First t1 is chopped with t2, then t2 is added to t1.
    
    for documentation of the intervaltree module refer to: 
        https://pypi.python.org/pypi/intervaltree
    
    Args:
        begin: float[0:1]: start of the interval 
        end: float[0:1]  end of the interval to be inserted. Note that begin
            can be larger than end (The layer is defined across the origin)
            
    Returns: 
        Projection: numpy array of the intervaltree structure of the form:
                  start end layer#
        np.array([[ 0.3  0.7  3. ]
                  [ 0.7  1.   2. ]
                  [ 0.   0.3  2. ]])                
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


def chop_interval_from_layup(layup,begin,end):
    ''' chop_interval_from_layup chops an intervaltree structure t1 from layup 
    to a given interval [begin, end] or [0,end][begin,1].
    
    for documentation of the intervaltree module refer to: 
        https://pypi.python.org/pypi/intervaltree
    
    Args:
        begin: float[0:1]: start of the interval 
        end: float[0:1]  end of the interval to which t1 is beeing chopped.
        Note that begin can be larger than end (The layer is defined across
        the origin)
            
    Returns: 
        Projection: numpy.array of the intervaltree structure of the form:
                  start end layer#
        np.array([[ 0.3  0.7  3. ]
                  [ 0.7  1.   2. ]
                  [ 0.   0.3  2. ]])                
    '''
    #get layup to intervaltree definition
    t1 = intervaltree.IntervalTree()
    for i,item in enumerate(layup):
        t1.addi(float(item[0]),float(item[1]),int(item[2]))
        
    #result_tree chop
    if begin>end:
        t1.chop(end,begin)
    else:
        t1.chop(0,begin)
        t1.chop(end,1)
        
    #transform back to python np.array
    tmp_pj = []
    for iv in t1:
        tmp_pj.append([iv.begin,iv.end,iv.data])
    Projection = np.asarray(tmp_pj)   
    return Projection


def cummulated_layup_boundaries(Layup):
    '''The cummulated_layup_boundaries fuctions generates a interval structure 
    from the Layup definitions so that for each layer the information of the 
    underlying layers that create and interval between zero and one is given.
    This information can be used to determine the offset origin for the layer. 
    the relevant_cummulated_layup_boundaries function chops this valuable 
    information to the interval of the Layup itself to determin the offset
    origin in the given interval.
        
    Args:
        Layup: float[0:1]: start of the interval 
            
    Returns: 
        Projection: numpy array of the intervaltree structure of the form:
                  start end layer#
        np.array([[ 0.3  0.7  3. ]
                  [ 0.7  1.   2. ]
                  [ 0.   0.3  2. ]])                
    '''
    
    nlayers = int(np.size(Layup,0))
    #clean up Layup array
    a = np.delete(Layup,[2,3,4],1)
    b = np.transpose(np.linspace(1,nlayers,nlayers))
    c = np.c_[a,b]  
    clean_layup = np.insert(c,0, np.array((0,1,0)),0)
    tmp = clean_layup[[0]]
    
    #Iterate over Layup
    projectionlist =  []
    projectionlist.append(tmp)
    for i in range(1,nlayers,1):
        begin   = float(clean_layup[i][0])
        end     = float(clean_layup[i][1])  
        tmp = insert_interval_in_layup(tmp,begin,end)
        tmp = tmp[np.lexsort(np.fliplr(tmp).T)]
        projectionlist.append(tmp)
        
    return projectionlist
        

def relevant_cummulated_layup_boundaries(Layup):
    '''The relevant_cummulated_layup_boundaries function chops this valuable 
    information to the interval of the Layup itself to determin the offset
    origin in the given interval.        
    
    Args:
        Layup: float[0:1]: start of the interval 
            
    Returns: 
        Projection: numpy array of the intervaltree structure of the form:
                  start end layer#
        np.array([[ 0.3  0.7  3. ]
                  [ 0.7  1.   2. ]
                  [ 0.   0.3  2. ]])                
    '''
    nlayers = int(np.size(Layup,0))
    #clean up Layup array
    a = np.delete(Layup,[2,3,4],1)
    b = np.transpose(np.linspace(1,nlayers,nlayers))
    c = np.c_[a,b]  
    clean_layup = np.insert(c,0, np.array((0,1,0)),0)
    tmp = clean_layup[[0]]
    #print clean_layup
    
    #Iterate over Layup
    relevant_projectionlist =  []
    for i in range(1,nlayers+1,1):
        begin   = float(clean_layup[i][0])
        end     = float(clean_layup[i][1])
        relevant_projectionlist.append(chop_interval_from_layup(tmp,begin,end))
        tmp = insert_interval_in_layup(tmp,begin,end)
        tmp = tmp[np.lexsort(np.fliplr(tmp).T)]
    
    return relevant_projectionlist
    
 
    
#==============================================================================       
if __name__ == '__main__':
    '''Executes the following code if the file is exceuted as the main file 
    directly. This is used in this context to test the fuctions from abobe. '''
    
    Layup = np.array([[  0.  ,   0.5  ,   0.23 , 45.  ,   1.  ],
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
    
    projectionlist=cummulated_layup_boundaries(Layup)
    relevant_projectionlist = relevant_cummulated_layup_boundaries(Layup)
        
    #=====================================================
    #               PLOT!!!!
    #=====================================================
    
    plt.figure(1)
    #============================= 
    plt.subplot(131)
    plt.xlim(0,1)
    plt.ylim(0,10)
    plt.title('Layup')    
    plt.xlabel('Coordinate')
    color=plt.cm.Paired(np.linspace(0,1,nlayers+4))
    
    for i,item in enumerate(clean_layup):
        timelines(i,item[0],item[1],color[int(item[2])])
        
    plt.subplot(132)
    plt.xlim(0,1)
    plt.ylim(0,10)
    plt.title('Cummulated Boundaries')    
    plt.xlabel('Coordinate')
    
    for i,idx in enumerate(projectionlist):
        for j,item in enumerate(idx):
            timelines(i+1,item[0],item[1],color[int(item[2])])
    
    plt.subplot(133)
    plt.xlim(0,1)
    plt.ylim(0,10)
    plt.title('Relevant Cummulated Boundaries')    
    plt.xlabel('Coordinate')
    
    for i,idx in enumerate(relevant_projectionlist):
        for j,item in enumerate(idx):
            timelines(i+1,item[0],item[1],color[int(item[2])])
    