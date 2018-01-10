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
    '''Truthfully, I currently don'tknow the exact itend of this funtion.'''
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

def insert_interval_in_layup(layup,begin,end,**kw):
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
    
    #KWARGS:
    if kw.get('value') !=  None:
        idx = kw.get('value')
    else:
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


def sort_layup_projection(Projection):
    '''The sort_layup_projection fuctions sorts the interval structure.
        
    Args:         
        Projection: list of numpy array of the intervaltree structure of the form:
                  start end layer#
        np.array([[ 0.3  0.5  3. ]
                  [ 0.7  1.   2. ]
                  [ 0.   0.3  2. ]])    
            
    Returns: 
        Projection: list of sorted numpy array of the intervaltree structure of the form:
                  start end layer#
        np.array([[ 0.7  1.   2. ]
                  [ 0.   0.3  2. ]
                  [ 0.3  0.5  3. ]])               
    '''
    sorted_Projection=[]
    for i,proj in enumerate(Projection):
        a =  proj[:,:2].flatten()
        unique, counts = np.unique(a, return_counts=True)
        iv_boundaries = np.sort(unique[np.where(counts == 1)])
        b = proj[proj[:,0].argsort()]
        #regular interval:
        if np.all(a>=iv_boundaries[0]) and np.all(a<=iv_boundaries[1]):
            tmp = b
        #regular interval:  
        else:
            splitter =  np.asscalar(np.where(b[:,0]==iv_boundaries[2])[0])
            c = [b[0:splitter,:],b[splitter::,:]]
            #print i,"Split:",b, '@:', splitter, '\n C:', c
            tmp = np.vstack((c[1],c[0]))
        
        sorted_Projection.append(tmp)  
    
    return sorted_Projection


def cummulated_layup_boundaries(Layup):
    '''The cummulated_layup_boundaries fuctions generates a interval structure 
    from the Layup definitions so that for each layer the information of the 
    underlying layers that create and interval between zero and one is given.
    This information can be used to determine the offset origin for the layer. 
    the relevant_cummulated_layup_boundaries function chops this valuable 
    information to the interval of the Layup itself to determin the offset
    origin in the given interval.
        
    Args:         
        Layup = np.array([[  0.  ,   0.5  ,   0.23 , 45.  ,   1.  ],
                  [  0.   ,  1.  ,   0.23 , 45.  ,   1.  ],
                  [  0.20 ,  0.532,  0.6  , 45.  ,   2.  ],
                  [  0.25 ,  0.3 ,  0.6  , 45.  ,   2.  ]],)
                  #Start[-]	End[-] thickness [mm] Orientation [deg] MatID
            
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
    for i in range(1,nlayers+1,1):
        begin   = float(clean_layup[i][0])
        end     = float(clean_layup[i][1])  
        tmp = insert_interval_in_layup(tmp,begin,end)
        tmp = tmp[np.lexsort(np.fliplr(tmp).T)]
        projectionlist.append(tmp)

    return sort_layup_projection(projectionlist)
        

def relevant_cummulated_layup_boundaries(Layup):
    '''The relevant_cummulated_layup_boundaries function chops this valuable 
    information to the interval of the Layup itself to determin the offset
    origin in the given interval.        
    
    Args:        
        Layup = np.array([[  0.  ,   0.5  ,   0.23 , 45.  ,   1.  ],
                  [  0.   ,  1.  ,   0.23 , 45.  ,   1.  ],
                  [  0.20 ,  0.532,  0.6  , 45.  ,   2.  ],
                  [  0.25 ,  0.3 ,  0.6  , 45.  ,   2.  ]],)
                  #Start[-]	End[-] thickness [mm] Orientation [deg] MatID
            
        Returns: Projection: list of numpy arrays of the intervaltree structure 
            of the form:
                          start end layer#
                  np.array([[ 0.3  0.7  3. ]
                          [ 0.7  1.   2. ]
                          [ 0.   0.3  2. ]])                
    '''
    if Layup.size == 0:
        return []    
    
    else:
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

        return sort_layup_projection(relevant_projectionlist)


def inverse_relevant_cummulated_layup_boundaries(Layup):
    ''' The function inverse_relevant_cummulated_layup_boundaries is helpful to
    determine the a_nodes during the meshing process.
    It basically uses the relevant_cummulated_layup_boundaries with a flipped
    Layup definition and assigns the correct layer number.
    
    Args:         
        Layup = np.array([[  0.  ,   0.5  ,   0.23 , 45.  ,   1.  ],
                  [  0.   ,  1.  ,   0.23 , 45.  ,   1.  ],
                  [  0.20 ,  0.532,  0.6  , 45.  ,   2.  ],
                  [  0.25 ,  0.3 ,  0.6  , 45.  ,   2.  ]],)
                  #Start[-]	End[-] thickness [mm] Orientation [deg] MatID
            
    Returns: Projection: list of numpy arrays of the intervaltree structure 
            of the form:
                      start end layer#
            np.array([[ 0.3  0.7  3. ]
                      [ 0.7  1.   2. ]
                      [ 0.   0.3  2. ]])                
    '''              
    flipped_Projection = relevant_cummulated_layup_boundaries(np.flipud(Layup))
    b = np.flipud(np.arange(1,len(Layup)+1))
    d = dict(enumerate(b))
    
    inverse_projection = []
    for iv_tree in flipped_Projection:
        part1 = iv_tree[:,:2]
        part2 = np.asarray([[d[l] for l in iv_tree[:,2]]])
        inverse_projection.append(np.concatenate((part1, part2.T), axis=1))

    inverse_projection = reversed(inverse_projection)
    return sort_layup_projection(inverse_projection)



def plot_layup_projection(Layup):
    ''' the function plot_layup_projection uses the functions from above and 
    plots the Layup as Timelines. Both for Layup and for the inverse Layup.
    
    the left column shows the Layup as it is defined. In both the regular and 
    the inverse representation a light blue bar is displayed to show the 
    Segment0 intervals respectivally unmeshed intervals.
    
    In the middle shows the relenvant part of the cummulated boundaries by 
    chopping them into the their relevant part . On the top the final Segment
    boundary is displayed as a cummulated interval between 0 and 1. 
    
    The right plot shows the inverse relevant cummulated boundaries that are 
    needed to determine the a_nodes during the meshing procedure. The color 
    grey marks the unmeshed intervals that need to me discretized by the the 
    equidistant function.
    
    Args:         
        Layup = np.array([[  0.  ,   0.5  ,   0.23 , 45.  ,   1.  ],
                  [  0.   ,  1.  ,   0.23 , 45.  ,   1.  ],
                  [  0.20 ,  0.532,  0.6  , 45.  ,   2.  ],
                  [  0.25 ,  0.3 ,  0.6  , 45.  ,   2.  ]],)
                  #Start[-]	End[-] thickness [mm] Orientation [deg] MatID
            
    Returns: None              
    '''
    relevant_projectionlist = relevant_cummulated_layup_boundaries(Layup)
    inverse_projection = inverse_relevant_cummulated_layup_boundaries(Layup)
    

    #=======PLOT=========
    nlayers = int(np.size(Layup,0))
    a = np.delete(Layup,[2,3,4],1)
    b = np.transpose(np.linspace(1,nlayers,nlayers))
    clean_layup = np.c_[a,b]

    plt.figure(3)
    
    #=======PLOT LAYUP=========  
    plt.subplot(131)
    plt.xlim(0,1)
    plt.ylim(-0.5,nlayers+2)
    plt.yticks([])
    plt.title('Layup')    
    plt.xlabel('Coordinate')

    color=plt.cm.Paired(np.linspace(0,1,nlayers+1))
    color[0] = [0, 0, 0, 1]
    color = np.vstack((color,np.array([0.7, 0.7, 0.7, 1])))
    
    plt.hlines(0, 0, 1, lw=4)
    plt.text(-0.05, 0, "Segment Boundary", fontsize=10, horizontalalignment='right')
    for i,item in enumerate(clean_layup):
        timelines(i+1,item[0],item[1],color[int(item[2])])
        string = "Layer "+str(i+1)
        plt.text(-0.05, i+1, string, fontsize=10, horizontalalignment='right')
    plt.text(-0.05, i+2, "Final Segment Boundary", fontsize=10, horizontalalignment='right')    
    
    
    #=======PLOT Relevant Cummulated Boundaries=========  
    plt.subplot(132)
    plt.xlim(0,1)
    plt.ylim(-0.5,nlayers+2)
    plt.yticks([])
    plt.title('Relevant Cummulated Boundaries\n(for layer creation)')    
    plt.xlabel('Coordinate')
    
    for i,idx in enumerate(relevant_projectionlist):
        for j,item in enumerate(idx):
            timelines(i+1,item[0],item[1],color[int(item[2])])
    
    final_ivlist=cummulated_layup_boundaries(Layup)[-1]
    for j,item in enumerate(final_ivlist):
        timelines(i+2,item[0],item[1],color[int(item[2])])

            
    #=======PLOT Inverse Relevant Cummulated Boundaries=========  
    plt.subplot(133)
    plt.xlim(0,1)
    plt.ylim(-0.5,nlayers+2)
    plt.yticks([])
    plt.title('Inverse Relevant Cummulated Boundaries\n(to determine a_nodes during meshing)')    
    plt.xlabel('Coordinate')
    
    for i,idx in enumerate(inverse_projection):
        for j,item in enumerate(idx):
            timelines(i+1,item[0],item[1],color[int(item[2])+1])
    
    plt.text(item[1], i+1, "Grey: unmeshed layer intervals", fontsize=8) 
    
    
    return None
    

#==============================================================================       
if __name__ == '__main__':
    '''Executes the following code if the file is exceuted as the main file 
    directly. This is used in this context to test the fuctions from abobe. '''
    
    plt.close('all')
    Layup = np.array([[  0.  ,   0.5  ,   0.23 , 45.  ,   1.  ],
                  [  0.   ,  1.  ,   0.23 , 45.  ,   1.  ],
                  [  0.3  ,  0.7  ,  0.6  , 45.  ,   2.  ],
                  [  0.05 ,  0.6 ,  0.6  , 45.  ,   2.  ],
                  [  0.85  , 0.2  ,  0.6  , 45.  ,   2.  ],
                  [  0.9 ,  0.4 ,  0.6  , 45.  ,   2.  ],
                  [  0.4 ,  0.6,  0.6  , 45.  ,   2.  ],
                  [  0.3 ,  0.7,  0.6  , 45.  ,   2.  ],
                  [  0.20 ,  0.532,  0.6  , 45.  ,   2.  ],
                  [  0.25 ,  0.3 ,  0.6  , 45.  ,   2.  ]],)
    

    plot_layup_projection(Layup)
    relevant_projectionlist = relevant_cummulated_layup_boundaries(Layup)
    inverse_projection = inverse_relevant_cummulated_layup_boundaries(Layup)