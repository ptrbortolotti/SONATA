import numpy as np
import matplotlib.pyplot as plt



def timelines(y, xstart, xstop, c='b'):
    """Plot timelines at y from xstart to xstop with given color."""   
    if xstart > xstop:
        plt.hlines(y,0,xstop, colors=c, lw=4)
        plt.hlines(y,xstart,1,colors=c,lw=4)
    else:
        plt.hlines(y, xstart, xstop,colors=c, lw=4)
    
    #plt.vlines(xstart, y+0.03, y-0.03,  lw=2)
    #plt.vlines(xstop, y+0.03, y-0.03,  lw=2)


Layup = np.array([[  0.  ,   1.  ,   0.23 , 45.  ,   1.  ],
                  [  0.   ,  1.  ,   0.23 , 45.  ,   1.  ],
                  [  0.3  ,  0.7  ,  0.6  , 45.  ,   2.  ],
                  [  0.05 ,  0.6 ,  0.6  , 45.  ,   2.  ],
                  [  0.85  , 0.2  ,  0.6  , 45.  ,   2.  ],
                  [  0.9 ,  0.4 ,  0.6  , 45.  ,   2.  ],
                  [  0.20 ,  0.532,  0.6  , 45.  ,   2.  ],
                  [  0.25 ,  0.3 ,  0.6  , 45.  ,   2.  ]],  dtype='float_')



#DETERMINE BOUNDING DEF FOR LAYER 9 
LayerNb = int(np.size(Layup,0))
Projection =  np.array([Layup[7][0], Layup[7][1], 7]) #Start with the layer below
for i in range(LayerNb-1,-1,-1):
    print i
    s1 = Layup[i][0]
    e1 = Layup[i][1]
    s2 = Layup[i-1][0]
    e2 = Layup[i-1][1]
    
    
    if (s1<e1 and s2<e2):
        #1. Possibility: No Overlap
        if s2>=e1:
            a = np.array([s2, e2, i-1])
            Projection = np.vstack((Projection, a))
            
        #2. Possibility: REAR Overlap
        elif (s2>=s1 and s2<=e1):
            a = np.array([e1, e2, i-1])
            Projection = np.vstack((Projection, a))
        
        #3. Possibility: FRONT Overlap
        elif (s2<=s1 and e2<=e1):
            a = np.array([s2, s1, i-1])
            Projection = np.vstack((Projection, a))
           
        #4. Possibility: Full Overlap
        elif (s2<=s1 and e2>=e1):
            a = np.array([s2, s1, i-1])
            b = np.array([e1, e2, i-1])
            Projection = np.vstack((Projection, a, b))

    elif (s1<e1 and s2>e2):    
        #1. Possibility: No Overlap 
        if s2>=e1 and e2<=s1:
            a = np.array([s2, 1, i-1])
            b = np.array([0, e2, i-1])
            Projection = np.vstack((Projection, a, b))
            
        
        
        
        
        
        
        
    start = Layup[i][0],
    end = Layup[i][1]
    a = np.array([Layup[i][0], Layup[i][1], i])

    








plt.figure(1)
color=iter(plt.cm.Paired(np.linspace(0,1,10)))
for i,item in enumerate(Layup):
    c=next(color)
    timelines(i+1,item[0],item[1],c)
    
plt.xlim(0,1)
plt.ylim(0,i+2)
plt.ylabel('Layer')    
plt.xlabel('Coordinate')
plt.show()


