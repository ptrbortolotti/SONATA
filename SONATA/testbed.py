import numpy as np
import matplotlib.pyplot as plt


def timelines(y, xstart, xstop, c='b'):
    """Plot timelines at y from xstart to xstop with given color."""   
    if xstart > xstop:
        plt.hlines(y,0,xstop, colors=c, lw=4)
        plt.hlines(y,xstart,1,colors=c,lw=4)
    else:
        plt.hlines(y, xstart, xstop,colors=c, lw=4)

Layup = np.array([[  0.  ,   1.  ,   0.23 , 45.  ,   1.  ],
                  [  0.   ,  1.  ,   0.23 , 45.  ,   1.  ],
                  [  0.3  ,  0.7  ,  0.6  , 45.  ,   2.  ],
                  [  0.05 ,  0.6 ,  0.6  , 45.  ,   2.  ],
                  [  0.85  , 0.2  ,  0.6  , 45.  ,   2.  ],
                  [  0.9 ,  0.4 ,  0.6  , 45.  ,   2.  ],
                  [  0.20 ,  0.532,  0.6  , 45.  ,   2.  ],
                  [  0.25 ,  0.3 ,  0.6  , 45.  ,   2.  ]],  dtype='float_')


Layup = np.array([[  0.1  ,   0.95  ,   0.23 , 45.  ,   1.  ],
                  [  0.5   ,  0.9  ,   0.23 , 45.  ,   1.  ],
                  [  0.0  ,   0.15  ,   0.23 , 45.  ,   1.  ],])

layup_intervals = []    
for i,item in enumerate(Layup):
    if item[1]<item[0]:
        tup = (0,float(item[1]),i+1)
        layup_intervals.append(tup)
        tup = (float(item[0]),1,i+1)
        layup_intervals.append(tup)
    else:
        tup = (float(item[0]),float(item[1]),i+1)
        layup_intervals.append(tup)
        


LayerNb = int(np.size(layup_intervals,0))
Projection = np.array([[layup_intervals[-1][0], layup_intervals[-1][1], LayerNb]]) #Start with the layer below
for i in range(LayerNb-1,0,-1):
    tmp_projection = np.array([Projection[-1]])
    
    for j,item in enumerate(Projection):
        s1 = float(tmp_projection[j][0])
        e1 = float(tmp_projection[j][1])
        s2 = float(layup_intervals[i-1][0])
        e2 = float(layup_intervals[i-1][1])
        
        #1. Possibility: No Overlap
        if s2>=e1 or e2<=s1:
            a = np.array([s2, e2, i])
            tmp_projection = np.vstack((tmp_projection, a))
            
        #2. Possibility: REAR Overlap
        elif (s2>=s1 and s2<=e1):
            a = np.array([e1, e2, i])
            tmp_projection = np.vstack((tmp_projection, a))
        
        #3. Possibility: FRONT Overlap
        elif (s2<=s1 and e2<=e1):
            a = np.array([s2, s1, i])
            tmp_projection = np.vstack((tmp_projection, a))
           
        #4. Possibility: Full Overlap
        elif (s2<=s1 and e2>=e1):
            a = np.array([s2, s1, i])
            b = np.array([e1, e2, i])
            tmp_projection = np.vstack((tmp_projection, a, b))
    
           
    
    #merged projection





#=====================================================
#               PLOT!!!!
#=====================================================

plt.figure(1)
#============================= 
plt.subplot(211)


color=plt.cm.Paired(np.linspace(0,1,10))

for i,item in enumerate(layup_intervals):
    timelines(i+1,item[0],item[1],color[int(item[2])])
 
for i,item in enumerate(Projection):
    timelines(0,item[0],item[1],color[int(item[2])])
    
plt.xlim(0,1)
plt.ylim(-0.1,i+2)
plt.ylabel('Layup')    
plt.xlabel('Coordinate')
plt.show()