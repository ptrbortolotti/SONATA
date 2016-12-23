import numpy as np
import matplotlib.pyplot as plt


def timelines(y, xstart, xstop, c='b'):
    """Plot timelines at y from xstart to xstop with given color."""   
    if xstart > xstop:
        plt.hlines(y,0,xstop, colors=c, lw=4)
        plt.hlines(y,xstart,1,colors=c,lw=4)
    else:
        plt.hlines(y, xstart, xstop,colors=c, lw=4)

def merge_intervals(intervals):
    """
    A simple algorithm can be used:
    1. Sort the intervals in increasing order
    2. Push the first interval on the stack
    3. Iterate through intervals and for each one compare current interval
       with the top of the stack and:
       A. If current interval does not overlap, push on to stack
       B. If current interval does overlap, merge both intervals in to one
          and push on to stack
    4. At the end return stack
    """
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []
    projection = []
    for i,higher in enumerate(sorted_by_lower_bound):
        if not merged:
            merged.append(higher)
            projection.append(higher)
        else:
            merged_lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                if not merged_lower[0] == higher[0]:
                    projection[-1] = (merged_lower[0],higher[0],merged_lower[2])
                projection.append((higher[0],higher[1],higher[2]))
                if merged_lower[1] > higher[1]:
                    projection.append((higher[1],merged_lower[1],merged_lower[2]))
                                     
                upper_bound = max(merged_lower[1], higher[1])
                merged[-1] = (merged_lower[0], upper_bound,higher[2])  # replace by merged interval
            else:
                merged.append(higher)
                projection.append(higher)
    return merged, projection

Layup = np.array([[  0.3  ,   0.6  ,   0.23 , 45.  ,   1.  ],
                  [  0.1   ,  0.4  ,   0.23 , 45.  ,   1.  ]])


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
        

merged_list,projection = merge_intervals(layup_intervals)


print("Original list of ranges: {}".format(layup_intervals))
#merged_list = merge_intervals(l)
print("List of ranges after merge_ranges: {}".format(merged_list))
print("projection: {}".format(projection))

sorted_by_lower_bound = sorted(layup_intervals, key=lambda tup: tup[0])






#=====================================================
#               PLOT!!!!
#=====================================================

plt.figure(1)
#============================= 
plt.subplot(211)
plt.xlim(0,1)
plt.ylim(0,i+2)
plt.ylabel('Layup')    
plt.xlabel('Coordinate')
color=plt.cm.Paired(np.linspace(0,1,6))

for i,item in enumerate(Layup):
    timelines(i+1,item[0],item[1],color[i])
 
for i,item in enumerate(projection):
    timelines(0,item[0],item[1],color[item[2]-1])
    
    
#=============================    
plt.subplot(212)
plt.xlim(0,1)
plt.ylim(0,i+2)
plt.ylabel('Sorted by Lower Bounds')    
plt.xlabel('Coordinate')
for i,item in enumerate(sorted_by_lower_bound):
    timelines(i+1,item[0],item[1],color[i])   

for i,item in enumerate(projection):
    timelines(0,item[0],item[1],color[item[2]-1])

plt.show()