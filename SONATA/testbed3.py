# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 15:37:44 2016

@author: TPflumm
"""
import numpy as np       
import shapely.geometry as shp
import matplotlib.pyplot as plt
from utils import calc_DCT_angles, P2Pdistance, isclose,unique_rows
from readinput import UIUCAirfoil2d,AirfoilDat2d

arrPts = AirfoilDat2d('hexagon.dat').T                      
arrPts = np.multiply(arrPts,100)
dist = 0.6
join_style=1
#OFFSET ALGORITHM
#join_style = 1#( 1:round,2:mitre,3:bevels)
side = 'left'    
res = 16 #resolution







if P2Pdistance(arrPts[0],arrPts[-1])<=1e-6:
    closed = True
    afpoly = shp.Polygon(arrPts)
    noffafpoly = afpoly.buffer(-dist)  # Inward offset
    data = np.array(noffafpoly.exterior)

else:
    closed = False
    line = shp.LineString(arrPts)
    offset = line.parallel_offset(dist,side,res,join_style)

    if isinstance(offset,shp.MultiLineString):    
        parts = hasattr(offset, 'geoms') and offset or [offset]
        for part in parts:
            if part.is_closed:
                x, y= part.xy
                data = np.vstack((x,y))
                data = data.T
            else:
                print "ERROR: \t A multilinestring has been created by shp_parallel_offset that is not closed!"
                
    elif isinstance(offset,shp.LineString):                 
        data = np.array(offset.coords)
    
    else:
        data = np.array(offset.coords)

shp_offset = data
#Interpolate Large spaces! 
seg_P2Plength = []
cumm_length = 0
Resolution = 100

for j in range(0,len(data)-1):
    seg_P2Plength.append(P2Pdistance(data[j],data[j+1]))
    cumm_length += P2Pdistance(data[j],data[j+1]) 
    
#Check if Refinement is necessary:
if len(seg_P2Plength) > 0 and max(seg_P2Plength) > cumm_length/Resolution :
    Refinement = True
else:
    Refinement = False

while Refinement == True:
    temp_data = []
    for i in range(0,len(data)-1):
        if P2Pdistance(data[i],data[i+1]) > (cumm_length/Resolution):
            p0 = data[i]
            p1 = data[i+1]  
            v1 = p1-p0
            p05 = p0+v1/2
            temp_data.append(p0)
            temp_data.append(p05)
        else:
            temp_data.append(data[i])
            
    temp_data.append(data[-1])        
    data = np.vstack(temp_data)        
     
    #Check if further Refinement is necessary
    seg_P2Plength = []
    cumm_length = 0
    for j in range(0,len(data)-1):
        seg_P2Plength.append(P2Pdistance(data[j],data[j+1]))
        cumm_length += P2Pdistance(data[j],data[j+1]) 
     
    if max(seg_P2Plength) > cumm_length/Resolution:
        Refinement = True
    else:
        Refinement = False   




##Make sure there are unique rows
#if np.allclose(data[0],data[-1]):
#    closed = True
#    
#else:
#    closed = False
#    before = len(data)
#    data = unique_rows(data.round(decimals=9))
#    after = len(data)
#    if after<before:
#        print 'non unique elements found',len(data)







fig = plt.figure(1)
ax = fig.add_subplot(111)
plt.clf()         
plt.plot(*arrPts.T, color='black', marker='.')

#    for i, item in enumerate(data):
#            plt.annotate(i, (item[0],item[1]), color='red')

plt.plot(*arrPts[0].T, color='green', marker='o')
plt.plot(*arrPts[-1].T, color='purple', marker='o')
plt.plot(*shp_offset.T, color='blue', marker='.')
plt.plot(*data.T, color='red', marker='.')
plt.plot(*data[0].T, color='green', marker='x')
plt.plot(*data[-1].T, color='purple', marker='x')
plt.axis('equal')  
plt.show()   
    

    
    
