# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 15:37:44 2016

@author: TPflumm
"""
import numpy as np       
import shapely.geometry as shp


def shp_parallel_offset(arrPts,dist):
#OFFSET ALGORITHM
    side = 'left'    
    line = shp.LineString(arrPts)
    res = 16 #resolution
    join_style = 1 #( 1:round,2:mitre,3:bevels)
    offset = line.parallel_offset(dist,side,res,join_style)
    
    if isinstance(offset,shp.MultiLineString):    
        parts = hasattr(offset, 'geoms') and offset or [offset]
        for part in parts:
            if part.is_closed:
                x, y= part.xy
                data = np.vstack((x,y))
                data = data.T
       
    elif isinstance(offset,shp.LineString):                 
        data = np.array(offset.coords)
    #plt.plot(*offlinepts.T, color='red', marker='.')
    
    else:
        data = np.array(offset.coords)
    print type(offset)
    
    return data

    
'''
Multiple Part Offset!!!

data = [] 
parts = hasattr(offset, 'geoms') and offset or [offset]
for part in parts:
    if part.is_closed:
        x, y= part.xy
        data = np.vstack((x,y))
        data = data.T

'''