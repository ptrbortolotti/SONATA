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
    afline = shp.LineString(arrPts)
    res = 16 #resolution
    join_style = 1 #( 1:round,2:mitre,3:bevels)
    offline = afline.parallel_offset(dist,side,res,join_style)
    offlinepts = np.array(offline.coords)
    #plt.plot(*offlinepts.T, color='red', marker='.')
    return offlinepts
