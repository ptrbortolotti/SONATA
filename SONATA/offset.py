# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 15:37:44 2016

@author: TPflumm
"""
import numpy as np       
import shapely.geometry as shp
import matplotlib.pyplot as plt
from utils import calc_DCT_angles, P2Pdistance


def shp_parallel_offset(arrPts,dist):
#OFFSET ALGORITHM
    side = 'left'    
    line = shp.LineString(arrPts)
    res = 16 #resolution
    join_style = 1#( 1:round,2:mitre,3:bevels)
    offset = line.parallel_offset(dist,side,res,join_style)
    #print type(offset)
    
    if isinstance(offset,shp.MultiLineString):    
        parts = hasattr(offset, 'geoms') and offset or [offset]
        for part in parts:
            if part.is_closed:
                x, y= part.xy
                data = np.vstack((x,y))
                data = data.T
                
    elif isinstance(offset,shp.LineString):                 
        data = np.array(offset.coords)
    
    
    else:
        data = np.array(offset.coords)
    
     #Interpolate Large spaces! 
    seg_P2Plength = []
    cumm_length = 0
    Resolution = 400
    
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

    return data


               
#    plt.figure(2)        
#    plt.plot(*arrPts.T, color='black', marker='.')
#    plt.plot(*data.T, color='red', marker='.')        
    

    
    
    
#    #Find corners and edges of original data   
#    DCT_angles1 = calc_DCT_angles(arrPts)
#    angular_deflection = 30
#    concave1 = []
#    for i in range(0,DCT_angles1.shape[0]):
#        if DCT_angles1[i] > (180 + angular_deflection): 
#            concave1.append(i)
#    NbConcave1 = np.size(concave1)
#    
#    for i,idx in enumerate(concave1):
#        ConcavePt = arrPts[idx]
#        v1 = arrPts[idx]-arrPts[idx-1]
#        n1 = np.array([-v1[1],v1[0]])
#        n1 = n1/np.linalg.norm(n1)
#        P1 = ConcavePt+n1*distt
#        v2 = arrPts[idx+1]-arrPts[idx]
#        n2 = np.array([-v2[1],v2[0]])
#        n2 = n2/np.linalg.norm(n2)
#        P2 = ConcavePt+n2*dist
#        
#        plt.plot(*P1.T, color='GREEN', marker='o')     
#        plt.plot(*P2.T, color='BLUE', marker='o')    
        
#    plt.axis('equal')
#    plt.show()    
#        