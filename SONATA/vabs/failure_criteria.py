# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 11:55:39 2018

@author: TPflumm

adapted from 
lamipy project - laminated composites calculations in Python.
failurecriteria.py - Module containing functions for calculating
					 safety factors according to different criteria.
                  Joao Paulo Bernhardt - September 2017   
                     
INPUTS:
mat -- instance of Material Class
stressM
strainM


OUTPUTS:
(sf, mode) -- safety factors(float), failure mode(string)

"""
import math


def tsaiwu_2D(mat, stressM, strainM):
    """ Calculates SF and mode according to Tsai-Wu criterion (layer-wise). """
    if mat.orth !=1:
        #print('WARNING: tsaiwu_2D criterion is not plausible for non orthotropic materials')
        pass

    f11 = 1/(mat.Xt*mat.Xc)
    f22 = 1/(mat.Yt*mat.Yc)
    f1 = 1/mat.Xt - 1/mat.Xc
    f2 = 1/mat.Yt - 1/mat.Yc
    f12 = -1/2*math.sqrt(f11*f22)
    f66 = 1/(mat.S21**2)

    a = f11*stressM.sigma11**2 + f22*stressM.sigma22**2 + f66*stressM.sigma12**2 + 2*f12*stressM.sigma11*stressM.sigma22
    b = f1*stressM.sigma11 + f2*stressM.sigma22;

    sf = (-b + (b**2 + 4*a)**(1/2))/(2*a)

    # Failure mode calculations  
    H1 = abs(f1*stressM.sigma11 + f11*stressM.sigma11**2)
    H2 = abs(f2*stressM.sigma22 + f22*stressM.sigma22**2)
    H6 = abs(f66*stressM.sigma12**2)
 
    if max(H1,H2,H6) == H1:
        mode = "fiber"        # fiber failure
    elif max(H1,H2,H6) == H2:
        mode = "matrix"        # matrix failure
    else:
        mode = "shear"        # shear failure
        
    return (sf, mode)       # Returns SF & mode    


def maxstress_2D(mat, stressM, strainM):
    """ Calc. SF and mode according to Max. Stress criterion (layer-wise). """

    Xt = mat_prop["Xt"]
    Xc = mat_prop["Xc"]
    Yt = mat_prop["Yt"]
    Yc = mat_prop["Yc"]
    S21 = mat_prop["S12"]

    # Verify for sig1
    if sig1 > 0:
        f_1 = (sig1/Xt)
    else:
        f_1 = (sig1/-Xc)

    # Verify for sig2
    if sig2 > 0:
         f_2 = (sig2/Yt)
    else:
        f_2 = (sig2/-Yc)

    # Verify for shear
    f_s = abs(tau)/S21

    f_max = max(f_1, f_2, f_s)

    # Find failure mode
    if f_max == f_1: 
    	mode = "fiber"
    elif f_max == f_2:
    	mode = "matrix"
    else:
    	mode = "shear"

    sf = 1/f_max

    # Result FS (1 / maximum of the 3 above)
    return [sf, mode]


def maxstrain_2D(mat, stressM, strainM):
    """ Calc. SF and mode according to Max. Strain criterion (layer-wise). """

    strainXt = mat_prop["Xt"] / mat_prop["E1"]
    strainXc = mat_prop["Xc"] / mat_prop["E1"]
    strainYt = mat_prop["Yt"] / mat_prop["E2"]
    strainYc = mat_prop["Yc"] / mat_prop["E2"]
    strainS21 = mat_prop["S12"] / mat_prop["G12"]

    # Verify for eps1
    if eps1 > 0:
        f_1 = (eps1/strainXt)
    else:
        f_1 = (eps1/-strainXc)

    # Verify for eps2
    if eps2 > 0:
         f_2 = (eps2/strainYt)
    else:
        f_2 = (eps2/-strainYc)

    # Verify for gamma
    f_s = abs(gamma)/strainS21

    f_max = max(f_1, f_2, f_s)

    # Find failure mode
    if f_max == f_1: 
    	mode = "fiber"
    elif f_max == f_2:
    	mode = "matrix"
    else:
    	mode = "shear"

    sf = 1/f_max

    # Result FS (1 / maximum of the 3 above)
    return [sf, mode]


def hashin_2D(mat, stressM, strainM):
    """ Calc. SF and mode according to Hashin criterion (layer-wise). """

    Xt = mat_prop["Xt"]
    Xc = mat_prop["Xc"]
    Yt = mat_prop["Yt"]
    Yc = mat_prop["Yc"]
    S21 = mat_prop["S12"]

    # Verify for sig1
    if sig1 >= 0:
        f_1 = (sig1/Xt)
    else:
        f_1 = -(sig1/Xc)

    # Verify for sig2
    if sig2 >= 0:
        f_2 = ((sig2/Yt)**2 + (tau/S21)**2)**0.5
    else:
        f_2 = ((sig2/Yc)**2 + (tau/S21)**2)**0.5

    f_max = max(f_1, f_2)

    if f_max == f_1:
    	mode = "fiber"
    else:
    	mode = "matrix"

    sf = 1/f_max

    # Result FS (1 / maximum of the 3 above)
    return [sf, mode]

def puck_2d(mat, stressM, strainM):
    pass


if __name__ == '__main__':
    pass
    
