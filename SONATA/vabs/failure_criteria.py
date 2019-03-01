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


Returns
----------
(sf, mode) -- safety factors(float), failure mode(string)
    

"""

def von_Mises(mat, stressM, strainM):
    """calculates the safety factor for the von Mises yield criterion with the
    material ultimate strenght"
    """
    sf = mat.UTS/stressM.sigma_vM
    mode = 'von_Mises'
    return (sf, mode)


def tsaiwu_2D(mat, stressM, strainM):
    """ Calculates SF and mode according to Tsai-Wu criterion (layer-wise). """

    Xt = mat.Xt
    Xc = mat.Xc
    Yt = mat.Yt
    Yc = mat.Yc
    S21 = mat.S21
    
    sig1 = stressM.sigma11
    sig2 = stressM.sigma22
    tau = stressM.sigma12

    f11 = 1/(Xt*Xc)
    f22 = 1/(Yt*Yc)
    f12 = -1/(2*(Xt*Xc*Yt*Yc)**(1/2))
    f66 = 1/(S21**2)
    f1 = 1/Xt - 1/Xc
    f2 = 1/Yt - 1/Yc

    a = f11*sig1**2 + f22*sig2**2 + f66*tau**2 + 2*f12*sig1*sig2
    b = f1*sig1 + f2*sig2;

    sf = (-b + (b**2 + 4*a)**(1/2))/(2*a)

    # Failure mode calculations  
    H1 = abs(f1*sig1 + f11*sig1**2)
    H2 = abs(f2*sig2 + f22*sig2**2)
    H6 = abs(f66*tau**2)
 
    if max(H1,H2,H6) == H1:
        mode = "fiber"        # fiber failure
    elif max(H1,H2,H6) == H2:
        mode = "matrix"        # matrix failure
    else:
        mode = "shear"        # shear failure
    
    # Returns SF & mode    
    return (sf, mode)


def maxstress_2D(mat, stressM, strainM):
    """ Calc. SF and mode according to Max. Stress criterion (layer-wise). """
    
    Xt = mat.Xt
    Xc = mat.Xc
    Yt = mat.Yt
    Yc = mat.Yc
    S21 = mat.S21
    
    sig1 = stressM.sigma11
    sig2 = stressM.sigma22
    tau = stressM.sigma12

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
    return (sf, mode)


def maxstrain_2D(mat, stressM, strainM):
    """ Calc. SF and mode according to Max. Strain criterion (layer-wise). """

    strainXt = mat.Xt / mat.E[0]
    strainXc = mat.Xc / mat.E[0]
    strainYt = mat.Yt / mat.E[1]
    strainYc = mat.Yc / mat.E[1]
    strainS21 = mat.S21 / mat.G[0]
    
    eps1 = strainM.epsilon11
    eps2 = strainM.epsilon22
    gamma = strainM.gamma12

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
    return (sf, mode)


def hashin_2D(mat, stressM, strainM):
    """ Calc. SF and mode according to Hashin criterion (layer-wise). """

    Xt = mat.Xt
    Xc = mat.Xc
    Yt = mat.Yt
    Yc = mat.Yc
    S21 = mat.S21
    
    sig1 = stressM.sigma11
    sig2 = stressM.sigma22
    tau = stressM.sigma12

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
    return (sf, mode)


if __name__ == '__main__':
    pass
    
