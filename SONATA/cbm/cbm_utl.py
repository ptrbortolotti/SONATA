#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 14:27:46 2019

@author: gu32kij
"""
# Third party modules
import numpy as np


def trsf_sixbysix(M, T):
    """
    Transform six-by-six compliance/stiffness matrix. 
    change of reference frame in engineering (or Voigt) notation.
    
    Parameters
    ----------
    M : np.ndarray
        6x6 Siffness or Mass Matrix
    T : np.ndarray
        Transformation Matrix
        
    Returns
    ----------
    res : np.ndarray
        Transformed 6x6 matrix
    """

    TS_1 = np.dot(np.dot(T.T, M[0:3, 0:3]), T)
    TS_2 = np.dot(np.dot(T.T, M[3:6, 0:3]), T)
    TS_3 = np.dot(np.dot(T.T, M[0:3, 3:6]), T)
    TS_4 = np.dot(np.dot(T.T, M[3:6, 3:6]), T)

    tmp_1 = np.vstack((TS_1, TS_2))
    tmp_2 = np.vstack((TS_3, TS_4))
    res = np.hstack((tmp_1, tmp_2))
    return res


def trsf_coords(coords, T):
    return np.dot(T, coords)


def dymore2sixbysix(arr):
    """
    Transforms the 29 elemet Dymore Beamproperties array back to 6x6 Mass and 
    6x6 TS Matrix, mu and eta
    
    Parameters
    --------
    arr : ndarray
        [Massterms(6) (m00, mEta2, mEta3, m33, m23, m22) 
        Stiffness(21) (k11, k12, k22, k13, k23, k33,... k16, k26, ...k66)
        Viscous Damping(1) mu, Curvilinear coordinate(1) eta]
        
    Returns
    --------
    tuple of (MM 6x6MassMatrix, TS 6x6 StiffnessMatrix, Viscous Damping mu, Curvilinear coordinate eta)

    """

    mu = arr[0]
    muXm2 = arr[1]
    muXm3 = arr[2]
    i33 = arr[3]
    i23 = arr[4]
    i22 = arr[5]
    
    MM = np.array([[mu  , 0   ,   0,      0, muXm3, -muXm2],
                   [0   , mu  ,   0, -muXm3,     0,      0],
                   [0   , 0  ,   mu,  muXm2,     0,      0],
                   [0   , -muXm3,   muXm2,  i22+i33, 0,  0],
                   [muXm3, 0  ,   0,  0,    i22,       i23],
                   [-muXm2, 0  ,  0,  0,     i23,     i33]])
    

    a = np.zeros((6, 6))
    for c, (i, j) in enumerate(zip(np.tril_indices(6)[1], np.tril_indices(6)[0])):
        # print(c,i,j)
        k = c + 6
        a[i, j] = arr[k]
    b = a.copy().T
    b[np.diag_indices(6)] = 0
    TS = a + b

    mu = arr[-2]
    eta = arr[-1]
    return (MM, TS, mu, eta)
