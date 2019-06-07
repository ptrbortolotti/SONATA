#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 09:48:14 2019

@author: gu32kij
"""
import numpy as np
from OCC.gp import gp_Pnt, gp_Pnt2d, gp_Trsf, gp_Vec2d, gp_Ax3, gp_Ax1, gp_Dir, gp_Ax2


def trsf_blfr_to_cbm(Ax2_blfr, Ax2_befr):
    """
    generates the transform from the SONATA reference frame to the local 2d cbm 
    frame
    
    Parameters
    --------
    Ax2_blfr : gp_Ax2
        Opencascade blade reference frame, 
        Describes a right-handed coordinate system in 3D space.
    Ax2_befr : gp_Ax2
        Opencascade locale beam reference frame, 
        Opencascade: Describes a right-handed coordinate system in 3D space.
    
    Returns
    --------
    Trsf : gp_Trsf
        Opencascade: non-persistent transformation in 3D space
    
    """
    Ax2_cbm = gp_Ax2(Ax2_befr.Location(), Ax2_befr.XDirection(), Ax2_befr.YDirection())
    #print('Ax2_cbm.Direction:', Ax2_cbm.XDirection().Coord(), Ax2_cbm.YDirection().Coord())
    trsf1 = gp_Trsf()
    trsf2 = gp_Trsf()

    trsf1.SetTransformation(gp_Ax3(Ax2_blfr),gp_Ax3(Ax2_befr))
    trsf2.SetTransformation(gp_Ax3(Ax2_befr),gp_Ax3(Ax2_cbm))

    Trsf = gp_Trsf()
    Trsf.Multiply(trsf2)
    Trsf.Multiply(trsf1)

    return Trsf


def trsf_cbm_to_blfr(Ax2_blfr, Ax2_befr):
    """
    generates the transform from the cbm local 2d frame to the blade frame
    
    Parameters
    --------
    Ax2_blfr : gp_Ax2
        Opencascade blade reference frame, 
        Describes a right-handed coordinate system in 3D space.
    Ax2_befr : gp_Ax2
        Opencascade locale beam reference frame, 
        Opencascade: Describes a right-handed coordinate system in 3D space.
    
    Returns
    --------
    Trsf : gp_Trsf
        Opencascade: non-persistent transformation in 3D space
    
    """
    Ax2_cbm = gp_Ax2(Ax2_befr.Location(), Ax2_befr.XDirection(), Ax2_befr.YDirection())
    #print('Ax2_cbm.Direction:', Ax2_cbm.XDirection().Coord(), Ax2_cbm.YDirection().Coord())

    trsf1 = gp_Trsf()
    trsf2 = gp_Trsf()

    trsf1.SetTransformation(gp_Ax3(Ax2_cbm),gp_Ax3(Ax2_befr))
    trsf2.SetTransformation(gp_Ax3(Ax2_befr),gp_Ax3(Ax2_blfr))

    Trsf = gp_Trsf()
    Trsf.Multiply(trsf2)
    Trsf.Multiply(trsf1)

    return Trsf

if __name__ == '__main__':
    Ax2_1 = gp_Ax2()
    Ax2_2 = gp_Ax2()
    
    trfs1 = trsf_cbm_to_blfr(Ax2_1, Ax2_2)
    trfs2 = trsf_blfr_to_cbm(Ax2_1, Ax2_2)