#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 16:10:35 2019

@author: gu32kij
"""
from OCC.gp import gp_Ax2, gp_Pnt, gp_Dir

class Component(gp_Ax2):
    '''Describes a right-handed coordinate system in 3D space. 
    Ax2 is an instance from the opencascade gp_Ax2 class.'''
    
    __slots__ = ('name')    
    def __init__(self, *args, **kwargs):
        super().__init__(*args)
        self.name = kwargs.get('name')
    
    def __repr__(self):
        """__repr__ is the built-in function used to compute the "official" 
        string reputation of an object, """
        return 'Component: '+ self.name
    
    def display_Ax2(self):
        pass
    
if __name__ == '__main__':
    C = Component(gp_Pnt(0,10,0),gp_Dir(1,0,0),name='TestComponent')
    