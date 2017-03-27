# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 13:30:28 2017

@author: TPflumm
"""

class Stress(object):
    #Object to represent the STRAIN Tensor within an Element    

    def __init__(self,Vec=None):
        if Vec is None:
            pass
        elif Vec is not None:
            self.Vector = Vec

    @property
    def sigma11(self):
        return self.Vector[0]
    
    @property
    def sigma12(self):
        return self.Vector[1]
    
    @property
    def sigma13(self):
        return self.Vector[2]
    
    @property
    def sigma22(self):
        return self.Vector[3]
    
    @property
    def sigma23(self):
        return self.Vector[4]
    
    @property
    def sigma33(self):
        return self.Vector[5]