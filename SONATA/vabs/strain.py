# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 13:25:49 2017

@author: TPflumm
"""

class Strain(object):
    #Object to represent the STRAIN Tensor within an Element    

    def __init__(self, Vec=None):
        if Vec is None:
            pass
        elif Vec is not None:
            self.Vector = Vec
        
    @property
    def epsilon11(self):    
        return self.Vector[0]
    
    @property
    def gamma12(self):          #2*epsilon12
        return self.Vector[1]
    
    @property
    def gamma13(self):          #2*epsilon13
        return self.Vector[2]
    
    @property
    def epsilon22(self):
        return self.Vector[3]
    
    @property
    def gamma23(self):          #2*epsilon23
        return self.Vector[4]
    
    @property
    def epsilon33(self):
        return self.Vector[5]
