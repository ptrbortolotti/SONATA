# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 13:11:55 2018

@author: TPflumm

"""
import numpy as np

    
class Material(object):
    __slots__ = ( 'id', 'name', 'orth', 'rho', 'E', 'nu', 'G', 'alpha', 'C') 
    
    def __init__(self, ID='NOID', name = 'NONAME', orth=None, rho=0, **kwargs):
        self.id = ID
        self.name = name
        self.orth = orth
        self.rho = rho/1000 #from g/cm3 to g/mm3
        
        #########ELASTIC PROPERTIES#####################
        #ISOTROPIC
        if orth == 0: 
            self.E = kwargs.get('E')
            self.nu = kwargs.get('nu')  
            self.alpha = kwargs.get('alpha')  

        #orthotropic material
        elif orth == 1:
            self.E = kwargs.get('E') 
            self.G = kwargs.get('G') 
            self.nu = kwargs.get('nu')  
            self.alpha = kwargs.get('alpha')  
             
        #general anisotropic material
        elif orth == 2:
            self.C = np.asarray(kwargs.get('C'))  
            self.alpha = kwargs.get('alpha')  


        #########Strength#####################
        if orth == 0:
            self.YS = None # Yield Strenth (Streckgrenze)
            self.UTS = None #Tensile Strenght (Zugfestigkeit)
        
        if orth == 1:
            self.Xt = None #Axial Tensile Strength of unidirectional composite MPa
            self.Xc = None #Axial Compression Strength of unidirectional composite MPa
            self.Yt = None #Transverse strenght of unidirectional composite
            self.Yc = None #Transverse strenght of unidirectional composite
            self.S21 = None # 
            self.S23 = None #

    def __repr__(self): 
        return  str('Material %s: %s' % (self.id, self.name))
     

class IsotropicMaterial(Material):
    __slots__ = ( 'E', 'nu', 'alpha', 'YS', 'UTS') 
    
    def __init__(self, **kw):
        kw['orth'] = 0
        Material.__init__(self, **kw)
        
        self.E = kw.get('E')
        self.nu = kw.get('nu')
        self.alpha = kw.get('alpha')
        
        self.YS = kw.get('YS') # Yield Strenth (Streckgrenze)
        self.UTS = kw.get('UTS') #Tensile Strenght (Zugfestigkeit)
        
    
class OrthotropicMaterial(Material):
    __slots__ = ( 'E', 'G', 'nu', 'alpha', 'Xt', 'Xc', 'Yt', 'Yc', 'S21', 'S23') 
    
    def __init__(self, **kw):
        kw['orth'] = 1
        Material.__init__(self, **kw)
        
        self.E = np.asarray(kw.get('E'))  
        self.G = np.asarray(kw.get('G')) 
        self.nu = np.asarray(kw.get('nu'))  
        self.alpha = kw.get('alpha')  
             
        self.Xt = kw.get('Xt')  #Axial Tensile Strength of unidirectional composite MPa
        self.Xc = kw.get('Xc')  #Axial Compression Strength of unidirectional composite MPa
        self.Yt = kw.get('Xc') #Transverse strenght of unidirectional composite
        self.Yc = kw.get('Yc')  #Transverse strenght of unidirectional composite
        self.S21 = kw.get('S21') # 
        self.S23 = kw.get('S23') # 


class AnisotropicMaterial(Material):
    __slots__ = ( 'E', 'G', 'nu', 'alpha', 'Xt', 'Xc', 'Yt', 'Yc', 'S21', 'S23')
    
    def __init__(self, **kw):
        kw['orth'] = 2
        Material.__init__(self, **kw)
        
        self.C = np.asarray(kw.get('C'))  
        self.alpha = np.asarray(kw.get('alpha'))  
        
        #TODO: Set Strenght Characteristics according for Anisotropic Material
        self.Xt = kw.get('Xt')  #Axial Tensile Strength of unidirectional composite MPa
        self.Xc = kw.get('Xc')  #Axial Compression Strength of unidirectional composite MPa
        self.Yt = kw.get('Xc') #Transverse strenght of unidirectional composite
        self.Yc = kw.get('Yc')  #Transverse strenght of unidirectional composite
        self.S21 = kw.get('S21') # 
        self.S23 = kw.get('S23') # 
        
        def __repr__(self): 
            return  str('AnisotropicMaterial %s: %s' % (self.id, self.name))


if __name__ == '__main__':
    a = IsotropicMaterial(MatID=1, name='iso_mat', rho=0.4, )
    b = OrthotropicMaterial(MatID=2, name='orth_mat', rho=0.5)
    c = AnisotropicMaterial(MatID=3, name='aniso_mat', rho=0.6)
        