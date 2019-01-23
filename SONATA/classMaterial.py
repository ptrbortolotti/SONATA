#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 09:38:33 2018

@author: Tobias Pflumm
"""
import os
import yaml
import numpy as np
from collections import OrderedDict

if __name__ == '__main__':
    os.chdir('..')
from SONATA.cbm.fileIO.read_yaml_input import clean_filestring

class Material(object):
    __slots__ = ( 'id', 'name', 'orth', 'rho') 
    
    def __init__(self, ID='NOID', name = 'NONAME', orth=None, rho=0, **kwargs):
        self.id = ID
        self.name = name
        self.orth = orth
        self.rho = rho/1000 #from g/cm3 to g/mm3
    
    def __repr__(self): 
        if self.orth==0:
            return  str('%s: IsotropicMaterial: %s' % (self.id, self.name))    
        elif self.orth == 1:
            return  str('%s: OrthotropicMaterial: %s' % (self.id, self.name))
        elif self.orth == 2:
            return  str('%s: OrthotropicMaterial: %s' % (self.id, self.name))
        else:
            return  str('%s: UndefinedMaterial: %s' % (self.id, self.name))
            

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
    __slots__ = ( 'C', 'alpha')
    
    def __init__(self, **kw):
        kw['orth'] = 2
        Material.__init__(self, **kw)
        
        self.C = np.asarray(kw.get('C'))  
        self.alpha = np.asarray(kw.get('alpha'))  
        
        #TODO: Set Strenght Characteristics according for Anisotropic Material


def read_IAE37_materials(yml):
    materials = OrderedDict()
    for i,mat in enumerate(yml):
        ID = i+1
        if mat.get('orth') == 0:
            materials[ID] = IsotropicMaterial(ID = ID, **mat)
        elif mat.get('orth') == 1:
            materials[ID] = OrthotropicMaterial(ID = ID, **mat)
        elif mat.get('orth') == 2:
             materials[ID] =AnisotropicMaterial(ID = ID, **mat)
    return materials


def read_yml_materials(fname):
    b_string = clean_filestring(fname,comments='#')
    mdb =  yaml.load(b_string)['Materials']
    
    materials = OrderedDict()
    for k,v in mdb.items():
        ID = int(k.split()[-1])
        if v['orth'] == 0:
            materials[ID] = IsotropicMaterial(ID = ID, **v)
        elif mdb[k]['orth'] == 1:
            materials[ID] = OrthotropicMaterial(ID = ID, **v)
        elif mdb[k]['orth'] == 2:
            materials[ID] = AnisotropicMaterial(ID = ID, **v)
    return materials
        

def find_material(materials, attr, value):
    return next((x for x in materials.values() if getattr(x,attr) == value), None)
    

if __name__ == '__main__':
    a = IsotropicMaterial(ID=1, name='iso_mat', rho=0.4, )
    b = OrthotropicMaterial(ID=2, name='orth_mat', rho=0.5)
    c = AnisotropicMaterial(ID=3, name='aniso_mat', rho=0.6)

    with open('jobs/PBortolotti/IEAonshoreWT.yaml', 'r') as myfile:
        inputs  = myfile.read()
    wt_data     = yaml.load(inputs)    
    materials1 = read_IAE37_materials(wt_data['materials'])
    
    materials2 = read_yml_materials('examples/mat_db.yml')