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
    """
    general material class
    
    Attributes
    ----------
    id : int
        material identifier
        
    name : str
        short material name
        
    description : str
        description of the material, e.g.: unidirektional ht-carbon fiber 
        composite with epoxy matrix (FVC of 60%)
        
    source : str
        source of the material properties e.g.: elasitc properties derived 
        from Schuermann, (p.184, 417) with a semi-empiric Puck approach
        
    orth : int
        orth is the flag to indicate whether the material is isotropic (0), 
        orthotropic (1) or general anisotropic (2) in consitency with VABS 
        Manual for Users (2011)
        
    rho : float 
        density in g/cm**3
    """
    
    __slots__ = ( 'id', 'name', 'description', 'source', 'orth', 'rho') 
    
    def __init__(self, ID='NOID', name = 'noname', description = 'nodescription', source = 'nosource', orth=None, rho=0, **kwargs):
        self.id = ID
        self.name = name
        self.description = description
        self.source = source
        self.orth = orth
        self.rho = rho
       
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
    """
    Isotropic Material

    
    Attributes
    ----------
    E : float
        in GPa; Young's modulus   
        
    nu : float
        nondimensional; Poisson's ratio         
    
    alpha : np.ndarray
        in 1/K; coefficient of thermal expansion in direction
  
    YS :  float
        in MPa; Yield Strenth (Streckgrenze)
    UTS : float
        in MPa; Ultimate Tensile Strenght (Zugfestigkeit)
    """
    
    __slots__ = ( 'E', 'nu', 'alpha', 'YS', 'UTS') 
    
    def __init__(self, **kw):
        kw['orth'] = 0
        Material.__init__(self, **kw)
        
        self.E = kw.get('E')
        self.nu = kw.get('nu')
        self.alpha = kw.get('alpha')
        
        self.YS = kw.get('YS')    
        self.UTS = kw.get('UTS')  
            
        
class OrthotropicMaterial(Material):
    """
    Orthotropic Material

    
    Attributes
    ----------
    E : np.ndarray
        in GPa; [E_1, E_2, E_3], with E_i: axial tensile modules in direction i 
        and E_2 and E_3 the transverse tensile modules respectively       
        
    G : np.ndarray
        in GPa; [G_12, G_13, G_23], with G_ij, is the shear modulus in
        direction j on the plane whose normal is in direction  i; for 
        transversal insotropic materials G_13 = G_12  
    
    nu : np.ndarray
        nondimensional; [nu12, nu_13, nu_23], nu_ij is the Poisson's ratio that 
        corresponds to a contraction in direction j when an extension is 
        applied in direction i.
    
    alpha_11 : np.ndarray
        in 1/K; [alpha_11, alpha_22, alpha_33], alpha_ii is the coefficient of 
        thermal expansion in direction ii
  
    Xt :  float
        in MPa; 0° tensile strenght
    Xc : float
        in MPa; 0° compressive strenght
    Yt : float
        in MPa; 90° tensile strenght
    Yc : float
        in MPa; 90° compressive strenght
    S21 :
        in MPa; in-/out of plane shear strength 

    """
    
    __slots__ = ( 'E', 'G', 'nu', 'alpha', 'Xt', 'Xc', 'Yt', 'Yc', 'S21', 'S23') 
    def __init__(self, **kw):
        kw['orth'] = 1
        Material.__init__(self, **kw)
        
        self.E = np.asarray(kw.get('E'))
        self.G = np.asarray(kw.get('G')) 
        self.nu = np.asarray(kw.get('nu'))  
        self.alpha = kw.get('alpha')  
        
        if all (k in kw for k in ('E_1','E_2','E_3')):
            self.E = np.array([kw.get('E_1'), kw.get('E_2'), kw.get('E_3')])
             
        if all (k in kw for k in ('G_12','G_13','G_23')):
            self.G = np.array([kw.get('G_12'), kw.get('G_13'), kw.get('G_23')])
            
        if all (k in kw for k in ('nu_12','nu_13','nu_23')):
            self.nu = np.array([kw.get('nu_12'), kw.get('nu_13'), kw.get('nu_23')])
            
        if all (k in kw for k in ('alpha_11','alpha_22','alpha_33')):
            self.alpha = np.array([kw.get('alpha_11'), kw.get('alpha_22'), kw.get('alpha_33')])

        self.Xt = kw.get('Xt')   #Axial Tensile Strength in [MPa]
        self.Xc = kw.get('Xc')   #Axial Compression Strength  [MPa]
        self.Yt = kw.get('Xc')   #Transverse Tensile strenght  [MPa]
        self.Yc = kw.get('Yc')   #Transverse  Compression strenght  [Mpa]
        self.S21 = kw.get('S21') #in-/out of plane shear strength [MPa]
        self.S23 = kw.get('S23') 


class AnisotropicMaterial(Material):
    __slots__ = ( 'C', 'alpha') 
    
    def __init__(self, **kw):
        kw['orth'] = 2
        Material.__init__(self, **kw)
        
        self.C = np.asarray(kw.get('C'))  
        self.alpha = np.asarray(kw.get('alpha'))  
        
        #TODO: Set Strenght Characteristics according for Anisotropic Material 
        #as second and forth order Strength Tensor Fij, and F

def read_IEA37_materials(yml):
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

    materials1 = read_yml_materials('examples/mat_db.yml')

    with open('jobs/PBortolotti/IEAonshoreWT.yaml', 'r') as myfile:
        inputs  = myfile.read()
    wt_data     = yaml.load(inputs)    
    materials2 = read_IEA37_materials(wt_data['materials'])
    

    with open('jobs/VariSpeed/UH-60A_adv.yml', 'r') as myfile:
        inputs  = myfile.read()
    data     = yaml.load(inputs)['materials']    
    materials3 = read_IEA37_materials(data)