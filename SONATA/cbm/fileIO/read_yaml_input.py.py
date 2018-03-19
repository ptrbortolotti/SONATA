# -*- coding: utf-8 -*-
"""
Created on Tue Apr 04 14:49:03 2017

@author: TPflumm
"""
import yaml
import os
import numpy as np

if __name__ == '__main__':
    os.chdir('C://TPflumm_local/work/SONATA')
from SONATA.cbm.fileIO.material import IsotropicMaterial, OrthotropicMaterial, AnisotropicMaterial
from SONATA.vabs.VABS_interface import VABS_config, export_cells_for_VABS, XSectionalProperties

def clean_filestring(fname, comments='#', skiprows=0):
    b = [x.partition(comments)[0].rstrip()+'\n' for it, x in enumerate(open(fname, 'r').readlines()) if it>=skiprows]
    return ''.join([s for s in b if s.strip("\r\n").strip()])        


def read_yaml_materialdb(fname):
    b_string = clean_filestring(fname,comments='#')
    mdb =  yaml.load(b_string)['Materials']
    
    MaterialLst=[]
    for i,d in mdb.items():
        d['ID'] = int(i.split()[-1])
        if d['orth'] == 0:
            MaterialLst.append(IsotropicMaterial(**d))
        elif mdb[i]['orth'] == 1:
            MaterialLst.append(OrthotropicMaterial(**d))
        elif mdb[i]['orth'] == 2:
            MaterialLst.append(AnisotropicMaterial(**d))
        
    return sorted(MaterialLst, key=lambda x: x.id)


def read_yaml_config(fname):
    pass

if __name__ == '__main__':
  
    os.chdir('C://TPflumm_local/work/SONATA')
    fname = 'examples/MaterialDB.yaml'      
    MaterialLst=read_yaml_materialdb(fname)

    fname = 'examples/sec_config.yaml'  
    b_string = clean_filestring(fname,comments='#')
    yDict =  yaml.load(b_string)
    
    '''TODO:
        Change Configuration Class with its attributes to something nicer!
        '''
        
    #read setup
    SETUP_mat_filename = yDict['Setup']['mat_filename']
    SETUP_NbOfWebs = yDict['Setup']['NbOfWebs']
    SETUP_BalanceWeight = yDict['Setup']['BalanceWeight']
    SETUP_input_type = yDict['Setup']['input_type']
    SETUP_datasource = yDict['Setup']['datasource']
    if SETUP_input_type == 3 or SETUP_input_type == 4:
        SETUP_radial_station = yDict['Setup']['radial_station']
    SETUP_scale_factor = yDict['Setup']['scale_factor']
    SETUP_Theta = yDict['Setup']['Theta']
    SETUP_mesh_resolution = yDict['Setup']['mesh_resolution']
    
    #read webs
    for i, d in yDict['Webs'].items():
        d['ID'] = int(i.split()[-1])
        
    #read Segments
    for i, d in yDict['Segments'].items():
        d['ID'] = int(i.split()[-1])