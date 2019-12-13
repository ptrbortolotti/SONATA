# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 14:53:09 2018

@author: TPflumm
"""
import os
from collections import OrderedDict
import numpy as np
import yaml

if __name__ == '__main__':
    os.chdir('/media/gu32kij/work/TPflumm/SONATA')
from SONATA.cbm.fileIO.read_yaml_input import clean_filestring
from SONATA.vabs.VABS_interface import VABS_config
from SONATA.classMaterial import read_IEA37_materials

class Configuration(object):
    '''
    Second Generation Configuration Class for the SONATA_CBM Disciplin  
    '''
    __slots__ = ( 'filename', 'setup', 'webs', 'segments', 'bw', 'flags', 'vabs_cfg') 
    def __init__(self, inputdata=None, materials = None, iae37 = False):
        self.setup, self.webs, self.segments, self.bw = {},{},{},{}
        self.filename = ''
        
        if isinstance(inputdata, str) and iae37 == False:
            self.filename = inputdata
            self._read_yaml_config(self.filename)
        
        elif isinstance(inputdata, dict) and iae37 == True:
            yml = inputdata
            self._read_IEA37(yml, materials)
            
        self.vabs_cfg = VABS_config()
        self.flags = {'mesh_core': True}


    def _read_IAE37(self, yml, materials):
        """ read the IAE37 style yml dictionary of a section and assign class 
        attributes to this configuration object
        
        Parameters:
        ----------
        yml : dict
           dictionary of the yaml style IAE37 section input    
        """
        #Setup:
        self.setup = {}
        self.setup['input_type'] = 5
        self.setup['datasource'] = None
        self.setup['material_db'] = None
        self.setup['radial_station'] = None
        self.setup['Theta'] = None  #GET FROM BLADE DEFINITION
        self.setup['scale_factor'] = 1 
        self.setup['BalanceWeight'] = False 
        self.setup['mesh_resolution'] = 280 #GET FROM function
       
        #Webs:
        self.webs = OrderedDict()
        for i,p1 in enumerate(sorted(yml.get('webs_positions'))):
            w = {}
            w['Pos1'] = p1
            w['Pos2'] = 1-p1
            self.webs[i+1] = w
        
        #Segments
        self.segments = OrderedDict()
        for i,s in enumerate(yml.get('segments')):
            d = {}
            key = s.get('id')
            if s.get('filler') == 'none': 
                d['CoreMaterial'] = 0
            else: 
                d['CoreMaterial'] = materials[s.get('filler')].id        
                    
            layerlst = s.get('layup')
            d['Layup'] = np.asarray([[l.get('start'), l.get('end'), l.get('thickness'), l.get('orientation'), materials[l.get('material_name')].id] for l in layerlst])
            d['Layup_names'] = [l.get('name') for l in layerlst]
            self.segments[key] = d
            
        #BalanceWeight
        if self.setup['BalanceWeight']:
            pass #TBD: Not yet needed!  
        
        
    def _read_yaml_config(self, fname):
        """ read the CBM config file and assign class attributes
        
        Parameters:
        ----------
        fname : str
            filename of the yaml style cbm input configuration    

        """
        b_string = clean_filestring(fname,comments='#')
        yDict =  yaml.load(b_string)
        self.setup = yDict['Setup']
        
        if 'Webs' in yDict.keys():
            D = {int(k.split()[-1]):v for (k, v) in yDict['Webs'].items()}    
            self.webs =  OrderedDict(sorted(D.items()))
        
        #print(self.setup['BalanceWeight'])
        if self.setup['BalanceWeight'] == True:
            self.bw = yDict['BalanceWeight']
                
        #read segments:
        D = {int(k.split()[-1]):v for (k, v) in yDict['Segments'].items()}

        for k in D:
            if D[k]['Layup'] != None:
                D[k]['Layup_names'] = np.asarray(D[k]['Layup'])[:,5].tolist()
                D[k]['Layup'] = np.asarray(D[k]['Layup'])[:,:5].astype(np.float)
            else:
                D[k]['Layup'] = np.empty((0,0))
                D[k]['Layup_names'] = np.empty((0,0))
                
        self.segments =  OrderedDict(sorted(D.items()))
        
        
        
        
    
if __name__ == '__main__':

    #classic configuration file: 
    os.chdir('/Users/rfeil/work/6_SONATA/SONATA')
    #fname = 'jobs/VariSpeed/advanced/sec_config.yml'
    fname = 'examples/sec_config.yml'
    cfg = Configuration(fname)


    #IAE37 Style configuration:
    with open('jobs/PBortolotti/IEAonshoreWT.yaml', 'r') as myfile:
        inputs  = myfile.read()
    
    yml = yaml.load(inputs)
    materials = {m.name:m for m in read_IAE37_materials(yml.get('materials'))}
    yml = yml.get('components').get('blade').get('2d_fem').get('sections')[1]
    wt_cfg = Configuration(yml, materials, iae37=True)       