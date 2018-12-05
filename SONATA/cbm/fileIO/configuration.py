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

class Configuration(object):
    '''
    Second Generation Configuration Class for the SONATA_CBM Disciplin  
    
    '''
    __slots__ = ( 'filename', 'setup', 'webs', 'segments', 'bw', 'flags', 'vabs_cfg') 
    def __init__(self, filename=None):
        if filename:
            self.setup, self.webs, self.segments, self.bw = {},{},{},{}
            self.read_yaml_config(filename)
            self.filename = filename
        else:
            self.filename = ''
            self.setup, self.webs, self.segments, self.bw = {},{},{},{}
            
        self.vabs_cfg = VABS_config()
        self.flags = {'mesh_core': True}

        
    def read_yaml_config(self, fname):
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
    os.chdir('/media/gu32kij/work/TPflumm/SONATA')
    #fname = 'jobs/VariSpeed/advanced/sec_config.yml'
    fname = 'examples/sec_config.yml'
    cfg = Configuration(fname)
