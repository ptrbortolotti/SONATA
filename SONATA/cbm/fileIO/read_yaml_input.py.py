# -*- coding: utf-8 -*-
"""
Created on Tue Apr 04 14:49:03 2017

@author: TPflumm
"""
import yaml


class SETUP(object):
    def __init__(self, mat_filename, NbOfWebs, BalanceWeight = False, input_type = 0):
        self.mat_filename = mat_filename
        self.NbOfWebs = NbOfWebs
        self.BalanceWeight = BalanceWeight
        self.input_type = input_type
    def __repr__(self):
        return "%s(name=%r, hp=%r, sp=%r)" % (
            self.__class__.__name__, self.mat_filename, self.NbOfWebs, self.BalanceWeight,  self.input_type)



if __name__ == '__main__':

    filename = 'yaml.input'
    a = ''
    with open(filename) as f:
        for line in f:
            line = line.partition('#')[0]
            line = line.rstrip()
            a += line 
            a += '\n'
         
    STR = ''.join([s for s in a.strip().splitlines(True) if s.strip("\r\n").strip()])    
    STR = STR.replace("\t", "    ") #replace tabs with 4 spaces
    
    
    cfg =  yaml.load(STR) 
    
    #for section in cfg:
        #print(section)
    print(cfg['Setup'])
    #print(cfg['Web'])
    #print(cfg['Seg_1'])