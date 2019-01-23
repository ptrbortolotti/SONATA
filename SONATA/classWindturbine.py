#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 09:37:42 2018

@author: Tobias Pflumm
"""
from jsonschema import validate
import yaml
import os

if __name__ == '__main__':
    os.chdir('..')

from SONATA.classPolar import Polar
from SONATA.classAirfoil import Airfoil
from SONATA.classMaterial import read_IAE37_materials

class Windturbine(object):

    __slots__ = ( 'name', 'materials', 'airfoils', 'components', )
    
    def __init__(self, yml=None, name = 'NONAME'):
        self.name = 'NONAME'
        self.components = None
        self.materials = None
        self.airfoils = None
        
        if isinstance(yml, dict): 
            self.read_IAE37(yml)
            
        if isinstance(name, str) and not 'NONAME': 
            self.name = name


    def __repr__(self):
        """__repr__ is the built-in function used to compute the "official" 
        string reputation of an object, """
        return 'Windturbine: '+ self.name
    
    
    def read_IAE37(self, yml, yml_schema):
        """
        reads the IAE Wind Task 37 style Airfoil dictionary and assigsn them to
        the class attributes
        """
        validate(yml, yml_schema)
        self.name = yml.get('name')
        self.airfoils = {af.get('name'):Airfoil(af) for af in yml.get('airfoils')}
        self.materials = read_IAE37_materials(yml.get('materials'))


    def plot_airfoilpolars():
        for af in wt.airfoils.values():
            af.plot_polars()




if __name__ == '__main__':   
    with open('jobs/PBortolotti/IEAonshoreWT.yaml', 'r') as myfile:
        inputs  = myfile.read()
    with open('jobs/PBortolotti/IEAturbine_schema.yaml', 'r') as myfile:
        schema  = myfile.read()
    
    yml = yaml.load(inputs)
    yml_schema = yaml.load(schema)
    wt = Windturbine()
    wt.read_IAE37(yml,yml_schema)