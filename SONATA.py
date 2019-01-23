#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SONATA: Multidiciplinary Rotor Blade Design Environment for Structural 
Optimization and Aeroelastic Analysis

The definition of the rotor blade topology is deliberately associated to the 
production of composite rotor blades. Thus, manufacturability is inherent from 
the geometric layup definition. Using orthogonal projection with corner-style 
differentiation the cross-section is discretized and can processed by the 
Variational Asymptotic Beam Sectional Analysis (VABS) afterwards. 

Date: 01/21/2019
@author: T.Pflumm, W.Garre, P.Bortolotti

Please Use the following Docstring styleguide:
https://numpydoc.readthedocs.io/en/latest/format.html

"""


    with open('jobs/PBortolotti/IEAonshoreWT.yaml', 'r') as myfile:
        inputs  = myfile.read()
    with open('jobs/PBortolotti/IEAturbine_schema.yaml', 'r') as myfile:
        schema  = myfile.read()
    validate(yaml.load(inputs), yaml.load(schema))
    wt_data     = yaml.load(inputs)
    
    
    Windturbine = Windturbine(readIAE37)