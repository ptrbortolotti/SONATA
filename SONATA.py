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
import yaml
from jsonschema import validate
from SONATA.classBlade import Blade
from SONATA.classAirfoil import Airfoil
from SONATA.classMaterial import read_IEA37_materials

#%% ====== WindTurbine ============== 
with open('./jobs/PBortolotti/IEAonshoreWT.yaml', 'r') as myfile:
    inputs  = myfile.read()
with open('jobs/PBortolotti/IEAontology_schema.yaml', 'r') as myfile:
    schema  = myfile.read()
validate(yaml.load(inputs), yaml.load(schema))    
yml = yaml.load(inputs)

airfoils = [Airfoil(af) for af in yml.get('airfoils')]
materials = read_IEA37_materials(yml.get('materials'))

job = Blade(name='IEAonshoreWT')
job.read_IEA37(yml.get('components').get('blade'), airfoils, materials, wt_flag=True)

job.blade_gen_section(mesh_flag = True)
job.blade_run_vabs()
job.blade_plot_sections()

job.blade_post_3dtopo(flag_lft = False, flag_topo = True)
