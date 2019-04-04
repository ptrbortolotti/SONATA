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

Date: 03/29/2019
@authors: T.Pflumm, W.Garre, P.Bortolotti, M. Morandini

Please Use the following Docstring styleguide:
https://numpydoc.readthedocs.io/en/latest/format.html

"""
import yaml
import numpy as np
from jsonschema import validate
import matplotlib.pyplot as plt

from SONATA.classBlade import Blade
from SONATA.classAirfoil import Airfoil
from SONATA.classMaterial import read_IEA37_materials


#with open('./jobs/PBortolotti/IEAonshoreWT.yaml', 'r') as myfile:
with open('./jobs/validation/boxbeam_config.yml', 'r') as myfile:
    inputs  = myfile.read()
yml = yaml.load(inputs)

airfoils = [Airfoil(af) for af in yml.get('airfoils')]
materials = read_IEA37_materials(yml.get('materials'))

job = Blade(name='IEAonshoreWT')
job.read_IEA37(yml.get('components').get('blade'), airfoils, materials, wt_flag=False)

job.blade_gen_section(mesh_flag = True, split_quads=True)
job.blade_run_vabs()
#job.blade_run_anbax()

#job.blade_plot_sections()
#job.blade_post_3dtopo(flag_lft = True, flag_topo = True, flag_mesh = True)
#anba = job.blade_exp_beam_props(solver='anbax')
vabs = job.blade_exp_beam_props(solver='vabs')
#anba_loc = job.blade_exp_beam_props(solver='anbax', cosy='local')
vabs_loc = job.blade_exp_beam_props(solver='vabs', cosy='local')

job.blade_plot_beam_props(ref = vabs_loc)