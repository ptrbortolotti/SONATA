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

from SONATA.classBlade import Blade
from SONATA.classAirfoil import Airfoil
from SONATA.classMaterial import read_IEA37_materials

with open('jobs/VariSpeed/UH-60A_adv.yml', 'r') as f:
     yml = yaml.load(f.read())

airfoils = [Airfoil(af) for af in yml.get('airfoils')]
materials = read_IEA37_materials(yml.get('materials'))

byml = yml.get('components').get('blade')
B2 = Blade(name='UH-60A_adv')
B2.read_IEA37(byml, airfoils, materials)     
# B2.post_3dtopo()

for key, cs in B2.sections.items():
   print('STATUS:\t Building Section at grid location %s' % (key))
   cs.cbm_gen_topo()
   cs.cbm_gen_mesh()
   cs.cbm_run_vabs()
   cs.cbm_post_2dmesh()