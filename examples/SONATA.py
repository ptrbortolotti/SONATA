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
import os
os.chdir('..')


import yaml
import numpy as np
from jsonschema import validate
import matplotlib.pyplot as plt
from SONATA.classBlade import Blade
from SONATA.classAirfoil import Airfoil
from SONATA.classMaterial import read_IEA37_materials
from SONATA.Pymore.app.marc_fanplot import frequency_analysis

plt.close('all')

job = Blade(filename='./jobs/MERIT/blade_config.yml')

job.blade_gen_section(mesh_flag = True, split_quads=False)

loads = {}
loads['F'] = np.array([[0.0, 72700, 0, 0], 
                       [1.0, 0, 0 , 0]]) 
loads['M'] = np.array([[0.0, -12, -636, -150], 
                       [1.0, 0, 0, 0 ]])     

job.blade_run_vabs(loads, rm_vabfiles=True)
#job.blade_run_anbax()

#job.blade_plot_sections(attribute='stressM.sigma11')
job.blade_plot_sections(savepath='merit.png')
job.blade_plot_sections(attribute='sf')
job.blade_post_3dtopo(flag_lft = True, flag_topo = True, flag_mesh = False)


beam = job.blade_exp_beam_props(solver='vabs', cosy='local', eta_offset=0.1491807)
dir_root = 'SONATA/Pymore/dym/mdl/03_rotormodel/02_UH60_rotor_fanplot_noaero/04_UH60_rotor_snglblade_noaero_pitchlink/'      
(freq, Omega, RPM_vec) = frequency_analysis(dir_root=dir_root)

