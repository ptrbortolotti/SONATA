#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np

from SONATA.classBlade import Blade
from SONATA.classAirfoil import Airfoil
from SONATA.classMaterial import read_IEA37_materials
#from Pymore.app.marc_fanplot import frequency_analysis

job = Blade(filename='UH-60A/UH-60A.yml')

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
dir_root = '/media/gu32kij/work/TPflumm/Pymore/Pymore/dym/mdl/03_rotormodel/02_UH60_rotor_fanplot_noaero/04_UH60_rotor_snglblade_noaero_pitchlink/'      
(freq, Omega, RPM_vec) = frequency_analysis(dir_root=dir_root)
plot_fandiagram(freq, Omega, RPM_vec)
