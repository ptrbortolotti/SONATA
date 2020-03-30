#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from SONATA.classBlade import Blade

job = Blade(filename='simple_blade.yml')
job.blade_gen_section(mesh_flag = True, split_quads=False)
job.blade_plot_sections()
job.blade_gen_loft(ruled=True, tolerance=1e-6, continuity=4, check_compatibility=True, filename="simple.iges")
job.blade_post_3dtopo(flag_lft = True, flag_topo = True, flag_mesh = False)

loads = { "F" : np.array([[0, 0, 0, 0], 
                         [1, 0, 0, 0]]),
          "M" : np.array([[0, 0, 200, 0],
                         [1, 0, 100, 0]]) }
 
job.blade_run_vabs(loads, vabs_path='C:\Program Files (x86)\VABS\VABSIII')
job.blade_plot_sections(attribute = "stressM.sigma11")
