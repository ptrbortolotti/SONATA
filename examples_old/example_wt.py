#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from SONATA.classBlade import Blade


flags = {"flag_ref_axes_wt":True,
         "flag_wt_ontology":True,}
    
job = Blade(filename="../jobs/RFeil/BAR0010n_curved_webs.yaml", flags=flags, stations = np.linspace(0.0, 1.0, 3))

#job = Blade(filename='UH-60A.yml')
#job.blade_gen_section(mesh_flag = True, split_quads=False)
#job.blade_plot_sections()
job.blade_gen_loft(ruled=False, tolerance=1e-6, max_degree=16, continuity=1, filename="wt.iges")
job.blade_post_3dtopo(flag_lft = True, flag_topo = True, flag_mesh = False)
#job.blade_run_vabs(vabs_path = 'C:\Program Files (x86)\VABS\VABSIII.exe')



