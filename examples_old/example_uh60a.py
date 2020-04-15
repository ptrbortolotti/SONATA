#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from SONATA.classBlade import Blade

job = Blade(filename='UH-60A.yml')
job.blade_gen_section(mesh_flag = True, split_quads=False)
job.blade_plot_sections()
#job.blade_gen_loft(ruled=True, tolerance=1e-6, continuity=4, check_compatibility=True, filename="uh60a.iges")
job.blade_post_3dtopo(flag_lft = True, flag_topo = True, flag_mesh = False)

job.blade_run_vabs(vabs_path = 'C:\Program Files (x86)\VABS\VABSIII.exe')