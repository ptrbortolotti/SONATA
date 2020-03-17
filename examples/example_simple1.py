#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from SONATA.classBlade import Blade

job = Blade(filename='UH-60A.yml')

job.blade_gen_section(mesh_flag = True, split_quads=False)
#job.blade_plot_sections()
job.blade_post_3dtopo(flag_lft = True, flag_topo = True, flag_mesh = False)