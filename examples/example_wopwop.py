#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

from SONATA.classBlade import Blade
from SONATA.cbm.topo.utils import PntLst_to_npArray

blade1 = Blade(name='UH-60A', filename = '../jobs/MonteCarlo/UH-60A_adv.yml')

MeshBSplineLst, MeshPnts, meshvecs = blade1.blade_gen_wopwop_mesh(np.linspace(0.1491807,1,40), 25, deformation=None)
#blade1.blade_post_3dtopo(flag_lft=True, flag_wopwop=True)
pntarray =  np.asarray([PntLst_to_npArray(cs) for cs in MeshPnts]) 
vecsarray =  np.asarray([PntLst_to_npArray(cs) for cs in meshvecs]) 

#blade1.blade_gen_section()
#blade1.blade_run_vabs(rm_vabfiles=False)
#blade1.blade_plot_sections(savepath = 'jobs/MERIT/test.jpg')
blade1.blade_post_3dtopo(flag_topo = True, flag_lft = True, flag_wopwop=True)