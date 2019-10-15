# -*- coding: utf-8 -*-
"""
Created on Mo Oct 07 11:18:28 2019

@author: Roland Feil
"""



import os
import numpy as np
from jsonschema import validate
import matplotlib.pyplot as plt
from SONATA.classBlade import Blade
from SONATA.classAirfoil import Airfoil
from SONATA.classMaterial import read_IEA37_materials




print('Current working directory is:', os.getcwd())
#os.chdir('..')

plt.close('all')

#job = Blade(name = 'UH-60A', filename = 'IEAontology_schema_new.yml') # initialize job with respective yaml input file
job = Blade(name = 'UH-60A', filename = '../../jobs/VariSpeed/UH-60A_adv.yaml') # initialize job with respective yaml input file


radial_stations = [0.1, 0.4]; # default: [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
#job = Blade(name = 'BAR_005a', filename = '/Users/rfeil/work/6_SONATA/SONATA/jobs/RFeil/IEAonshoreWT_BAR_005a.yaml', wt_flag = True, stations = radial_stations) # initialize job with respective yaml input file


job.blade_gen_section() # generate blade section(s) - build & mesh segments

job.blade_run_vabs()
#job.blade_run_anbax()




# === PLOTS ===
# --- plot blade sections ---
# ---------- FLAGS ----------
# plotTheta11 - ?
# theta_3 - show Theta_3 angles; Theta_3: rotational angles of individual cells in space in relative to the y-axis
# plotDisplacement - ToDO
# savepath - provide path for saving figures (save only if provided)
# Example: job.blade_plot_sections(attribute='theta_3', plotTheta11=True, savepath = '/Users/rfeil/work/6_SONATA/SONATA/jobs/RFeil/figures')
job.blade_plot_sections()



# --- plot 3D topology of the blade ---
# ---------- FLAGS ----------
# flag_wf - wireframe; i.e. maps all airfoil wires, scaled and rotated at every grid point
# flag_lft - loft flag to generate blade surface
# flag_topo - topology of selected blade sections
# flag_mesh - ToDo
#job.blade_post_3dtopo(flag_lft = True, flag_topo = True, flag_mesh = False)
job.blade_post_3dtopo(flag_wf = True, flag_lft = False, flag_topo = False)

#beam = job.blade_exp_beam_props(solver='vabs', cosy='local', eta_offset=0.1491807)


# EOF