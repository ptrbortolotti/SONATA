# -*- coding: utf-8 -*-
"""
Created on Mo Oct 07 11:18:28 2019

@author: Roland Feil
"""

import os
import yaml
import numpy as np
from scipy.interpolate import interp1d
from jsonschema import validate
from SONATA.classBlade import Blade
from SONATA.classAirfoil import Airfoil
from SONATA.classMaterial import read_IEA37_materials
import matplotlib.pyplot as plt

from jobs.RFeil.write_sonata2beamdyn import write_beamdyn_axis, write_beamdyn_prop

# ==============
# Main
# ==============

print('Current working directory is:', os.getcwd())
# os.chdir('..')

plt.close('all')  # close existing plots

folder_str = '/Users/rfeil/work/6_SONATA/SONATA/jobs/RFeil/'
job_str = 'BAR006.yaml'  # 'IEAonshoreWT_BAR_005a.yaml'
# job_str = 'IEAonshoreWT.yaml'
filename_str = (folder_str + job_str)

# Define flags
flag_wt_ontology = True  # if true, use ontology definition of wind turbines for yaml files
flag_ref_axes_wt = True  # if true, rotate reference axes from wind definition to comply with SONATA (rotorcraft # definition)
# create flag dictionary
flags_dict = {"flag_wt_ontology": flag_wt_ontology, "flag_ref_axes_wt": flag_ref_axes_wt}

# User defined radial stations; default: [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
radial_stations = [0.4]
# job = Blade(name = 'BAR_005a', filename = filename_str, folder = folder_str, wt_flag = True, s2bd_flag = True, stations = radial_stations) # initialize job with respective yaml input file

# job = Blade(name='BAR_005a', filename=filename_str, wt_flag=True,
#             stations=radial_stations)  # initialize job with respective yaml input file


job = Blade(name='BAR_005a', filename=filename_str, flags=flags_dict,
            stations=radial_stations)  # initialize job with respective yaml input file

job.blade_gen_section()  # generate blade section(s) - build & mesh segments

job.blade_run_vabs()
# job.blade_run_anbax()


beam_stiff = np.zeros([len(radial_stations), 6, 6])
beam_inertia = np.zeros([len(radial_stations), 6, 6])

# get stiffness and inertia martrices out
for i in range(len(job.beam_properties)):
    for j in range(5):
        # receive 6x6 timoshenko stiffness matrix
        beam_stiff[i, j, :] = job.beam_properties[i, 1].TS[j, :]
        # receive 6x6 mass matrix
        beam_stiff[i, j, :] = job.beam_properties[i, 1].TS[j, :]

# write BeamDyn input files
write_beamdyn_axis(folder_str, job.yml.get('name'), job.yml.get('components').get('blade'))
write_beamdyn_prop(folder_str, job.yml.get('name'), radial_stations, beam_stiff, beam_inertia)



job.blade_plot_sections()


job.blade_post_3dtopo(flag_wf = True, flag_lft = False, flag_topo = False) #  ToDO flag_lft

# EOF
