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

# ===== Provide Path Directory & Yaml Filename ===== #
folder_str = '/Users/rfeil/work/6_SONATA/SONATA/jobs/RFeil/'
job_str = 'BAR007_YAML_out_NPTS50.yaml'  # 'IEAonshoreWT_BAR_005a.yaml'
job_str = 'example_rectangular_beam.yaml'
# job_str = 'IEAonshoreWT.yaml'
filename_str = folder_str + job_str

# ===== Define flags ===== #
# output flags
flag_write_BeamDyn = False
# numerical flags
flag_wt_ontology = False  # if true, use ontology definition of wind turbines for yaml files
flag_ref_axes_wt = True  # if true, rotate reference axes from wind definition to comply with SONATA (rotorcraft # definition)
flags_dict = {"flag_wt_ontology": flag_wt_ontology, "flag_ref_axes_wt": flag_ref_axes_wt}  # create flag dictionary

# ===== User defined radial stations ===== #
# default (if default should be applied unhinge input radial_stations from fob = Blade): [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
# radial_stations = [0., 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8]  # ToDO: fix yaml for r/R = 0.0 and 1.0 (also 0.3 and others?)
radial_stations = [0., 0.5, 1.0]

# ===== Execute SONATA Blade Component Object ===== #
# job = Blade(name='BAR_005a', filename=filename_str, folder=folder_str, wt_flag=True, s2bd_flag = True, stations = radial_stations) # initialize job with respective yaml input file
# job = Blade(name='BAR_005a', filename=filename_str, wt_flag=True, stations=radial_stations)  # initialize job with respective yaml input file

job = Blade(name='BAR_005a', filename=filename_str, flags=flags_dict, stations=radial_stations)  # initialize job with respective yaml input file

# ===== Build & mesh segments ===== #
job.blade_gen_section()  # generate blade section(s) - build & mesh segments

# ===== VABS / ANBAX ===== #
job.blade_run_vabs()
# job.blade_run_anbax()


# ===== BeamDyn Export ===== #
beam_stiff = np.zeros([len(radial_stations), 6, 6])
beam_inertia = np.zeros([len(radial_stations), 6, 6])

# get stiffness and inertia martrices out
for i in range(len(job.beam_properties)):
    for j in range(5):
        beam_stiff[i, j, :] = job.beam_properties[i, 1].TS[j, :]  # receive 6x6 timoshenko stiffness matrix
        beam_stiff[i, j, :] = job.beam_properties[i, 1].TS[j, :]  # receive 6x6 mass matrix

# write BeamDyn input files
if flag_write_BeamDyn:
    print('STATUS:\t Write BeamDyn input files')
    write_beamdyn_axis(folder_str, job.yml.get('name'), job.yml.get('components').get('blade'))
    write_beamdyn_prop(folder_str, job.yml.get('name'), radial_stations, beam_stiff, beam_inertia)



# ===== PLOTS ===== #
# flag_plotTheta11 = True  # description ?
attribute_str = 'MatID'  # default: MatID; others: theta_3 (ply orientation angles of individual cell alignment in space in relative to the y-axis)
flag_plotDisplacement = False  # description ? ToDO
flag_plotTheta11 = True  # description ? ToDo # material orientation angle (eg +-45 deg)
# saves figures in folder_str/figures if savepath is provided:
# job.blade_plot_sections(attribute=attribute_str, plotTheta11=flag_plotTheta11, plotDisplacement=False, savepath=folder_str)


flag_wf = True      # plot wire-frame
flag_lft = False    # plot loft of blade surface (flag_wf=True obligatory)
flag_topo = True   # plot mesh topology
job.blade_post_3dtopo(flag_wf=flag_wf, flag_lft=flag_lft, flag_topo=flag_topo) #  ToDO flag_lft

# EOF
