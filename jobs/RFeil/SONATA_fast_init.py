# -*- coding: utf-8 -*-
"""
Created on Mo Oct 07 11:18:28 2019

@author: Roland Feil
"""

import os
from builtins import len, range


import yaml
import numpy as np
from scipy.interpolate import interp1d
from jsonschema import validate
from SONATA.classBlade import Blade
from SONATA.classAirfoil import Airfoil
from SONATA.classMaterial import read_IEA37_materials
import matplotlib.pyplot as plt

from jobs.RFeil.write_sonata2beamdyn import write_beamdyn_axis, write_beamdyn_prop

from jobs.RFeil.plot_vabs_anbax_verification import plot_vabs_anbax, export_beam_struct_properties

# ==============
# Main
# ==============

print('Current working directory is:', os.getcwd())
# os.chdir('..')

plt.close('all')  # close existing plots

# ===== Provide Path Directory & Yaml Filename ===== #
folder_str = '/Users/rfeil/work/6_SONATA/SONATA/jobs/RFeil/'

# job_str = 'BAR007.yaml'  # 'IEAonshoreWT_BAR_005a.yaml'
# job_str = 'BAR009.yaml'
job_str = 'BAR032.yaml'

# job_str = 'yaml_examples/example_rectangular_beam_ht_ontology.yaml'  # apply ht ontology for rectangular beam example: flag_wt_ontology: False; flag_ref_axes_wt = False
# job_str = 'yaml_examples/example_circular_beam_ht_ontology.yaml'
# job_str = 'yaml_examples/example_circular_beam_wt_ontology.yaml'

job_name = 'SONATA_job'  # can also be more specified, e.g. 'BAR_005a', 'example_rectangular_beam_ht_ontology' or else
filename_str = folder_str + job_str

# ===== Define flags ===== #
# --- output flags ---
flag_write_BeamDyn      = False      # write BeamDyn input files for follow-up OpenFAST analysis
# --- numerical flags ---
flag_wt_ontology        = True      # if true, use ontology definition of wind turbines for yaml files
flag_ref_axes_wt        = True      # if true, rotate reference axes from wind definition to comply with SONATA (rotorcraft # definition)

# --- plotting flags ---
# For plots within blade_plot_sections
attribute_str           = 'MatID'   # default: MatID; others: theta_3 (ply orientation angles of individual cell alignment in space in relative to the y-axis)
flag_plotDisplacement   = False     # description ? ToDO
flag_plotTheta11        = True      # description ? ToDo # material orientation angle (eg +-45 deg)
# For plots within blade_post_3dtopo
flag_wf                 = True      # plot wire-frame
flag_lft                = True     # plot loft of blade surface (flag_wf=True obligatory); Note: create loft with grid refinement without too many radial_stations
flag_topo               = True      # plot mesh topology

# create flag dictionary
flags_dict = {"flag_wt_ontology": flag_wt_ontology, "flag_ref_axes_wt": flag_ref_axes_wt, \
              "attribute_str": attribute_str, \
              "flag_plotDisplacement": flag_plotDisplacement, "flag_plotTheta11": flag_plotTheta11, \
              "flag_wf": flag_wf, "flag_lft": flag_lft, "flag_topo": flag_topo}


# === Add additional radial stations === #
# important for resolving lofting issues
if flags_dict['flag_lft']:
    npts = 100
    radial_stations_sine = []
    for n in range(npts):
        radial_stations_sine.append(np.sin(n/npts*np.pi/2))
else:
    radial_stations_sine = []

# ===== User defined radial stations ===== #
# Define the radial stations for cross sectional analysis (only used for flag_wt_ontology = True -> otherwise, sections from yaml file are used!)
# radial_stations = [0., 0.1, 0.2, 0.301, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]  # BAR209 ToDO: fix yaml for r/R = 0.3 and 0.9
# radial_stations = [0., 0.1, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]  # ToDO: fix yaml for r/R = 0.3 and 0.9
radial_stations = [0.4]

# ===== Execute SONATA Blade Component Object ===== #
# name          - job name of current task
# filename      - string combining the defined folder directory and the job name
# flags         - communicates flag dictionary (defined above)
# stations      - input of radial stations for cross sectional analysis
# stations_sine - input of radial stations for refinement (only and automatically applied when lofing flag flag_lft = True)
# job = Blade(name=job_name, filename=filename_str, flags=flags_dict, stations=radial_stations, stations_sine=radial_stations_sine)  # initialize job with respective yaml input file

# ===== Build & mesh segments ===== #
job.blade_gen_section(mesh_flag = True, split_quads=True)  # vabs working with False; anbax working with True
# job.blade_gen_section()  # generate blade section(s) - build & mesh segments

# ===== VABS / ANBAX ===== #
job.blade_run_vabs()
# job.blade_run_anbax()


# The following code is only used for VABS/anbax verification studies
flag_verify_vabs_anbax = False

if flag_verify_vabs_anbax:
    # ------------------------ #
    # VABS/ANBAX verification

    job.blade_run_vabs()  # run VABS
    vabs_beam_stiff = np.zeros([len(radial_stations), 6, 6])
    vabs_beam_inertia = np.zeros([len(radial_stations), 6, 6])
    vabs_beam_section_mass = np.zeros([len(radial_stations), 1])
    for i in range(len(job.beam_properties)):
        vabs_beam_section_mass[i] = job.beam_properties[i, 1].m00  # receive mass per unit span
        for j in range(6):
            vabs_beam_stiff[i, j, :] = np.array(job.beam_properties[i, 1].TS[j, :])  # receive 6x6 timoshenko stiffness matrix
            vabs_beam_inertia[i, j, :] = np.array(job.beam_properties[i, 1].MM[j, :])  # receive 6x6 mass matrix
    # export_beam_struct_properties(folder_str, job_str, radial_stations, solver='vabs', beam_stiff=vabs_beam_stiff, beam_inertia=vabs_beam_inertia, beam_mass_per_length=vabs_beam_section_mass)  # export beam structural properties to csv file

    job.beam_properties = None  # clear var before initializing anbax

    job.blade_run_anbax()  # run anbax
    anbax_beam_stiff = np.zeros([len(radial_stations), 6, 6])
    anbax_beam_inertia = np.zeros([len(radial_stations), 6, 6])
    anbax_beam_section_mass = np.zeros([len(radial_stations), 1])
    for i in range(len(job.beam_properties)):
        anbax_beam_section_mass[i] = job.beam_properties[i, 1].m00  # receive mass per unit span
        for j in range(5):
            anbax_beam_stiff[i, j, :] = np.array(job.beam_properties[i, 1].TS[j, :])  # receive 6x6 timoshenko stiffness matrix
            anbax_beam_inertia[i, j, :] = np.array(job.beam_properties[i, 1].MM[j, :])  # receive 6x6 mass matrix

    plot_vabs_anbax(radial_stations, vabs_beam_stiff, vabs_beam_inertia, anbax_beam_stiff, anbax_beam_inertia, job_str)  # plot 6x6 stiffness and mass matrices (vabs/anbax comparison)
    # ------------------------ #


# ===== BeamDyn Export ===== #
# initialize
beam_stiff = np.zeros([len(radial_stations), 6, 6])  # adapt length init for flag_wt_ontology=False !!!
beam_inertia = np.zeros([len(radial_stations), 6, 6])  # adapt length init for flag_wt_ontology=False !!!

# write BeamDyn input files
if flag_write_BeamDyn:
    # get stiffness and inertia martrices out
    for i in range(len(job.beam_properties)):
        for j in range(5):
            beam_stiff[i, j, :] = job.beam_properties[i, 1].TS[j, :]  # receive 6x6 timoshenko stiffness matrix
            beam_stiff[i, j, :] = job.beam_properties[i, 1].MM[j, :]  # receive 6x6 mass matrix

    print('STATUS:\t Write BeamDyn input files')
    write_beamdyn_axis(folder_str, job.yml.get('name'), job.yml.get('components').get('blade'))
    write_beamdyn_prop(folder_str, job.yml.get('name'), radial_stations, beam_stiff, beam_inertia)



# ===== PLOTS ===== #
# saves figures in folder_str/figures if savepath is provided:
job.blade_plot_sections(attribute=attribute_str, plotTheta11=flag_plotTheta11, plotDisplacement=False, savepath=folder_str)

job.blade_post_3dtopo(flag_wf=flags_dict['flag_wf'], flag_lft=flags_dict['flag_lft'], flag_topo=flags_dict['flag_topo'])

# EOF
