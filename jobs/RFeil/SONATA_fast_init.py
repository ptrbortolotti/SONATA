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
job_str = 'BAR009.yaml'

# job_str = 'yaml_examples/example_rectangular_beam_ht_ontology.yaml'  # apply ht ontology for rectangular beam example: flag_wt_ontology: False; flag_ref_axes_wt = False
# job_str = 'yaml_examples/example_circular_beam_ht_ontology.yaml'
# job_str = 'yaml_examples/example_circular_beam_wt_ontology.yaml'

job_name = 'SONATA_job'  # can also be more specified, e.g. 'BAR_005a', 'example_rectangular_beam_ht_ontology' or else
filename_str = folder_str + job_str

# ===== Define flags ===== #
# --- output flags ---
flag_write_BeamDyn      = True      # write BeamDyn input files for follow-up OpenFAST analysis
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
    # sine functional characteristics of Npt stations
    # 100 points in a sine distribution
    radial_stations_sine = [0.00000000,0.01570732,0.03141076,0.04710645,0.06279052,0.07845910,0.09410831,0.10973431,0.12533323,0.14090123,0.15643447,0.17192910,0.18738131,0.20278730,0.21814324,0.23344536,0.24868989,0.26387305,0.27899111,0.29404033,0.30901699,0.32391742,0.33873792,0.35347484,0.36812455,0.38268343,0.39714789,0.41151436,0.42577929,0.43993917,0.45399050,0.46792981,0.48175367,0.49545867,0.50904142,0.52249856,0.53582679,0.54902282,0.56208338,0.57500525,0.58778525,0.60042023,0.61290705,0.62524266,0.63742399,0.64944805,0.66131187,0.67301251,0.68454711,0.69591280,0.70710678,0.71812630,0.72896863,0.73963109,0.75011107,0.76040597,0.77051324,0.78043041,0.79015501,0.79968466,0.80901699,0.81814972,0.82708057,0.83580736,0.84432793,0.85264016,0.86074203,0.86863151,0.87630668,0.88376563,0.89100652,0.89802758,0.90482705,0.91140328,0.91775463,0.92387953,0.92977649,0.93544403,0.94088077,0.94608536,0.95105652,0.95579301,0.96029369,0.96455742,0.96858316,0.97236992,0.97591676,0.97922281,0.98228725,0.98510933,0.98768834,0.99002366,0.99211470,0.99396096,0.99556196,0.99691733,0.99802673,0.99888987,0.99950656,0.99987663,1.00000000]
    # 200 points in a sine distribution
    # radial_stations_sine = [0.00000000,0.00785390,0.01570732,0.02355976,0.03141076,0.03925982,0.04710645,0.05495018,0.06279052,0.07062699,0.07845910,0.08628637,0.09410831,0.10192446,0.10973431,0.11753740,0.12533323,0.13312134,0.14090123,0.14867243,0.15643447,0.16418685,0.17192910,0.17966075,0.18738131,0.19509032,0.20278730,0.21047176,0.21814324,0.22580127,0.23344536,0.24107506,0.24868989,0.25628937,0.26387305,0.27144045,0.27899111,0.28652455,0.29404033,0.30153796,0.30901699,0.31647697,0.32391742,0.33133789,0.33873792,0.34611706,0.35347484,0.36081083,0.36812455,0.37541557,0.38268343,0.38992769,0.39714789,0.40434360,0.41151436,0.41865974,0.42577929,0.43287258,0.43993917,0.44697862,0.45399050,0.46097437,0.46792981,0.47485639,0.48175367,0.48862124,0.49545867,0.50226553,0.50904142,0.51578590,0.52249856,0.52917900,0.53582679,0.54244154,0.54902282,0.55557023,0.56208338,0.56856185,0.57500525,0.58141318,0.58778525,0.59412106,0.60042023,0.60668235,0.61290705,0.61909395,0.62524266,0.63135280,0.63742399,0.64345586,0.64944805,0.65540017,0.66131187,0.66718277,0.67301251,0.67880075,0.68454711,0.69025124,0.69591280,0.70153143,0.70710678,0.71263852,0.71812630,0.72356978,0.72896863,0.73432251,0.73963109,0.74489406,0.75011107,0.75528181,0.76040597,0.76548321,0.77051324,0.77549574,0.78043041,0.78531693,0.79015501,0.79494435,0.79968466,0.80437564,0.80901699,0.81360845,0.81814972,0.82264052,0.82708057,0.83146961,0.83580736,0.84009355,0.84432793,0.84851021,0.85264016,0.85671752,0.86074203,0.86471344,0.86863151,0.87249601,0.87630668,0.88006330,0.88376563,0.88741345,0.89100652,0.89454464,0.89802758,0.90145512,0.90482705,0.90814317,0.91140328,0.91460716,0.91775463,0.92084548,0.92387953,0.92685660,0.92977649,0.93263902,0.93544403,0.93819134,0.94088077,0.94351216,0.94608536,0.94860019,0.95105652,0.95345417,0.95579301,0.95807290,0.96029369,0.96245524,0.96455742,0.96660010,0.96858316,0.97050647,0.97236992,0.97417339,0.97591676,0.97759994,0.97922281,0.98078528,0.98228725,0.98372863,0.98510933,0.98642926,0.98768834,0.98888650,0.99002366,0.99109975,0.99211470,0.99306846,0.99396096,0.99479214,0.99556196,0.99627038,0.99691733,0.99750280,0.99802673,0.99848910,0.99888987,0.99922904,0.99950656,0.99972243,0.99987663,0.99996916,1.00000000]
else:
    radial_stations_sine = []

# ===== User defined radial stations ===== #
# Define the radial stations for cross sectional analysis (only used for flag_wt_ontology = True -> otherwise, sections from yaml file are used!)
# radial_stations = [0., 0.1, 0.2, 0.301, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]  # ToDO: fix yaml for r/R = 0.3 and 0.9
radial_stations = [0.4]

# ===== Execute SONATA Blade Component Object ===== #
# name          - job name of current task
# filename      - string combining the defined folder directory and the job name
# flags         - communicates flag dictionary (defined above)
# stations      - input of radial stations for cross sectional analysis
# stations_sine - input of radial stations for refinement (only and automatically applied when lofing flag flag_lft = True)
job = Blade(name=job_name, filename=filename_str, flags=flags_dict, stations=radial_stations, stations_sine=radial_stations_sine)  # initialize job with respective yaml input file

# ===== Build & mesh segments ===== #
job.blade_gen_section()  # generate blade section(s) - build & mesh segments

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
            vabs_beam_stiff[i, j, :] = job.beam_properties[i, 1].TS[j, :]  # receive 6x6 timoshenko stiffness matrix
            vabs_beam_inertia[i, j, :] = job.beam_properties[i, 1].MM[j, :]  # receive 6x6 mass matrix
    # export_beam_struct_properties(folder_str, job_str, radial_stations, solver='vabs', beam_stiff=vabs_beam_stiff, beam_inertia=vabs_beam_inertia, beam_mass_per_length=vabs_beam_section_mass)  # export beam structural properties to csv file

    job.beam_properties = None  # clear var before initializing anbax

    job.blade_run_anbax()  # run anbax
    anbax_beam_stiff = np.zeros([len(radial_stations), 6, 6])
    anbax_beam_inertia = np.zeros([len(radial_stations), 6, 6])
    anbax_beam_section_mass = np.zeros([len(radial_stations), 1])
    for i in range(len(job.beam_properties)):
        anbax_beam_section_mass[i] = job.beam_properties[i, 1].m00  # receive mass per unit span
        for j in range(5):
            anbax_beam_stiff[i, j, :] = job.beam_properties[i, 1].TS[j, :]  # receive 6x6 timoshenko stiffness matrix
            anbax_beam_inertia[i, j, :] = job.beam_properties[i, 1].MM[j, :]  # receive 6x6 mass matrix

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
