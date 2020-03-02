# -*- coding: utf-8 -*-
"""
Created on Mo Oct 07 11:18:28 2019

@author: Roland Feil
"""

import os, sys
sys.path.insert(0,'/Users/rfeil/work/6_SONATA/SONATA')  # import sys path to import 'SONATA' & 'job' modules

import subprocess
import matplotlib.pyplot as plt

from SONATA.classBlade import Blade
from SONATA.utl.beam_struct_eval import beam_struct_eval
import SONATA.utl_openfast.fast_out_utilities as wtc_utilities

from SONATA.airconics_blade_cad.airfoil_cst import TurbineCAD
import SONATA.airconics_blade_cad.airconics.liftingSurface as liftingSurface


# ==============
# Main
# ==============

print('Current working directory is:', os.getcwd())
# os.chdir('..')

plt.close('all')  # close existing plots

# ===== Provide Path Directory & Yaml Filename ===== #
# folder_str = '/Users/rfeil/work/6_SONATA/SONATA/jobs/RFeil/'
folder_str = '/Users/rfeil/work/6_SONATA/SONATA/jobs/RFeil/HFM/'
# folder_str = '/Users/rfeil/work/6_SONATA/SONATA/jobs/RFeil/yaml_examples/'

# job_str = 'BAR0010n.yaml' # baseline
# job_str = 'BAR0010n_inflatable.yaml'
# job_str = 'BAR0013s_adapted.yaml'
# job_str = 'IEA-15-240-RWT_V1.yaml'
# job_str = 'box_beam_HT_layup1.yaml'
job_str = 'box_beam_HT_circular_tube.yaml'


# job_str = 'example_circular_beam_ht_ontology.yaml'  # apply ht ontology for circular beam example: flag_wt_ontology: False; flag_ref_axes_wt = False
# job_str = 'example_circular_beam_wt_ontology.yaml'
# job_str = 'example_circular_beam_wt_ontology_varied_materials.yaml'
# job_str = 'example_rectangular_beam_ht_ontology_anbax_verify.yaml'  # apply ht ontology for rectangular beam example: flag_wt_ontology: False; flag_ref_axes_wt = False
# job_str = 'example_rectangular_beam_ht_ontology_anbax_verify_II_plyangles.yaml'
# job_str = 'example_egg_beam_wt_ontology_varied_materials.yaml'

job_name = 'SONATA_job'  # can also be more specified, e.g. 'BAR_005a', 'example_rectangular_beam_ht_ontology' or else
filename_str = folder_str + job_str

# ===== Define flags ===== #
# --- numerical flags ---
flag_wt_ontology        = False      # if true, use ontology definition of wind turbines for yaml files
flag_ref_axes_wt        = False      # if true, rotate reference axes from wind definition to comply with SONATA (rotorcraft # definition)

# --- plotting flags ---
# Define mesh resolution, i.e. the number of points along the profile that is used for out-to-inboard meshing of a 2D blade cross section
mesh_resolution = 200
# For plots within blade_plot_sections
attribute_str           = 'MatID'   # default: MatID; others: theta_3 (ply orientation angles of individual cell alignment in space in relative to the y-axis)
flag_plotDisplacement   = False     # description ? ToDO
flag_plotTheta11        = False      # description ? ToDo # material orientation angle (eg +-45 deg)
# For plots within blade_post_3dtopo
flag_wf                 = True      # plot wire-frame
flag_lft                = True     # plot lofted shape of blade surface (flag_wf=True obligatory); Note: create loft with grid refinement without too many radial_stations; also exports step file of lofted shape
flag_topo               = True      # plot mesh topology

flag_airconics_CAD      = False      # create *.iges file

# create flag dictionary
flags_dict = {"flag_wt_ontology": flag_wt_ontology, "flag_ref_axes_wt": flag_ref_axes_wt, \
              "attribute_str": attribute_str, \
              "flag_plotDisplacement": flag_plotDisplacement, "flag_plotTheta11": flag_plotTheta11, \
              "flag_wf": flag_wf, "flag_lft": flag_lft, "flag_topo": flag_topo, "mesh_resolution": mesh_resolution}



# === Add additional radial stations === #
# important for resolving lofting issues
if flags_dict['flag_lft']:
    # npts = 40  # 200
    radial_stations_add = []
    # for n in range(npts):
    #     radial_stations_add.append(np.sin(n/npts*np.pi/2))
    # radial_stations_add = np.round(np.arange(0, 1.0, 0.0001), 2)
else:
    radial_stations_add = []

# ===== User defined radial stations ===== #
# Define the radial stations for cross sectional analysis (only used for flag_wt_ontology = True -> otherwise, sections from yaml file are used!)

# BAR0010
# radial_stations = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
# radial_stations = [0., 0.25, 0.5, 0.75, 1.0]

# IEA-15-240_RWT
# radial_stations = [0., 0.1, 0.2, 0.3, 0.4, 0.48, 0.6, 0.7, 0.8, 0.9, 1.0]
# radial_stations = [0., 0.1, 0.2, 0.3, 0.4, 0.501, 0.6, 0.7, 0.8, 0.9, 1.0]

# radial_stations = [0.1, 0.2, 0.3, 0.4, 0.5]
radial_stations = [0.0, 1.0]


# ===== Execute SONATA Blade Component Object ===== #
# name          - job name of current task
# filename      - string combining the defined folder directory and the job name
# flags         - communicates flag dictionary (defined above)
# stations      - input of radial stations for cross sectional analysis
# stations_sine - input of radial stations for refinement (only and automatically applied when lofing flag flag_lft = True)
job = Blade(name=job_name, filename=filename_str, flags=flags_dict, stations=radial_stations, stations_add=radial_stations_add)  # initialize job with respective yaml input file

# ===== Build & mesh segments ===== #
# job.blade_gen_section()  # generate blade section(s) - build & mesh segments
job.blade_gen_section(topo_flag=True, mesh_flag = True, split_quads=True)  # split_quads: vabs working with False & True; anbax working with True (only)


# ===== VABS / ANBAX individually ===== #
# job.blade_run_vabs()
# job.blade_run_anbax()

# ===== VABS / ANBAX combination/verification ===== #
# Define flags
flag_run_vabs = True
flag_run_anbax = False
flag_verify_vabs_anbax = False              # needs flag_run_vabs & flag_run_anbax set to True to be valid!
flag_DeamDyn_def_transform = False           # transform from SONATA to BeamDyn coordinate system
flag_plot_vabs_struct_characteristics = False     # plots and saves figures of structural characteristics of respective analysis for the range of defined radial stations
flag_csv_export = True                      # export csv files with structural data
flag_write_BeamDyn = False                   # write BeamDyn input files for follow-up OpenFAST analysis (requires flag_DeamDyn_def_transform = True)
# Update flags dictionary
flags_dict['flag_run_vabs'] = flag_run_vabs
flags_dict['flag_run_anbax'] = flag_run_anbax
flags_dict['flag_verify_vabs_anbax'] = flag_verify_vabs_anbax
flags_dict['flag_DeamDyn_def_transform'] = flag_DeamDyn_def_transform
flags_dict['flag_plot_vabs_struct_characteristics'] = flag_plot_vabs_struct_characteristics
flags_dict['flag_csv_export'] = flag_csv_export
flags_dict['flag_write_BeamDyn'] = flag_write_BeamDyn

# run evalutation
beam_struct_eval(flags_dict, radial_stations, job, folder_str, job_str)

# ===== PLOTS ===== #
# saves figures in folder_str/figures if savepath is provided:
job.blade_plot_sections(attribute=attribute_str, plotTheta11=flag_plotTheta11, plotDisplacement=False, savepath=folder_str)

# job.blade_post_3dtopo(flag_wf=flags_dict['flag_wf'], flag_lft=flags_dict['flag_lft'], flag_topo=flags_dict['flag_topo'])





# ===== WRITE *.iges CAD File ===== #
if flag_airconics_CAD == True & flag_wt_ontology == True:
    # So far, this toolbox relates to the yaml/BeamDyn definition; i.e. only works with flag_wt_ontology = True
    CAD_shape = TurbineCAD(filename_str)
    CAD_shape.create_blade_cst()

    blade = liftingSurface.LiftingSurface(CAD_shape.airfoil_func, NSegments=100)
    print("Finished creating loft")
    status = blade.Write(job_str[:-5] + '.iges')
    print(status)

    # optional plot that shows wires
    # fig, ax = plt.subplots()
    # # ax = fig.gca(projection='3d')
    # for eta in np.arange(0.0, 1.01, 0.05):
    #     CAD_shape.blade_af(eta, ax=ax)
    # plt.show()







# # ===== RUN BeamDyn or OpenFAST ===== #
# BeamDyn_driver = '/Users/rfeil/work/2_OpenFAST/openfast_BDloads_output/build/modules-local/beamdyn/beamdyn_driver'
# FAST_driver = '/Users/rfeil/work/2_OpenFAST/openfast_BDloads_output/build/glue-codes/openfast/openfast'
# analysis_dir = '/Users/rfeil/work/6_SONATA/SONATA/jobs/RFeil/00_analysis/analysis'
# input_file = 'bd_driver.inp'

# # Mac
# # olddir = os.getcwd()
# os.chdir(analysis_dir)
# # FNULL = open(os.devnull, 'w')
# print('Run BeamDyn')
# cmd = [BeamDyn_driver, input_file]
# stdout = subprocess.run(cmd, stdout=subprocess.PIPE).stdout.decode('utf-8')
# print(stdout)
# print('Wohoo - BeamDyn successfully executed!')

# # ===== Read output ===== #
# FAST_IO = wtc_utilities.FAST_IO()
# allinfo, alldata = FAST_IO.load_output([analysis_dir + '/' + input_file[:-4] + '.out'])


# # ===== PLOT ===== #
# cases = {}
# cases['RootForces'] = ['RootFxr', 'RootFyr', 'RootFzr']
# cases['RootMoments'] = ['RootMxr', 'RootMyr', 'RootMzr']
# cases['TipTransDefl'] = ['TipTDxr', 'TipTDyr', 'TipTDzr']
# cases['TipAngDefl'] = ['TipRDxr', 'TipRDyr', 'TipRDzr']
# flag_save_figs = False
# FAST_IO.plot_fast_out(cases, allinfo, alldata, flag_save_figs=flag_save_figs)

# EOF
