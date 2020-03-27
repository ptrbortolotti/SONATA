# -*- coding: utf-8 -*-
"""
Created on Mo Oct 07 11:18:28 2019

@author: Roland Feil
"""

import os, sys
sys.path.insert(0,'/Users/rfeil/work/6_SONATA/SONATA')  # import sys path to import 'SONATA' & 'job' modules

import time
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import csv

from openmdao.api import Group, Problem, ScipyOptimizer, IndepVarComp, ScipyOptimizeDriver, ExplicitComponent, ExecComp  #, pyOptSparseDriver
from openmdao.drivers.genetic_algorithm_driver import SimpleGADriver
# from concurrent import futures
from openmdao.drivers.pyoptsparse_driver import pyOptSparseDriver

from SONATA.classBlade import Blade
from SONATA.utl_openfast.utls_sonata2beamdyn import convert_structdef_SONATA_to_beamdyn


print('Current working directory is:', os.getcwd())



class StructOpt(ExecComp):

    def setup(self):

        # inputs
        self.add_input('OptVar1')
        self.add_input('OptVar2')
        self.add_input('OptVar3')
        # outputs
        self.add_output('Stiffness_obj')

    def compute(self, inputs, outputs):
        # OptVar = inputs['OptVar1']  # to-be-varied variable
        # OptVar = [inputs['OptVar1'], inputs['OptVar2']]  # to-be-varied variable
        OptVar = [inputs['OptVar1'], inputs['OptVar2'], inputs['OptVar3']]  # to-be-varied variable
        print('NEW OptVar1 = ' + str(float(inputs['OptVar1'])))
        print('NEW OptVar2 = ' + str(float(inputs['OptVar2'])))
        print('NEW OptVar3 = ' + str(float(inputs['OptVar3'])))

        # run job
        job = Blade(name=job_name, filename=filename_str, flags=flags_dict, stations=radial_stations, stations_add=radial_stations_add, flag_opt=flag_opt, opt_vars=OptVar)
        job.blade_gen_section(mesh_flag=True, split_quads=True)
        job.blade_run_vabs()

        beam_prop = job.beam_properties
        beam_prop_MM = beam_prop[0][1].MM
        beam_prop_TS = beam_prop[0][1].TS

        print('*******************')
        print('STIFFNESS:')
        print('GJ = ' + str(beam_prop_TS[3, 3]))  # GJ
        print('EIflap = ' + str(beam_prop_TS[4, 4]))  # EIFLAP
        print('EIlag = ' + str(beam_prop_TS[5, 5]))  # EILAG
        print('*******************')


        # outputs['EIflap'] = beam_prop['beam_stiff'][0,4,4]
        GJ_diff = np.abs(beam_prop_TS[3, 3] - beam_prop_init_TS[3, 3])/beam_prop_init_TS[3, 3]
        EIflap_diff = np.abs(beam_prop_TS[4, 4] - beam_prop_init_TS[4, 4])/beam_prop_init_TS[4, 4]
        EIlag_diff = np.abs(beam_prop_TS[5, 5] - beam_prop_init_TS[5, 5])/beam_prop_init_TS[5, 5]

        # outputs['Stiffness_obj'] = np.sqrt(GJ_diff**2 + EIflap_diff**2 + EIlag_diff**2)
        outputs['Stiffness_obj'] = np.sqrt(EIflap_diff**2 + EIlag_diff**2)
        # outputs['Stiffness_obj'] = GJ_diff
        # outputs['Stiffness_obj'] = EIflap_diff
        # outputs['Stiffness_obj'] = EIlag_diff

        print('OptVar1 = ' + str(float(OptVar[0])))
        print('OptVar2 = ' + str(float(OptVar[1])))
        print('OptVar3 = ' + str(float(OptVar[2])))
        print('Delta_GJ = ' + str(float(GJ_diff)))
        print('Delta_EIflap = ' + str(float(EIflap_diff)))
        print('Delta_EIlag = ' + str(float(EIlag_diff)))
        print('Stiffness OBJECTIVE = ' + str(float(outputs['Stiffness_obj'])))

        # job.blade_plot_sections(attribute='MatID', plotTheta11=False, plotDisplacement=False, savepath=folder_str, opt_var=str(float(OptVar1)))

        # write EIflap value of each iteration to file
        with open('inflatable_blade_studies/opt_temp.csv', 'a') as opt_csvfile:
            opt_csvfile_writer = csv.writer(opt_csvfile, delimiter=',')
                                        # OptVar1                   init_value                    iter_value              Objective - Delta
            opt_csvfile_writer.writerow([str(float(OptVar[0])), str(float(OptVar[1])), str(float(OptVar[2])),
                                         str(beam_prop_init_TS[3, 3]), str(beam_prop_TS[3, 3]), str(float(GJ_diff)),        # Torsion
                                         str(beam_prop_init_TS[4, 4]), str(beam_prop_TS[4, 4]), str(float(EIflap_diff)),    # Flap
                                         str(beam_prop_init_TS[5, 5]), str(beam_prop_TS[5, 5]), str(float(EIlag_diff)), str(float(outputs['Stiffness_obj']))])  # Lag




if __name__ == "__main__":

    start_time = time.time()

    job_name = 'SONATA_job'
    folder_str = '/Users/rfeil/work/6_SONATA/SONATA/jobs/RFeil/inflatable_blade_studies/'

    job_str_init = 'BAR0010n_HT_baseline_0a.yaml'  # Reference job for optimization objectives
    # job_str_init = 'BAR0010n_HT_baseline_kevlar_10.yaml'  # Reference job for optimization objectives
    filename_str_init = folder_str + job_str_init

    job_str = 'BAR0010n_HT_baseline_kevlar_13_carbon.yaml'  # To be optimized job
    # job_str = 'BAR0010n_HT_baseline_kevlar_11.yaml'  # To be optimized job
    filename_str = folder_str + job_str


    # ======================== #
    # ===== Define flags ===== #
    # ======================== #
    flag_wt_ontology = False
    flag_ref_axes_wt = False
    mesh_resolution = 200
    attribute_str = 'MatID'
    flag_plotDisplacement = False     # description ? ToDO
    flag_plotTheta11 = False
    flag_wf = True      # plot wire-frame
    flag_lft = True
    flag_topo = True      # plot mesh topology

    # create flag dictionary
    flags_dict = {"flag_wt_ontology": flag_wt_ontology, "flag_ref_axes_wt": flag_ref_axes_wt,
                "attribute_str": attribute_str,
                "flag_plotDisplacement": flag_plotDisplacement, "flag_plotTheta11": flag_plotTheta11,
                "flag_wf": flag_wf, "flag_lft": flag_lft, "flag_topo": flag_topo, "mesh_resolution": mesh_resolution}

    radial_stations_add = []
    radial_stations = [0.76458948]


    # ============================= #
    # ===== Compute Reference ===== #
    # ============================= #
    # flag_opt = False         # no optimization during initalization
    # Initialize Reference job for scaling of optimization objectives
    job_init = Blade(name=job_name, filename=filename_str_init, flags=flags_dict, stations=radial_stations, stations_add=radial_stations_add)
    job_init.blade_gen_section(topo_flag=True, mesh_flag = True, split_quads=True)
    job_init.blade_run_vabs()
    # job_init.blade_plot_sections(attribute='MatID', plotTheta11=False, plotDisplacement=False, savepath=folder_str, opt_var='init')

    beam_prop_init = job_init.beam_properties  # use of SONATA/VABS coordinate system

    beam_prop_init_MM = beam_prop_init[0][1].MM
    beam_prop_init_TS = beam_prop_init[0][1].TS

    print('*******************')
    print('REFERENCE STIFFNESS:')
    print('GJ = ' + str(beam_prop_init_TS[3, 3]))  # GJ
    print('EIflap = ' + str(beam_prop_init_TS[4, 4]))  # EIFLAP
    print('EIlag = ' + str(beam_prop_init_TS[5, 5]))  # EILAG
    print('*******************')



    # ======================== #
    # ===== Optimization ===== #
    # ======================== #

    flag_opt = True         # activate optimization flag after model initialization
    try:
        os.remove('inflatable_blade_studies/opt_temp.csv')  # remove *.csv file before running another optimization
    except:
        print(' ')

    # set up Optimization Problem
    p = Problem()

    # indeps = p.model.add_subsystem('indeps', IndepVarComp(), promotes=['OptVar1'])       # One design variables
    # indeps = p.model.add_subsystem('indeps', IndepVarComp(), promotes=['OptVar1', 'OptVar2'])     # Two design variables
    indeps = p.model.add_subsystem('indeps', IndepVarComp(), promotes=['OptVar1', 'OptVar2', 'OptVar3'])     # Two design variables
    # indeps.add_output('OptVar1', 0.0005)
    indeps.add_output('OptVar1', 0.02) #0.025)
    indeps.add_output('OptVar2', 0.31)
    indeps.add_output('OptVar3', 0.0) #4.5)
    # p.model.add_subsystem('structopt', StructOpt(), promotes_inputs=['OptVar1'], promotes_outputs=['Stiffness_obj'])
    # p.model.add_subsystem('structopt', StructOpt(), promotes_inputs=['OptVar1', 'OptVar2'], promotes_outputs=['Stiffness_obj'])
    p.model.add_subsystem('structopt', StructOpt(), promotes_inputs=['OptVar1', 'OptVar2', 'OptVar3'], promotes_outputs=['Stiffness_obj'])

    # --------------------------------- #
    # SNOPT optimizer from pyOptSParse
    # --------------------------------- #
    p.driver = pyOptSparseDriver()  # ScipyOptimizeDriver()
    p.driver.options['optimizer'] = "SNOPT"
    p.driver.opt_settings['Major feasibility tolerance'] = 1e-3  #1e-7
    p.driver.opt_settings['Major iterations limit'] = 20
    p.driver.opt_settings['Summary file'] = 'SNOPT_Summary_file.txt'
    p.driver.opt_settings['Print file'] = 'SNOPT_Print_file.txt'
    p.driver.opt_settings['Major step limit'] = 0.1 # 1.0e-02
    # stepsize = 1.0e-03  # stepsize in finite difference (fd) method


    # --------------------------------- #
    # Define Design Variables, Objecties and Constraints
    # --------------------------------- #
    # p.model.add_design_var('OptVar1', lower=0.00031, upper=0.02)  # Kevlar ply thickness: 0.01 in = 0.000254 m
    p.model.add_design_var('OptVar1', lower=0.001, upper=0.2)  #  filament wound thickness
    p.model.add_design_var('OptVar2', lower=0.05, upper=0.5)   #  ellipsis stretching parameter in m; i.e. curvature value
    p.model.add_design_var('OptVar3', lower=-20.0, upper=20.0, scaler=1e-3/90)   #  ply orientation angle of the materials
    p.model.add_objective('Stiffness_obj', scaler = 1)


    p.setup()
    p.model.approx_totals(method='fd') #, step=stepsize)
    p.run_driver() 


    print("--- Computational time: %s seconds ---" % (time.time() - start_time))


    # # ===== PLOTS ===== #
    # # saves figures in folder_str/figures if savepath is provided:
    # job.blade_plot_sections(attribute=attribute_str, plotTheta11=flag_plotTheta11, plotDisplacement=False, savepath=folder_str)
