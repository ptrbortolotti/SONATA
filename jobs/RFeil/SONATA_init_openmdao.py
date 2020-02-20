"""
openMDAO wrapper for SONATA
"""
import os
from datetime import datetime
import numpy as np
from openmdao.api import Group, Problem, ScipyOptimizer, IndepVarComp, ScipyOptimizeDriver, ExplicitComponent, ExecComp  #, pyOptSparseDriver
from openmdao.drivers.genetic_algorithm_driver import SimpleGADriver
# from concurrent import futures
from openmdao.drivers.pyoptsparse_driver import pyOptSparseDriver

import sys
sys.path.insert(0,'/Users/rfeil/work/6_SONATA/SONATA')  # import sys path to import 'SONATA' & 'job' modules

from SONATA.classBlade import Blade
from jobs.RFeil.utls.utls_sonata2beamdyn import convert_structdef_SONATA_to_beamdyn

os.chdir('/Users/rfeil/work/6_SONATA/SONATA/jobs/RFeil')  # go back to chosen working directory
print('Current working directory is:', os.getcwd())

import csv

class StructOpt(ExecComp):

    def setup(self):
        # self.counter = 0
        # self.startTime = datetime.now()
        # print(self.startTime)

        self.add_input('pos_web0')
        self.add_output('EIflap')

    def compute(self, inputs, outputs):
        OptVar1 = inputs['pos_web0']  # to-be-varied variable

        # run job
        job = Blade(name=job_name, filename=filename_str, flags=flags_dict, stations=radial_stations, stations_add=radial_stations_add, flag_opt=flag_opt, opt_vars=OptVar1)
        job.blade_gen_section(mesh_flag=True, split_quads=True)
        job.blade_run_vabs()
        beam_prop = convert_structdef_SONATA_to_beamdyn(radial_stations, job.beam_properties)  # convert to BeamDyn coord sys def


        test_value = 4406212.6351

        # outputs['EIflap'] = np.abs(test_value * (1 + OptVar1) - 1.06*test_value)
        # outputs['EIflap'] = -0.5 * OptVar1**2


        outputs['EIflap'] = beam_prop['beam_stiff'][0,4,4]


        print('OptVar1 = ' + str(OptVar1))
        print('EIflap = ' + str(float(outputs['EIflap'])))

        job.blade_plot_sections(attribute='MatID', plotTheta11=False, plotDisplacement=False, savepath=folder_str, opt_var=str(float(OptVar1)))


        # write EIflap value of each iteration to file
        with open('figures/opt_temp.csv', 'a') as opt_csvfile:
            opt_csvfile_writer = csv.writer(opt_csvfile, delimiter=',')
            opt_csvfile_writer.writerow([str(float(OptVar1)), str(float(outputs['EIflap']))])




if __name__ == "__main__":


    job_name = 'SONATA_job'
    folder_str = '/Users/rfeil/work/6_SONATA/SONATA/jobs/RFeil/'
    job_str = 'BAR0010n.yaml'  # baseline
    filename_str = folder_str + job_str

    # ===== Define flags ===== #
    flag_wt_ontology = True
    flag_ref_axes_wt = True
    mesh_resolution = 200
    attribute_str = 'MatID'
    flag_plotDisplacement = False     # description ? ToDO
    flag_plotTheta11 = True
    flag_wf = True      # plot wire-frame
    flag_lft = True
    flag_topo = True      # plot mesh topology

    # create flag dictionary
    flags_dict = {"flag_wt_ontology": flag_wt_ontology, "flag_ref_axes_wt": flag_ref_axes_wt,
                "attribute_str": attribute_str,
                "flag_plotDisplacement": flag_plotDisplacement, "flag_plotTheta11": flag_plotTheta11,
                "flag_wf": flag_wf, "flag_lft": flag_lft, "flag_topo": flag_topo, "mesh_resolution": mesh_resolution}

    radial_stations_add = []
    radial_stations = [0.95]


    # flag_opt = False         # no optimization during initalization
    # Initialize Reference job for scaling of optimization objectives
    job_init = Blade(name=job_name, filename=filename_str, flags=flags_dict, stations=radial_stations, stations_add=radial_stations_add)
    job_init.blade_gen_section(mesh_flag=True, split_quads=True)
    job_init.blade_run_vabs()
    job_init.blade_plot_sections(attribute='MatID', plotTheta11=False, plotDisplacement=False, savepath=folder_str, opt_var='init')

    beam_prop_init = convert_structdef_SONATA_to_beamdyn(radial_stations, job_init.beam_properties)  # convert to BeamDyn coord sys def

    print('EIflap_init = ' + str(beam_prop_init['beam_stiff'][0][4, 4]))

    flag_opt = True         # activate optimization flag after model initialization
    # set up Optimization Problem
    p = Problem()

    indeps = p.model.add_subsystem('indeps', IndepVarComp(), promotes=['pos_web0'])
    indeps.add_output('pos_web0', -0.01)
    p.model.add_subsystem('structopt', StructOpt(), promotes_inputs=['pos_web0'], promotes_outputs=['EIflap']) 

    p.driver = pyOptSparseDriver()  # ScipyOptimizeDriver()
    # p.driver.options['optimizer'] = 'SLSQP'  # 'SLSQP' 'SNOPT' 'CONMIN'
    # p.driver.options['debug_print'] = ['desvars', 'objs']
    # p.driver.opt_settings['MAXIT'] = 10    # Maximum Iterations
    # p.driver.opt_settings['ACC'] = 1e-5     # Convergence Accurancy
    # p.driver.opt_settings['IPRINT'] = 1     # Output Level (<0 - None, 0 - Screen, 1 - File)
    # p.driver.opt_settings['IOUT'] = 6       # Output Unit Number
    # p.driver.opt_settings['IFILE'] = 'SLSQP.out'    # Output File Name

    p.driver.options['optimizer'] = "SNOPT"
    p.driver.opt_settings['Major feasibility tolerance'] = 1e-7
    p.driver.opt_settings['Major iterations limit'] = 10
    p.driver.opt_settings['Summary file'] = 'SNOPT_Summary_file.txt'
    p.driver.opt_settings['Print file'] = 'SNOPT_Print_file.txt'
    p.driver.opt_settings['Major step limit'] = 0.01


    p.model.add_design_var('pos_web0', lower=-0.02, upper=0.15)
    p.model.add_objective('EIflap', scaler = -1)  #/beam_prop_init['beam_stiff'][0][4,4])  #, adder = -beam_prop_init['beam_stiff'][0][4,4])  # 1./23511118.672)  # 1./2.445646288e+09
    # p.model.add_objective('EIflap', scaler = 1./23511118.672)  # 1./2.445646288e+09 

    p.setup()
    p.model.approx_totals()
    p.run_driver() 

    p.model.list_inputs(values = True, hierarchical=False)
    p.model.list_outputs(values = True, hierarchical=False)


    # # ===== PLOTS ===== #
    # # saves figures in folder_str/figures if savepath is provided:
    # job.blade_plot_sections(attribute=attribute_str, plotTheta11=flag_plotTheta11, plotDisplacement=False, savepath=folder_str)