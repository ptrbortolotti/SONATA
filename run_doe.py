#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 12:34:05 2019

@author: gu32kij
"""#
import os
from mpi4py import MPI
import numpy as np
from openmdao.api import Problem, IndepVarComp
from openmdao.test_suite.components.paraboloid import Paraboloid

from openmdao.api import DOEDriver, FullFactorialGenerator, LatinHypercubeGenerator, ListGenerator, UniformGenerator
from openmdao.api import SqliteRecorder

from SONATA.openmdao_utl.doe_utl import doe_sampler
from SONATA.openmdao_utl.ExplComp_blade import ExplComp_Blade
from SONATA.openmdao_utl.ExplComp_frequency_analysis import ExplComp_Frequency_Analysis
from SONATA.openmdao_utl.ExplComp_hover import ExplComp_Hover_Analysis
from SONATA.openmdao_utl.ExplComp_flight import ExplComp_Flight_Analysis
from SONATA.openmdao_utl.doe_utl import filename_generator, directory_creator

#=====================Determine recorder filename=====================
#dirName = directory_creator('/scratch/gu32kij', string='doe_testcases', MPI=True)
dirName = '/scratch/gu32kij/2019-07-12_test'
recorder_fname = filename_generator(directory = dirName, string='cases', file_extension= '.sql')

#=====================prepare Samples=====================
nsamples = 100
desvar_dict = {}
#Material 1
desvar_dict.update({  'm1_E'  : {'means':np.array([139e9, 12.615e9])},
                      'm1_G_lt': {'means':np.array([5.892e9])},
                      'm1_rho': {'means':np.array([1.536e3])}})

desvar_dict.update({  'm3_E'  : {'means':np.array([177.760e9, 12.615e9])},
                      'm3_G_lt': {'means':np.array([5.892e9])},
                      'm3_rho': {'means':np.array([1.572e3])}})

desvar_dict['m1_E']['stdvs'] = [0.07, 0.04]*desvar_dict['m1_E']['means']
desvar_dict['m1_G_lt']['stdvs'] = 0.12*desvar_dict['m1_G_lt']['means']
desvar_dict['m1_rho']['stdvs'] = 0.05*desvar_dict['m1_rho']['means']

desvar_dict['m3_E']['stdvs'] = [0.07, 0.04]*desvar_dict['m3_E']['means']
desvar_dict['m3_G_lt']['stdvs'] = 0.12*desvar_dict['m3_G_lt']['means']
desvar_dict['m3_rho']['stdvs'] = 0.05*desvar_dict['m3_rho']['means']

#Thickness of box spar layers
#Web 1 and 2 Position

samples = doe_sampler(desvar_dict, nsamples=nsamples)

#=====================Define Problem=====================
prob = Problem()
model = prob.model

ivc = IndepVarComp()
ivc.add_output('m1_E', val=np.array([139.360e9, 12.615e9]), desc='material 1 - elastic modulus')
ivc.add_output('m1_G_lt', val=np.array([5.892e9]), desc='material 1 - shear modulus')
ivc.add_output('m1_rho', val=1.536e3, desc='material 1 - density')
ivc.add_output('m3_E', val=np.array([177.760e9, 12.615e9]), desc='material 3 - elastic modulus')
ivc.add_output('m3_G_lt', val=np.array([5.892e9]), desc='material 3 - shear modulus')
ivc.add_output('m3_rho', val=1.572e3, desc='material 3 - density')

model.add_subsystem('ivc', ivc, promotes_outputs=['m1_E', 'm1_G_lt', 'm1_rho', 'm3_E', 'm3_G_lt', 'm3_rho'])
model.add_subsystem('blade', ExplComp_Blade(), promotes_inputs=['m1_E', 'm1_G_lt', 'm1_rho', 'm3_E', 'm3_G_lt', 'm3_rho'], promotes_outputs=['BeamProps'])
model.add_subsystem('marc_freqanalysis', ExplComp_Frequency_Analysis(), promotes_inputs=['BeamProps'], promotes_outputs=['freq', 'Omega', 'RPM_vec'])
model.add_subsystem('marc_hover', ExplComp_Hover_Analysis(), promotes_inputs=['BeamProps'], promotes_outputs=['steady_mean_elastic_tip_response', 'steady_stdvs_elastic_tip_response', 'steady_mean_bearing_response', 'steady_stdvs_bearing_response' ])
model.add_subsystem('marc_flight', ExplComp_Flight_Analysis(), promotes_inputs=['BeamProps'], promotes_outputs=['mean_elastic_tip_response', 'stdvs_elastic_tip_response', 'vibratory_hubloads'])

model.add_design_var('m1_E')
model.add_design_var('m1_G_lt')
model.add_design_var('m1_rho')
model.add_design_var('m3_E')
model.add_design_var('m3_G_lt')
model.add_design_var('m3_rho')

model.add_objective('BeamProps')
model.add_objective('freq')
model.add_objective('Omega')
model.add_objective('mean_elastic_tip_response')
model.add_objective('stdvs_elastic_tip_response')
model.add_objective('vibratory_hubloads')
model.add_objective('steady_mean_elastic_tip_response')
model.add_objective('steady_stdvs_elastic_tip_response')
model.add_objective('steady_mean_bearing_response')
model.add_objective('steady_stdvs_bearing_response')

prob.driver = DOEDriver(ListGenerator(samples))
prob.driver.options['run_parallel'] =  True
prob.driver.options['procs_per_model'] =  1
prob.driver.add_recorder(SqliteRecorder(recorder_fname))

#=====================run driver=====================
prob.setup()
prob.run_driver()
#prob.run_model()
prob.cleanup()