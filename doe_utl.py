#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 12:34:05 2019

@author: gu32kij
"""
from mpi4py import MPI

from datetime import datetime
from openmdao.api import Problem, IndepVarComp
from openmdao.test_suite.components.paraboloid import Paraboloid

from openmdao.api import DOEDriver, FullFactorialGenerator, LatinHypercubeGenerator, ListGenerator, UniformGenerator
from openmdao.api import SqliteRecorder, CaseReader

from SONATA.openmdao_utl.doe_sampler import doe_sampler, plot_samples

#%% ========================= M A I N =========================================
#if __name__ == '__main__':
samples = doe_sampler(1000, des_vars = ['x','y'], means=[1,0], stdvs= [0.5,1])

s = datetime.now().isoformat(sep='_',timespec='minutes')
jobid =  s.replace(':','').replace('.','')
dbname =  '/scratch/doe_cases_'+jobid+'.sql'

prob = Problem()
model = prob.model
    
model.add_subsystem('p1', IndepVarComp('x', 0.0), promotes=['x'])
model.add_subsystem('p2', IndepVarComp('y', 0.0), promotes=['y'])
model.add_subsystem('comp', Paraboloid(), promotes=['x', 'y', 'f_xy'])

model.add_design_var('x', lower=0.0, upper=1.0)
model.add_design_var('y', lower=0.0, upper=1.0)
model.add_objective('f_xy')

#prob.driver = DOEDriver(ListGenerator(samples))
prob.driver = DOEDriver(ListGenerator(samples))
prob.driver.options['run_parallel'] =  True
prob.driver.options['procs_per_model'] =  1
prob.driver.add_recorder(SqliteRecorder(dbname))

prob.setup()
prob.run_driver()
prob.cleanup()

    #cr = CaseReader(dbname)
    #cases = cr.list_cases('driver')
    
    #plot_samples(samples, upper_tri = False)
    #(outputs, values) = read_doe_data([dbname])