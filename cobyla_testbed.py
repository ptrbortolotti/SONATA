#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 09:08:33 2018

@author: gu32kij
"""

from openmdao.api import Problem, ScipyOptimizeDriver, ExecComp, IndepVarComp

# We'll use the component that was defined in the last tutorial
from openmdao.test_suite.components.paraboloid import Paraboloid

# build the model
prob = Problem()
indeps = prob.model.add_subsystem('indeps', IndepVarComp())
indeps.add_output('x', 3.0)
indeps.add_output('y', -4.0)

prob.model.add_subsystem('parab', Paraboloid())

# define the component whose output will be constrained
prob.model.add_subsystem('const', ExecComp('g = x + y'))

prob.model.connect('indeps.x', ['parab.x', 'const.x'])
prob.model.connect('indeps.y', ['parab.y', 'const.y'])

# setup the optimization
#prob.driver = ScipyOptimizeDriver()
#prob.driver.options['optimizer'] = 'COBYLA'
#prob.driver.options['debug_print'] = ['desvars','objs']
#prob.driver.options['maxiter'] = 5
#prob.driver.options['tol'] = 1e-4
#prob.driver.opt_settings['rhobeg'] = 0.1

# setup the optimization
prob.driver = ScipyOptimizeDriver()
prob.driver.options['optimizer'] = 'SLSQP'
prob.driver.options['debug_print'] = ['desvars','objs']
prob.driver.options['maxiter'] = 150
prob.driver.options['tol'] = 1e-4
#prob.driver.opt_settings['eps'] = 0.1

prob.model.add_design_var('indeps.x', lower=-50, upper=50, ref0 = -50, ref=50)
prob.model.add_design_var('indeps.y', lower=-50, upper=50, ref0 = -50, ref=50)
prob.model.add_objective('parab.f_xy')

# to add the constraint to the model
prob.model.add_constraint('const.g', lower=0, upper=10.)
# prob.model.add_constraint('const.g', equals=0.)

prob.setup()
prob.run_driver()

# minimum value
print(prob['parab.f_xy'])

# location of the minimum
print(prob['indeps.x'])
print(prob['indeps.y'])