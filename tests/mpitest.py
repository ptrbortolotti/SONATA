#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 18:02:35 2019

@author: gu32kij
"""

from openmdao.api import Problem, Group, IndepVarComp, SimpleGADriver
from openmdao.test_suite.components.branin import Branin

prob = Problem()
model = prob.model = Group()

model.add_subsystem('p1', IndepVarComp('xC', 7.5))
model.add_subsystem('p2', IndepVarComp('xI', 0.0))
model.add_subsystem('comp', Branin())

model.connect('p2.xI', 'comp.x0')
model.connect('p1.xC', 'comp.x1')

model.add_design_var('p2.xI', lower=-5.0, upper=10.0)
model.add_design_var('p1.xC', lower=0.0, upper=15.0)
model.add_objective('comp.f')

prob.driver = SimpleGADriver()
prob.driver.options['bits'] = {'p1.xC': 8}
prob.driver.options['max_gen'] = 10
prob.driver.options['run_parallel'] = True

prob.setup(check=False)
prob.run_driver()

# Optimal solution
print('comp.f', prob['comp.f'])