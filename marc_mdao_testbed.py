""" Definition of the Paraboloid component, which evaluates the equation
(x-3)^2 + xy + (y+4)^2 = 3
"""
from __future__ import division, print_function
from SONATA.Pymore.marc_explcomp import MARC_ExplComp
from CBM_dummy_explcomp import CBM_dummy
from openmdao.core.explicitcomponent import ExplicitComponent
from openmdao.api import Problem, ScipyOptimizeDriver, ExecComp, IndepVarComp, Group, SimpleGADriver

import numpy as np

if __name__ == "__main__":
#    from openmdao.core.problem import Problem
#    from openmdao.core.group import Group
#    from openmdao.core.indepvarcomp import IndepVarComp
    import utl.coef as coef

    beamProp0 = coef.refBeamProp()

    model = Group()
    ivc = IndepVarComp()
#    ivc.add_output('m', 0.52705078)
#    ivc.add_output('m', 26.17265625)
    ivc.add_output('m0', 0.0)
#    ivc.add_output('m1', 0.0)
#    ivc.add_output('m2', 0.0)
#    ivc.add_output('m3', 0.0)
    
    ivc.add_output('k0', 1.0)
    ivc.add_output('k1', 1.0)
    ivc.add_output('k2', 1.0)
    ivc.add_output('k3', 1.0)
    
    model.add_subsystem('des_vars', ivc)
#    model.add_subsystem('CBM_dummy', CBM_dummy())
    model.add_subsystem('MARC', MARC_ExplComp())

#    model.connect('des_vars.m', 'CBM_dummy.m')
#    model.connect('CBM_dummy.beamProp', 'MARC.beamProp')
    model.connect('des_vars.m0', 'MARC.massProp0')
#    model.connect('des_vars.m1', 'MARC.massProp1')
#    model.connect('des_vars.m2', 'MARC.massProp2')
#    model.connect('des_vars.m3', 'MARC.massProp3')
    
    model.connect('des_vars.k0', 'MARC.beamProp0')
    model.connect('des_vars.k1', 'MARC.beamProp1')
    model.connect('des_vars.k2', 'MARC.beamProp2')
    model.connect('des_vars.k3', 'MARC.beamProp3')

    prob = Problem(model) 
    
    print (' ** mdao **: prepare optimization')
    # setup the optimization

    prob.driver = SimpleGADriver()
    prob.driver.options['max_gen'] = 40
    prob.driver.options['pop_size'] = 80
    prob.driver.options['bits'] = {'des_vars.m0' : 16}
#    prob.driver.options['bits'] = {'des_vars.m1' : 16}
#    prob.driver.options['bits'] = {'des_vars.m2' : 16}
#    prob.driver.options['bits'] = {'des_vars.m3' : 16}
    
    prob.driver.options['bits'] = {'des_vars.k0' : 16}
    prob.driver.options['bits'] = {'des_vars.k1' : 16}
    prob.driver.options['bits'] = {'des_vars.k2' : 16}
    prob.driver.options['bits'] = {'des_vars.k3' : 16}

#    prob.driver = ScipyOptimizeDriver()    
#    [Nelder-Mead, Powell, CG, BFGS, Newton-CG, L-BFGS-B, TNC, COBYLA, SLSQP] 	
#    prob.driver.options['optimizer'] = 'Nelder-Mead'
#    prob.driver.options['maxiter'] = 100
#    prob.driver.options['tol'] = 1e-3
#    prob.driver.options['disp'] = True
    
    prob.model.add_design_var('des_vars.k0', lower=0.5, upper=2.0)
    prob.model.add_design_var('des_vars.k1', lower=0.5, upper=2.0)
    prob.model.add_design_var('des_vars.k2', lower=0.5, upper=2.0)
    prob.model.add_design_var('des_vars.k3', lower=0.5, upper=2.0)
    
    prob.model.add_design_var('des_vars.m0', lower=0.0, upper=1.0)
#    prob.model.add_design_var('des_vars.m1', lower=0.0, upper=1.0)
#    prob.model.add_design_var('des_vars.m2', lower=0.0, upper=1.0)
#    prob.model.add_design_var('des_vars.m3', lower=0.0, upper=1.0)
    
    prob.model.add_objective('MARC.measureFreq')
    
    # to add the constraint to the model
    # prob.model.add_constraint('const.g', lower=0, upper=10.)
    # prob.model.add_constraint('const.g', equals=0.)
    
    print (' ** mdao **: run optimization driver: ')# + prob.driver.options['optimizer'])
    prob.setup()
#    prob.run_model()
    prob.run_driver()
    print(prob['MARC.measureFreq'])
#    print(prob['des_vars.m0'], prob['des_vars.m1'], prob['des_vars.m2'], prob['des_vars.m3'], prob['des_vars.k0'], prob['des_vars.k1'], prob['des_vars.k2'], prob['des_vars.k3'])
    print(prob['des_vars.m0'], prob['des_vars.k0'], prob['des_vars.k1'], prob['des_vars.k2'], prob['des_vars.k3'])