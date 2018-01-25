# -*- coding: utf-8 -*-
"""
This is currently the main excecution script for SONATA.

SONATA is a preprocessor for parametric analysis and design of composite beam 
cross-sections in a multidisciplinary rotor design environment. A helicopter 
rotor blade represents a classical aeroelastic problem, where the aerodynamic 
behavior, the structural elasticity and vibrational dynamics have to be studied 
simultaneously.  While a geometric definition of a rotorblade with CAD tools is 
simple, the transfer to a meshed cross-sectional representation may prohibit 
automated design optimization. Consequently, most researches have developed 
individual parametric mesh generators for the cross-sectional analysis, that 
reduces their structural model to few design variables in the process. 
SONATA represents such a preprocessor. SONATA is written in python and is using
for a lot of operations the Opencascade (CAD) kernel with its python wrapper 
(pythonocc).

In the future, the SONATA execution script shall inlude an openmdao structure 
which can call the then unlying functionalities of SONATA.CBM (Composite Beam Model)

Date: 01/02/2017
@author: TPflumm
"""
import matplotlib.pyplot as plt
from datetime import datetime
from multiprocessing import Pool
from shutil import copyfile

from SONATA.fileIO.configuration import Configuration
from SONATA.fileIO.readinput import read_material_input
from SONATA.fileIO.hiddenprints import HiddenPrints

from SONATA.mdao.cbm_explcomp import CBM_ExplComp

from openmdao.api import Problem, ScipyOptimizer, IndepVarComp, ExplicitComponent, SimpleGADriver
    


def f(tpl):
    print tpl
    (radial_station, src, cnt) = tpl
    dst = src[:-6]+str(cnt)+src[-6:]
    copyfile(src, dst)
    config = Configuration(dst)
        
    p = Problem()
    #Generate independentVariableComponent
    ivc = p.model.add_subsystem('ivc', IndepVarComp())
    ivc.add_output('wp1', 0.38)
    ivc.add_output('wp2', 0.62)
    ivc.add_output('spar_lt', 2.2, units='mm')
    ivc.add_output('skin_lt', 1.2, units='mm')
    ivc.add_output('rho_1', 0.05)
    
    #Generate Group of two Components
    p.model.add_subsystem('cbm_comp', CBM_ExplComp(config))
    p.model.connect('ivc.wp1', 'cbm_comp.WEB_Pos1')
    p.model.connect('ivc.wp2', 'cbm_comp.WEB_Pos2')
    p.model.connect('ivc.spar_lt', 'cbm_comp.Spar_layer_thickness')
    p.model.connect('ivc.skin_lt', 'cbm_comp.Skin_layer_thickness')
    p.model.connect('ivc.rho_1', 'cbm_comp.Core1_density')
    
    p.model.add_design_var('ivc.wp1', lower=0.22, upper=0.38)
    p.model.add_design_var('ivc.wp2', lower=0.6, upper=0.8)
    p.model.add_design_var('ivc.spar_lt', lower=0.3, upper=2.2)
    p.model.add_design_var('ivc.skin_lt', lower=0.3, upper=1.6)
    p.model.add_design_var('ivc.rho_1', lower=0.05, upper=19.25)
    
    p.model.add_objective('cbm_comp.obj1')
    
    p.model.add_constraint('cbm_comp.Xm2', lower=-20.0, upper=2.0)
    p.model.add_constraint('cbm_comp.GJ', lower=1.2e11, upper=1.6e11)
    p.model.add_constraint('cbm_comp.EI2', lower=1.3e11, upper=1.4e11)
    p.model.add_constraint('cbm_comp.EI3', lower=3.2e12, upper=3.6e12)

    #Setup the Problem
    p.driver = ScipyOptimizer()
    p.driver.options['optimizer'] = 'COBYLA'
    p.driver.options['disp'] = True
    p.driver.options['tol'] = 1e-4
    p.driver.options['maxiter'] = 2
    p.driver.opt_settings['rhobeg'] = 0.1 
 
    
#    p.driver = SimpleGADriver()
#    p.driver.options['bits'] = {'ivc.wp1' : 8}
#    p.driver.options['bits'] = {'ivc.wp2' : 8}
##    p.driver.options['bits'] = {'ivc.spar_lt' : 4}
##    p.driver.options['bits'] = {'ivc.skin_lt' : 4}
#    p.driver.options['max_gen'] = 20

    

    #p.driver.options['optimizer'] = 'COBYLA'
    p.setup()
    p.run_driver()
    print p['cbm_comp.obj1'], (p['ivc.wp1'], p['ivc.wp2'], p['ivc.spar_lt'], p['ivc.skin_lt'])
    return p.model.cbm_comp.job
    
if __name__ == '__main__':
    
    
    pool = Pool(processes=2)
    filename = 'jobs/VariSpeed/01_simple/sec_config.input'
    
    test = [[2500,filename,1],[3000,filename,2]]
    
    job = pool.map(f,test)
    
