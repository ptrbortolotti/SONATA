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

from SONATA.fileIO.configuration import Configuration
from SONATA.fileIO.readinput import read_material_input
from SONATA.cbm import CBM

from openmdao.api import ExplicitComponent

class CBM_Component(ExplicitComponent):
    """
    A simple CBM Component that computes the the composite beam properties.
    """

    def setup(self):
        
        # Inputs
        self.add_input('WEB_Pos1', val=0.38, desc='first Position (between 0 and 1)')
        self.add_input('WEB_Pos2', val=0.62, desc='second Position (between 0 and 1)')
        self.add_input('Skin_layer_thickness', val=2.5, units='mm', desc='first Position (between 0 and 1) ')
        self.add_input('Spar_layer_thickness', val=2.2, units='mm', desc='first Position (between 0 and 1) ')  
        self.add_input('Core1_density', val=2.5, units='g/cm^3', desc='first Position (between 0 and 1) ')
        self.add_input('Core2_density', val=2.2, units='g/cm^3', desc='first Position (between 0 and 1) ')  

        # Outputs
        self.add_output('Xm2', units='mm', desc='x coordinate of the mass center of gravity')
        self.add_output('Xs2', units='mm', desc='x coordinate of the shear center')
        self.add_output('EA', units='N')
        self.add_output('GJ', units='N')
        self.add_output('EI2', units='N')
        self.add_output('EI3', units='N')
        
        self.add_output('obj1', units='mm', desc='objective function1')   
        
        # Finite difference all partials.
        self.declare_partials('*', '*', method='fd')
        
    def compute(self, inputs, outputs):
                
        filename = 'jobs/VariSpeed/sec_config.input'
        config = Configuration(filename)
        MaterialLst = read_material_input(config.SETUP_mat_filename)

        config.BW_XPos = inputs['BW_XPos'][0]
        config.

        job1 = CBM(config,MaterialLst)
        job1.cbm_gen_topo()
        job1.cbm_gen_mesh()
        job1.cbm_run_vabs(filename)
        
        outputs['Xm2'] = job1.BeamProperties.Xm2
        outputs['Xs2'] = job1.BeamProperties.Xs2 
        outputs['EA'] = job1.BeamProperties.CS[0][0]
        outputs['GJ'] = job1.BeamProperties.CS[1][1]
        outputs['EI2'] = job1.BeamProperties.CS[2][2]
        outputs['EI3'] = job1.BeamProperties.CS[3][3]


if __name__ == '__main__':
    from openmdao.api import Problem, ScipyOptimizer, IndepVarComp
    
    p = Problem()
    #Generate independentVariableComponent
    ivc = p.model.add_subsystem('ivc', IndepVarComp())
    ivc.add_output('x', -33.5, units='mm')
    #ivc.add_output('y', 0.0, units='mm')
    
    #Generate Group of two Components
    p.model.add_subsystem('cbm_comp', CBM_Component())
    p.model.connect('ivc.x', 'cbm_comp.BW_XPos')
    #p.model.connect('ivc.y', 'cbm_comp.BW_YPos')
    
    p.model.add_design_var('ivc.x', lower=-34, upper=0)
    #p.model.add_design_var('ivc.y', lower=-2, upper=2)
    p.model.add_objective('cbm_comp.obj1')
    
    #Setup the Problem

    p.driver = ScipyOptimizer()
    p.driver.options['optimizer'] = 'COBYLA'
    p.driver.options['disp'] = True
    p.driver.options['tol'] = 1e-3
    p.driver.options['maxiter'] = 30


    #p.driver.options['optimizer'] = 'COBYLA'
    p.setup()
    p.run_driver()
    print p['cbm_comp.obj1'], p['ivc.x'], 

    
