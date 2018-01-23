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

from SONATA.fileIO.configuration import Configuration
from SONATA.fileIO.readinput import read_material_input
from SONATA.cbm import CBM
from SONATA.fileIO.hiddenprints import HiddenPrints

from openmdao.api import Problem, ScipyOptimizer, IndepVarComp, ExplicitComponent, SimpleGADriver


class CBM_Component(ExplicitComponent):
    """
    A simple CBM Component that computes the the composite beam properties.
    """
    
    def setup(self):
        self.counter = 0
        self.startTime = datetime.now()
        # Inputs
        self.add_input('WEB_Pos1', val=0.38, desc='first Position (between 0 and 1)')
        self.add_input('WEB_Pos2', val=0.62, desc='second Position (between 0 and 1)')
        self.add_input('Skin_layer_thickness', val=1.2, units='mm')
        self.add_input('Spar_layer_thickness', val=2.2, units='mm')  
        self.add_input('Core1_density', val=0.05)
        self.add_input('Core2_density', val=0.05)  

        # Outputs
#        self.add_output('Xm2', units='mm', desc='x coordinate of the mass center of gravity')
#        self.add_output('Xs2', units='mm', desc='x coordinate of the shear center')
#        self.add_output('EA', units='N')
#        self.add_output('GJ', units='Nmm^2')
#        self.add_output('EI2', units='Nmm^2')
        #self.add_output('EI3')
        self.add_output('obj1', desc='objective function')   
        self.add_output('Xm2', desc='x location of Center of Gravity')
        self.add_output('Xm3', desc='y location of Center of Gravity')

        self.add_output('EA', desc='Axial Stiffness')
        self.add_output('GJ', desc='Torsional Stiffness')
        self.add_output('EI2', desc='Flapping Bending Stiffness')
        self.add_output('EI3', desc='Lagging Bending Stiffness')
        #self.add_output('obj2', desc='objective function')  
        
        # Finite difference all partials.
        #self.declare_partials('*', '*', method='fd')
        
    def compute(self, inputs, outputs):
        print '----------------------------------------------------------------'
        print  self.counter, datetime.now()-self.startTime,
        print '%.2f  ' %inputs['WEB_Pos1'][0], 
        print '%.2f  ' %inputs['WEB_Pos2'][0], 
        print '%.2f  ' %inputs['Spar_layer_thickness'][0], 
        print '%.2f  ' %inputs['Skin_layer_thickness'][0], 
        print '%.2f  ' %inputs['Core1_density'][0]
                
        filename = 'jobs/VariSpeed/01_simple/sec_config.input'
        config = Configuration(filename)
        MaterialLst = read_material_input(config.SETUP_mat_filename)

        #config.SETUP_radial_station = 2500
        config.WEB_Pos1[0] = inputs['WEB_Pos1'][0]
        config.WEB_Pos2[0] = inputs['WEB_Pos2'][0]
        config.SEG_Layup[0][0][2] = inputs['Skin_layer_thickness'][0]
        config.SEG_Layup[0][1][2] = inputs['Skin_layer_thickness'][0]
        
        config.SEG_Layup[1][0][2] = inputs['Spar_layer_thickness'][0]
        config.SEG_Layup[1][1][2] = inputs['Spar_layer_thickness'][0]
        config.SEG_Layup[1][2][2] = inputs['Spar_layer_thickness'][0]
        config.SEG_Layup[1][3][2] = inputs['Spar_layer_thickness'][0]
        
        MaterialLst[11].rho = inputs['Core1_density'][0]
        MaterialLst[12].rho = inputs['Core2_density'][0]
       
        
        
        with HiddenPrints():
            self.job = CBM(config,MaterialLst)
            self.job.cbm_gen_topo()
            self.job.cbm_gen_mesh()
            self.job.cbm_run_vabs()
        
        
        outputs['obj1'] = self.job.BeamProperties.MpUS
        
        outputs['Xm2'] = self.job.BeamProperties.Xm2
        outputs['Xm3'] = self.job.BeamProperties.Xm3
        outputs['EA']   = self.job.BeamProperties.CS[0][0]
        outputs['GJ']   = self.job.BeamProperties.CS[1][1]
        outputs['EI2']  = self.job.BeamProperties.CS[2][2]
        outputs['EI3']  = self.job.BeamProperties.CS[3][3]
        
         print outputs   
        self.counter += 1

    def post_cbm(self):
        return self.job.cbm_post_2dmesh()
    
if __name__ == '__main__':

    p = Problem()
    #Generate independentVariableComponent
    ivc = p.model.add_subsystem('ivc', IndepVarComp())
    ivc.add_output('wp1', 0.38)
    ivc.add_output('wp2', 0.62)
    ivc.add_output('spar_lt', 2.2, units='mm')
    ivc.add_output('skin_lt', 1.2, units='mm')
    ivc.add_output('rho_1', 0.05)

    #Generate Group of two Components
    p.model.add_subsystem('cbm_comp', CBM_Component())
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
    p.driver.options['maxiter'] = 1000
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
    job = p.model.cbm_comp.post_cbm()
    
