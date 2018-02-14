# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 13:54:37 2018

@author: TPflumm
"""
import numpy as np
from datetime import datetime
from openmdao.api import ExplicitComponent

from SONATA.cbm.fileIO.hiddenprints import HiddenPrints
from SONATA.cbm.sonata_cbm import CBM

class CBM_ExplComp(ExplicitComponent):
    """
    A simple CBM Component that computes the the composite beam properties.
    """
    def __init__(self, config):
        super().__init__()
        self.config = config
        self.ref_dct = {}

    
    def setup(self):
        self.counter = 0
        self.startTime = datetime.now()

        self.set_input()
        self.set_output()

        # Finite difference all partials.
        #self.declare_partials('*', '*', method='fd')
    
    
    def set_input(self):
        self.add_input('WEB_Pos1', val=0.38, desc='first Position (between 0 and 1)')
        self.add_input('WEB_Pos2', val=0.62, desc='second Position (between 0 and 1)')
        self.add_input('Skin_layer_thickness', val=1.2, units='mm')
        self.add_input('Spar_layer_thickness', val=2.2, units='mm')  
        self.add_input('Core1_density', val=0.05)
        self.add_input('Core2_density', val=0.05)  


    def connect_input_to_config(self,inputs):
        self.job.config.WEB_Pos1[0] = inputs['WEB_Pos1'][0]
        self.job.config.WEB_Pos2[0] = inputs['WEB_Pos2'][0]
        self.job.config.SEG_Layup[0][0][2] = inputs['Skin_layer_thickness'][0]
        self.job.config.SEG_Layup[0][1][2] = inputs['Skin_layer_thickness'][0]
        
        self.job.config.SEG_Layup[1][0][2] = inputs['Spar_layer_thickness'][0]
        self.job.config.SEG_Layup[1][1][2] = inputs['Spar_layer_thickness'][0]
        self.job.config.SEG_Layup[1][2][2] = inputs['Spar_layer_thickness'][0]
        self.job.config.SEG_Layup[1][3][2] = inputs['Spar_layer_thickness'][0]
        
        #self.job.MaterialLst[11].rho = inputs['Core1_density'][0]
        #self.job.MaterialLst[12].rho = inputs['Core2_density'][0]


    def set_output(self):
        self.add_output('obj', desc='objective_function')   
        #self.add_output('obj2', desc='objective function')  
        self.add_output('MpUS', desc='Mass per unit span (kg/m)')   
        self.add_output('Xm2', desc='x location of Center of Gravity')
        self.add_output('Xm3', desc='y location of Center of Gravity')
        self.add_output('Xs2', desc='x location of Shear Center')
        self.add_output('Xs3', desc='y location of Shear Center')
        self.add_output('EA', desc='Axial Stiffness')
        self.add_output('GJ', desc='Torsional Stiffness')
        self.add_output('EI2', desc='Flapping Bending Stiffness')
        self.add_output('EI3', desc='Lagging Bending Stiffness')


    def connect_output_from_job(self, outputs):
        outputs['MpUS'] = self.job.BeamProperties.MpUS
        outputs['Xm2'] = self.job.BeamProperties.Xm2
        outputs['Xm3'] = self.job.BeamProperties.Xm3
        outputs['EA']   = self.job.BeamProperties.CS[0][0]
        outputs['GJ']   = self.job.BeamProperties.CS[1][1]
        outputs['EI2']  = self.job.BeamProperties.CS[2][2]
        outputs['EI3']  = self.job.BeamProperties.CS[3][3]
        outputs['EI3']  = self.job.BeamProperties.MpUS
        outputs['obj'] = self.compute_objective()
        
    
    def set_references(self,ref_dct):
        self.ref_dct = ref_dct


    def compute_objective(self):
        o1 = abs(self.job.BeamProperties.CS[2][2]*1e-6 - self.ref_dct['bending_stiffnesses'][0]) / self.ref_dct['bending_stiffnesses'][0]
        o2 = abs(self.job.BeamProperties.CS[3][3]*1e-6 - self.ref_dct['bending_stiffnesses'][1]) / self.ref_dct['bending_stiffnesses'][1]
        o3 = abs(self.job.BeamProperties.CS[1][1]*1e-6 - self.ref_dct['torsional_stiffness']) / self.ref_dct['torsional_stiffness']
        o4 = abs(self.job.BeamProperties.CS[0][0] - self.ref_dct['axial_stiffness']) / self.ref_dct['axial_stiffness']
        o5 = abs(self.job.BeamProperties.MpUS - self.ref_dct['mass_per_unit_span']) / self.ref_dct['mass_per_unit_span']
              
        obj_arr = np.hstack((o1**2,o2**2,o3**2,o4**2,o5**2))
        return np.average(obj_arr)


    def compute(self, inputs, outputs):
        elapsed_t = datetime.now() - self.startTime
        m, s = divmod(elapsed_t.seconds, 60)
        h, m = divmod(m, 60)
        if self.counter == 0:
            print('--time----wp1---wp2---spar_lt-skin_lt--obj---------------------------------')
        print(self.counter, end=' ')
        print('%02d:%02d:%02d ' % (h,m,s), end=' ')
        print('%.2f, ' %inputs['WEB_Pos1'][0], end=' ') 
        print('%.2f, ' %inputs['WEB_Pos2'][0], end=' ') 
        print('%.2f, ' %inputs['Spar_layer_thickness'][0], end=' ') 
        print('%.2f, ' %inputs['Skin_layer_thickness'][0], end=' ')
        #print 'rho_1%.2f  ' %inputs['Core1_density'][0]
        #SETUP A CBM JOB:
        self.job = None
        self.job = CBM(self.config)
        self.connect_input_to_config(inputs)
        with HiddenPrints():
            self.job.cbm_gen_topo()
            self.job.cbm_gen_mesh()
            self.job.cbm_run_vabs()

        self.connect_output_from_job(outputs)
        print(outputs['obj'])
        self.counter += 1

    def post_cbm(self):
        return self.job.cbm_post_2dmesh()
