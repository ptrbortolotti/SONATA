# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 13:54:37 2018

@author: TPflumm
"""
from datetime import datetime
from openmdao.api import ExplicitComponent

from SONATA.fileIO.hiddenprints import HiddenPrints
from SONATA.cbm import CBM

class CBM_ExplComp(ExplicitComponent):
    """
    A simple CBM Component that computes the the composite beam properties.
    """
    def __init__(self, config):
        super(CBM_ExplComp, self).__init__()
        self.config = config

    
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
        self.add_output('MpUS', desc='Mass per unit span (kg/m)')   
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
        elapsed_t = datetime.now() - self.startTime
        m, s = divmod(elapsed_t.seconds, 60)
        h, m = divmod(m, 60)
        print '----------------------------------------------------------------'
        print  self.counter,
        print '%02d:%02d:%02d' % (h,m,s),
        print '%.2f  ' %inputs['WEB_Pos1'][0], 
        print '%.2f  ' %inputs['WEB_Pos2'][0], 
        print '%.2f  ' %inputs['Spar_layer_thickness'][0], 
        print '%.2f  ' %inputs['Skin_layer_thickness'][0], 
        print '%.2f  ' %inputs['Core1_density'][0]
                
        #SETUP A CBM JOB:
        self.job = CBM(self.config)
        
        #config.SETUP_radial_station = 2500
        self.job.config.WEB_Pos1[0] = inputs['WEB_Pos1'][0]
        self.job.config.WEB_Pos2[0] = inputs['WEB_Pos2'][0]
        self.job.config.SEG_Layup[0][0][2] = inputs['Skin_layer_thickness'][0]
        self.job.config.SEG_Layup[0][1][2] = inputs['Skin_layer_thickness'][0]
        
        self.job.config.SEG_Layup[1][0][2] = inputs['Spar_layer_thickness'][0]
        self.job.config.SEG_Layup[1][1][2] = inputs['Spar_layer_thickness'][0]
        self.job.config.SEG_Layup[1][2][2] = inputs['Spar_layer_thickness'][0]
        self.job.config.SEG_Layup[1][3][2] = inputs['Spar_layer_thickness'][0]
        
        self.job.MaterialLst[11].rho = inputs['Core1_density'][0]
        self.job.MaterialLst[12].rho = inputs['Core2_density'][0]
               
        with HiddenPrints():
            self.job.cbm_gen_topo()
            self.job.cbm_gen_mesh()
            self.job.cbm_run_vabs()
        
        
        outputs['MpUS'] = self.job.BeamProperties.MpUS
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
