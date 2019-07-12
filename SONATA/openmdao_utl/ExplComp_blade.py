#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 13:34:08 2019

@author: gu32kij
"""
import numpy as np
from datetime import datetime
import sys

from openmdao.api import ExplicitComponent
from SONATA.cbm.fileIO.hiddenprints import HiddenPrints
from SONATA.classBlade import Blade
#from SONATA.classAirfoil import Airfoil
#from SONATA.classMaterial import read_IEA37_materials
import SONATA.Pymore.utl.coef as coef

class ExplComp_Blade(ExplicitComponent):
    
    """
    A simple Blade ExplicitComponent that computes the the composite beam properties of the blade
    """
    def __init__(self):
        super().__init__()

    def setup(self):
        self.counter = 0
        self.startTime = datetime.now()
        self.set_input()
        self.set_output()
        self.set_partials()
        
    def set_input(self):        
        self.add_input('m1_E', val=np.array([139.360e9, 12.615e9, 12.615e9]), desc='material 1 - elastic modulus')  
        self.add_input('m1_rho', val=1.536e3, desc='material 1 - density')
 
#        self.add_input('m3_E', val=np.array([139.360e9, 12.615e9, 12.615e9]))  
#        self.add_input('m3_rho', val=1.572e3)
#
#        self.add_input('m3_E', val=2.05)
#        self.add_input('m3_rho', val=2.05)
#       
#        self.add_input('t_box1', val=1.35)
#        self.add_input('t_box2', val=1.35)
#        self.add_input('t_box3', val=1.45)
#        self.add_input('t_box4', val=0.50) 
#
#        self.add_input('t_skin', val=0.25)
#        self.add_input('t_erosion', val=0.82)
#        self.add_input('s_w1', val=0.43)
#        self.add_input('s_w2', val=0.57)
#        
    def set_output(self):
        self.add_output('BeamProps', val=np.zeros((11,29)), desc='Massterms(6), Stiffness(21), damping(1) and coordinate(1)')   
        
    def set_partials(self):
         self.declare_partials('*', '*', method='fd', step=0.05) #finite differences all partials
#        self.declare_partials('BeamProps', 'rho_mat3', method='fd', step = 2e-2)
#        self.declare_partials('BeamProps', 't_sparcap1', method='fd', step = 1e-2)
#        self.declare_partials('BeamProps', 't_sparcap2', method='fd', step = 1e-2)
#        self.declare_partials('BeamProps', 't_sparcap3', method='fd', step = 1e-2)
#        self.declare_partials('BeamProps', 't_sparcap4', method='fd', step = 1e-2)
#        self.declare_partials('BeamProps', 's_w1', method='fd', step = 1e-2)
#        self.declare_partials('BeamProps', 's_w2', method='fd', step = 1e-2)
#        self.declare_partials('BeamProps', 's_spar2', method='fd', step = 1e-2)
        
    def compute(self, inputs, outputs):
#        elapsed_t = datetime.now() - self.startTime
#        m, s = divmod(elapsed_t.seconds, 60)
#        h, m = divmod(m, 60)
#        if self.counter == 0:
#            print('--time',  end=' ')  
#            for k in inputs:
#                print('--%s, ' %k, end=' ')  
#            print('')
#        print(('%2i' % self.counter), end=' ')
#        print('%02d:%02d:%02d [' % (h,m,s), end=' ')
        with HiddenPrints():
            self.job = Blade(name='UH-60A', filename = 'jobs/MonteCarlo/UH-60A_adv.yml')
        self.connect_input_to_config(inputs)

        #SETUP SONATA-CBM JOB:        
        try:
            with HiddenPrints():
                self.job.blade_gen_section(mesh_flag = True, split_quads=False)
                print(inputs['m1_E'], self.job.materials[1].E)
            self.job.blade_run_vabs(ramdisk=True)
            beam = self.job.blade_exp_beam_props(solver='vabs', cosy='local', eta_offset=0.1)
            beamProps = coef.join_beam_props(beam, coef.refBeamProp())
            outputs['BeamProps'] = beamProps
                                       
        except KeyboardInterrupt:
            raise Exception    
        
        except:
           print('] [Unexpected error:', sys.exc_info()[0], ']')
        self.counter += 1   
        
        
    def connect_input_to_config(self, inputs):
        
        #Materialprops:
        self.job.materials[1].E = inputs['m1_E']
        self.job.materials[1].rho = inputs['m1_rho']
        
        #Architecture:
#        self.job.config.webs[1]['Pos1'] = inputs['s_w1'][0]
#        self.job.config.webs[1]['Pos2'] = 1-self.job.config.webs[1]['Pos1']

#        self.job.config.webs[2]['Pos1'] = inputs['s_w2'][0]
#        self.job.config.webs[2]['Pos2'] = 1-self.job.config.webs[2]['Pos1']
#        
#        self.job.config.segments[0]['Layup'][6][0] = inputs['s_spar2'][0]
#        self.job.config.segments[0]['Layup'][6][1] = 1-inputs['s_spar2'][0]
        
        #Segment 0 :
        #self.job.config.segments[0]['Layup'][0][2] = inputs['t_erosion'][0]
#        self.job.config.segments[0]['Layup'][1][2] = inputs['t_overwrap'][0]
#        self.job.config.segments[0]['Layup'][2][2] = inputs['t_overwrap'][0]
#        self.job.config.segments[0]['Layup'][3][2] = inputs['t_overwrap'][0]
#        self.job.config.segments[0]['Layup'][4][2] = inputs['t_overwrap'][0]
        
        #Segment 1:
        #self.job.config.segments[1]['Layup'][0][2] = inputs['t_spar1'][0]
        
        #Segment 2:
#        self.job.config.segments[2]['Layup'][0][2] = inputs['t_sparcap1'][0]
#        self.job.config.segments[2]['Layup'][1][2] = inputs['t_sparcap2'][0]
#        self.job.config.segments[2]['Layup'][2][2] = inputs['t_sparcap3'][0]
#        self.job.config.segments[2]['Layup'][3][2] = inputs['t_sparcap4'][0]
#        
        #self.job.config.segments[2]['Layup'][3][4] = inputs['o_spar_cap1']
        
        #Segment 3:
#        self.job.MaterialLst[10].rho = inputs['rho_mat11'][0]/1000
#        self.job.MaterialLst[2].rho = inputs['rho_mat3'][0]/1000

    def post_cbm(self):
        return self.job.blade_plot_sections(plotTheta11=False)