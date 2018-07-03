################################################################################
#                                                                              #
#                                                                              #
#                        Author:  Willem Garre                                 #
#                        File:    marc_explcomp.py                             #
#                        Date:    24.01.2018                                   #
#                        Version: V1.0                                         #
#                                                                              #
#                                                                              #
################################################################################
# Preamble
from __future__ import division, print_function
from openmdao.core.explicitcomponent import ExplicitComponent
#import os

import numpy as np
#from datetime import datetime
import SONATA.Pymore.utl.optimization as obj

from SONATA.Pymore.marc.marc import MARC
from SONATA.Pymore.utl.plot import fan_plot

import SONATA.Pymore.utl.coef as coef

#%%
class MARC_ExplComp(ExplicitComponent):
    """
    Evaluates the fanplot of the MARC rotor system (UH60A as default)
    Or: the consecutive measure for placement of rotor frequencies
    
    Later: measures from time-domain calculations (dynamic analysis) in addition
    """
    def initialize(self):
        # self.job = MARC()
        dir_root = 'SONATA/Pymore/dym/mdl/03_rotormodel/05_UH60_rotor_optimization/01_UH60_rotor_snglblade_static/'
    #    dir_root = '/work/WGarre/02_Projects/01_pymore/11_dymore/00_mdl/03_rotormodel/03_UH60_rotor_hover/05_UH60_rotor_fourblade_3Dinfl_baseline/'

        self.job = MARC(dir_root, 'rotor_assembly.dym')
        
        self.classEval = 0
        
        #-------------------------------------#
        #   prepare permutations              #
        #-------------------------------------#
        nbOfLoc = 5
        self.RPM_vec = np.linspace(4.3*2*np.pi*0.7, 4.3*2*np.pi*1.1, nbOfLoc)

        #-------------------------------------#
        #   location to store results         #
        #-------------------------------------#        
        self.result_dir = 'SONATA/Pymore/rlt/'

    def setup(self):        
        #%%===================================#
        #   Setup                             #
        #=====================================#        
#        self.add_input('beamProp', val=np.ones([46,29]))
        self.add_input('massProp0', val=0.0)
#        self.add_input('massProp1', val=0.0)
#        self.add_input('massProp2', val=0.0)
#        self.add_input('massProp3', val=0.0)
        
        self.add_input('beamProp0', val=1.0)
        self.add_input('beamProp1', val=1.0)
        self.add_input('beamProp2', val=1.0)
        self.add_input('beamProp3', val=1.0)
        self.add_input('BeamPropSec', val=np.zeros((2,29)), desc='Massterms(6), Stiffness(21), damping(1) and coordinate(1)') 

        self.add_output('obj', val=1.0e05)

        # Finite difference all partials.
        # THIS IS USEFUL TO CALCULATE THE TRIM SOLUTION IN TIME DOMAIN
        # self.declare_partials('*', '*', method='fd')

    def compute(self, inputs, outputs):        
        #%%===================================#
        #   Compute                           #
        #=====================================#
        self.classEval = self.classEval + 1
#        print (' ** MARC **: computation initiated, fun evals: ' + str(self.classEval))   
#        beamP0 = inputs['beamProp0']
#        beamP1 = inputs['beamProp1']
#        beamP2 = inputs['beamProp2']
#        beamP3 = inputs['beamProp3']
#        
#        massP0 = inputs['massProp0']
#        massP1 = inputs['massProp0']
#        massP2 = inputs['massProp0']
#        massP3 = inputs['massProp0']
        
        #-------------------------------------#
        #   Run Simulations                   #
        #-------------------------------------#
#        starttime = datetime.now()
        
        #-------------------------------------#
        #   Change Beam Properties            #
        #-------------------------------------#        
        # Assemble Beam Properties
        beamProp = inputs['BeamPropSec']
#        beamProp = coef.refBeamProp()
        self.job.marc_set_beamProp('BLADE_BP_CG01', beamProp)

        #-------------------------------------#
        #   Change Mass Properties            #
        #-------------------------------------#            
#        self.job.marc_set_massProp('BLADE_MP_D01', massP0)
#        self.job.marc_set_massProp('BLADE_MP_E01', massP1)
#        self.job.marc_set_massProp('BLADE_MP_F01', massP2)
#        self.job.marc_set_massProp('BLADE_MP_G01', massP3)
        
        #-------------------------------------#
        #   Run Simulations                   #
        #-------------------------------------#  
        self.job.fanplot(self.RPM_vec, self.result_dir)
        
#        val_filename = '/work/WGarre/02_Projects/03_varispeed/01_UH60_rotor_validation/02_UH60_rotor_validation_fanplot/Fanplot_Bowen_Davies_Diss.csv'
#        fan_plot(np.real(freq), val_filename, RPM_vec, 'rlt/')        
        #%%===================================#
        #   Process Results                   #
        #=====================================#
        objFun = obj.gradPlacement(np.real(self.job.analysis.freq), self.RPM_vec)
#        print (' ** MARC **: Time Analysis: ' + str(datetime.now() - starttime)     )  
#        print (' ** MARC **: objective = ' +str(objFun) +', input = ' + str(massP0) +', ' + str(massP1) +', ' + str(massP2) +', ' + str(massP3) +
#                                                                 ', ' + str(beamP0) +', ' + str(beamP1) +', ' + str(beamP2) +', ' + str(beamP3) )
        # get dymore results
#        outputs['measureFreq'] = freq[1:,:]/(2*np.pi)
        outputs['obj'] = objFun
        
        
    def compute_partials(self, inputs, partials):
        #%%===================================#
        #   Compute Partials                  #
        #=====================================#
        pass