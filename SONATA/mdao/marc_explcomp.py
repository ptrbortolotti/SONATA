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

from openmdao.core.explicitcomponent import ExplicitComponent
import numpy as np
from datetime import datetime
from Pymore import marc

_marc = __import__('15_marc.marc', globals(), locals(), ['MARC'], -1)
MARC = _marc.MARC

#%%
class MARC_ExplComp(ExplicitComponent):
    
    def __init_(self):
        self.job = MARC()
    
    def setup(self):
        #%%===================================#
        #   Setup                             #
        #=====================================#
        self.add_input('x', val=0.0)
        self.add_input('y', val=0.0)

        self.add_output('f', val=0.0)

        # Finite difference all partials.
        self.declare_partials('*', '*', method='fd')

    def compute(self, inputs, outputs):
        #%%===================================#
        #   Compute                           #
        #=====================================#
        #-------------------------------------#
        #   prepare permutations              #
        #-------------------------------------#
        nbOfEig = self.job.analysis.sta_get_eigNb()
        RPM_vec = np.linspace(15.0, 4.3*2*np.pi, nbOfEig)
        
        #-------------------------------------#
        #   initialize storage                #
        #-------------------------------------#
        freq = np.empty([nbOfEig])
        
        #-------------------------------------#
        #   Run Simulations                   #
        #-------------------------------------#
        starttime = datetime.now()
        
        #-------------------------------------#
        #   Change Beam Properties            #
        #-------------------------------------#        
        beamProp = inputs['beamProp']
        self.job.analysis.sta_set_beamProp(np.array([beamProp]))
        
        #-------------------------------------#
        #   Run Simulations                   #
        #-------------------------------------#  
        for i in RPM_vec:
        
            # set dymore input
            self.job.analysis.sta_set_rigRot(np.array([i, 0.0, 0.0]))
        
            # run dymore single job
            self.job.marc_init()
            self.job.marc_step()
            
            # get dymore results
            result = self.job.analysis.sta_get_eig()
        
            freq = np.vstack((freq, result))
        
        print('\n *** Time Analysis: ' + str(datetime.now() - starttime))       
        
        # get dymore results
        outputs['f_xy'] = freq[1:,:]/(2*np.pi)