#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 14:12:19 2018

@author: gu32kij
"""

import numpy as np
import copy 
from openmdao.api import Group, IndepVarComp

from jobs.VariSpeed.cbm_explcomp import CBM_ExplComp
from jobs.VariSpeed.marc_explcomp import MARC_ExplComp

class MDAO_Group(Group):
    '''Generate Group of two Components'''
    def __init__(self, config, **kw):
        super().__init__()
        self.config = config
        #self.ref_config = copy.deepcopy(config) 
        #self.ref_dct = {}
        self.kw = kw
    
    def setup(self):
        
        ivc = self.add_subsystem('ivc', IndepVarComp(), promotes=['*'])
        ivc.add_output('s_w1', 0.44)
        ivc.add_output('s_w2', 0.3)
        ivc.add_output('t_erosion', 0.91)
        ivc.add_output('t_overwrap',0.25)  
        ivc.add_output('t_spar1', 3.00)
        ivc.add_output('rho_mat3', 0.05)  
        ivc.add_output('t_sparcap1',2.5)
        ivc.add_output('t_sparcap2', 1.833)
        ivc.add_output('t_sparcap3', 1.833)
        ivc.add_output('t_sparcap4', 1.842) 
        ivc.add_output('rho_mat11', 0.05)  
        #ivc.add_output('rho_1', 0.05)

        promo_inputlst = ['s_w1','s_w2','t_erosion','t_overwrap','t_spar1', 'rho_mat3', 't_sparcap1', 't_sparcap2', 't_sparcap3', 't_sparcap4', 'rho_mat11']
        self.add_subsystem('cbm_comp', CBM_ExplComp(self.config, self.kw), promotes_inputs=promo_inputlst, promotes_outputs=['BeamPropSec'])
        
        #Generate MARC Subsystem
        self.add_subsystem('marc_comp', MARC_ExplComp(), promotes_inputs=['BeamPropSec'])