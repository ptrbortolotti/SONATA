#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 16:00:30 2019

@author: gu32kij
"""
import numpy as np
from datetime import datetime
import sys

from openmdao.api import ExplicitComponent
from SONATA.cbm.fileIO.hiddenprints import HiddenPrints
from SONATA.Pymore.app.marc_fanplot import frequency_analysis

class ExplComp_Frequency_Analysis(ExplicitComponent):
    """
    A simple Frequency Analysis ExplicitComponent that computes the the fanplot 
    for the rotor with the input beam properties of the blade
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
        self.add_input('BeamProps', val=np.zeros((11,29)), desc='Massterms(6), Stiffness(21), damping(1) and coordinate(1)')
        
    def set_output(self):
        self.add_output('freq', val=np.zeros((11,14)))   
        self.add_output('Omega') 
        self.add_output('RPM_vec', val=np.zeros((11)))
        
    def set_partials(self):
         self.declare_partials('*', '*', method='fd', step=0.05) #finite differences all partials

    def compute(self, inputs, outputs):       
        dir_root = 'SONATA/Pymore/dym/mdl/03_rotormodel/02_UH60_rotor_fanplot_noaero/04_UH60_rotor_snglblade_noaero_pitchlink/'      
        with HiddenPrints():
            (freq, Omega, RPM_vec) = frequency_analysis(dir_root=dir_root, beamProp = inputs['BeamProps'])
        outputs['freq'] = freq
        outputs['Omega'] = Omega
        outputs['RPM_vec'] = RPM_vec                                  
