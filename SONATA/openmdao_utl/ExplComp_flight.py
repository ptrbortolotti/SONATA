#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 16:00:30 2019

@author: gu32kij
"""
import yaml
import numpy as np
from datetime import datetime
import sys

from openmdao.api import ExplicitComponent
from SONATA.cbm.fileIO.hiddenprints import HiddenPrints

from SONATA.openmdao_utl.doe_utl import filename_generator
import SONATA.Pymore.utl.coef as coef
from SONATA.Pymore.app.marc_flight_analysis import flight_analysis
from SONATA.Pymore.app.marc_reference_utl import plot_ref_data
from SONATA.Pymore.app.marc_flight_analysis_utl import interp_dui, plot_dui, \
                    plot_sensors, extract_blade_data, plot_rotor_polarcontour


class ExplComp_Flight_Analysis(ExplicitComponent):
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
        self.add_output('CT_ff', val=np.zeros((465)))  
        
    def set_partials(self):
        self.declare_partials('*', '*', method='fd', step=0.05) #finite differences all partials

    def compute(self, inputs, outputs):       
        with open('SONATA/Pymore/app/uh60a_8513test.yml', 'r') as f:
             yml = yaml.load(f.read(), Loader = yaml.FullLoader).get('pmfa')
            
        dui = interp_dui(yml.get('dui'))    
        path = yml.get('path')
        sensors = yml.get('sensors')
        
        BeamProps = {'BLADE_BP_AB01': inputs['BeamProps']}
        savepath = filename_generator(directory = '/scratch/gu32kij', string='uh60a_c8513', file_extension= '.pkl')
        data = flight_analysis(path, dui, beamprops = BeamProps, sensors=sensors, savepath=savepath)
        
        data['Psi'] = np.deg2rad(-data['SHAFT_SS_POS'][:,3]+180)
        data['CT'] = data['SHAFT_SS_HUBFORCES'][:,2]*2.4
        outputs['CT_ff'] = data['CT'][-465:]