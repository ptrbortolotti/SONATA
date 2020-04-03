# -*- coding: utf-8 -*-
"""
Created on Thursday Oct 10 10:52:28 2019

@author: Roland Feil
"""
# ============================================= #
"""
Runs OpenFAST or BeamDyn

Outputs: 
    - Plots -

"""
# ============================================= #


import os
import numpy as np
import subprocess
import SONATA.utl_openfast.fast_out_utilities as wtc_utilities


def utl_run_fast():

    # Mac
    # olddir = os.getcwd()
    os.chdir(analysis_dir)
    # FNULL = open(os.devnull, 'w')
    print('Run BeamDyn')
    cmd = [BeamDyn_driver, input_file]
    stdout = subprocess.run(cmd, stdout=subprocess.PIPE).stdout.decode('utf-8')
    print(stdout)
    print('Wohoo - BeamDyn successfully executed!')

    # ===== Read output ===== #
    FAST_IO = wtc_utilities.FAST_IO()
    allinfo, alldata = FAST_IO.load_output([analysis_dir + '/' + input_file[:-4] + '.out'])


    # ===== PLOT ===== #
    cases = {}
    cases['RootForces'] = ['RootFxr', 'RootFyr', 'RootFzr']
    cases['RootMoments'] = ['RootMxr', 'RootMyr', 'RootMzr']
    cases['TipTransDefl'] = ['TipTDxr', 'TipTDyr', 'TipTDzr']
    cases['TipAngDefl'] = ['TipRDxr', 'TipRDyr', 'TipRDzr']
    flag_save_figs = False
    FAST_IO.plot_fast_out(cases, allinfo, alldata, flag_save_figs=flag_save_figs)


if __name__ == '__main__':

    BeamDyn_driver = '/Users/rfeil/work/2_OpenFAST/openfast_BDloads_output/build/modules-local/beamdyn/beamdyn_driver'
    FAST_driver = '/Users/rfeil/work/2_OpenFAST/openfast_BDloads_output/build/glue-codes/openfast/openfast'
    analysis_dir = '/Users/rfeil/work/6_SONATA/SONATA/jobs/RFeil/00_analysis/analysis'
    input_file = 'bd_driver.inp'

    utl_run_fast(BeamDyn_driver, FAST_driver, analysis_dir, input_file)