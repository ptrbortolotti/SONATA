#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 16:01:55 2018

@author: gu32kij
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import os
from concurrent import futures
from time import sleep, time
import copy
from multiprocessing import Pool
import functools
from tqdm import tqdm


if __name__ == '__main__':  
    os.chdir('../..')
from SONATA.cbm.fileIO.configuration import Configuration
from SONATA.cbm.sonata_cbm import CBM
from SONATA.cbm.fileIO.hiddenprints import HiddenPrints

from SONATA.cbm.fileIO.dymore_utils import read_dymore_beam_properties, interp1d_dymore_beam_properties
import SONATA.Pymore.utl.coef as coef
from SONATA.utl.plot import plot_fandiagram, plot_histogram_2Ddata, plot_beam_properties

from marc_fanplot import calc_fandiagram

def run_monte_carle_sample(job, sample, hide=True):
    if hide==True:
        with HiddenPrints():
            job.MaterialLst[8].E = sample
            job.cbm_run_vabs()
            
            x_offset = 0.81786984
            x0 = coef.refBeamProp()[0]
            job.config.setup['radial_station']=2000
            x1 = job.cbm_set_DymoreMK(x_offset) 
            job.config.setup['radial_station']=7500
            x2 = job.cbm_set_DymoreMK(x_offset)
            x3 = coef.refBeamProp()[-1]
            bp = np.vstack((x0,x1,x2,x3))
            (freq, Omega, RPM_vec, beamProp, Kmat) = calc_fandiagram(beamProp = bp)

    #RUN DYMORE CRUISE FLIGHT
       
    return (job.BeamProperties, freq, Omega, RPM_vec, beamProp)


if __name__ == '__main__':    
    
    filename = 'jobs/VariSpeed/uh60a_cbm_advanced/sec_config.yml'
    config = Configuration(filename)
    job = CBM(config)
    job.cbm_gen_topo()
    job.cbm_gen_mesh()
    #job.cbm_review_mesh()
    job.cbm_run_vabs()
    ref_cs = copy.copy(job.BeamProperties)
    #job.cbm_post_2dmesh(title='NoTitle')
    #job.cbm_post_3dtopo()
        
    #VARIATE THE Modulus of Material 1.
    N = 150 #number of Samples
    n = np.linspace(0,N-1,N)
    cv = 0.05 #coefficient of variation
    mu_E9 = np.array([165012,   8798,   8798])
    sigma = cv*mu_E9 # mean and standard deviation
    samples = np.random.normal(mu_E9, sigma, (N,3))
    
    t0 = time()
    #Multiprocessing 
    with Pool(processes=7) as pool:
        temp = [tup for tup in tqdm(pool.imap(functools.partial(run_monte_carle_sample, job), samples), total=len(samples))]
        
    #ThreadPoolExecutor
#    with futures.ThreadPoolExecutor(max_workers=7) as e:
#        BeamPropLst = [bp for bp in tqdm(e.map(functools.partial(foo, job), samples), total=len(samples))]
    
    #Single-Core
#    temp = [tup for tup in map(functools.partial(foo, job), samples)]
        
    BeamPropLst =[t[0] for t in temp]
    FreqArr = np.asarray([t[1] for t in temp])
    OmegaArr = np.asarray([t[2] for t in temp])
    RPM_vec = np.asarray([t[3] for t in temp])
    bp = np.asarray([t[4] for t in temp])
    t1 = time()
    t = t1-t0
    print('MP Time: %f' % t)
    

#%% ================= P O S T  -  P R O C E S S I N G =========================  
    plt.close('all')
    
    #%FAN PLOT POSTPROCESSING
    f_mean = np.mean(FreqArr.real, axis=0)
    f_std = np.std(FreqArr.real, axis=0)
    plot_fandiagram(f_mean,  OmegaArr[0], RPM_vec[0], sigma = f_std)
    
    #PLOT 2D Histogram
    tmp = np.asarray([b.CS for b in BeamPropLst])
    plot_histogram_2Ddata(tmp, ref_cs.CS, title='4x4 Stiffness Matrix')
    
    tmp = np.asarray([b.TS for b in BeamPropLst])
    plot_histogram_2Ddata(tmp, ref_cs.TS, title='6x6 Stiffness Matrix')

    tmp = np.asarray([b.MM for b in BeamPropLst])
    plot_histogram_2Ddata(tmp, ref_cs.MM, title='6x6 Mass Matrix')

    #PLOT Beam Properties
    data = np.mean(bp, axis=0) #Dymore 29, Beam Properties! to standartdize the procedure! 
    sigma = np.std(bp, axis=0)
    ref = coef.refBeamProp()
    
    plot_beam_properties(data, sigma, ref)