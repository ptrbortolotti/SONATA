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
from SONATA.cbm.classCBM import CBM#
from SONATA.cbm.classCBMConfig import CBMConfig
from SONATA.cbm.fileIO.hiddenprints import HiddenPrints

from SONATA.cbm.fileIO.dymore_utils import read_dymore_beam_properties, interp1d_dymore_beam_properties
import SONATA.Pymore.utl.coef as coef
from SONATA.utl.plot import plot_fandiagram, plot_histogram_2Ddata, plot_beam_properties

from SONATA.utl.marc_fanplot import calc_fandiagram

def run_MC_sample(fname, sample, hide=True):
    if hide==True:

        with HiddenPrints():
            #create job instance!
            config = CBMConfig(fname)
            job = CBM(config)
            job.cbm_gen_topo()
            job.cbm_gen_mesh()
            
            job.materials[9].E = sample

            job.cbm_run_vabs(ramdisk=True)
            
            x_offset = 0.81786984
            x0 = coef.refBeamProp()[0]
            job.config.setup['radial_station']=2000
            x1 = job.cbm_set_DymoreMK(x_offset)
            job.config.setup['radial_station']=7500
            x2 = job.cbm_set_DymoreMK(x_offset)
            x3 = coef.refBeamProp()[-1]
            bp = np.vstack((x0,x1,x2,x3))
            (freq, Omega, RPM_vec) = calc_fandiagram(beamProp = bp)
            # (bladetip_displacement) = calc_hover(bp1)
            # (calc_forwardflight)
            

        print('Sample finished:', sample)
        print(x1[6:9])      

    #RUN DYMORE CRUISE FLIGHT
    return (job.BeamProperties, freq, Omega, RPM_vec, bp)


#def CBMsampler(fname, nsamples=1000, style='latin_hypercube'):
#    from pyDOE import lhs
#    from scipy.stats.distributions import norm
#    
#    '''[('MaterialLst[0].E', distr = 'normal', mu='ref', cv=0.05), 
#        ('config.webs[1].Pos1', distr = 'normal', mu='ref', cv=0.05),
#        ('config.webs[2].Pos2', distr = 'exp', mu=0.23, cv=0.05)]'''
#    
#    
#    design = lhs(4, samples=10)
#    from scipy.stats.distributions import norm
#    means = [1, 2, 3, 4]
#    stdvs = [0.1, 0.5, 1, 0.25]
#    for i in xrange(4):
#        design[:, i] = norm(loc=means[i], scale=stdvs[i]).ppf(design[:, i])
#    
#    
#    
#    #identify value for P and run function to sample it.
#    for pTup in pLst:        
#           
##    cbmLst = []
##    if 'MaterialLst' in pstr:
##        getattr(cbm, 'MaterialLst')
##        #Select 
##    
##    if 'config' in pstr:
##        getattr(config, 'MaterialLst')
##      
#    #Sampling
#    return cbmLst 

if __name__ == '__main__':    
    

    #Define a Sampling Function that Generates a Sampling list of CBM instances. That can be passed to the Monte Carlo Simulation!
    
    
    fname = 'jobs/VariSpeed/uh60a_cbm_advanced/sec_config.yml'
    #fname = 'jobs/VariSpeed/uh60a_cbm_simple/sec_config.yml'
    config = CBMConfig(fname)
    job = CBM(config)
    job.cbm_gen_topo()
    job.cbm_gen_mesh()
    #job.cbm_review_mesh()
    job.cbm_run_vabs()
    ref_cs = copy.copy(job.BeamProperties)
    #job.cbm_post_2dmesh(title='NoTitle')
    #job.cbm_post_3dtopo()
    
    
    #VARIATE THE Modulus of Material 1.
    N = 250   #number of Samples
    n = np.linspace(0,N-1,N)
    cv = 0.05 #coefficient of variation
    mu_E9 = job.materials[9].E
    #mu_E1 = job.MaterialLst[1].E
    #mu_rho = job.MaterialLst[11].rho
    
    ID = np.array([np.linspace(0,N-1,N)]).T
    samples = np.random.normal(mu_E9, cv*mu_E9, (N,3))
    #sample2 = np.random.normal(mu_E0, cv*mu_E0, (N,3))
    #sample3 = np.random.normal(mu_rho, cv*mu_rho, (N,1))
    #samples = np.hstack((ID,sample1))                                        

    t0 = time()
#    #Multiprocessing 
    with Pool(processes=7) as pool:
        temp = [tup for tup in pool.imap(functools.partial(run_MC_sample, fname), samples)]
        
    #ThreadPoolExecutor
#    with futures.ThreadPoolExecutor(max_workers=8) as e:
#        temp = [tup for tup in tqdm(e.map(functools.partial(run_monte_carle_sample, config), samples), total=len(samples))]
    
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
    from matplotlib2tikz import save as tikz_save
    plt.close('all')
    
    job.cbm_post_2dmesh(title='UH-60A Simple Config')
    #%FAN PLOT POSTPROCESSING
    f_mean = np.mean(FreqArr.real, axis=0)
    f_std = np.std(FreqArr.real, axis=0)
    plot_fandiagram(f_mean,  OmegaArr[0], RPM_vec[0], sigma = f_std, ), #ref_str = ['l1','f1','f2','f3','l2','t1','f4']
    tikz_save('/media/gu32kij/HTMWTUM/Oeffentlich/Publikationen/ERF/2019/TPflumm, WGarre/abstract/img/fanplot.tikz', figureheight='\\figureheight', figurewidth='\\figurewidth' )
    #PLOT 2D Histogram
    tmp = np.asarray([b.CS for b in BeamPropLst])
    plot_histogram_2Ddata(tmp, ref = ref_cs.CS, title='4x4 Stiffness Matrix')
#    tikz_save('/media/gu32kij/HTMWTUM/Oeffentlich/Publikationen/ERF/2019/TPflumm, WGarre/abstract/img/CS.tikz', figureheight='\\figureheight', figurewidth='\\figurewidth' )
#    tmp = np.asarray([b.TS for b in BeamPropLst])
#    plot_histogram_2Ddata(tmp, ref = ref_cs.TS, title='6x6 Stiffness Matrix')
#
#    tmp = np.asarray([b.MM for b in BeamPropLst])
#    plot_histogram_2Ddata(tmp, ref = ref_cs.MM, title='6x6 Mass Matrix')

    #PLOT Beam Properties
    data = np.mean(bp, axis=0) #Dymore 29, Beam Properties! to standartdize the procedure! 
    sigma = np.std(bp, axis=0)
    ref = coef.refBeamProp()
    
    plot_beam_properties(data, sigma, ref)