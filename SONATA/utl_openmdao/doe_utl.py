#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 14:56:31 2019

@author: gu32kij
"""
import os
import glob
import numpy as np
from pyDOE2 import lhs
from datetime import date

from scipy.stats.distributions import norm, norminvgauss, pareto
import matplotlib.pyplot as plt

#from openmdao.api import CaseReader
#from SONATA.utl.plot import plot_histogram_2Ddata

def doe_sampler(ivc_dct, nsamples = 100, pdf = 'norm', design='lhs'):
    """
    generates samples for the openmdao ListGenerator(DOEGenerator) 
    
    Returns
    --------
    samples : list of (name, value) tuples for the design variables.
    """
    
    
    N = 0
    for key,values in ivc_dct.items():
        N += len(values.get('means'))
    
    design = lhs(N, nsamples)
    n = 0
    for i,(k,v) in enumerate(ivc_dct.items()):
        n1 = n+len(v.get('means'))
        ivc_dct[k]['samples'] = norm(loc=v['means'], scale=v['stdvs']).ppf(design[:, n:n1])
        n = n1
        
    samples = []
    for i in range(nsamples):
        tmp = []
        for k,v in ivc_dct.items():
            tmp.append((k,v['samples'][i]))
        samples.append(tmp)
        
    return samples


def plot_samples(samples, title='design variables', **kwargs):
    """
    plots the openmdao samples for the DOEListGenerator instance
    
    Parameters
    ---------
    samples : list of name, value tuples for the design variables.

    """
    data = np.asarray(samples)[:,:,1:].astype(float)
    a = np.asarray(samples)[0,:,0]
    b = np.expand_dims(a, 0)
    labels = np.expand_dims(b, 2)
    imax = data.shape[1]
    jmax = data.shape[2]
    fig, ax = plt.subplots(imax,jmax)
    fig.suptitle(title, fontsize=16)
    fig.subplots_adjust(wspace=0.25, hspace=0.25)         

    for i in range(imax):
        for j in range(jmax):
            if imax==1 and jmax==1:
                axh = ax
            elif imax>1 and jmax==1:
                axh = ax[i]
            else:
                axh = ax[i][j]
            
            arr = data[:,i,j]
            count, bins, ignored =  axh.hist(arr, 30, density=True, alpha=0.5, label='histogram')
            sigma = arr.std()
            mu = arr.mean()
            if sigma > 0:
                axh.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) *
                               np.exp( - (bins - mu)**2 / (2 * sigma**2) ),
                         linewidth=2, color='r', linestyle='--', alpha=0.5, label='ML approx.')       
                string = r'$\sigma$ = %.2f $\%%$' % sigma
                axh.text(mu,0,string)
                axh.set_xlim(-20,20)
            axh.set_ylabel(labels[0,i,j])
    plt.legend()
    plt.show()


def filename_generator(directory = '/scratch/gu32kij', string = 'doe_cases', file_extension='.sql'):
    """
    generates a filename based on isoformat date and adds an consecutive 
    integer at the end if if does exist.
    """
    day = date.today().isoformat()
    i=0
    if directory != '':
        tmp = directory + '/' + string + '_' + str(day) + '_' + str(i) + file_extension + '*'
    else:
        tmp =string + '_' + str(day) + '_' + str(i) + file_extension + '*'
        
    while len(glob.glob(tmp)) != 0:
        i += 1
        if directory != '':
            tmp = directory + '/' + string + '_' + str(day) + '_' + str(i) + file_extension + '*'
        else:
            tmp =string + '_' + str(day) + '_' + str(i) + file_extension + '*'
        
        
    if directory != '':
        fname = directory + '/' + string + '_' + str(day) + '_' + str(i) + file_extension
    else:
        fname =string + '_' + str(day) + '_' + str(i) + file_extension
    return fname

def directory_creator(parent_directory = '/scratch/gu32kij', string = 'doe_cases', MPI=False):
    """
    generates a filename based on isoformat date and adds an consecutive 
    integer at the end if if does exist.
    """
    day = date.today().isoformat()
    i=0
    if parent_directory != '':
        tmp = parent_directory + '/' + string + '_' + str(day) + '_' + str(i)
    else:
        tmp =string + '_' + str(day) + '_' + str(i)
        
        
    while os.path.exists(tmp):
        i += 1
        if parent_directory != '':
            tmp = parent_directory + '/' + string + '_' + str(day) + '_' + str(i)
        else:
            tmp =string + '_' + str(day) + '_' + str(i)
        
        
    if parent_directory != '':
        dirname = parent_directory + '/' + string + '_' + str(day) + '_' + str(i)
    else:
        dirname =string + '_' + str(day) + '_' + str(i)
    
    os.mkdir(dirname)
    return dirname


if __name__ == '__main__':
    
    desvar_dict = {'m1_E': {'means':np.array([0,0,0]), 'stdvs':np.array([1,1,1])},
        'm1_rho': {'means':np.array([0]), 'stdvs':np.array([1])}}
    
    
    samples = doe_sampler(desvar_dict, nsamples=100)
    plot_samples(samples, upper_tri = False)