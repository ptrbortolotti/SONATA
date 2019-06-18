#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 14:56:31 2019

@author: gu32kij
"""
import numpy as np
from pyDOE2 import lhs
from scipy.stats.distributions import norm, norminvgauss, pareto
import matplotlib.pyplot as plt
from SONATA.utl.plot import plot_histogram_2Ddata

def doe_sampler(samples = 100, des_vars = ['x', 'y', 'z'],  means=[0, 0, 0], stdvs=[1.0, 1.0, 1.0], pdf = 'norm', design='lhs'):
    """
    
    
    Returns
    --------
    samples : list of name, value tuples for the design variables.
    """
    
    
    N = len(des_vars)
    design = lhs(N, samples)
    
    for i in range(N):
        design[:, i] = norm(loc=means[i], scale=stdvs[i]).ppf(design[:, i])
    
    #print(design)
    samples = []
    for s in design:
        tmp = [tuple(a) for a in zip(des_vars, list(s))]
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

if __name__ == '__main__':
    samples = doe_sampler(500, des_vars = ['x','y'], means=[1,0], stdvs= [0.5,1])
    plot_samples(samples, upper_tri = False)