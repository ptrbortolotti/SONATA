#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 10:23:01 2019

@author: gu32kij
"""
import numpy as np
from datetime import date
from collections import OrderedDict
import subprocess
import glob
from openmdao.api import SqliteRecorder, CaseReader
import matplotlib.pyplot as plt

from SONATA.utl.plot import plot_fandiagram, plot_beam_properties, plot_2dhist, plot_beam_properties_histogram, plot_histogram_2Ddata
from SONATA.cbm.cbm_utl import dymore2sixbysix
from SONATA.Pymore.app.marc_fanplot import frequency_analysis
from SONATA.openmdao_utl.doe_utl import filename_generator

#Prepare Samples


#%%============================================================================
#    E X E C U T E  -  D O E 
#==============================================================================

recalc = False
if recalc:
    #Prepare database filenames
    dbnamestring = filename_generator(directory = '/scratch/gu32kij', string='doe_cases', file_extension= '.sql')

    #RUN DOE
    s = 'mpirun -n 6 python run_doe.py'
    stdout = subprocess.run(s.split(), stdout=subprocess.PIPE).stdout.decode('utf-8')
    
else: 
    dbnamestring = '/scratch/gu32kij/cases_2019-07-20_3.sql*'
    
#%%============================================================================
#    P O S T  -  P R O C E S S I N G
#==============================================================================
plt.close('all')

files = glob.glob(dbnamestring)

cases = OrderedDict()
for f in files: 
    cr = CaseReader(f)
    cases.update(cr.get_cases('driver'))

mc_m1_E = np.asarray([c.outputs['m1_E'] for c in cases])
mc_m1_rho = np.asarray([c.outputs['m1_rho'] for c in cases])
mc_freq = np.asarray([c.outputs['freq'] for c in cases])
mc_BeamProps = np.asarray([c.outputs['BeamProps'] for c in cases])
mc_mean_elastic_tip_response = np.asarray([c.outputs['mean_elastic_tip_response'] for c in cases])
mc_stdvs_elastic_tip_response = np.asarray([c.outputs['stdvs_elastic_tip_response'] for c in cases])
mc_vibratory_hubloads = np.asarray([c.outputs['vibratory_hubloads'] for c in cases])
#mc_mean_bearing_response = np.asarray([c.outputs['mean_bearing_response'] for c in cases])
#mc_stdvs_bearing_response = np.asarray([c.outputs['stdvs_bearing_response'] for c in cases])


mean_BeamProps = np.mean(mc_BeamProps, axis=0)
stdvs_BeamProps = np.std(mc_BeamProps, axis=0)

# ============ Plot Beam-Properties ============================
#plot_beam_properties(mean_BeamProps, sigma=stdvs_BeamProps, description=False)
#
## ============ Plot 2D Histograms ============================
## - Plot 2D histogram of input design variables
plot_2dhist(mc_m1_E, np.array([r'$E_{11}$', r'$E_{22}$', r'$E_{33}$']))
plot_2dhist(mc_m1_rho, np.array([r'$\rho_{m_1}$']))
#
## - Plot 2D histogram of selected section 6x6 matrices
#plot_beam_properties_histogram(mc_BeamProps, station_idx=[8], bins  = 20)
## - Plot stdvs as scatterplot of selected section 6x6 matrices
#plot_beam_properties_histogram(mc_BeamProps, station_idx=[8], ptype='scatter')

#============  P Y M O R E  -   F R E Q .  A N A L Y S I S  ===================
Omega = 4.3*2*np.pi
RPM_vec = np.linspace(0.2*Omega, 1.2*Omega, 11)
mean_freq = np.mean(mc_freq, axis=0)
stdvs_freq = np.std(mc_freq, axis=0)
plot_fandiagram(mean_freq, Omega, RPM_vec, sigma=stdvs_freq)

#plot_2dhist(mc_freq[:,8,[2,6]], ['2nd flap frequency', '4th flap frequency', ], ptype = 'scatter')
#plot_2dhist(mc_freq[:,8,[2,6]], ['2nd flap frequency', '4th flap frequency', ], ptype = 'scatter')

#============  P Y M O R E  -  H O V E R   A N A L Y S I S  ===================
#plot_2dhist(mc_mean_elastic_tip_response, ['Elastic flap response', 'Elastic lag response', 'Elastic torsion response'])
#plot_2dhist(mc_mean_bearing_response, ['Bearing flap response', 'Bearing lag response', 'Bearing torsion response'])


#============  P Y M O R E  -  F L I G H T   A N A L Y S I S  ===================
test = np.reshape(mc_vibratory_hubloads, (28, 2,3))
plot_2dhist(test, np.array([['Fx/T', 'Fy/T', 'Fz/T'],[ 'Mx/Q', 'My/Q', 'Mz/Q']]))

#data['Psi'] = np.deg2rad(-data['SHAFT_SS_POS'][:,3]+180)
#psi = data['Psi'][-samples:]
  
xint = np.linspace(0,2*np.pi,181) #TODO!! GET COORECT DATA FROM ANALYIS!
fig, ax = plt.subplots(1,3)
title = 'Elastic Blade Tip-Response: '
fig.suptitle(title, fontsize=16)
fig.subplots_adjust(wspace=0.4, hspace=0.25)    

axh = ax[0]

for i, item in enumerate(mc_mean_elastic_tip_response):
    ymean = mc_mean_elastic_tip_response[i,0]
    ysigma = mc_stdvs_elastic_tip_response[i,0]
    axh.plot(xint, ymean, 'k')
    axh.fill_between(xint, ymean-ysigma, ymean+ysigma, color='k', alpha=0.4, linestyle='--', antialiased=True) 
    axh.set_ylabel(r'Elastic flap response [m]')
    axh.set_xlabel(r'Azimuth, $\Psi$')
    axh.set_xticks([0,np.pi/2,np.pi,3/2*np.pi,2*np.pi])
    axh.set_xticklabels(['0', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
    axh.set_xlim(0,2*np.pi)

axh = ax[1]
for i, item in enumerate(mc_mean_elastic_tip_response):
    ymean = mc_mean_elastic_tip_response[i,1]
    ysigma = mc_stdvs_elastic_tip_response[i,1]
    axh.plot(xint, ymean, 'k')
    axh.fill_between(xint, ymean-ysigma, ymean+ysigma, color='k', alpha=0.4, linestyle='--', antialiased=True) 
    axh.set_ylabel(r'Elastic lag response [m]')
    axh.set_xlabel(r'Azimuth, $\Psi$')
    axh.set_xticks([0,np.pi/2,np.pi,3/2*np.pi,2*np.pi])
    axh.set_xticklabels(['0', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
    axh.set_xlim(0,2*np.pi)

axh = ax[2]
for i, item in enumerate(mc_mean_elastic_tip_response):
    ymean = mc_mean_elastic_tip_response[i,2]
    ysigma = mc_stdvs_elastic_tip_response[i,2]
    axh.plot(xint, ymean, 'k')
    axh.fill_between(xint, ymean-ysigma, ymean+ysigma, color='k', alpha=0.4, linestyle='--', antialiased=True) 
    axh.set_ylabel(r'Elastic torsion response [deg]')
    axh.set_xlabel(r'Azimuth, $\Psi$')
    axh.set_xticks([0,np.pi/2,np.pi,3/2*np.pi,2*np.pi])
    axh.set_xticklabels(['0', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
    axh.set_xlim(0,2*np.pi)