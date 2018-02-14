# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 11:18:28 2018

@author: TPflumm
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import BlendedGenericTransform
from openmdao.api import Problem, ScipyOptimizer, IndepVarComp 
from openmdao.api import ScipyOptimizer, SimpleGADriver
print(os.getcwd())
os.chdir('C://TPflumm_local/work/SONATA')

from SONATA.cbm.fileIO.configuration import Configuration
from SONATA.cbm.fileIO.dymore_utils import read_dymore_beam_properties, interp1d_dymore_beam_properties
from SONATA.cbm.fileIO.readinput import read_material_input
from SONATA.cbm.sonata_cbm import CBM
from SONATA.vabs.VABS_interface import VABS_config, export_cells_for_VABS, XSectionalProperties
from SONATA.mdao.cbm_explcomp import CBM_ExplComp

#==============================================================================
#%%      Interpolate Values from dct_dym for optimization
#==============================================================================

#READ DYMORE BEAM PROPERTIES:
folder = 'jobs/VariSpeed/uh60a_blade/'
filename = folder + 'dymore_uh60a_rotor_blade.dat'
dct_dym = read_dymore_beam_properties(filename, x_offset = 0.81786984)
x = 2500 
dct_interp = interp1d_dymore_beam_properties(dct_dym,x)

#READ DAVIS UH-60A BLADE PROPERTIES:
dct_davis = {}
dct_davis['torsional_stiffness'] = np.loadtxt(folder + 'torsional_stiffness.dat')
dct_davis['torsional_inertia'] = np.loadtxt(folder + 'torsional_inertia.dat')
dct_davis['flapping_stiffness'] = np.loadtxt(folder + 'flapping_stiffness.dat')
dct_davis['edgewise_stiffness'] = np.loadtxt(folder + 'edgewise_stiffness.dat')
dct_davis['edgewise_inertia'] = np.loadtxt(folder + 'edgewise_inertia.dat')
dct_davis['mass'] = np.loadtxt(folder + 'mass.dat')
dct_davis['cg'] = np.loadtxt(folder + 'cg.dat')


#=============================================================================
#%%      SONATA - CBM
#==============================================================================
filename = 'jobs/VariSpeed/02_advanced/sec_config.input'
config = Configuration(filename)
config.SETUP_radial_station = x
flag_opt = True

if flag_opt:   
    
    p = Problem()
    #Generate independentVariableComponent
    ivc = p.model.add_subsystem('ivc', IndepVarComp())
    ivc.add_output('wp1', 0.32)
    ivc.add_output('wp2', 0.68)
    ivc.add_output('spar_lt', 2.2, units='mm')
    ivc.add_output('skin_lt', 1.2, units='mm')
    #ivc.add_output('rho_1', 0.05)
    
    #Generate Group of two Components
    p.model.add_subsystem('cbm_comp', CBM_ExplComp(config))
    p.model.cbm_comp.set_references(dct_interp) 
    p.model.connect('ivc.wp1', 'cbm_comp.WEB_Pos1')
    p.model.connect('ivc.wp2', 'cbm_comp.WEB_Pos2')
    p.model.connect('ivc.spar_lt', 'cbm_comp.Spar_layer_thickness')
    p.model.connect('ivc.skin_lt', 'cbm_comp.Skin_layer_thickness')
    #p.model.connect('ivc.rho_1', 'cbm_comp.Core1_density')
    
    p.model.add_design_var('ivc.wp1', lower=0.25, upper=0.45)
    p.model.add_design_var('ivc.wp2', lower=0.6, upper=0.8)
    p.model.add_design_var('ivc.spar_lt', lower=0.3, upper=2.2)
    p.model.add_design_var('ivc.skin_lt', lower=0.3, upper=1.6)
    #p.model.add_design_var('ivc.rho_1', lower=0.05, upper=19.25)
    
    lo = 0.8
    up = 1.2
    p.model.add_objective('cbm_comp.obj')
    #p.model.add_constraint('cbm_comp.Xm2', lower=dct_interp['centre_of_mass_location'][0]*lo, upper=dct_interp['centre_of_mass_location'][0]*lo)
    #p.model.add_constraint('cbm_comp.EI2', lower=dct_interp['bending_stiffnesses'][0]*lo, upper=dct_interp['bending_stiffnesses'][0]*lo)
    #p.model.add_constraint('cbm_comp.EI3', lower=dct_interp['bending_stiffnesses'][1]*lo, upper=dct_interp['bending_stiffnesses'][1]*lo)
    
    #Setup the Problem
#    p.driver = ScipyOptimizer()
#    p.driver.options['optimizer'] = 'COBYLA'
#    p.driver.options['disp'] = True
#    p.driver.options['tol'] = 1e-3
#    p.driver.options['maxiter'] = 1000
#    p.driver.opt_settings['rhobeg'] = 0.05 
    
    p.driver= SimpleGADriver()
    p.driver.options['bits'] = {'ivc.skin_lt' : 4}
    p.driver.options['bits'] = {'ivc.spar_lt' : 4}
    p.driver.options['bits'] = {'ivc.wp1' : 4}
    p.driver.options['bits'] = {'ivc.wp2' : 4}
    p.driver.options['pop_size'] = 5
    p.driver.options['max_gen'] = 5


    p.setup()
    p.run_driver()
    #print p['cbm_comp.MpUS'], (p['ivc.wp1'], p['ivc.wp2'], p['ivc.spar_lt'], p['ivc.skin_lt'])
    job_opt = p.model.cbm_comp.job
    job_opt.cbm_post_2dmesh(title = 'Optimization')
    

ref_config = Configuration(filename)
ref_config.SETUP_radial_station = x
job = CBM(ref_config)
job.cbm_gen_topo()
job.cbm_gen_mesh()
job.cbm_run_vabs()
job.cbm_post_2dmesh(title = 'Reference')


#==============================================================================
#%%      P L O T
#==============================================================================
plt.rc('text', usetex=False)
f, axarr = plt.subplots(5,2, sharex=True)    


axarr[0,0].plot(dct_davis['mass'][:,0],dct_davis['mass'][:,1],'r:', label='from S.J. Davis (1981, Sikorsky Aircraft Division)')
axarr[0,0].plot(dct_dym['x'],dct_dym['mass_per_unit_span'],'--', label='from DYMORE UH-60A (Yeo)')
axarr[0,0].plot(dct_interp['x'],dct_interp['mass_per_unit_span'],'gx', label='lin. interp. from DYMORE')
axarr[0,0].plot(job.config.SETUP_radial_station,job.BeamProperties.MpUS,'ko', label='SONATA CBM w. VABS')
if flag_opt:
    axarr[0,0].plot(job_opt.config.SETUP_radial_station,job_opt.BeamProperties.MpUS,'P', label='SONATA CBM OPT w. VABS', color='orange')
axarr[0,0].set_ylim([5,40])
axarr[0,0].set_ylabel(r'$m_{00}$ [kg/m]')
axarr[0,0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.5), ncol=2)

axarr[1,0].plot(dct_davis['cg'][:,0],dct_davis['cg'][:,1],'r:')
axarr[1,0].plot(dct_dym['x'],dct_dym['centre_of_mass_location'][:,0]*1000,'--')
axarr[1,0].plot(dct_interp['x'],dct_interp['centre_of_mass_location'][0]*1000,'gx')
axarr[1,0].plot(job.config.SETUP_radial_station,job.BeamProperties.Xm2,'ko')
if flag_opt:
    axarr[1,0].plot(job_opt.config.SETUP_radial_station,job_opt.BeamProperties.Xm2,'P', color='orange')
axarr[1,0].set_ylim([-100,100])
axarr[1,0].set_ylabel(r' $X_{m2}$ [mm]')


axarr[2,0].plot(dct_davis['torsional_inertia'][:,0],dct_davis['torsional_inertia'][:,1],'r:')
axarr[2,0].plot(dct_dym['x'],dct_dym['moments_of_inertia'][:,0],'--')
axarr[2,0].plot(dct_interp['x'],dct_interp['moments_of_inertia'][0],'gx')
axarr[2,0].plot(job.config.SETUP_radial_station,job.BeamProperties.MMatMC[3,3]*1e-6,'ko')
if flag_opt:
    axarr[2,0].plot(job_opt.config.SETUP_radial_station,job_opt.BeamProperties.MMatMC[3,3]*1e-6,'P', color='orange')
axarr[2,0].set_ylabel(r'$m_{11}$ [kg-m]')


axarr[3,0].plot(dct_dym['x'],dct_dym['moments_of_inertia'][:,1],'--')
axarr[3,0].plot(dct_interp['x'],dct_interp['moments_of_inertia'][1],'gx')
axarr[3,0].plot(job.config.SETUP_radial_station,job.BeamProperties.MMatMC[4,4]*1e-6,'ko')
if flag_opt:
    axarr[3,0].plot(job_opt.config.SETUP_radial_station,job_opt.BeamProperties.MMatMC[4,4]*1e-6,'P', color='orange')
axarr[3,0].set_ylabel(r'$m_{22}$ [kg-m]')


axarr[4,0].plot(dct_davis['edgewise_inertia'][:,0],dct_davis['edgewise_inertia'][:,1],'r:')
axarr[4,0].plot(dct_dym['x'],dct_dym['moments_of_inertia'][:,2],'--')
axarr[4,0].plot(dct_interp['x'],dct_interp['moments_of_inertia'][2],'gx')
axarr[4,0].plot(job.config.SETUP_radial_station,job.BeamProperties.MMatMC[5,5]*1e-6,'ko')
if flag_opt:
    axarr[4,0].plot(job_opt.config.SETUP_radial_station,job_opt.BeamProperties.MMatMC[5,5]*1e-6,'P', color='orange')
axarr[4,0].set_ylabel(r'$m_{33}$ [kg-m]')
axarr[4,0].set_xlabel('Radius [mm]')

#------------------------------------

axarr[0,1].plot(dct_dym['x'],dct_dym['shear_centre_location'][:,0]*1000,'--')
axarr[0,1].plot(dct_interp['x'],dct_interp['shear_centre_location'][0]*1000,'gx')
axarr[0,1].plot(job.config.SETUP_radial_station,job.BeamProperties.Xs2,'ko')
if flag_opt:
    axarr[0,1].plot(job_opt.config.SETUP_radial_station,job_opt.BeamProperties.Xs2,'P', color='orange')
axarr[0,1].set_ylim([-100,100])
axarr[0,1].set_ylabel(r'$X_{s2}$ [mm]')


axarr[1,1].plot(dct_dym['x'],dct_dym['axial_stiffness'],'--')
axarr[1,1].plot(dct_interp['x'],dct_interp['axial_stiffness'],'gx')
axarr[1,1].plot(job.config.SETUP_radial_station,job.BeamProperties.CS[0,0],'ko')
if flag_opt:
    axarr[1,1].plot(job_opt.config.SETUP_radial_station,job_opt.BeamProperties.CS[0,0],'P', color='orange')
axarr[1,1].set_ylabel(r'$EA \; [N]$')
axarr[1,1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))


axarr[2,1].plot(dct_davis['torsional_stiffness'][:,0],dct_davis['torsional_stiffness'][:,1],'r:')
axarr[2,1].plot(dct_dym['x'],dct_dym['torsional_stiffness'],'--')
axarr[2,1].plot(dct_interp['x'],dct_interp['torsional_stiffness'],'gx')
axarr[2,1].plot(job.config.SETUP_radial_station, job.BeamProperties.CS[1,1]*1e-6,'ko')
if flag_opt:
    axarr[2,1].plot(job_opt.config.SETUP_radial_station, job_opt.BeamProperties.CS[1,1]*1e-6,'P', color='orange')
axarr[2,1].set_ylabel(r'$GJ \; [Nm^2]$')
axarr[2,1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))


axarr[3,1].plot(dct_davis['flapping_stiffness'][:,0],dct_davis['flapping_stiffness'][:,1],'r:')
axarr[3,1].plot(dct_dym['x'],dct_dym['bending_stiffnesses'][:,0],'--') 
axarr[3,1].plot(dct_interp['x'],dct_interp['bending_stiffnesses'][0],'gx') 
axarr[3,1].plot(job.config.SETUP_radial_station, job.BeamProperties.CS[2,2]*1e-6,'ko') 
if flag_opt:
    axarr[3,1].plot(job_opt.config.SETUP_radial_station,job_opt.BeamProperties.CS[2,2]*1e-6,'P', color='orange') 
axarr[3,1].set_ylabel(r'$EI_{2} \; [Nm^2]$')
axarr[3,1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))


axarr[4,1].plot(dct_davis['edgewise_stiffness'][:,0],dct_davis['edgewise_stiffness'][:,1],'r:')
axarr[4,1].plot(dct_dym['x'],dct_dym['bending_stiffnesses'][:,1],'--')
axarr[4,1].plot(dct_interp['x'],dct_interp['bending_stiffnesses'][1],'gx')
axarr[4,1].plot(job.config.SETUP_radial_station,job.BeamProperties.CS[3,3]*1e-6,'ko') 
if flag_opt:
    axarr[4,1].plot(job_opt.config.SETUP_radial_station,job_opt.BeamProperties.CS[3,3]*1e-6,'P', color='orange') 
axarr[4,1].set_ylabel(r'$EI_3 \; [Nm^2]$')
axarr[4,1].set_xlabel('Radius [mm]')
axarr[4,1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))