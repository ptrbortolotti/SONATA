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
from openmdao.drivers.genetic_algorithm_driver import SimpleGADriver

#Sonata Modules: Make Sure to be in the SONATA working directory!
os.chdir('../..')  #print(os.getcwd())

from SONATA.cbm.fileIO.configuration import Configuration
from SONATA.cbm.fileIO.dymore_utils import read_dymore_beam_properties, interp1d_dymore_beam_properties
from SONATA.cbm.fileIO.readinput import read_material_input
from SONATA.cbm.sonata_cbm import CBM
from SONATA.vabs.VABS_interface import VABS_config, export_cells_for_VABS, XSectionalProperties

#job specific modules!
from jobs.VariSpeed.mdao_group import MDAO_Group

__spec__ = None

#==============================================================================
#%%      Interpolate Values from dct_dym for optimization
#==============================================================================
#READ DYMORE BEAM PROPERTIES:
folder = 'SONATA/Pymore/dym/mdl/03_rotormodel/05_UH60_rotor_optimization/01_UH60_rotor_snglblade_static/'
filename = folder + 'rotor_blade.dat'
dct_dym = read_dymore_beam_properties(filename, x_offset = 0.81786984)

#READ DAVIS UH-60A BLADE PROPERTIES:
folder = 'jobs/VariSpeed/uh60a_data_blade/uh60a_blade/'
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
filename = 'jobs/VariSpeed/uh60a_cbm_advanced/sec_config.yml'

config = Configuration(filename)
config.setup['radial_station'] = 2500
config.setup['BalanceWeight'] = False
dct_interp = interp1d_dymore_beam_properties(dct_dym,config.setup['radial_station'])

flag_ref = False
if flag_ref:
    job = CBM(config)
    job.cbm_gen_topo()
    job.cbm_gen_mesh()
    job.cbm_run_vabs()
    job.cbm_post_2dmesh(title = 'Reference')
    job.cbm_set_DymoreMK(x_offset = 0.81786984)

#=============================================================================
#%%      SONATA - Pymore
#==============================================================================
flag_opt = True
if flag_opt:   
    p = Problem()
    p.model = MDAO_Group(config, ref_dct = dct_interp)

    p.model.add_design_var('s_w1', lower=0.35, upper=0.44, ref=0.45, ref0 = 0.44)
    p.model.add_design_var('s_w2', lower=0.2,  upper=0.31, ref=0.31, ref0 = 0.30)
    p.model.add_design_var('t_sparcap1', lower=0.4, upper=2.7, ref=2.7, ref0 = 0.4)
    p.model.add_design_var('t_sparcap2', lower=0.4, upper=2.7, ref=2.7, ref0 = 0.4)
    p.model.add_design_var('t_sparcap3', lower=0.4, upper=2.7, ref=2.7, ref0 = 0.4)
    p.model.add_design_var('t_sparcap4', lower=0.4, upper=2.7, ref=2.7, ref0 = 0.4)
    p.model.add_design_var('rho_mat11', lower=0.05, upper=19.25,   ref=19.25, ref0 = 0)
    
    p.model.add_objective('cbm_comp.obj')
    #p.model.add_objective('marc_comp.obj')
    
    p.driver= SimpleGADriver()
    p.set_solver_print(level=0)
    p.driver.options['debug_print'] = ['desvars','objs','totals']
    p.driver.options['bits'] = {'s_w1' : 8}
    p.driver.options['bits'] = {'s_w2' : 8}
    p.driver.options['bits'] = {'t_sparcap1' : 8}
    p.driver.options['bits'] = {'t_sparcap2' : 8}
    p.driver.options['bits'] = {'t_sparcap3' : 8}
    p.driver.options['bits'] = {'t_sparcap4' : 8}
    p.driver.options['bits'] = {'rho_mat11' : 8}

    p.driver.options['pop_size'] = 80
    p.driver.options['max_gen'] = 10
    p.driver.options['run_parallel'] = False

    p.setup()
    p.run_driver()
    #p.run_model()
    #print p['cbm_comp.MpUS'], (p['ivc.wp1'], p['ivc.wp2'], p['ivc.spar_lt'], p['ivc.skin_lt'])
    job_opt = p.model.cbm_comp.job
    job_opt.cbm_post_2dmesh(title = 'Optimization')
    val_fname = 'jobs/VariSpeed/uh60a_data_blade/Fanplot_Bowen_Davies_Diss.csv'
    p.model.marc_comp.job.fanplot_show(p.model.marc_comp.RPM_vec, p.model.marc_comp.result_dir,val_fname=val_fname)

#==============================================================================
#%%      P L O T
#==============================================================================
plt.rc('text', usetex=False)
f, axarr = plt.subplots(4,2, sharex=True)    

axarr[0,0].plot(dct_davis['mass'][:,0],dct_davis['mass'][:,1],'r:', label='from S.J. Davis (1981, Sikorsky Aircraft Division)')
axarr[0,0].plot(dct_dym['x'],dct_dym['mass_per_unit_span'],'--', label='from DYMORE UH-60A (Yeo)')
#axarr[0,0].plot(dct_interp['x'],dct_interp['mass_per_unit_span'],'gx', label='lin. interp. from DYMORE')
axarr[0,0].plot(job.config.setup['radial_station'],job.BeamProperties.MpUS,'ko', label='SONATA CBM (VABS)')
if flag_opt:
    axarr[0,0].plot(job_opt.config.setup['radial_station'],job_opt.BeamProperties.MpUS,'P', label='SONATA CBM OPT w. VABS', color='orange')
axarr[0,0].set_ylim([5,40])
axarr[0,0].set_ylabel(r'$m_{00}$ [kg/m]')
axarr[0,0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.5), ncol=2)

axarr[1,0].plot(dct_davis['cg'][:,0],dct_davis['cg'][:,1],'r:')
axarr[1,0].plot(dct_dym['x'],dct_dym['centre_of_mass_location'][:,0]*1000,'--')
#axarr[1,0].plot(dct_interp['x'],dct_interp['centre_of_mass_location'][0]*1000,'gx')
axarr[1,0].plot(job.config.setup['radial_station'],job.BeamProperties.Xm2,'ko')
if flag_opt:
    axarr[1,0].plot(job_opt.config.setup['radial_station'],job_opt.BeamProperties.Xm2,'P', color='orange')
axarr[1,0].set_ylim([-100,100])
axarr[1,0].set_ylabel(r' $X_{m2}$ [mm]')


axarr[2,0].plot(dct_davis['torsional_inertia'][:,0],dct_davis['torsional_inertia'][:,1],'r:')
axarr[2,0].plot(dct_dym['x'],dct_dym['moments_of_inertia'][:,0],'--')
#axarr[2,0].plot(dct_interp['x'],dct_interp['moments_of_inertia'x_offset = 0.81786984][0],'gx')
axarr[2,0].plot(job.config.setup['radial_station'],job.BeamProperties.MMatMC[3,3]*1e-6,'ko')
if flag_opt:
    axarr[2,0].plot(job_opt.config.setup['radial_station'],job_opt.BeamProperties.MMatMC[3,3]*1e-6,'P', color='orange')
axarr[2,0].set_ylabel(r'$m_{11}$ [kg-m]')


#axarr[3,0].plot(dct_dym['x'],dct_dym['moments_of_inertia'][:,1],'--')
#axarr[3,0].plot(dct_interp['x'],dct_interp['moments_of_inertia'][1],'gx')
#axarr[3,0].plot(job.config.SETUP_radial_station,job.BeamProperties.MMatMC[4,4]*1e-6,'ko')
#if flag_opt:
#    axarr[3,0].plot(job_opt.config.SETUP_radial_station,job_opt.BeamProperties.MMatMC[4,4]*1e-6,'P', color='orange')
#axarr[3,0].set_ylabel(r'$m_{22}$ [kg-m]')


axarr[3,0].plot(dct_davis['edgewise_inertia'][:,0],dct_davis['edgewise_inertia'][:,1],'r:')
axarr[3,0].plot(dct_dym['x'],dct_dym['moments_of_inertia'][:,2],'--')
#axarr[3,0].plot(dct_interp['x'],dct_interp['moments_of_inertia'][2],'gx')
axarr[3,0].plot(job.config.setup['radial_station'],job.BeamProperties.MMatMC[5,5]*1e-6,'ko')
if flag_opt:
    axarr[3,0].plot(job_opt.config.setup['radial_station'],job_opt.BeamProperties.MMatMC[5,5]*1e-6,'P', color='orange')
axarr[3,0].set_ylabel(r'$m_{33}$ [kg-m]')
axarr[3,0].set_xlabel('Radius [mm]')

#------------------------------------

#axarr[0,1].plot(dct_dym['x'],dct_dym['shear_centre_location'][:,0]*1000,'--')
#axarr[0,1].plot(dct_interp['x'],dct_interp['shear_centre_location'][0]*1000,'gx')
#axarr[0,1].plot(job.config.SETUP_radial_station,job.BeamProperx_offset = 0.81786984ties.Xs2,'ko')
#if flag_opt:
#    axarr[0,1].plot(job_opt.config.SETUP_radial_station,job_opt.BeamProperties.Xs2,'P', color='orange')
#axarr[0,1].set_ylim([-100,100])
#axarr[0,1].set_ylabel(r'$X_{s2}$ [mm]')


axarr[0,1].plot(dct_dym['x'],dct_dym['axial_stiffness'],'--')
#axarr[0,1].plot(dct_interp['x'],dct_interp['axial_stiffness'],'gx')
axarr[0,1].plot(job.config.setup['radial_station'],job.BeamProperties.CS[0,0],'ko')
if flag_opt:
    axarr[0,1].plot(job_opt.config.setup['radial_station'],job_opt.BeamProperties.CS[0,0],'P', color='orange')
axarr[0,1].set_ylabel(r'$EA \; [N]$')
axarr[0,1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))


axarr[1,1].plot(dct_davis['torsional_stiffness'][:,0],dct_davis['torsional_stiffness'][:,1],'r:')
axarr[1,1].plot(dct_dym['x'],dct_dym['torsional_stiffness'],'--')
#axarr[1,1].plot(dct_interp['x'],dct_interp['torsional_stiffness'],'gx')
axarr[1,1].plot(job.config.setup['radial_station'], job.BeamProperties.CS[1,1]*1e-6,'ko')
if flag_opt:
    axarr[1,1].plot(job_opt.config.setup['radial_station'], job_opt.BeamProperties.CS[1,1]*1e-6,'P', color='orange')
axarr[1,1].set_ylabel(r'$GJ \; [Nm^2]$')
axarr[1,1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
x_offset = 0.81786984

axarr[2,1].plot(dct_davis['flapping_stiffness'][:,0],dct_davis['flapping_stiffness'][:,1],'r:')
axarr[2,1].plot(dct_dym['x'],dct_dym['bending_stiffnesses'],'--') 
#axarr[2,1].plot(dct_interp['x'],dct_interp['bending_stiffnesses'][0],'gx') 
axarr[2,1].plot(job.config.setup['radial_station'], job.BeamProperties.CS[2,2]*1e-6,'ko') 
if flag_opt:
    axarr[2,1].plot(job_opt.config.setup['radial_station'],job_opt.BeamProperties.CS[2,2]*1e-6,'P', color='orange') 
axarr[2,1].set_ylabel(r'$EI_{2} \; [Nm^2]$')
axarr[2,1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))


axarr[3,1].plot(dct_davis['edgewise_stiffness'][:,0],dct_davis['edgewise_stiffness'][:,1],'r:')
axarr[3,1].plot(dct_dym['x'],dct_dym['bending_stiffnesses'][:,1],'--')
#axarr[3,1].plot(dct_interp['x'],dct_interp['bending_stiffnesses'][1],'gx')
axarr[3,1].plot(job.config.setup['radial_station'],job.BeamProperties.CS[3,3]*1e-6,'ko') 
if flag_opt:
    axarr[3,1].plot(job_opt.config.setup['radial_station'],job_opt.BeamProperties.CS[3,3]*1e-6,'P', color='orange') 
axarr[3,1].set_ylabel(r'$EI_3 \; [Nm^2]$')
axarr[3,1].set_xlabel('Radius [mm]')
axarr[3,1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))

#from matplotlib2tikz import save as tikz_save
#tikz_save('UH60A_beam.tikz', figureheight='\\figureheight', figurewidth='\\figurewidth' )

