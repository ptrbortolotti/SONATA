# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 11:18:28 2018

@author: TPflumm
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import BlendedGenericTransform
from openmdao.api import Problem, ScipyOptimizer, IndepVarComp, ScipyOptimizeDriver 
from openmdao.drivers.genetic_algorithm_driver import SimpleGADriver
from concurrent import futures

#Sonata Modules: Make Sure to be in the SONATA working directory!
os.chdir('../..')  #print(os.getcwd())

from SONATA.cbm.fileIO.configuration import Configuration
from SONATA.cbm.fileIO.dymore_utils import read_dymore_beam_properties, interp1d_dymore_beam_properties
from SONATA.cbm.fileIO.readinput import read_material_input
from SONATA.cbm.sonata_cbm import CBM
from SONATA.vabs.VABS_interface import VABS_config, export_cells_for_VABS, XSectionalProperties

#job specific modules!
from jobs.VariSpeed.sonata_group import Sonata_Group
plt.close('all')
__spec__ = None

#==============================================================================
#%%      Interpolate Values from dct_dym for optimization
#==============================================================================
#READ DYMORE BEAM PROPERTIES:
folder = 'SONATA/Pymore/dym/mdl/03_rotormodel/05_UH60_rotor_optimization/01_UH60_rotor_snglblade_static/'
filename = folder + 'rotor_blade.dat'
x_offset = 0.81786984
dct_dym = read_dymore_beam_properties(filename, x_offset = x_offset)

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

def run_cbm_optimization(radial_station, dct_dym, flag_ref=True, flag_opt=False, solver='slsqp'):
    filename = 'jobs/VariSpeed/uh60a_cbm_advanced/sec_config_R%s.yml' % (radial_station)
    config = Configuration(filename)
    config.setup['radial_station'] = radial_station
    config.setup['BalanceWeight'] = True
    dct_interp = interp1d_dymore_beam_properties(dct_dym,config.setup['radial_station'])
    
    #=============================================================================
    #%%      SONATA - Pymore
    #==============================================================================
    job=None
    #flag_ref = True
    if flag_ref:
        job = CBM(config)
        job.cbm_gen_topo()
        job.cbm_gen_mesh()
        job.cbm_run_vabs(rm_vabfiles=True)
    
    job_opt=None
    #flag_opt = False
    #solver = 'slsqp'
    if flag_opt:   
        p = Problem()
        p.model = Sonata_Group(config, ref_dct = dct_interp)
    
        p.model.add_design_var('rho_mat3', lower=0.05, upper=19.25, ref0 = 0, ref=19.25)
        #p.model.add_design_var('rho_mat11', lower=0.05, upper=19.25, ref0 = 0, ref=19.25)
        p.model.add_design_var('t_sparcap1', lower=0.35, upper=2.7, ref0 = 0.4, ref=2.7)
        p.model.add_design_var('t_sparcap2', lower=0.35, upper=2.7, ref0 = 0.4, ref=2.7)
        p.model.add_design_var('t_sparcap3', lower=0.35, upper=2.7, ref0 = 0.4, ref=2.7)
        p.model.add_design_var('t_sparcap4', lower=0.35, upper=2.7, ref0 = 0.4, ref=2.7)
        p.model.add_design_var('s_w1', lower=0.35, upper=0.425, ref0 = 0.35, ref=0.42)
        p.model.add_design_var('s_w2', lower=0.2,  upper=0.32, ref0 = 0.2, ref=0.31)
        p.model.add_design_var('s_spar2', lower=0.445,  upper=0.476, ref0 = 0.445, ref=0.476)
        
        p.model.add_objective('cbm_comp.obj')
    
        #p.model.add_objective('marc_comp.obj')
        
        
        if solver == 'simpleGA':
            p.driver= SimpleGADriver()
            p.set_solver_print(level=2)
            p.driver.options['debug_print'] = ['desvars','objs','totals']
            p.driver.options['bits'] = {'s_w1' : 8}
            p.driver.options['bits'] = {'s_w2' : 8}
            p.driver.options['bits'] = {'t_sparcap1' : 8}
    #        p.driver.options['bits'] = {'t_sparcap2' : 8}
    #        p.driver.options['bits'] = {'t_sparcap3' : 8}
    #        p.driver.options['bits'] = {'t_sparcap4' : 8}
            p.driver.options['bits'] = {'rho_mat3' : 8}
        
            p.driver.options['pop_size'] = 10
            p.driver.options['max_gen'] = 10
            p.driver.options['run_parallel'] = False
            
        elif solver == 'slsqp':
            p.driver = ScipyOptimizeDriver()
            p.driver.options['optimizer'] = 'SLSQP'
            p.driver.options['debug_print'] = ['desvars','objs']
            p.driver.options['maxiter'] = 100
            p.driver.options['tol'] = 2e-5
            p.driver.opt_settings['eps'] = 1.0
            
        elif solver == 'newton-cg':
            p.driver = ScipyOptimizeDriver()
            p.driver.options['optimizer'] = 'Newton-CG'
            p.driver.options['debug_print'] = ['desvars','objs']
            p.driver.options['maxiter'] = 100
            p.driver.options['tol'] = 1e-5
            p.driver.opt_settings['eps'] = 0.1
            
        elif solver == 'cobyla':
            p.driver = ScipyOptimizeDriver()
            p.driver.options['optimizer'] = 'COBYLA'
            p.driver.options['debug_print'] = ['desvars','objs']
            p.driver.options['maxiter'] = 300
            p.driver.options['tol'] = 1e-5
            p.driver.opt_settings['rhobeg'] = 0.8
        
        p.setup()
        #p.run_model()
        p.run_driver()
    
    if flag_opt:    
        job_opt = p.model.cbm_comp.job
    return (job, job_opt)

radial_station = [2000,7500]
with futures.ProcessPoolExecutor(max_workers=6) as e:
    fs = {n:e.submit(run_cbm_optimization, n, dct_dym) for n in radial_station}
    print('Alle Aufgaben gestartet')
print('Alle Aufgaben erledigt.')

#==================================================================================
#%%      P Y M O R E
#==================================================================================
from SONATA.Pymore.marc.marc import MARC
a = fs[2000].result()[0].cbm_set_DymoreMK(x_offset)
a[-1] = +0.000e+00
b = fs[2000].result()[0].cbm_set_DymoreMK(x_offset)
c = fs[7500].result()[0].cbm_set_DymoreMK(x_offset)
d = fs[7500].result()[0].cbm_set_DymoreMK(x_offset)
d[-1] = +7.361e+00
beamProp = np.vstack((a,b,c,d))

nbOfLoc = 11
Omega = 4.3*2*np.pi #in rad/sec
RPM_vec = np.linspace(0.2*Omega, 1.15*Omega, nbOfLoc)

dir_root = 'SONATA/Pymore/dym/mdl/03_rotormodel/05_UH60_rotor_optimization/01_UH60_rotor_snglblade_static/'
result_dir = 'SONATA/Pymore/rlt/'

job_pym = MARC(dir_root, 'rotor_assembly.dym')
job_pym.marc_set_beamProp('BLADE_BP_CD01', beamProp)
job_pym.fanplot(RPM_vec, result_dir)

#==================================================================================
#%%      P L O T


#==========================FAN-PLOT==========================================
res = np.real(job_pym.analysis.freq)

plt.figure()
plt.subplot(121)
plt.grid(True)

#plot rotor-harmonics 
x =  np.linspace(0, 1.2, 20)
y =  x*Omega/(2*np.pi)
for i in range(1,9):
    color = '#333333'
    plt.plot(x,i*y,'--',color='grey')
    string = r'$%i\Omega$' % (i)
    plt.text(x[-1]+.01, i*y[-1], string)

#read and plot reference data:
fname = 'jobs/VariSpeed/uh60a_data_blade/fanplot_uh60a_bowen-davies-PhD.txt'
ref_data = np.loadtxt(fname,skiprows=1,delimiter=',')
ref_str = open(fname).readline().replace('\n','').split(',')
x = ref_data[:,0]
for i,d in enumerate(ref_data.T):
    s=ref_str[i]
    if 'f' in s:
        colorhex = 'blue'
        plt.plot(x, d, '--',color=colorhex)
    elif 'l' in s:
        colorhex = 'red'
        plt.plot(x, d,'--', color=colorhex)
    elif 't' in s:
        colorhex = 'green'
        plt.plot(x,d,'--',color=colorhex)
        
#plot dymore frequencies:
x = RPM_vec/Omega
ref_str = ['l1','f1','f2','f3','f4','t1','f5']
for i,d in enumerate(res[:,:len(ref_str)].T):
    s=ref_str[i]
    plt.plot(x,d,'b')
    if 'f' in s:
        colorhex = 'blue'
        plt.plot(x, d, 'o-', color=colorhex)
        string = r'%s flap' % (s[-1])
        plt.text(x[-1]+.01, d[-1], string, color=colorhex)
    elif 'l' in s:
        colorhex = 'red'
        plt.plot(x, d, 'o-', color=colorhex)
        string = r'%s lead-lag' % (s[-1])
        plt.text(x[-1]+.01, d[-1], string, color=colorhex)
    elif 't' in s:
        colorhex = 'green'
        plt.plot(x, d, 'o-', color=colorhex)
        string = r'%s torsion' % (s[-1])
        plt.text(x[-1]+.01, d[-1], string, color=colorhex)
        
plt.ylim((0,40))
plt.xlim((0,1.2))
plt.title('Fan-Plot')
plt.xlabel(r'Main Rotor Speed $\Omega$ [1/Rev]')
plt.ylabel(r'Eigenfrequencies $\omega$ [Hz]')


#==========================CROSS-SECTIONS==========================================

if flag_ref:
    job.cbm_post_2dmesh(title = 'Reference')
    
if flag_opt:
        val_fname = 'jobs/VariSpeed/uh60a_data_blade/Fanplot_Bowen_Davies_Diss.csv'
        #p.model.marc_comp.job.fanplot_show(p.model.marc_comp.RPM_vec, p.model.marc_comp.result_dir,val_fname=val_fname)

        job_opt.cbm_post_2dmesh(title = 'Optimization')
    
    
#==========================BEAM-PROPERTIES=======================================    
plt.rc('text', usetex=False)
f, axarr = plt.subplots(3,2, sharex=True)    

#---------------m00------------------------------------------------------------
axarr[0,0].plot(dct_davis['mass'][:,0],dct_davis['mass'][:,1],'r:', label='from S.J. Davis (1981, Sikorsky Aircraft Division)')
axarr[0,0].plot(dct_dym['x'],dct_dym['mass_per_unit_span'],'--', label='from DYMORE UH-60A (Yeo)')
#axarr[0,0].plot(dct_interp['x'],dct_interp['mass_per_unit_span'],'gx', label='lin. interp. from DYMORE')
axarr[0,0].plot(job.config.setup['radial_station'],job.BeamProperties.MpUS,'o', color='grey', label='SONATA CBM (VABS)')
if flag_opt:
    axarr[0,0].plot(job_opt.config.setup['radial_station'],job_opt.BeamProperties.MpUS,'k^', label='SONATA CBM OPT w. VABS')
axarr[0,0].set_ylim([5,40])
axarr[0,0].set_ylabel(r'$m_{00}$ [kg/m]')
axarr[0,0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.5), ncol=2)

#---------------Xm2------------------------------------------------------------
axarr[0,1].plot(dct_davis['cg'][:,0],dct_davis['cg'][:,1],'r:')
axarr[0,1].plot(dct_dym['x'],dct_dym['centre_of_mass_location'][:,0]*1000,'--')
#axarr[1,0].plot(dct_interp['x'],dct_interp['centre_of_mass_location'][0]*1000,'gx')
axarr[0,1].plot(job.config.setup['radial_station'],job.BeamProperties.Xm2,'o',color='grey')
if flag_opt:
    axarr[0,1].plot(job_opt.config.setup['radial_station'],job_opt.BeamProperties.Xm2,'k^')
axarr[0,1].set_ylim([-100,100])
axarr[0,1].set_ylabel(r' $X_{m2}$ [mm]')

#---------------EA-------------------------------------------------------------
axarr[1,0].plot(dct_dym['x'],dct_dym['axial_stiffness'],'--')
#axarr[0,1].plot(dct_interp['x'],dct_interp['axial_stiffness'],'gx')
axarr[1,0].plot(job.config.setup['radial_station'],job.BeamProperties.CS[0,0],'o', color='grey')
if flag_opt:
    axarr[1,0].plot(job_opt.config.setup['radial_station'],job_opt.BeamProperties.CS[0,0],'k^')
axarr[1,0].set_ylabel(r'$EA \; [N]$')
axarr[1,0].ticklabel_format(style='sci', axis='y', scilimits=(0,0))

#---------------GJ-------------------------------------------------------------
axarr[1,1].plot(dct_davis['torsional_stiffness'][:,0],dct_davis['torsional_stiffness'][:,1],'r:')
axarr[1,1].plot(dct_dym['x'],dct_dym['torsional_stiffness'],'--')
#axarr[1,1].plot(dct_interp['x'],dct_interp['torsional_stiffness'],'gx')
axarr[1,1].plot(job.config.setup['radial_station'], job.BeamProperties.CS[1,1]*1e-6,'o', color='grey')
if flag_opt:
    axarr[1,1].plot(job_opt.config.setup['radial_station'], job_opt.BeamProperties.CS[1,1]*1e-6,'k^')
axarr[1,1].set_ylabel(r'$GJ \; [Nm^2]$')
axarr[1,1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
x_offset = 0.81786984

#---------------EI2------------------------------------------------------------
axarr[2,0].plot(dct_davis['flapping_stiffness'][:,0],dct_davis['flapping_stiffness'][:,1],'r:')
axarr[2,0].plot(dct_dym['x'],dct_dym['bending_stiffnesses'][:,0],'--') 
#axarr[2,1].plot(dct_interp['x'],dct_interp['bending_stiffnesses'][0],'gx') 
axarr[2,0].plot(job.config.setup['radial_station'], job.BeamProperties.CS[2,2]*1e-6,'o', color='grey') 
if flag_opt:
    axarr[2,0].plot(job_opt.config.setup['radial_station'],job_opt.BeamProperties.CS[2,2]*1e-6,'k^') 
axarr[2,0].set_ylabel(r'$EI_{2} \; [Nm^2]$')
axarr[2,0].ticklabel_format(style='sci', axis='y', scilimits=(0,0))

#---------------EI3------------------------------------------------------------
axarr[2,1].plot(dct_davis['edgewise_stiffness'][:,0],dct_davis['edgewise_stiffness'][:,1],'r:')
axarr[2,1].plot(dct_dym['x'],dct_dym['bending_stiffnesses'][:,1],'--')
#axarr[3,1].plot(dct_interp['x'],dct_interp['bending_stiffnesses'][1],'gx')
axarr[2,1].plot(job.config.setup['radial_station'],job.BeamProperties.CS[3,3]*1e-6,'o', color='grey') 
if flag_opt:
    axarr[2,1].plot(job_opt.config.setup['radial_station'],job_opt.BeamProperties.CS[3,3]*1e-6,'k^') 
axarr[2,1].set_ylabel(r'$EI_3 \; [Nm^2]$')
axarr[2,1].set_xlabel('Radius [mm]')
axarr[2,1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))

#---------------m11------------------------------------------------------------
#axarr[2,0].plot(dct_davis['torsional_inertia'][:,0],dct_davis['torsional_inertia'][:,1],'r:')
#axarr[2,0].plot(dct_dym['x'],dct_dym['moments_of_inertia'][:,0],'--')
##axarr[2,0].plot(dct_interp['x'],dct_interp['moments_of_inertia'x_offset = 0.81786984][0],'gx')
#axarr[2,0].plot(job.config.setup['radial_station'],job.BeamProperties.MMatMC[3,3]*1e-6,'o', color='grey')
#if flag_opt:
#    axarr[2,0].plot(job_opt.config.setup['radial_station'],job_opt.BeamProperties.MMatMC[3,3]*1e-6,'k^')
#axarr[2,0].set_ylabel(r'$m_{11}$ [kg-m]')

#---------------m22------------------------------------------------------------
#axarr[3,0].plot(dct_dym['x'],dct_dym['moments_of_inertia'][:,1],'--')
#axarr[3,0].plot(dct_interp['x'],dct_interp['moments_of_inertia'][1],'gx')
#axarr[3,0].plot(job.config.SETUP_radial_station,job.BeamProperties.MMatMC[4,4]*1e-6,'ko')
#if flag_opt:
#    axarr[3,0].plot(job_opt.config.SETUP_radial_station,job_opt.BeamProperties.MMatMC[4,4]*1e-6,'P', color='orange')
#axarr[3,0].set_ylabel(r'$m_{22}$ [kg-m]')

#---------------m33------------------------------------------------------------
#axarr[3,0].plot(dct_davis['edgewise_inertia'][:,0],dct_davis['edgewise_inertia'][:,1],'r:')
#axarr[3,0].plot(dct_dym['x'],dct_dym['moments_of_inertia'][:,2],'--')
##axarr[3,0].plot(dct_interp['x'],dct_interp['moments_of_inertia'][2],'gx')
#axarr[3,0].plot(job.config.setup['radial_station'],job.BeamProperties.MMatMC[5,5]*1e-6,'o', color='grey')
#if flag_opt:
#    axarr[3,0].plot(job_opt.config.setup['radial_station'],job_opt.BeamProperties.MMatMC[5,5]*1e-6,'k^')
#axarr[3,0].set_ylabel(r'$m_{33}$ [kg-m]')
#axarr[3,0].set_xlabel('Radius [mm]')

#---------------Xs2------------------------------------------------------------
#axarr[0,1].plot(dct_dym['x'],dct_dym['shear_centre_location'][:,0]*1000,'--')
#axarr[0,1].plot(dct_interp['x'],dct_interp['shear_centre_location'][0]*1000,'gx')
#axarr[0,1].plot(job.config.SETUP_radial_station,job.BeamProperx_offset = 0.81786984ties.Xs2,'ko')
#if flag_opt:
#    axarr[0,1].plot(job_opt.config.SETUP_radial_station,job_opt.BeamProperties.Xs2,'P', color='orange')
#axarr[0,1].set_ylim([-100,100])
#axarr[0,1].set_ylabel(r'$X_{s2}$ [mm]')


#from matplotlib2tikz import save as tikz_save
#tikz_save('UH60A_beam.tikz', figureheight='\\figureheight', figurewidth='\\figurewidth' )