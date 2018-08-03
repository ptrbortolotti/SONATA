# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 11:18:28 2018

@author: TPflumm
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.transforms import BlendedGenericTransform
from openmdao.api import Problem, ScipyOptimizer, IndepVarComp, ScipyOptimizeDriver 
from openmdao.drivers.genetic_algorithm_driver import SimpleGADriver
from concurrent import futures
from matplotlib2tikz import save as tikz_save

#Sonata Modules: Make Sure to be in the SONATA working directory!
os.chdir('../..')  #print(os.getcwd())

from SONATA.cbm.fileIO.configuration import Configuration
from SONATA.cbm.fileIO.dymore_utils import read_dymore_beam_properties, interp1d_dymore_beam_properties
from SONATA.cbm.fileIO.readinput import read_material_input
from SONATA.cbm.sonata_cbm import CBM
from SONATA.vabs.VABS_interface import VABS_config, export_cells_for_VABS, XSectionalProperties
import SONATA.Pymore.utl.coef as coef

#job specific modules!
from jobs.VariSpeed.sonata_group import Sonata_Group

plt.close('all')
__spec__ = None

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
    job_opt=None
    
    #flag_ref = True
    if flag_ref:
        job = CBM(config)
        job.cbm_gen_topo()
        job.cbm_gen_mesh()
        job.cbm_run_vabs(rm_vabfiles=True)
    
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

#%%#######RUN CBM OPTIMIZATION#############################################################
#==============================================================================
#      Interpolate Values from dct_dym for optimization
#==============================================================================
#READ DYMORE BEAM PROPERTIES:
folder = 'SONATA/Pymore/dym/mdl/03_rotormodel/05_UH60_rotor_optimization/01_UH60_rotor_snglblade_static/'
filename = folder + 'rotor_blade.dat'
x_offset = 0.81786984
dct_dym = read_dymore_beam_properties(filename, x_offset = x_offset)

flag_ref=True 
flag_opt=True 
solver='slsqp'
radial_station = [2000,7500]
with futures.ProcessPoolExecutor(max_workers=6) as e:
    fs = {n:e.submit(run_cbm_optimization, n, dct_dym, flag_ref, flag_opt, solver) for n in radial_station}
    print('Alle Aufgaben gestartet')
print('Alle Aufgaben erledigt.')

#==================================================================================
#%%      P Y M O R E
#==================================================================================
from SONATA.Pymore.marc.marc import MARC
from collections import OrderedDict
dct_cbm_job_ref = OrderedDict()
dct_cbm_job_opt = OrderedDict()

for k in fs:
    res = fs[k].result()
    if res[0] != None: dct_cbm_job_ref[k] = res[0]
    if res[1] != None: dct_cbm_job_opt[k] = res[1]

#TODO: Use a proper interpolation scheme and include start and end properties 
tmp = []
if flag_opt==False:
    for k in dct_cbm_job_ref:
        tmp.append(dct_cbm_job_ref[k].cbm_set_DymoreMK(x_offset))
else:
    for k in dct_cbm_job_opt:
        tmp.append(dct_cbm_job_opt[k].cbm_set_DymoreMK(x_offset))

tmp.insert(0,coef.refBeamProp()[0])
tmp.append(coef.refBeamProp()[-1])

beamProp = np.asarray(tmp)
#dct_dym['mass_per_unit_span'][0]

nbOfLoc = 11
Omega = 4.3*2*np.pi #in rad/sec
RPM_vec = np.linspace(0.2*Omega, 1.2*Omega, nbOfLoc)

dir_root = 'SONATA/Pymore/dym/mdl/03_rotormodel/05_UH60_rotor_optimization/01_UH60_rotor_snglblade_static/'
result_dir = 'SONATA/Pymore/rlt/'

job_pym = MARC(dir_root, 'rotor_assembly.dym')
job_pym.marc_set_beamProp('BLADE_BP_CD01', beamProp)
job_pym.fanplot(RPM_vec, result_dir)

#==================================================================================
#%%      P L O T
#==========================FAN-PLOT==========================================
plt.close('all')
res = np.real(job_pym.analysis.freq)
plt.figure()
#plt.subplot(121)

plt.grid(True)

#plot rotor-harmonics 
x =  np.linspace(0, 1.2, 20)
y =  x*Omega/(2*np.pi)
for i in range(1,9):
    color = '#333333'
    plt.plot(x,i*y,'--',color='grey')
    string = r'$%i\Omega$' % (i)
    plt.text(x[-1]-.06, i*y[-1]+.5, string, color='grey')

#read and plot reference data:
fname = 'jobs/VariSpeed/uh60a_data_blade/fanplot_uh60a_bowen-davies-PhD.txt'
ref_data = np.loadtxt(fname,skiprows=1,delimiter=',')
ref_str = open(fname).readline().replace('\n','').split(',')
x = ref_data[:,0]
for i,d in enumerate(ref_data.T):
    s=ref_str[i]
    if 'f' in s:
        colorhex = 'blue'
        plt.plot(x, d, ':',color=colorhex)
    elif 'l' in s:
        colorhex = 'red'
        plt.plot(x, d,':', color=colorhex)
    elif 't' in s:
        colorhex = 'green'
        plt.plot(x,d,':',color=colorhex)
        
#plot dymore frequencies:
x = RPM_vec/Omega
ref_str = ['l1','f1','f2','f3','l2','t1','f4']
D = {'1':'s','2':'^','3':'o','4':'d'}
ms = 3
for i,d in enumerate(res[:,:len(ref_str)].T):
    s=ref_str[i]
    plt.plot(x,d,'b')
    m = D[s[-1]] 
    if 'f' in s:
        colorhex = 'blue'
        plt.plot(x, d, '-', color=colorhex, marker=m, markersize = ms)
        string = r'%s flap' % (s[-1])
        plt.text(x[-1]+.01, d[-1], string, color=colorhex)
    elif 'l' in s:
        colorhex = 'red'
        plt.plot(x, d, 'o-', color=colorhex, marker=m, markersize = ms)
        string = r'%s lead-lag' % (s[-1])
        plt.text(x[-1]+.01, d[-1], string, color=colorhex)
    elif 't' in s:
        colorhex = 'green'
        plt.plot(x, d, 'o-', color=colorhex, marker=m, markersize = ms)
        string = r'%s torsion' % (s[-1])
        plt.text(x[-1]+.01, d[-1], string, color=colorhex)
        
line1 = mlines.Line2D([], [], color='black', linestyle='-', marker='o', label='Eigenfrequencies')
line2 = mlines.Line2D([], [], color='black', linestyle=':', label='UH-60A Reference')

plt.ylim((0,45))
plt.xlim((0,1.2))
plt.title('Fan-Plot')
plt.xlabel(r'Rotor Rotational Speed, $\Omega / \Omega_{ref}$')
plt.ylabel(r'Eigenfrequencies, $\omega$ [Hz]')
plt.legend(handles=[line1,line2])
plt.show()
#tikz_save('/media/gu32kij/HTMWTUM/Oeffentlich/Publikationen/ERF/2018/TPflumm, WGarre/paper/img/UH60A_fanplot.tikz', figureheight='\\figureheight', figurewidth='\\figurewidth' )

#==========================EIGEN-MODES==========================================
eigv = job_pym.analysis.eigv

blade_len = 7.36082856
r_attachment = 0.81786984
r_hinge = 0.378
station = 8
blade_with_att = blade_len + r_attachment - r_hinge
R = blade_len + r_attachment
#        fan_plot(np.real(freq), val_fname, RPM_vec, result_dir)
#sim_plot(freq, RPM_vec)

plt.figure()
plt.subplot(311)
i = 1
#            plt.plot(eigVec[0,70*i:70*(i+1)], 'k--')
pos = (np.hstack((0.0,np.linspace(3.39,42.39,40.0)))*blade_with_att/42.39 + 0.378)
IDs = np.hstack((70*(i+1), np.linspace(30+70*i,70*(i+1)-1,40,dtype=int)))


for x in range(eigv[:,IDs,:].shape[0]):
    for z in range(eigv[:,IDs,:].shape[2]):
        eigv[x,0,z] = np.interp(r_attachment, pos, eigv[x,IDs,z])

pos = (np.hstack((0.0,np.linspace(2.39,42.39,41.0)))*blade_with_att/42.39 + 0.378)/R
IDs = np.hstack((70*(i+1), 0, np.linspace(30+70*i,70*(i+1)-1,40,dtype=int)))

scale_vals = np.amax(eigv[:,IDs,:], axis = 1)
scale_vals_min = np.amin(eigv[:,IDs,:], axis = 1)

for index, x in np.ndenumerate(scale_vals):
    if np.absolute(scale_vals_min[index]) > x:
        scale_vals[index] = scale_vals_min[index]

#scale_vals = np.ones(scale_vals.shape)

plt.plot(pos,eigv[0,IDs,station]/scale_vals[0,station], color='red', marker='s', markersize=1.5, label='1. lead-lag')
#plt.plot(pos,eigv[1,IDs,station]/scale_vals[1,station], 'k-.^', label='2.mode: lead-lag')
#plt.plot(pos,eigv[2,IDs,station]/scale_vals[2,station], 'k:',   label='3.mode: lead-lag')
plt.plot(pos,eigv[3,IDs,station]/scale_vals[3,station],  color='blue', marker='o', markersize=1.5, label='3. flap')
plt.plot(pos,eigv[4,IDs,station]/scale_vals[4,station],  color='red', marker='^', markersize=1.5, label='2. lead-lag')
#plt.plot(pos,eigv[5,IDs,station]/scale_vals[5,station], 'k--+', label='6.mode: lead-lag')
#        plt.plot(pos,eigv[6,IDs,station]/scale_vals[6,station], 'k-',   label='7.mode: lead-lag')

plt.grid()
#plt.legend()
plt.ylim([-1.05,1.05])
plt.ylabel('Normalized Lead-Lag')

plt.subplot(312)
i = 2
#            plt.plot(eigVec[0,70*i:70*(i+1)], 'k--')
#            plt.plot(eigVec[0,30+70*i:70*(i+1)], 'r--', label='1.mode: flap')
pos = (np.hstack((0.0,np.linspace(3.39,42.39,40.0)))*blade_with_att/42.39 + 0.378)
IDs = np.hstack((70*(i+1), np.linspace(30+70*i,70*(i+1)-1,40,dtype=int)))
for x in range(eigv[:,IDs,:].shape[0]):
    for z in range(eigv[:,IDs,:].shape[2]):
        eigv[x,0,z] = np.interp(r_attachment, pos, eigv[x,IDs,z])

pos = (np.hstack((0.0,np.linspace(2.39,42.39,41.0)))*blade_with_att/42.39 + 0.378)/R
IDs = np.hstack((70*(i+1), 0, np.linspace(30+70*i,70*(i+1)-1,40,dtype=int)))

scale_vals = np.amax(eigv[:,IDs,:], axis = 1)
scale_vals_min = np.amin(eigv[:,IDs,:], axis = 1)

for index, x in np.ndenumerate(scale_vals):
    if np.absolute(scale_vals_min[index]) > x:
        scale_vals[index] = scale_vals_min[index]

#scale_vals = np.ones(scale_vals.shape)

#plt.plot(pos,eigv[0,IDs,station]/scale_vals[0,station], 'k-.',  label='1.mode: flap')       
plt.plot(pos,eigv[1,IDs,station]/scale_vals[1,station], color='blue', marker='s', markersize=1.5, label='1. flap')
plt.plot(pos,eigv[2,IDs,station]/scale_vals[2,station], color='blue', marker='^', markersize=1.5, label='2. flap')
plt.plot(pos,eigv[3,IDs,station]/scale_vals[3,station], color='blue', marker='o', markersize=1.5, label='3. flap')
plt.plot(pos,eigv[4,IDs,station]/scale_vals[4,station], color='red', marker='^', markersize=1.5, label='2. lead-lag')
#plt.plot(pos,eigv[5,IDs,station]/scale_vals[4,station], 'k--+', label='6.mode: flap')
plt.plot(pos,eigv[6,IDs,station]/scale_vals[6,station], color='blue', marker='d', markersize=1.5, label='4. flap')
plt.grid()
#plt.legend(loc='lower center', ncol=5)
plt.ylim([-1.05,1.05])
plt.ylabel('Normalized Flap')

plt.subplot(313)
i = 3
pos = (np.hstack((0.0,np.linspace(3.39,42.39,40.0)))*blade_with_att/42.39 + 0.378)
IDs = np.hstack((70*(i+1), np.linspace(30+70*i,70*(i+1)-1,40,dtype=int)))
for x in range(eigv[:,IDs,:].shape[0]):
    for z in range(eigv[:,IDs,:].shape[2]):
        eigv[x,0,z] = 0.0

pos = (np.hstack((0.0,np.linspace(2.39,42.39,41.0)))*blade_with_att/42.39 + 0.378)/R
IDs = np.hstack((70*(i+1), 0, np.linspace(30+70*i,70*(i+1)-1,40,dtype=int)))


scale_vals = np.amax(eigv[:,IDs,:], axis = 1)
scale_vals_min = np.amin(eigv[:,IDs,:], axis = 1)

for index, x in np.ndenumerate(scale_vals):
    if np.absolute(scale_vals_min[index]) > x:
        scale_vals[index] = scale_vals_min[index]

#scale_vals = np.ones(scale_vals.shape)

#plt.plot(pos,eigv[0,IDs,station]/scale_vals[0,station], 'k-.',  label='1.mode: torsion')
#plt.plot(pos,eigv[1,IDs,station]/scale_vals[1,station], 'k-.^', label='2.mode: torsion')
#plt.plot(pos,eigv[2,IDs,station]/scale_vals[2,station], 'k:',   label='3.mode: torsion')
#plt.plot(pos,eigv[3,IDs,station]/scale_vals[3,station], 'k--',  label='4.mode: torsion')
#plt.plot(pos,eigv[4,IDs,station]/scale_vals[4,station], 'k--*', label='5.mode: torsion')
plt.plot(pos,eigv[5,IDs,station]/scale_vals[5,station], color='green', marker='s', markersize=1.5)
#plt.plot(pos,eigv[6,IDs,station]/scale_vals[6,station], 'k-',   label='7.mode: torsion')

line1 = mlines.Line2D([], [], color='red', linestyle='-', marker='s', label='1 lead-lag')
line2 = mlines.Line2D([], [], color='blue', linestyle='-',marker='s', label='1 flap')
line3 = mlines.Line2D([], [], color='blue', linestyle='-',marker='^', label='2 flap')
line4 = mlines.Line2D([], [], color='blue', linestyle='-',marker='o', label='3 flap')
line5 = mlines.Line2D([], [], color='red', linestyle='-',marker='^',  label='2 lead-lag')
line6 = mlines.Line2D([], [], color='green',linestyle='-',marker='s', label='1 torsion')
line7 = mlines.Line2D([], [], color='blue', linestyle='-', marker='d', label='4 flap')

plt.grid()
#plt.legend()
plt.subplots_adjust(bottom=0.3)
plt.legend(handles=[line1,line2,line3,line4,line5,line6,line7],loc='lower center', ncol=4, bbox_to_anchor=(0.5, -0.7),)
plt.ylim([-1.05,1.05])
plt.xlabel('Radial Station, r/R')
plt.ylabel('Normalized Torsion')
#plt.subplots_adjust(wspace=0.3, left=0.1, right=0.9)
#tikz_save('/media/gu32kij/HTMWTUM/Oeffentlich/Publikationen/ERF/2018/TPflumm, WGarre/paper/img/UH60A_eigenmodes.tikz', figureheight='\\figureheight', figurewidth='\\figurewidth' )

#%%==========================CROSS-SECTIONS==========================================
dest_folder = '/media/gu32kij/HTMWTUM/Oeffentlich/Publikationen/ERF/2018/TPflumm, WGarre/paper/img/'
for k in dct_cbm_job_ref:
    title = 'Reference at R=%i' % (k)
    dct_cbm_job_ref[k].cbm_post_2dmesh(title = title)
    savepath = dest_folder+'Reference_R%i.jpg' % (k)
    f = plt.gcf()  # f = figure(n) if you know the figure number
    f.set_size_inches(50,25)
    plt.tight_layout()
    plt.savefig(savepath, dpi=200)
    
for k in dct_cbm_job_opt:
    title = 'Optimization Result at R=%i' % (k)
    dct_cbm_job_opt[k].cbm_post_2dmesh(title = title)
    plt.tight_layout()
    savepath = dest_folder+'Optimization_%i.jpg' % (k)
    f = plt.gcf()  # f = figure(n) if you know the figure number
    f.set_size_inches(50,25)
    plt.savefig(savepath, dpi=200)
    
    
#%%==========================BEAM-PROPERTIES=======================================    
#plt.rc('text', usetex=False)
f, axarr = plt.subplots(3,2, sharex=True)    

#---------------m00------------------------------------------------------------
#axarr[0,0].plot(dct_davis['mass'][:,0],dct_davis['mass'][:,1],'r:', label='from S.J. Davis (1981, Sikorsky Aircraft Division)')
axarr[0,0].plot(dct_dym['x'],dct_dym['mass_per_unit_span'],'--', label='UH-60A Reference')
#axarr[0,0].plot(dct_interp['x'],dct_interp['mass_per_unit_span'],'gx', label='lin. interp. from DYMORE')

tmp_lst = []
for k in dct_cbm_job_ref:
    job = dct_cbm_job_ref[k]
    tmp_lst.append([job.config.setup['radial_station'],job.BeamProperties.MpUS])
tmp_arr = np.asarray(tmp_lst)
axarr[0,0].plot(tmp_arr[:,0],tmp_arr[:,1],'o', color='grey', label='Initial Cross-Sections')

tmp_lst = [[dct_dym['x'][0][0],dct_dym['mass_per_unit_span'][0][0]]]
for k in dct_cbm_job_opt:
    job_opt = dct_cbm_job_opt[k]
    tmp_lst.append([job_opt.config.setup['radial_station'],job_opt.BeamProperties.MpUS])
tmp_lst.append([dct_dym['x'][-1][0],dct_dym['mass_per_unit_span'][-1][0]])    
tmp_arr = np.asarray(tmp_lst)

axarr[0,0].plot(tmp_arr[:,0],tmp_arr[:,1],'k^:', markerfacecolor='none', label='New Composite Beam')
axarr[0,0].plot(tmp_arr[1:-1,0],tmp_arr[1:-1,1],'k^', label='Optimized Cross-Sections')

axarr[0,0].set_ylim([0,50])
axarr[0,0].set_ylabel(r'$m_{00}$ [kg/m]')


#---------------Xm2------------------------------------------------------------
#axarr[0,1].plot(dct_davis['cg'][:,0],dct_davis['cg'][:,1],'r:')
axarr[0,1].plot(dct_dym['x'],dct_dym['centre_of_mass_location'][:,0]*1000,'--')
axarr[0,1].plot(dct_dym['x'][0],dct_dym['centre_of_mass_location'][0,0]*1000,'k^',markerfacecolor='none')
#axarr[1,0].plot(dct_interp['x'],dct_interp['centre_of_mass_location'][0]*1000,'gx')
for k in dct_cbm_job_ref:
    job = dct_cbm_job_ref[k]
    axarr[0,1].plot(job.config.setup['radial_station'],job.BeamProperties.Xm2,'o',color='grey')
    
tmp_lst = [[dct_dym['x'][0][0],dct_dym['centre_of_mass_location'][0,0]*1000]]   
for k in dct_cbm_job_opt:
    job_opt = dct_cbm_job_opt[k]
    tmp_lst.append([job_opt.config.setup['radial_station'],job_opt.BeamProperties.Xm2])
tmp_lst.append([dct_dym['x'][-1][0], dct_dym['centre_of_mass_location'][-1,0]*1000])    
tmp_arr = np.asarray(tmp_lst)
axarr[0,1].plot(tmp_arr[:,0],tmp_arr[:,1],'k^:', markerfacecolor='none')
axarr[0,1].plot(tmp_arr[1:-1,0],tmp_arr[1:-1,1],'k^')
axarr[0,1].set_ylim([-120,120])
axarr[0,1].set_ylabel(r' $X_{m2}$ [mm]')

#---------------EA-------------------------------------------------------------
axarr[1,0].plot(dct_dym['x'],dct_dym['axial_stiffness'],'--')
#axarr[0,1].plot(dct_interp['x'],dct_interp['axial_stiffness'],'gx')
for k in dct_cbm_job_ref:
    job = dct_cbm_job_ref[k]
    axarr[1,0].plot(job.config.setup['radial_station'],job.BeamProperties.CS[0,0],'o', color='grey')

    
tmp_lst = [[dct_dym['x'][0][0],dct_dym['axial_stiffness'][0][0]]]       
for k in dct_cbm_job_opt:
    job_opt = dct_cbm_job_opt[k]
    tmp_lst.append([job_opt.config.setup['radial_station'],job_opt.BeamProperties.CS[0,0]])
tmp_lst.append([dct_dym['x'][-1][0],dct_dym['axial_stiffness'][-1][0]])      
tmp_arr = np.asarray(tmp_lst)
axarr[1,0].plot(tmp_arr[:,0],tmp_arr[:,1],'k^:', markerfacecolor='none')
axarr[1,0].plot(tmp_arr[1:-1,0],tmp_arr[1:-1,1],'k^') 

axarr[1,0].set_ylabel(r'$EA \; [N]$')
axarr[1,0].set_ylim((0,1.4e9))
axarr[1,0].ticklabel_format(style='sci', axis='y', scilimits=(0,0))

#---------------GJ-------------------------------------------------------------
#axarr[1,1].plot(dct_davis['torsional_stiffness'][:,0],dct_davis['torsional_stiffness'][:,1],'r:')
axarr[1,1].plot(dct_dym['x'],dct_dym['torsional_stiffness'],'--')
#axarr[1,1].plot(dct_interp['x'],dct_interp['torsional_stiffness'],'gx')

for k in dct_cbm_job_ref:
    job = dct_cbm_job_ref[k]
    axarr[1,1].plot(job.config.setup['radial_station'], job.BeamProperties.CS[1,1]*1e-6,'o', color='grey')
    
tmp_lst = [[dct_dym['x'][0][0],dct_dym['torsional_stiffness'][0][0]]]       
for k in dct_cbm_job_opt:
    job_opt = dct_cbm_job_opt[k]
    tmp_lst.append([job_opt.config.setup['radial_station'],job_opt.BeamProperties.CS[1,1]*1e-6])
tmp_lst.append([dct_dym['x'][-1][0],dct_dym['torsional_stiffness'][-1][0]])      
tmp_arr = np.asarray(tmp_lst)
axarr[1,1].plot(tmp_arr[:,0],tmp_arr[:,1],'k^:', markerfacecolor='none')
axarr[1,1].plot(tmp_arr[1:-1,0],tmp_arr[1:-1,1],'k^')  

tmp_arr = np.asarray(tmp_lst)
axarr[1,1].set_ylabel(r'$GJ \; [Nm^2]$')
axarr[1,1].set_ylim((0.5e5,2.25e5))
axarr[1,1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
x_offset = 0.81786984

#---------------EI2------------------------------------------------------------
#axarr[2,0].plot(dct_davis['flapping_stiffness'][:,0],dct_davis['flapping_stiffness'][:,1],'r:')
axarr[2,0].plot(dct_dym['x'],dct_dym['bending_stiffnesses'][:,0],'--') 
#axarr[2,1].plot(dct_interp['x'],dct_interp['bending_stiffnesses'][0],'gx') 
for k in dct_cbm_job_ref:
    job = dct_cbm_job_ref[k]
    axarr[2,0].plot(job.config.setup['radial_station'], job.BeamProperties.CS[2,2]*1e-6,'o', color='grey') 
    
tmp_lst = [[dct_dym['x'][0][0],dct_dym['bending_stiffnesses'][:,0][0]]]       
for k in dct_cbm_job_opt:
    job_opt = dct_cbm_job_opt[k]
    tmp_lst.append([job_opt.config.setup['radial_station'],job_opt.BeamProperties.CS[2,2]*1e-6])
tmp_lst.append([dct_dym['x'][-1][0],dct_dym['bending_stiffnesses'][:,0][-1]])      
tmp_arr = np.asarray(tmp_lst)
axarr[2,0].plot(tmp_arr[:,0],tmp_arr[:,1],'k^:', markerfacecolor='none')    
axarr[2,0].plot(tmp_arr[1:-1,0],tmp_arr[1:-1,1],'k^')  

axarr[2,0].set_ylabel(r'$EI_{2} \; [Nm^2]$')
axarr[2,0].set_ylim((0,3e5))
axarr[2,0].ticklabel_format(style='sci', axis='y', scilimits=(0,0))

#---------------EI3------------------------------------------------------------
#axarr[2,1].plot(dct_davis['edgewise_stiffness'][:,0],dct_davis['edgewise_stiffness'][:,1],'r:')
axarr[2,1].plot(dct_dym['x'],dct_dym['bending_stiffnesses'][:,1],'--')
#axarr[3,1].plot(dct_interp['x'],dct_interp['bending_stiffnesses'][1],'gx')
for k in dct_cbm_job_ref:
    job = dct_cbm_job_ref[k]
    axarr[2,1].plot(job.config.setup['radial_station'],job.BeamProperties.CS[3,3]*1e-6,'o', color='grey') 
    
tmp_lst = [[dct_dym['x'][0][0],dct_dym['bending_stiffnesses'][:,1][0]]]     
for k in dct_cbm_job_opt:
    job_opt = dct_cbm_job_opt[k]
    tmp_lst.append([job_opt.config.setup['radial_station'],job_opt.BeamProperties.CS[3,3]*1e-6])
tmp_lst.append([dct_dym['x'][-1][0],dct_dym['bending_stiffnesses'][:,1][-1]])      
tmp_arr = np.asarray(tmp_lst)
axarr[2,1].plot(tmp_arr[:,0],tmp_arr[:,1],'k^:', markerfacecolor='none')      
axarr[2,1].plot(tmp_arr[1:-1,0],tmp_arr[1:-1,1],'k^')  
    
axarr[2,1].set_ylabel(r'$EI_3 \; [Nm^2]$')
axarr[2,1].set_ylim((0,3.5e6))
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

axarr[0,0].legend(loc='upper center', bbox_to_anchor=(0.1, 1.5), ncol=2)
f.subplots_adjust(bottom=0.3, wspace=0.5)
plt.show()
#from matplotlib2tikz import save as tikz_save
#tikz_save('/media/gu32kij/HTMWTUM/Oeffentlich/Publikationen/ERF/2018/TPflumm, WGarre/paper/img/UH60A_beam.tikz', figureheight='\\figureheight', figurewidth='\\figurewidth' )