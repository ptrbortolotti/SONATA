# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 10:42:32 2018

@author: TPflumm
"""

from SONATA.Pymore.marc.marc import MARC
import numpy as np
import matplotlib.pyplot as plt


nbOfLoc = 11
Omega = 4.3*2*np.pi #in rad/sec
RPM_vec = np.linspace(0.2*Omega, 1.15*Omega, nbOfLoc)


dir_root = 'SONATA/Pymore/dym/mdl/03_rotormodel/05_UH60_rotor_optimization/01_UH60_rotor_snglblade_static/'
result_dir = 'SONATA/Pymore/rlt/'

job_pym = MARC(dir_root, 'rotor_assembly.dym')
#job_pym.marc_set_beamProp('BLADE_BP_CD01', beamProp)
job_pym.fanplot(RPM_vec, result_dir)
#job_pym.fanplot_show(RPM_vec, result_dir)

#%% plot fan-plot

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
