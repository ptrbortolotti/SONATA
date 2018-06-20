# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 14:36:09 2018

@author: TPflumm
"""
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

def read_dymore_beam_properties(filename, x_offset = 0.81786984):
    #READ FILE AND CLEAN UP COMMENTS, EMPTY LINES WITH SPACES AND NEWLINES
    a = ''
    with open(filename) as f:
        for line in f:
            line = line.partition('#')[0]
            line = line.rstrip()
            a += line 
            a += '\n'
         
    STR = ''.join([s for s in a.strip().splitlines(True) if s.strip("\r\n").strip()])
    STR = STR.split('STRUCTURAL_DEFINITION')[1]
    STR = STR.split('AERODYNAMIC_DEFINITION')[0]
    
    keys =   ['@CURVILINEAR_COORDINATE      {',\
              '@AXIAL_STIFFNESS         {',\
              '@BENDING_STIFFNESSES     {',\
              '@TORSIONAL_STIFFNESS     {',\
              '@SHEARING_STIFFNESSES    {',\
              '@MASS_PER_UNIT_SPAN      {',\
              '@MOMENTS_OF_INERTIA      {',\
              '@CENTRE_OF_MASS_LOCATION {',\
              '@SHEAR_CENTRE_LOCATION   {',\
              '@CENTROID_LOCATION       {']
              
    dct = {}
    for k in keys:
        regex = k+'(.+?)\}+?'
        tmp_lst = re.findall(regex, STR)
        tmp_lst = [s.split(',') for s in tmp_lst]
        tmp_lst = [list(map(float,s)) for s in tmp_lst]
        arr = np.asarray(tmp_lst)
        dct[k[1:-1].strip().lower()] = arr
    
    dct['x'] = dct['curvilinear_coordinate'] + x_offset
    dct['x'] = dct['x']*1000 #in mm
    return dct


def interp1d_dymore_beam_properties(dct_dym, x):
    '''Interpolate Values from dct_dym at x for optimization'''
    dct_interp = {}
    dct_interp['x'] = x
    for k in dct_dym:
        if k != 'x':
            ynew_tmp = []
            for i, item in enumerate(dct_dym[k].T):
                f = interpolate.interp1d(dct_dym['x'][:,0], item) 
                ynew_tmp.append(f(dct_interp['x']))         
            dct_interp[k] = np.asarray(ynew_tmp)
    return dct_interp



if __name__ == '__main__':
    filename = 'C://TPflumm_local/work/SONATA/jobs/VariSpeed/uh60a_blade/dymore_uh60a_rotor_blade.dat'
    dct = read_dymore_beam_properties(filename,x_offset = 0.81786984)
    plt.rc('text', usetex=False)
    
    f, axarr = plt.subplots(4,2, sharex=True)    
    
    axarr[0,0].plot(dct['x'],dct['mass_per_unit_span'],'o--', label='experiments')
    axarr[0,0].set_ylabel('mass per unit span [kg/m]')
    
    axarr[1,0].plot(dct['x'],dct['centre_of_mass_location'][:,0],'o--')
    axarr[1,0].plot(dct['x'],dct['centre_of_mass_location'][:,1],'o--')  
    axarr[1,0].set_ylabel('CG location [m]')
    
    axarr[3,0].set_xlabel('Radius [m]')
    
    
    axarr[0,1].plot(dct['x'],dct['axial_stiffness'],'o--')
    axarr[0,1].set_ylabel('EA [Nm^2]')
    
    axarr[1,1].plot(dct['x'],dct['torsional_stiffness'],'o--')
    axarr[1,1].set_ylabel('GJ [Nm^2]')
    
    axarr[2,1].plot(dct['x'],dct['bending_stiffnesses'][:,1],'o--')
    axarr[2,1].set_ylabel('EI_2_2')
    
    axarr[3,1].plot(dct['x'],dct['bending_stiffnesses'][:,0],'o--')  
    axarr[3,1].set_ylabel('EI_3_3')
    axarr[3,1].set_xlabel('Radius [m]')
