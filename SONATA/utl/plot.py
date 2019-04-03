#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 14:49:51 2018

@author: gu32kij
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.lines as mlines
from matplotlib.widgets import TextBox
#from matplotlib2tikz import save as tikz_save
#plt.rc('text', usetex=True)


def plot_histogram_2Ddata(data, **kwargs):
    """
    simple plotting procedure to illustrate multiple datasets of a 2D array in 
    in a histogram.
    
    Parameters
    ----------
    data : ndarray
        array.shape = (sample Nb, i, j)
        
    ref : ndarray, optional
    upper_tri : bool
    title : str
    ylabel : list of str
        
    """ 
    
    #Define function Args and Defaults for kwargs
    shape = data.shape[1:]
    ref  = None
    upper_tri = True
    title='No Title'
    ylabel = [[('k%s%s' % (i+1,j+1)) for j in range(shape[0])] for i in range(shape[1])]
    
    if 'ref' in kwargs:
        if isinstance(kwargs['ref'], (np.ndarray)) and kwargs['ref'].shape == data.shape[1:]:
            ref = kwargs['ref']
        else:
            print ('ref must provide a ndarray with shape %s' % str(data.shape[1:]))
        
    if 'upper_tri' in kwargs:
        if type(kwargs['upper_tri']) == bool:
            upper_tri = kwargs['upper_tri']
        else:
            print('upper_tri must provide a bool value')
        
    if 'title' in kwargs:
        if type(kwargs['title']) == str:
            title = kwargs['title']
        else:
            print('title must provide as string')
            
    if 'ylabel' in kwargs:
        #list of strings
        ylabel = kwargs['ylabel']
        
    #init figure
    shape = data.shape[1:]
    fig, ax = plt.subplots(shape[0],shape[1], sharex=False)
    fig.suptitle(title, fontsize=14)
    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    
    #upper triangle 
    if upper_tri:
        ut = np.zeros((shape[0],shape[1]))            
        ut[np.triu_indices(shape[0])] = 1      
    else:
        ut = np.ones((shape[0],shape[1]))            
    
    for i in range(shape[0]):
        for j in range(shape[1]):
            
            if ut[i,j] == 1 and ref[i,j] != 0:     
          
                if isinstance(ref, (np.ndarray)): 
                    if ref[i,j] != 0:
                        arr = ( data[:,i,j] - ref[i,j] ) / ref[i,j] * 100
                        sigma = arr.std()
                        mu = arr.mean()
                        string = r'$\sigma$ = %.2f $\%%$' % sigma
                        ax[i][j].text(mu,0,string)
                        ax[i][j].set_xlim(-20,20)

                    else:
                        arr = data[:,i,j]

                    count, bins, ignored = ax[i][j].hist(arr, 30, density=True, alpha=0.5, label='histogram')
                    #print(count,bins,ignored)
                    sigma = arr.std()
                    mu = arr.mean()
                    if sigma > 0:
                        ax[i][j].plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) *
                                       np.exp( - (bins - mu)**2 / (2 * sigma**2) ),
                                 linewidth=2, color='r', linestyle='--', alpha=0.5, label='ML approx.')                
                else:
                    arr = data[:,i,j]   
                    count, bins, ignored = ax[i][j].hist(arr, 30, density=True, alpha=0.5, label='histogram')
  
                
            elif ut[i,j] == 0:
                fig.delaxes(ax[i][j])
                
            if i==j:
                ax[i][j].set_xlabel(r'$\%%$ deviation from baseline value')
                if j==shape[1]-1:
                    ax[i][j].legend()
            
            ax[i][j].set_ylabel(ylabel[i][j])
        
    return None


def plot_beam_properties(data, sigma=None, ref=None, x_offset = 0, description = True):
    '''
    generates a plot of the beamproperties using the 
    
    Parameters
    ----------
    Data : np.ndarray 
        Massterms(6), Stiffness(21), damping(1), and coordinate(1) stacked 
        horizontally for the use in DYMORE/PYMORE/MARC  
    sigma: np.ndarray
        standart deviation of Data
    
    ref: np.ndarray
        of reference to select = 'all', 'Massterms', 'Stiffness'
    
    x_offset : float 
        x offset if curvilinear coordinate doen't start at the rotational 
        center of the rotordefault UH-60A value = 0.8178698
            
    '''
      
    #Inertia Properties:
    ylabel = [(r'm_{00}','kg/m'), (r'm_{00}X_{m2}','kg'), (r'm_{00}X_{m3}','kg'), (r'm_{33}','kg m'), (r'm_{23}','kg m'), (r'm_{22}','kg m')]
    x = data[:,-1] + x_offset
    if isinstance(ref, (np.ndarray)):
        ref_x = ref[:,-1]  + x_offset
    fig1, ax1 = plt.subplots(2, 3, sharex=True)
    fig1.subplots_adjust(wspace=0.3, hspace=0.3)
    fig1.suptitle('Inertia Properties', fontsize=14)

    c = 0
    for i in range(2):
        for j in range(3):
            if isinstance(ref, (np.ndarray)):
                ax1[i][j].plot(ref_x,ref[:,c],'-xb')
            if sigma:
                ax1[i][j].fill_between(x, data[:,c]-sigma[:,c], data[:,c]+sigma[:,c], alpha=0.25, edgecolor='r',  linestyle=':', facecolor='r', antialiased=True,)
            ax1[i][j].plot(x,data[:,c],'--k.')
            ax1[i][j].ticklabel_format(axis='y',style='sci')
            ax1[i][j].set_xlabel(r'$r \quad [1/R]$')
            tmp = r'$%s \quad [%s]$' % (ylabel[c][0],ylabel[c][1])
            ax1[i][j].set_ylabel(tmp)

            c += 1
    
    #6x6 Stiffness Properties:
    TSunits = np.array([['N','N','N','Nm','Nm','Nm'],
                        ['N','N','N','Nm','Nm','Nm'],
                        ['N','N','N','Nm','Nm','Nm'],
                        ['Nm','Nm','Nm','Nm^2','Nm^2','Nm^2'],
                        ['Nm','Nm','Nm','Nm^2','Nm^2','Nm^2'],
                        ['Nm','Nm','Nm','Nm^2','Nm^2','Nm^2']])
    
    ut = np.zeros((6,6))            
    ut[np.triu_indices(6)] = 1          
    
    fig2, ax2 = plt.subplots(6, 6)
    fig2.subplots_adjust(wspace=0.5, hspace=0.5)
    fig2.suptitle('6x6 Stiffness Matrix', fontsize=14)
    
    c = 6
    for j in range(6):
        for i in range(6):
            if ut[i,j] == 1:
                if isinstance(ref, (np.ndarray)):
                    ax2[i][j].plot(ref_x,ref[:,c],'-xb')
                if sigma:
                    ax2[i][j].fill_between(x, data[:,c]-sigma[:,c], data[:,c]+sigma[:,c], alpha=0.25, edgecolor='r',  linestyle=':', facecolor='r', antialiased=True,)
                ax2[i][j].plot(x,data[:,c],'--k.')
                ylabel = r'$k_{%s%s} \quad [%s]$' % (i+1,j+1,TSunits[i,j])
                ax2[i][j].set_ylabel(ylabel)
                #ax2[i][j].set_ylim(0)
                #ax2[i][j].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
                
                if i==j:
                    ax2[i][j].set_xlabel(r'$r \quad [1/R]$')
                
                c += 1
            else:
                fig2.delaxes(ax2[i][j])
                #fig2.text(0,0,r'$A=B\cdotC$')
                #fig2.text(0,0,r'$\left(\begin{matrix}k11&k12&k13&k14&k15&k16\\k12&k22&k23&k24&k25&k26\\k13&k23&k33&k34&k35&k36\\k14&k24&k34&k44&k45&k46\\k15&k25&k35&k45&k55&k56\\k16&k26&k36&k46&k56&k66\end{matrix}\right)$')


    desc = r'\begin{minipage}[b]{10cm} '\
        r'\underline{\textbf{Description:}} \\' \
        r'The 6x6 sectional stiffness matrix, TS (Timoshenko Stiffness Matrix) ' \
        r'(1-extension; 2,3-shear, 4-twist; 5,6-bending) relates the sectional axial ' \
        r'strain, $\epsilon_1$, transverse shearing strains, $\epsilon_2$ and $\epsilon_3$, '\
        r'twisting curvatures, $\kappa_1$ and two bending curvatures, $\kappa_2$ and $\kappa_3$, '\
        r'to the axial force, $F_1$, transverse shear forces, $F_2$ and $F_3$, twisting '\
        r'moment, $M_1$, and two bending moments, $M_2$ and $M_3$. The relationship between '\
        r'these sectional strains and sectional stress resultants takes the form of a'\
        r'symmetric, 6x6 matrix: \\' \
        r'$$ \left( \begin{matrix} F_{1} \\ F_{2} \\ F_{3} \\ M_{1} \\ M_{2} \\ M_{3} \end{matrix} \right) = ' \
        r'\left( \begin{matrix} k_{11} & k_{12} & k_{13} & k_{14} & k_{15} & k_{16} \\ ' \
                                r'k_{12} & k_{22} & k_{23} & k_{24} & k_{25} & k_{26} \\ ' \
                                r'k_{13} & k_{23} & k_{33} & k_{34} & k_{35} & k_{36} \\ ' \
                                r'k_{14} & k_{24} & k_{34} & k_{44} & k_{45} & k_{46} \\ ' \
                                r'k_{15} & k_{25} & k_{35} & k_{45} & k_{55} & k_{56} \\ ' \
                                r'k_{16} & k_{26} & k_{36} & k_{46} & k_{56} & k_{66} \end{matrix} \right) '\
          r'\cdot \left( \begin{matrix} \epsilon_{1} \\ \epsilon_{2} \\ \epsilon_{3} \\ \kappa_{1} \\ \kappa_{2} \\ \kappa_{3} \end{matrix} \right) $$ \end{minipage} '

    if description == True:
        from matplotlib import rcParams
        rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
        plt.figtext(0.05, 0.05, desc, usetex=True, wrap=True,  bbox=dict(ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)))        
    plt.show()



def plot_fandiagram(res, Omega, RPM_vec, **kwargs):
    """
    plotting procedure to generate a Fan-Diagram (rotor eigenfrequencies vs.
    rotor rotational speed).
    
    Parameters
    ----------
    res : np.ndarray 
        result
    
    Omega : float
        reference Rotational Speed to normalize frequencies in rad/sec
    
    RPM_Vec : np.ndarray, of 
    sigma : np.ndarray of the standart deviation corresponding to the fanplot.
    ref_fname : filename of the reference data.
    ref_str : list of strings to indentify the eigenmodes 
                (f1= first flap, l2= second lag, t1=first torsion)

    """
    
    #Define function Args and Defaults for kwargs
    sigma  = None
    ref_fname = None
    ref_str = ['x1','x1','x1','x1','x1','x1','x1']
    title='Fan Diagram'
    
    if 'sigma' in kwargs:
        sigma = kwargs['sigma']
    if 'ref_fname' in kwargs:
        ref_fname = kwargs['ref_fname']
    if 'ref_str' in kwargs:
        ref_str = kwargs['ref_str']
    if 'title' in kwargs:
        title = kwargs['title']


    #init figure
    plt.figure()
    plt.grid(True)
    legend_lines = []
    
    #plot rotor-harmonics 
    x =  np.linspace(0, 1.2, 20)
    y =  x*Omega/(2*np.pi)
    for i in range(1,9):
        #color = '#333333'
        plt.plot(x,i*y,'--',color='grey')
        string = r'$%i\Omega$' % (i)
        plt.text(x[-1]-.06, i*y[-1]+.5, string, color='grey')
    
    
    #read and plot reference data:
    if ref_fname != None:
        #fname = 'jobs/VariSpeed/uh60a_data_blade/fanplot_uh60a_bowen-davies-PhD.txt'
        ref_data = np.loadtxt(ref_fname,skiprows=1,delimiter=',')
        ref_str2 = open(ref_fname).readline().replace('\n','').split(',')
        x = ref_data[:,0]
        for i,d in enumerate(ref_data.T):
            s=ref_str2[i]
            if 'f' in s:
                colorhex = 'blue'
                plt.plot(x, d, ':',color=colorhex)
            elif 'l' in s:
                colorhex = 'red'
                plt.plot(x, d,':', color=colorhex)
            elif 't' in s:
                colorhex = 'green'
                plt.plot(x,d,':',color=colorhex)
        legend_lines.append(mlines.Line2D([], [], color='black', linestyle=':', label='Reference')) 
    

    #plot standart deviation and fill between -sigma and +sigma
    res = np.real(res)
    x = RPM_vec/Omega
    if isinstance(sigma, (np.ndarray)):
        for i,d in enumerate(res[:,:len(ref_str)].T):
            #print(d.shape, sigma[:,i].shape)
            plt.fill_between(x, d-sigma[:,i], d+sigma[:,i], alpha=0.25, edgecolor='r',  linestyle=':', facecolor='r', antialiased=True,)
        legend_lines.append(mlines.Line2D([], [], color='red', linestyle=':', label=r'standard deviation $\pm \sigma$'))    
         
        
    #plot dymore frequencies:
    x = RPM_vec/Omega
    #ref_str = ['l1','f1','f2','f3','l2','t1','f4']
    D = {'1':'s','2':'^','3':'o','4':'d'}
    ms = 3
    for i,d in enumerate(res[:,:len(ref_str)].T):
        s=ref_str[i]
        #plt.plot(x,d,'b')
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
            
        else:
            colorhex = 'black'
            plt.plot(x, d, 'o-', color=colorhex, marker=m, markersize = ms)
            
    legend_lines.append(mlines.Line2D([], [], color='black', linestyle='-', marker='o', label='mean eigenfrequencies'))


    plt.ylim((0,45))
    plt.xlim((0,1.2))
    plt.title(title)
    plt.xlabel(r'Rotor Rotational Speed, $\Omega / \Omega_{ref}$')
    plt.ylabel(r'Eigenfrequencies, $\omega$ [Hz]')
    plt.legend(handles=legend_lines)
    plt.show()

    return None


def plot_eigenmodes(eigv, string='BLADE_BP_CG01', **kwargs):
    
    '''
    TBD:
    For the general Illustration of the eigenmodes
    - MARC / DYMORE Interface:
     getter function for pos and IDs of the body to extract coordinates and 
     entry location in eigv matrix
        
    '''
        #==========================EIGEN-MODES==========================================
    blade_len = 7.36082856
    r_attachment = 0.81786984
    r_hinge = 0.378
    station = 8
    blade_with_att = blade_len + r_attachment - r_hinge
    R = blade_len + r_attachment
    
    
    plt.figure()
    plt.subplot(311)
    i = 1
    
    #            plt.plot(eigVec[0,70*i:70*(i+1)], 'k--')
    pos = (np.hstack((0.0,np.linspace(3.39,42.39,40.0)))*blade_with_att/42.39 + 0.378)
    IDs = np.hstack((70*(i+1), np.linspace(30+70*i,70*(i+1)-1,40,dtype=int)))
    
    print(pos, IDs)

    for x in range(eigv[:,IDs,:].shape[0]):
        for z in range(eigv[:,IDs,:].shape[2]):
            eigv[x,0,z] = np.interp(r_attachment, pos, eigv[x,IDs,z])
    
    pos = (np.hstack((0.0,np.linspace(2.39,42.39,41.0)))*blade_with_att/42.39 + 0.378)/R
    IDs = np.hstack((70*(i+1), 0, np.linspace(30+70*i,70*(i+1)-1,40,dtype=int)))
    
    
    print(pos)
    
    scale_vals = np.amax(eigv[:,IDs,:], axis = 1)
    scale_vals_min = np.amin(eigv[:,IDs,:], axis = 1)
    
    for index, x in np.ndenumerate(scale_vals):
        if np.absolute(scale_vals_min[index]) > x:
            scale_vals[index] = scale_vals_min[index]
    
    #scale_vals = np.ones(scale_vals.shape)
    
    plt.plot(pos,eigv[0,IDs,station]/scale_vals[0,station], color='red', marker='s', markersize=1.5, label='1. lead-lag')
    plt.plot(pos,eigv[1,IDs,station]/scale_vals[1,station], 'k-.^', label='2.mode: lead-lag')
    plt.plot(pos,eigv[2,IDs,station]/scale_vals[2,station], 'k:',   label='3.mode: lead-lag')
    plt.plot(pos,eigv[3,IDs,station]/scale_vals[3,station],  color='blue', marker='o', markersize=1.5, label='3. flap')
    plt.plot(pos,eigv[4,IDs,station]/scale_vals[4,station],  color='red', marker='^', markersize=1.5, label='2. lead-lag')
    plt.plot(pos,eigv[5,IDs,station]/scale_vals[5,station], 'k--+', label='6.mode: lead-lag')
    plt.plot(pos,eigv[6,IDs,station]/scale_vals[6,station], 'k-',   label='7.mode: lead-lag')
    
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
   


def sim_plot(fq_dym, RPM_vec, target_dir):

    #plt.figure(figsize=(12,8))
    plt.figure(figsize=(12,8))
    array_per = RPM_vec
    plt.plot(array_per, fq_dym[:,1],'k*')
    plt.plot(array_per, fq_dym[:,2],'k*')
    plt.plot(array_per, fq_dym[:,3],'k*')
    plt.plot(array_per, fq_dym[:,6],'k*')
    plt.plot(array_per, fq_dym[:,7],'k*')
    plt.plot(array_per, fq_dym[:,10],'k*')
    plt.plot(array_per, fq_dym[:,0],'r*')
    plt.plot(array_per, fq_dym[:,4],'r*')
    plt.plot(array_per, fq_dym[:,5],'g*')
    plt.plot(array_per, fq_dym[:,9],'k*')

    plt.xlabel(r'RPM / RPM_{ref}')
    plt.ylabel(r'Freq / RPM_{ref}')
    plt.legend()
    #tikz_save(target_dir+'F2.tikz',figurewidth='\\figurewidth', figureheight='\\figureheight')
    plt.savefig(target_dir+'F2.png')
    #plt.plot(np.array([0, 100]),np.array([0, 100]),'k-')
