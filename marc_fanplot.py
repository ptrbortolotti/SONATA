################################################################################
#                                                                              #
#                                                                              #
#                        Author:  Willem Garre                                 #
#                        File:    PymoreRun.py                                 #
#                        Date:    04.04.2017                                   #
#                        Version: V0.2                                         #
#                                                                              #
#                                                                              #
################################################################################
#%%===================================#
#   Import Required Modules           #
#=====================================#
import numpy as np
import os

os.chdir('../..')

from SONATA.Pymore.marc.marc import MARC
import SONATA.Pymore.utl.optimization as obj
import SONATA.Pymore.utl.coef as coef

def fanplot_show(job, RPM_vec, result_dir, **kw):
    #%% plots the cambell-diagramm or fan-plot
    '''    
    arguments:
        res -- results array with the rotational speed and the first two flatwise and 
        first two edgewise eigenfrequencies
        
    '''
    Omega = 4.3 * 2 * np.pi
    print("Omega hard coded for UH-60A!")
    res = self.analysis.freq
    
    plt.figure()
    #plt.subplot(121)
    plt.grid(True)
    
    #plot rotor-harmonics 
    x =  np.linspace(0, 1.2, 20)
    y =  x*Omega/(2*np.pi)
    for i in range(1,9):
        plt.plot(x,i*y,'k--')
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
    for i,d in enumerate(res[:len(ref_str)].T):
        print(i)
        s=ref_str[i]
        plt.plot(x,d,'b')
        if 'f' in s:
            colorhex = 'blue'
            plt.plot(x, d, 'o-', color=colorhex)
            string = r'%s flap' % (s[-1])
            plt.text(x[-1]+.01, d[-1], string, color=colorhex)
        elif 'l' in s:
            colorhex = 'red'
            print(colorhex)
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
        
        
# old fanplot routine below
# delete as soon as new routine above is tested!
#        if 'val_fname' in kw:
#            val_fname = kw['val_fname']
#        else:
#            val_fname = '../11_varispeed/01_UH60_rotor_validation/02_UH60_rotor_validation_fanplot/Fanplot_Bowen_Davies_Diss.csv'
#        #%%===================================#
#        #   Plot Results                      #
#        #=====================================#
#        freq = self.analysis.freq
#        eigv = self.analysis.eigv
#        
#        blade_len = 7.36082856
#        r_attachment = 0.81786984
#        r_hinge = 0.378
#        
#        blade_with_att = blade_len + r_attachment - r_hinge
#        
#        fan_plot(np.real(freq), val_fname, RPM_vec, result_dir)
#        #sim_plot(freq, RPM_vec)
#        
#        plt.figure(figsize=(16,18))
#        plt.subplot(311)
#        i = 1
#        #            plt.plot(eigVec[0,70*i:70*(i+1)], 'k--')
#        pos = np.hstack((0.0,np.linspace(3.39,42.39,40.0)))*blade_with_att/42.39 + 0.378
#        IDs = np.hstack((70*(i+1), np.linspace(30+70*i,70*(i+1)-1,40,dtype=int)))
#        
#        for x in range(eigv[:,IDs,:].shape[0]):
#            for z in range(eigv[:,IDs,:].shape[2]):
#                eigv[x,0,z] = np.interp(r_attachment, pos, eigv[x,IDs,z])
#
#        pos = np.hstack((0.0,np.linspace(2.39,42.39,41.0)))*blade_with_att/42.39 + 0.378
#        IDs = np.hstack((70*(i+1), 0, np.linspace(30+70*i,70*(i+1)-1,40,dtype=int)))
#
#        scale_vals = np.amax(eigv[:,IDs,:], axis = 1)
#        scale_vals_min = np.amin(eigv[:,IDs,:], axis = 1)
#        
#        for index, x in np.ndenumerate(scale_vals):
#            if np.absolute(scale_vals_min[index]) > x:
#                scale_vals[index] = scale_vals_min[index]
#        
##        scale_vals = np.ones(scale_vals.shape)
#        
#        plt.plot(pos,eigv[0,IDs,3]/scale_vals[0,3], 'k-.',  label='1.mode: lead-lag')
##        plt.plot(pos,eigv[1,IDs,3]/scale_vals[1,3], 'k-.*', label='2.mode: lead-lag')
##        plt.plot(pos,eigv[2,IDs,3]/scale_vals[2,3], 'k:',   label='3.mode: lead-lag')
#        plt.plot(pos,eigv[3,IDs,3]/scale_vals[3,3], 'k--',  label='4.mode: lead-lag')
#        plt.plot(pos,eigv[4,IDs,3]/scale_vals[4,3], 'k--*', label='5.mode: lead-lag')
##        plt.plot(pos,eigv[5,IDs,3]/scale_vals[5,3], 'k--+', label='6.mode: lead-lag')
##        plt.plot(pos,eigv[6,IDs,3]/scale_vals[6,3], 'k-',   label='7.mode: lead-lag')
#
#        plt.grid()
#        plt.legend()
#        plt.xlabel('radial station [m]')
#        
#        plt.subplot(312)
#        i = 2
#        #            plt.plot(eigVec[0,70*i:70*(i+1)], 'k--')
#        #            plt.plot(eigVec[0,30+70*i:70*(i+1)], 'r--', label='1.mode: flap')
#        pos = np.hstack((0.0,np.linspace(3.39,42.39,40.0)))*blade_with_att/42.39 + 0.378
#        IDs = np.hstack((70*(i+1), np.linspace(30+70*i,70*(i+1)-1,40,dtype=int)))
#        
#        for x in range(eigv[:,IDs,:].shape[0]):
#            for z in range(eigv[:,IDs,:].shape[2]):
#                eigv[x,0,z] = np.interp(r_attachment, pos, eigv[x,IDs,z])
#
#        pos = np.hstack((0.0,np.linspace(2.39,42.39,41.0)))*blade_with_att/42.39 + 0.378
#        IDs = np.hstack((70*(i+1), 0, np.linspace(30+70*i,70*(i+1)-1,40,dtype=int)))
#       
#        scale_vals = np.amax(eigv[:,IDs,:], axis = 1)
#        scale_vals_min = np.amin(eigv[:,IDs,:], axis = 1)
#        
#        for index, x in np.ndenumerate(scale_vals):
#            if np.absolute(scale_vals_min[index]) > x:
#                scale_vals[index] = scale_vals_min[index]
#        
##        scale_vals = np.ones(scale_vals.shape)
#
##        plt.plot(pos,eigv[0,IDs,3]/scale_vals[0,3], 'k-.',  label='1.mode: flap')       
#        plt.plot(pos,eigv[1,IDs,3]/scale_vals[1,3], 'k-.*', label='2.mode: flap')
#        plt.plot(pos,eigv[2,IDs,3]/scale_vals[2,3], 'k:',   label='3.mode: flap')
#        plt.plot(pos,eigv[3,IDs,3]/scale_vals[3,3], 'k--',  label='4.mode: flap')
#        plt.plot(pos,eigv[4,IDs,3]/scale_vals[4,3], 'k--*', label='5.mode: flap')
#        plt.plot(pos,eigv[5,IDs,3]/scale_vals[4,3], 'k--+', label='6.mode: flap')
#        plt.plot(pos,eigv[6,IDs,3]/scale_vals[6,3], 'k-',   label='7.mode: flap')
#        plt.grid()
#        plt.legend()
#        plt.xlabel('radial station [m]')
#        
#        plt.subplot(313)
#        i = 3
#        pos = np.hstack((0.0,np.linspace(3.39,42.39,40.0)))*blade_with_att/42.39 + 0.378
#        IDs = np.hstack((70*(i+1), np.linspace(30+70*i,70*(i+1)-1,40,dtype=int)))
#        
#        for x in range(eigv[:,IDs,:].shape[0]):
#            for z in range(eigv[:,IDs,:].shape[2]):
#                eigv[x,0,z] = 0.0
#
#        pos = np.hstack((0.0,np.linspace(2.39,42.39,41.0)))*blade_with_att/42.39 + 0.378
#        IDs = np.hstack((70*(i+1), 0, np.linspace(30+70*i,70*(i+1)-1,40,dtype=int)))
#        
#        scale_vals = np.amax(eigv[:,IDs,:], axis = 1)
#        scale_vals_min = np.amin(eigv[:,IDs,:], axis = 1)
#        
#        for index, x in np.ndenumerate(scale_vals):
#            if np.absolute(scale_vals_min[index]) > x:
#                scale_vals[index] = scale_vals_min[index]
#
##        scale_vals = np.ones(scale_vals.shape)
#
##        plt.plot(pos,eigv[0,IDs,3]/scale_vals[0,3], 'k-.',  label='1.mode: torsion')
##        plt.plot(pos,eigv[1,IDs,3]/scale_vals[1,3], 'k-.*', label='2.mode: torsion')
##        plt.plot(pos,eigv[2,IDs,3]/scale_vals[2,3], 'k:',   label='3.mode: torsion')
##        plt.plot(pos,eigv[3,IDs,3]/scale_vals[3,3], 'k--',  label='4.mode: torsion')
#        plt.plot(pos,eigv[4,IDs,3]/scale_vals[4,3], 'k--*', label='5.mode: torsion')
#        plt.plot(pos,eigv[5,IDs,3]/scale_vals[5,3], 'k--+', label='6.mode: torsion')
##        plt.plot(pos,eigv[6,IDs,3]/scale_vals[6,3], 'k-',   label='7.mode: torsion')
#
#        plt.grid()
#        plt.legend()
#        plt.xlabel('radial station [m]')
#        #            plt.plot(eigVec[4,70*i:70*(i+1)], 'r-.')
#        #            plt.plot(eigVec[5,70*i:70*(i+1)], 'b-.')
#        
##        i = 2
##        k = 6
##        df = eigv[k,30+70*i:70*(i+1),:]-eigv[k,29+70*i:70*(i+1)-1,:]
##        
##        dxq = -np.divide(eigv[k,29+70*i:70*(i+1)-1,:], df)
##        ds  = np.sign(eigv[k,30+70*i:70*(i+1),:])-np.sign(eigv[k,29+70*i:70*(i+1)-1,:])
##        print(ds)
#        
##        kkk = np.linspace(0,df.shape[0]-1,df.shape[0])
##        x1 = np.tile(kkk[:,None],[1,df.shape[1]])
##        
##        xq = x1 + dxq
##        
##        xq = xq / 39.0 * 7.36082856
##        
##        xq = xq[np.where(ds!=0)]
#        #print(np.where(ds!=0))
##        xq = np.reshape(xq,(int(xq.size/df.shape[1]),df.shape[1]))
##        print(xq)
#        
##        plt.figure(figsize=(6,8))
#        #print( np.amin(xq[lb:ub,:] , axis=0))
##        print(np.mean(xq , axis=1))
##        plt.plot(xq, '*')

if __name__ == "__main__":
#%%===================================#
#   Simulation Settings               #
#=====================================#
recalc = True
result_dir = 'SONATA/Pymore/rlt/'


if recalc:
    #%%===================================#
    #   Prepare Simulations               #
    #=====================================#    
    #-------------------------------------#
    #   Create MARC instance              #
    #-------------------------------------#f
#    dir_root = 'dym/mdl/03_rotormodel/05_UH60_rotor_optimization/03_UH60_rotor_snglblade_static_locBladeMass/'
    dir_root = 'SONATA/Pymore/dym/mdl/03_rotormodel/05_UH60_rotor_optimization/01_UH60_rotor_snglblade_static/'
#    dir_root = '/work/WGarre/02_Projects/01_pymore/11_dymore/00_mdl/03_rotormodel/03_UH60_rotor_hover/05_UH60_rotor_fourblade_3Dinfl_baseline/'

    job = MARC(dir_root, 'rotor_assembly.dym')
#    job = MARC()
    
    #-------------------------------------#
    #   prepare permutations              #
    #-------------------------------------#
    nbOfEig = job.analysis.sta_get_eigNb()
    nbOfNod = job.analysis.sta_get_nodNb()
    nbOfLoc = 5
    RPM_vec = np.linspace(4.3*2*np.pi*0.7, 4.3*2*np.pi*1.1, nbOfLoc)

    beamProp = coef.refBeamProp()
#    beamProp = coef.testBeamProp()
    # in case of more than one blade but only one table of beam properties
    # all blades are updated simultaneously! 
#    beamProp = coef.optBaseLineBeamProp()
#    k = 0.5
#    beamProp[0:12,20] = beamProp[0:12,20] * k
#    beamProp[12:24,20] = beamProp[12:24,20] * k
#    beamProp[24:36,20] = beamProp[24:36,20] * k
#    beamProp[36:,20] = beamProp[36:,20] * k
#    job.marc_set_beamProp('BLADE_BP_CD01', beamProp)
    job.marc_set_beamProp('BLADE_BP_CG01', beamProp)

    #=====================================#
    #   Calculate Fanplot                 #
    #=====================================#   
    #-------------------------------------#
    #   initialize storage                #
    #-------------------------------------#
    nbOfEig = job.analysis.sta_get_eigNb()
    nbOfNod = job.analysis.sta_get_nodNb()
    freq = np.empty([nbOfEig])
    eigv = np.expand_dims(np.ones([nbOfEig, nbOfNod*6]), 2)

    #=====================================#
    #   Run Fanplot                       #
    #=====================================#
#   starttime = datetime.now()
    for RPM in RPM_vec:
        
        job.analysis.sta_eigenAnalysis(RPM)
        eigFreq, eigVec = job.analysis.sta_get_eigenSol()
        
        # store result
        freq = np.vstack((freq, eigFreq))
        eigv = np.concatenate((eigv, eigVec), axis=2)

    # scale and skip first entry    
    job.analysis.freq = freq[1:,:]/(2*np.pi)
    job.analysis.eigv = eigv[:,:,1:]

    #=====================================#
    #   Save Results                      #
    #=====================================#
    np.save(result_dir+'freq', freq, allow_pickle=True, fix_imports=True)
    np.save(result_dir+'RPM', RPM_vec, allow_pickle=True, fix_imports=True)
    
    job.fanplot_show(RPM_vec, result_dir)

else:
    #%%===================================#
    #   Load Results                      #
    #=====================================#
    freq = np.load(result_dir+'freq.npy', mmap_mode=None, allow_pickle=True, fix_imports=True, encoding='ASCII')
    RPM_vec = np.load(result_dir+'RPM.npy', mmap_mode=None, allow_pickle=True, fix_imports=True, encoding='ASCII')

    job.fanplot_show()

#%%===================================#
#   Process Results                   #
#=====================================#
#objFun = obj.gradPlacement(np.real(freq), RPM_vec)
