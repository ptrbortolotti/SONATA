import numpy as np

from SONATA.Pymore.marc.marc import MARC
import SONATA.Pymore.utl.optimization as obj
import SONATA.Pymore.utl.coef as coef
from SONATA.Pymore.utl.plot import plot_fandiagram

if __name__ == "__main__":
    #-----------Simulation Setting---------------------------------
    recalc = True
    result_dir = 'SONATA/Pymore/rlt/'
    Omega = 4.3*2*np.pi  # [rad/sec]
    
    if recalc:
        #-----------Create Marc Instance---------------------------------
    #    dir_root = 'dym/mdl/03_rotormodel/05_UH60_rotor_optimization/03_UH60_rotor_snglblade_static_locBladeMass/'
    #    dir_root = '/work/WGarre/02_Projects/01_pymore/11_dymore/00_mdl/03_rotormodel/03_UH60_rotor_hover/05_UH60_rotor_fourblade_3Dinfl_baseline/'
        dir_root = 'SONATA/Pymore/dym/mdl/03_rotormodel/05_UH60_rotor_optimization/01_UH60_rotor_snglblade_static/'
        job = MARC(dir_root, 'rotor_assembly.dym')

        #-----------prepare permutations---------------------------------
        RPM_vec = np.linspace(0.2*Omega, 1.2*Omega, 11)

        # in case of more than one blade but only one table of beam properties all blades are updated simultaneously! 
        beamProp = coef.refBeamProp()
        job.marc_set_beamProp('BLADE_BP_CG01', beamProp)

        #-----------Calculate Fanplot---------------------------------    
        for i, rpm in enumerate(RPM_vec):
            job.analysis.sta_iterate(rpm)
            eigFreq, eigVec = job.analysis.sta_get_eigenSol()

            if i == 0:
                freq = eigFreq/(2*np.pi) #in [Hz]
                eigv = eigVec
            
            else:
                freq = np.vstack((freq, eigFreq/(2*np.pi)))   
                eigv = np.concatenate((eigv, eigVec), axis=2)
    
        job.analysis.freq = freq
        job.analysis.eigv = eigv

        #-----------Plot--------------------------------
        ref_fname = 'jobs/VariSpeed/uh60a_data_blade/fanplot_uh60a_bowen-davies-PhD.txt'
        ref_str = ['l1','f1','f2','f3','l2','t1','f4']
        plot_fandiagram(job.analysis.freq, Omega, RPM_vec, ref_fname=ref_fname)
        
        #-----------Save Results--------------------------------
        np.save(result_dir+'freq', freq, allow_pickle=True, fix_imports=True)
        np.save(result_dir+'RPM', RPM_vec, allow_pickle=True, fix_imports=True)

    else:
        #-----------Load Results----------------------                   
        freq = np.load(result_dir+'freq.npy', mmap_mode=None, allow_pickle=True, fix_imports=True, encoding='ASCII')
        RPM_vec = np.load(result_dir+'RPM.npy', mmap_mode=None, allow_pickle=True, fix_imports=True, encoding='ASCII')
        ref_fname = 'jobs/VariSpeed/uh60a_data_blade/fanplot_uh60a_bowen-davies-PhD.txt'
        plot_fandiagram(freq, Omega, RPM_vec, ref_fname=ref_fname)