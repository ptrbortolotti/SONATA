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

from SONATA.Pymore.marc.marc import MARC
import SONATA.Pymore.utl.optimization as obj
import SONATA.Pymore.utl.coef as coef

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
    dir_root = 'SONATA/Pymore/dym/mdl/04_rotortestrig/'
#    dir_root = '/work/WGarre/02_Projects/01_pymore/11_dymore/00_mdl/03_rotormodel/03_UH60_rotor_hover/05_UH60_rotor_fourblade_3Dinfl_baseline/'

    job = MARC(dir_root, 'Assembly.dym')
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
    job.fanplot(RPM_vec, result_dir)
    
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
