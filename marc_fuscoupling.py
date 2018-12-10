################################################################################
#                                                                              #
#                                                                              #
#                        Author:  Willem Garre                                 #
#                        File:    uh60a_testbed.py                             #
#                        Date:    08.01.2018                                   #
#                        Version: V1.0                                         #
#                                                                              #
#                                                                              #
################################################################################
#%%===================================#
#   Import Required Modules           #
#=====================================#
import numpy as np
import matplotlib.pyplot as plt

from datetime import datetime

from SONATA.Pymore.marc.marc import MARC

import SONATA.Pymore.utl.coef as coef
import SONATA.Pymore.utl.slope as slope

from SONATA.Pymore.linSys import LinSys
from SONATA.Pymore.runge import IterCtrl, rk4

def plot_results():
    #%%===================================#
    #   Load Reference                    #
    #=====================================#
    vel = np.linspace(0.0, 75.0, T.shape[0])
    plt.plot(vel, T)
#    plt.plot(u_res[:24000,1], M_hub[:24000,0])
    
    plt.figure(figsize=(8,6))
    plt.plot(vel, M_hub[:,2])
    #plt.plot(u_res[:24000,1], M_hub1[:24000,1])
    #plt.plot(u_res[:24000,1], M_hub2[:24000,1])
    #plt.plot(u_res[:24000,1], M_hub3[:24000,1])
    #plt.plot(u_res[:24000,1], M_hub4[:24000,1])
    
#    plt.plot(vel, M_hub1[:,1])
#    plt.plot(vel, M_hub2[:,1])
#    plt.plot(vel, M_hub3[:,1])
#    plt.plot(vel, M_hub4[:,1])

if __name__ == "__main__":
    #%%===================================#
    #   Simulation Settings               #
    #=====================================#
    recalc = True
    result_dir = 'SONATA/Pymore/rlt/'
    identifier = 'fuscoupling'
    N_b = 4.0
    c_ref = 0.53
    R_ref = 8.179
    P2CPs, M2CPs, T2CTs, P2CP, M2CP, T2CT = coef.dim2coef(N_b, c_ref, R_ref, 4.3*np.pi*2, 1.225)
    sigma_ref = coef.sigma(N_b, c_ref, R_ref)

    rpm = 4.3*2*np.pi

    deltat = 0.001
    tmin = 0.0
    tmax = 2.0
    

    # create objects
    fuselage = LinSys('SONATA/Pymore/fuselage_studi.txt')
    process = IterCtrl(deltat, tmin, tmax)
    
    x_rk4 = np.zeros((fuselage.dof, 1))
    
    if recalc:
        #%%===================================#
        #   Prepare Simulations               #
        #=====================================#
        #-------------------------------------#
        #   prepare permutations              #
        #-------------------------------------#
     
        psi = np.array([7.492924700099042e+002])

        steps = (tmax-tmin) / deltat + 1
        
        PsiMin = np.interp(tmin, np.array([0.0, 1.0e+003]), np.array([7.492924700099042e+002, 2.7766989290882124e+004]))
        
        PsiMax = np.interp(tmax, np.array([0.0, 1.0e+003]), np.array([7.492924700099042e+002, 2.7766989290882124e+004]))
        
        Psi = np.linspace(PsiMin, PsiMax, steps)
        
        #-------------------------------------#
        #   Create MARC instance              #
        #-------------------------------------#
    #    dir_root = 'SONATA/Pymore/dym/mdl/03_rotormodel/05_UH60_rotor_optimization/02_UH60_rotor_fourblade_3Dinfl_testmdl/'
#        dir_root = 'SONATA/Pymore/dym/mdl/03_rotormodel/03_UH60_rotor_hover/05_UH60_rotor_fourblade_3Dinfl_baseline/'
        dir_root = 'SONATA/Pymore/dym/mdl/03_rotormodel/04_UH60_rotor_cruise/02_UH60_rotor_fourblade_3Dinfl_trim'
        job = MARC(dir_root, 'rotor_assembly.dym')
        
    #    job = MARC()
        #-------------------------------------#
        #   initialize storage                #
        #-------------------------------------#
        u_stack = np.empty([4])
        y_stack = np.empty([3])
        x_stack = v_stack = f_stack = np.empty([24])
        
        #%%===================================#
        #   Run Simulations                   #
        #=====================================#
        starttime = datetime.now()
    
        for i, val in enumerate(Psi):
            if i==0:
                Omega = (Psi[i]-psi[0]) / deltat
            else:
                Omega = (Psi[i]-Psi[i-1]) / deltat
            
            job.marc_set_tableData('SHAFT_TF_A03A04_PD',   val)
            job.analysis.dyn_timeStep(deltat)
#            T, Q, P, f, x, v = job.analysis.dyn_get_hubloads(job)
            
            shaftLoads = np.asarray(job.marc_get_sensorData('ASSEMBLY_SS_FOURBLADE_FORCES_SHAFT_INERTIAL', 6))
            
            loads1 = np.hstack((shaftLoads, np.zeros(18)))
            # TIME STEP OF FUSELAGE
            # calculate generalized forces and moments
            if i == 0: 
                loads0 = loads1
            
            Pi  = np.matmul(fuselage.ev, loads0)
            Pi1 = np.matmul(fuselage.ev, loads1)
        
            # do rk4 step
            rk4(fuselage, process, Pi, Pi1)
        
            # update load vector
            loads0 = loads1
        
            # recover real states
            x_rk4 = np.hstack((x_rk4, np.matmul(np.transpose(fuselage.ev), np.transpose([fuselage.xi[:fuselage.dim]]))))
            
            # write values to dymore
            job.marc_set_tableData('SHAFT_TF_A03A04', x_rk4[2,-1]*0.1)
            job.marc_set_tableData('SHAFT_TF_A06A07', x_rk4[5,-1]*0.1)
            job.marc_set_tableData('SHAFT_TF_A08A09', x_rk4[0,-1]*0.1)
            job.marc_set_tableData('SHAFT_TF_A10A11', x_rk4[1,-1]*0.1)
            job.marc_set_tableData('SHAFT_TF_A13A14', x_rk4[3,-1]*0.1)
            job.marc_set_tableData('SHAFT_TF_A15A16', x_rk4[4,-1]*0.1)
            
            # STORE DYMORE MODEL OUTPUT
#            u_stack = np.vstack((u_stack, np.array([tim[4],Col[i],Psi[i],Omega])))
#            y_stack = np.vstack((y_stack, np.array([T,Q,P])))
#            x_stack = np.vstack((x_stack, x))
#            v_stack = np.vstack((v_stack, v))
            f_stack = np.vstack((f_stack, loads1))       
            
    
        print('\n *** Time Dynamic Analysis: ' + str(datetime.now() - starttime))
            
        #%%===================================#
        #   Save Results                      #
        #=====================================#
#        np.save(result_dir+'u_'+identifier, u_stack[1:,:], allow_pickle=True, fix_imports=True)
#        np.save(result_dir+'y_'+identifier, y_stack[1:,:], allow_pickle=True, fix_imports=True)
#        np.save(result_dir+'x_'+identifier, x_stack[1:,:], allow_pickle=True, fix_imports=True)
#        np.save(result_dir+'v_'+identifier, v_stack[1:,:], allow_pickle=True, fix_imports=True)
        np.save(result_dir+'f_'+identifier, f_stack[1:,:], allow_pickle=True, fix_imports=True)
#        np.save(result_dir+'s_'+identifier, Sty  , allow_pickle=True, fix_imports=True)   
    
    else:
        #%%===================================#
        #   Load Results                      #
        #=====================================#
        u_res = np.load(result_dir+'u_'+identifier+'.npy', mmap_mode=None, allow_pickle=True, fix_imports=True, encoding='ASCII')
        y_res = np.load(result_dir+'y_'+identifier+'.npy', mmap_mode=None, allow_pickle=True, fix_imports=True, encoding='ASCII')
        x_res = np.load(result_dir+'x_'+identifier+'.npy', mmap_mode=None, allow_pickle=True, fix_imports=True, encoding='ASCII')
        v_res = np.load(result_dir+'v_'+identifier+'.npy', mmap_mode=None, allow_pickle=True, fix_imports=True, encoding='ASCII')
        f_res = np.load(result_dir+'f_'+identifier+'.npy', mmap_mode=None, allow_pickle=True, fix_imports=True, encoding='ASCII')
        Sty   = np.load(result_dir+'s_'+identifier+'.npy', mmap_mode=None, allow_pickle=True, fix_imports=True, encoding='ASCII')
        
        
        
        M_hub1 = np.cross(f_res[:,0:3], x_res[:,0:3]) - f_res[:,3:6]
        M_hub2 = np.cross(f_res[:,6:9], x_res[:,6:9]) - f_res[:,9:12]
        M_hub3 = np.cross(f_res[:,12:15], x_res[:,12:15]) - f_res[:,15:18]
        M_hub4 = np.cross(f_res[:,18:21], x_res[:,18:21]) - f_res[:,21:24]
        M_hub = M_hub1 + M_hub2 + M_hub3 + M_hub4
        Q = M_hub[:,2]
        P = N_b * np.einsum('ij,ij->i', -f_res[:,0:3], v_res[:,0:3]) + N_b * np.einsum('ij,ij->i', -f_res[:,3:6], v_res[:,3:6])
        T = f_res[:,2] + f_res[:,8] + f_res[:,14] + f_res[:,20]  
 
    #%%===================================#
    #   Plot Results                      #
    #=====================================#         
    plot_results()
