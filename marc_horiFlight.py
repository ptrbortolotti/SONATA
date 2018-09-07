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
    recalc = False
    result_dir = 'SONATA/Pymore/rlt/'
    identifier = 'vertFlight'
    N_b = 4.0
    c_ref = 0.53
    R_ref = 8.179
    P2CPs, M2CPs, T2CTs, P2CP, M2CP, T2CT = coef.dim2coef(N_b, c_ref, R_ref, 4.3*np.pi*2, 1.225)
    sigma_ref = coef.sigma(N_b, c_ref, R_ref)
    
    if recalc:
        #%%===================================#
        #   Prepare Simulations               #
        #=====================================#
        #-------------------------------------#
        #   prepare permutations              #
        #-------------------------------------#
        #        0,    1,      2,      3,      4
        # col = [col0, colmin, colmax, colstp]
        # rpm = [rpm0, rpmmin, rpmmax, rpmref, rpmstp]
        # tim = [t0,   tcol,   trpm,   tstdy,  tstpsize]
        # psi = [psi0]
     
        col = np.array([0.02,  0.019, 0.021, 3])
    #    col = np.array([0.02,  0.015, 0.025, 3])
    #    col = np.array([0.02,  0.005, 0.04, 15])
    #    rpm = np.array([100.0, 70.0,  110.0, 4.3*2*np.pi, 9])
        rpm = np.array([100.0, 100.0, 100.0, 4.3*2*np.pi, 2])
    #    tim = np.array([0.0,   3.0,   3.0,   1.5,         0.001])
        tim = np.array([0.0,   1.0,   1.0,   0.5,         0.001])
        psi = np.array([172.0+tim[4]*rpm[3]])
        
        Col, RPM, Psi, Sty = slope.rpm_col(col, rpm, tim, psi)
        LatCyc = LngCyc = 0.0
        
        vx = np.linspace(0.0, -75.0, Col.shape[0])
        Vinf = np.vstack((vx, np.zeros((2, vx.shape[0]))))
        
        
        #-------------------------------------#
        #   Create MARC instance              #
        #-------------------------------------#
    #    dir_root = 'SONATA/Pymore/dym/mdl/03_rotormodel/05_UH60_rotor_optimization/02_UH60_rotor_fourblade_3Dinfl_testmdl/'
        dir_root = 'SONATA/Pymore/dym/mdl/03_rotormodel/03_UH60_rotor_hover/05_UH60_rotor_fourblade_3Dinfl_baseline/'
    
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
                Omega = (Psi[i]-psi[0]) / tim[4]
            else:
                Omega = (Psi[i]-Psi[i-1]) / tim[4]
                
            job.marc_time_step(Omega, tim[4], Psi[i], Col[i], LatCyc[i], LngCyc[i], Vinf[i])
            T, Q, P, f, x, v = job.marc_get_hubloads()
            
            # STORE DYMORE MODEL OUTPUT
            u_stack = np.vstack((u_stack, np.array([tim[4],Col[i],Psi[i],Omega])))
            y_stack = np.vstack((y_stack, np.array([T,Q,P])))
            x_stack = np.vstack((x_stack, x))
            v_stack = np.vstack((v_stack, v))
            f_stack = np.vstack((f_stack, f))       
            
    
        print('\n *** Time Dynamic Analysis: ' + str(datetime.now() - starttime))
            
        #%%===================================#
        #   Save Results                      #
        #=====================================#
        np.save(result_dir+'u_'+identifier, u_stack[1:,:], allow_pickle=True, fix_imports=True)
        np.save(result_dir+'y_'+identifier, y_stack[1:,:], allow_pickle=True, fix_imports=True)
        np.save(result_dir+'x_'+identifier, x_stack[1:,:], allow_pickle=True, fix_imports=True)
        np.save(result_dir+'v_'+identifier, v_stack[1:,:], allow_pickle=True, fix_imports=True)
        np.save(result_dir+'f_'+identifier, f_stack[1:,:], allow_pickle=True, fix_imports=True)
        np.save(result_dir+'s_'+identifier, Sty  , allow_pickle=True, fix_imports=True)   
    
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