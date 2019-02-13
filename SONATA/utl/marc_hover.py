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
import SONATA.Pymore.utl.read as read

#%%===================================#
#   Simulation Settings               #
#=====================================#
recalc = True
result_dir = 'SONATA/Pymore/rlt/'
identifier = 'hover'
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
    col = np.array([0.02,  0.015, 0.025, 3])
#    col = np.array([0.02,  0.005, 0.04, 15])
#    rpm = np.array([100.0, 70.0,  110.0, 4.3*2*np.pi, 9])
    rpm = np.array([100.0, 100.0, 100.0, 4.3*2*np.pi, 2])

    #rpm_ref = 4.3*2*np.pi 
    
    
#    tim = np.array([0.0,   3.0,   3.0,   1.5,         0.001])
    tim = np.array([0.0,   1.0,   1.0,   0.5,         0.001])
    psi = np.array([172.0+tim[4]*rpm[3]])
    
    Col, RPM, Psi, Sty = slope.rpm_col(col, rpm, tim, psi)
    #Sty: 
    #-------------------------------------#
    #   Create MARC instance              #
    #-------------------------------------#
#    job = MARC()  
    dir_root = 'SONATA/Pymore/dym/mdl/03_rotormodel/03_UH60_rotor_hover/05_UH60_rotor_fourblade_3Dinfl_baseline/'
    job = MARC(dir_root, 'rotor_assembly.dym')        
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
            
        #print(Omega, tim[4], Psi[i], Col[i])
        job.analysis.dyn_timeStep(Omega, tim[4], Psi[i], Col[i], 0.0, 0.0, np.zeros(3))
        T, Q, P, f, x, v = job.analysis.dyn_get_hubloads(job)
        print(T)
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
    
    M_hub = np.cross(f_res[:,0:3], x_res[:,0:3]) - f_res[:,3:6]
    Q = N_b * M_hub[:,2]
    P = N_b * np.einsum('ij,ij->i', -f_res[:,0:3], v_res[:,0:3]) + N_b * np.einsum('ij,ij->i', -f_res[:,3:6], v_res[:,3:6])
    T = N_b * f_res[:,2]  
    

#%%===================================#
#   Load Reference                    #
#=====================================#
val_filename = '/work/WGarre/02_Projects/03_varispeed/01_UH60_rotor_validation/03_UH60_rotor_validation_performance/CPs_CTs_UH60_WT.csv'
polar_hover  = read.Polar_Hover(val_filename)

CT_ref = np.asarray(polar_hover.CT)
CP_ref = np.asarray(polar_hover.CP)
    
    
#%%===================================#
#   Plot Results                      #
#=====================================#    
startID = 40
endID = 39

plt.plot(u_res[:24000,1], T[:24000])

Col = u_res[:,1]
RPM = u_res[:,3]

Sty = np.arange(499,np.size(Col,0),4500)

ID = np.append(np.arange(4500*41+2999,np.size(Col,0)-1*4500,4500),np.arange(1499,4500*39,4500))
#ID = np.append(np.arange(startID,np.size(Sty)-1,1,dtype=int), np.arange(0,endID,1,dtype=int))
plt.figure(figsize=(16,12))
plt.plot(Col)
plt.plot(ID,Col[ID], '*r')
plt.show()

plt.figure(figsize=(16,12))
plt.plot(RPM)
plt.plot(ID,RPM[ID], '*r')
plt.show()
#ID = np.delete(ID, np.arange(9, np.size(ID)-1, 10))
ID = np.concatenate((ID[0:15],np.flip(ID[15:30],0),ID[30:45],np.flip(ID[45:60],0),ID[60:75],np.flip(ID[75:90],0),ID[90:105],np.flip(ID[105:120],0),ID[120:135]))

RPMj = RPM[ID]
RPMj = RPMj.reshape((9,15))
Colj = Col[ID]
Colj = Colj.reshape((9,15))
Tj = T[ID]
Tj = Tj.reshape((9,15))
Qj = Q[ID]
Qj = Qj.reshape((9,15))
#Qjref = np.tile(Qj[6,:],(9,1))
Pj = P[ID]
Pj = Pj.reshape((9,15))
#Pjref = np.tile(Pj[6,:],(9,1))

plt.figure(figsize=(16,12))
plt.plot(u_res[:,1])
plt.plot(Sty,u_res[Sty,1], 'r*')

plt.figure(figsize=(16,12))
plt.plot(u_res[:,3]/np.max(u_res[:,3]))
plt.plot(u_res[:,1]/np.max(u_res[:,1]))

plt.figure(figsize=(14,16))
#plt.plot(y[:,0]*T2CT, y[:,2]*P2CP,'b-.', label = 'pymore')
#plt.plot(T*T2CT, P*P2CP,'g--', label = 'pymore P')
#plt.plot(T*T2CT, Q*M2CP,'b--', label = 'pymore P')
#plt.plot(Tj*T2CT, Qj*M2CP,'r--', label = 'pymore Q')
#plt.plot(np.transpose(Tj[6,:])*T2CT, np.transpose(Qj[6,:])*M2CP,'r--', label = 'pymore Q')
plt.plot(np.transpose(Tj)*T2CT, np.transpose(Pj)*P2CP,'k--', label = 'pymore P') # lines of const rotor speed
Pp = np.copy(Pj)
for i in range(np.size(Tj[:,0])):
    Pp[i-1,:] = np.divide((Pj[-3,:]-np.interp(Tj[-3,:],Tj[i-1,:],Pj[i-1,:])),Pj[-3,:])*100
    for j in range(np.size(Tj[-3,:])):
        if Tj[i-1,j-1] < Tj[-3,0] or Tj[i-1,j-1] > Tj[-3,-1] or Tj[-3,j-1] > Tj[i-1,-1] or Tj[-3,j-1] < Tj[i-1,0]:
            Pp[i-1,j-1] = np.nan
#plt.plot(np.transpose(T[ID])*T2CT, np.transpose(Q[ID])*M2CP,'r*', label = 'pymore Q')
plt.plot(np.transpose(T[ID])*T2CT, np.transpose(P[ID])*P2CP,'r+', label = 'pymore P') # lines of const rotor speed

#plt.plot(Tj[:,1:8]*T2CT, Pj[:,1:8]*P2CP,'k--', label = 'pymore Q') # lines of const collective
#plt.plot(Tj[i,:]*T2CT, Pj[i,:]*P2CP,'k--', label = 'pymore Q')
#plt.plot(np.transpose(Tj)*T2CT, np.transpose(Pj)*P2CP,'k--', label = 'pymore Q')
plt.plot(CT_ref*sigma_ref,CP_ref*sigma_ref,'k*', label = 'hover wind tunnel test')
plt.legend()

plt.xlabel(r'C_T')
plt.ylabel(r'C_P')
plt.legend()
plt.grid()

#tikz_save(result_dir+'hover.tikz',figurewidth='\\figurewidth', figureheight='\\figureheight')
plt.savefig(result_dir+'hover.png')

plt.figure(figsize=(14,16))
plt.plot(np.transpose(Tj)*T2CTs, np.transpose(Pj),'k--', label = 'pymore P') # lines of const rotor speed
plt.grid()
#plt.figure(figsize=(14,6))
#plt.plot(RPMj, 'k-')
#plt.figure(figsize=(14,6))
#plt.plot(colj[:,1], 'k-')
plt.figure(figsize=(14,16))
plt.plot(np.transpose(Tj)*T2CTs, np.transpose(Pp),'k--', label = 'pymore Pp') # lines of const rotor speed
