#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 13:02:57 2018

@author: gu32kij
"""
import os
import platform
import subprocess
import shutil

if __name__ == '__main__':   #Make Sure to be in the SONATA working directory!
    os.chdir('../..')  #print(os.getcwd())
    
from SONATA.cbm.display.display_mesh import plot_cells
from SONATA.cbm.mesh.mesh_utils import sort_and_reassignID
from SONATA.cbm.fileIO.read_yaml_input import read_yaml_materialdb
from SONATA.cbm.fileIO.nastran_utils import read_nastran_bulkdata

from SONATA.vabs.VABS_interface import VABS_config, export_cells_for_VABS, XSectionalProperties
from SONATA.vabs.strain import Strain
from SONATA.vabs.stress import Stress


MaterialLst = read_yaml_materialdb('jobs/JBlaut/mat_db.yml')
filename = 'jobs/JBlaut/Querschnitt.nas'
nodes,mesh = read_nastran_bulkdata(filename)

#CONVERT MATERIAL DATAENTRIES: 
#Mat 2: UD-Roving,  -
#Mat 3: +-45 Gewebe ->Y Mat1: theta3 = +45
#Mat 11: Schaum/Gleitlage 
#Mat 10: +-0/90 Gewebe -> Mat1, 
            
for c in mesh:
    if c.MatID == 3:
        c.MatID = 1
        c.theta_3 = 45

    elif c.MatID == 10:
        c.MatID = 1
        
    elif c.MatID == 11:
        c.MatID = 3   

#RUN VABS
mesh,nodes = sort_and_reassignID(mesh)
vabs_filename = filename.replace('.nas', '.vab')
print('STATUS:\t Running VABS for Constitutive modeling:')

#Copy licensefile to 
src = os.getcwd()+'/SONATA/vabs/licenses/license.'+platform.node().lower()
dst = os.getcwd()+'/license'
if os.path.isfile(src):
    shutil.copyfile(src, dst)
else:
    print('no license file found at',src)

if platform.system() == 'Linux':
    executable = 'SONATA/vabs/bin/VABSIII'
elif platform.system == 'Windows':
    executable = 'SONATA/vabs/bin/VABSIII.exe'
    
vabs_cfg = VABS_config()
vabs_cfg.recover_flag = 1
vabs_cfg.M = [2000e3,0,0]

if vabs_cfg.recover_flag == 1:
    vabs_cfg.recover_flag=0
    export_cells_for_VABS(mesh, nodes, vabs_filename, vabs_cfg, MaterialLst)
    
    stdout = subprocess.run([executable,vabs_filename], stdout=subprocess.PIPE).stdout.decode('utf-8')
    vabs_cfg.recover_flag=1
    export_cells_for_VABS(mesh,nodes,vabs_filename, vabs_cfg, MaterialLst)
    print('STATUS:\t Running VABS for 3D Recovery:')
    stdout = subprocess.run([executable,vabs_filename], stdout=subprocess.PIPE).stdout.decode('utf-8')
else:
    export_cells_for_VABS(mesh, nodes ,vabs_filename, vabs_cfg, MaterialLst)
    stdout = subprocess.run([executable,vabs_filename], stdout=subprocess.PIPE).stdout.decode('utf-8')

stdout = stdout.replace('\r\n\r\n','\n\t   -')
stdout = stdout.replace('\r\n','\n\t   -')
stdout = stdout.replace('\n\n','\n\t   -')
stdout = stdout[:-2]
        
if ' VABS finished successfully' in stdout:
    stdout = 'STATUS:\t VABS Calculations Completed: \n\t   -' + stdout
else:
    stdout = 'ERROR:\t VABS Calculations Incomplete!: \n\t   -' + stdout
print(stdout) 

#VABS Postprocessing:
BeamProperties = XSectionalProperties(vabs_filename+'.K')
if vabs_cfg.recover_flag == 1:
    BeamProperties.read_all_VABS_Results()
    #ASSIGN Stress and strains to elements:
    for i,c in enumerate(mesh):
        c.strain = Strain(BeamProperties.ELE[i][1:7])
        c.stress = Stress(BeamProperties.ELE[i][7:13])
        c.strainM = Strain(BeamProperties.ELE[i][13:19])
        c.stressM = Stress(BeamProperties.ELE[i][19:25])
    
    #ASSIGN Displacement U to nodes:
    for i,n in enumerate(nodes):
        n.displacement = BeamProperties.U[i][3:6]

#PLOT!
plot_cells(mesh,nodes,'MatID')
plot_cells(mesh,nodes,'theta_3', BeamProperties, title='None', plotTheta11=False, plotDisplacement=False)
plot_cells(mesh,nodes,'stress.sigma12', BeamProperties, title='None', plotTheta11=False, plotDisplacement=False)
