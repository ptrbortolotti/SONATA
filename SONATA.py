# -*- coding: utf-8 -*-
"""
This is currently the main excecution script for SONATA.

SONATA is a preprocessor for parametric analysis and design of composite beam 
cross-sections in a multidisciplinary rotor design environment. A helicopter 
rotor blade represents a classical aeroelastic problem, where the aerodynamic 
behavior, the structural elasticity and vibrational dynamics have to be studied 
simultaneously.  While a geometric definition of a rotorblade with CAD tools is 
simple, the transfer to a meshed cross-sectional representation may prohibit 
automated design optimization. Consequently, most researches have developed 
individual parametric mesh generators for the cross-sectional analysis, that 
reduces their structural model to few design variables in the process. 
SONATA represents such a preprocessor. SONATA is written in python and is using
for a lot of operations the Opencascade (CAD) kernel with its python wrapper 
(pythonocc).

Because it is currently still under extensive development the code isn't as 
clean as it could be but allows easy debugging, printing and 3D displaying. In 
the future, the SONATA execution script shall inlude an openmdao structure 
which can call the then unlying functionalities. 

Date: 01/02/2017
@author: TPflumm
"""

import matplotlib.pyplot as plt

from SONATA.fileIO.configuration import Configuration
from SONATA.fileIO.readinput import read_material_input
from SONATA.cbm import CBM

plt.close('all')    
filename = 'jobs/EPiet/AREA_R1000/sec_config.input'

config = Configuration(filename)
#config.SETUP_radial_station = 160
#config.flag_mesh_core = True
#config.VABS.recover_flag = 1
#config.VABS.M = [0,2000e3,0]

MaterialLst = read_material_input(config.SETUP_mat_filename)

#TODO: Comment the CBM Class and memeber functions properly!
#TODO: include optionflags and Vabs_setup in Configuration
#TODO: include yaml style input files

job1 = CBM(config,MaterialLst)
job1.cbm_gen_topo()
job1.cbm_display_config()
job1.cbm_gen_mesh()
job1.cbm_review_mesh()
#job1.cbm_run_vabs(filename)
job1.cbm_post_2dmesh()
#job1.cbm_post_3dtopo()
#job1.cbm_post_3dmesh()


#job2 = CBM(config,MaterialLst)
#job2.config.BW_XPos = job1.config.BW_XPos+5
#job2.config.SEG_Layup[0][4][2] = job1.config.SEG_Layup[0][4][2]+0.3
#job2.cbm_gen_topo()
#job2.cbm_gen_mesh()
##job2.cbm_post_3dmesh()
#job1.cbm_run_vabs(filename)
