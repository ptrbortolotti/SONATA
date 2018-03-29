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

from SONATA.cbm.fileIO.configuration import Configuration
from SONATA.cbm.sonata_cbm import CBM


plt.close('all')    
#TODO: Comment the CBM Class and memeber functions properly!
#TODO: include optionflags and Vabs_setup in Configuration

filename = 'examples/sec_config.yml'
filename = 'jobs/VariSpeed/advanced/sec_config.yml'
config = Configuration(filename)

job = CBM(config)
#job.cbm_save()
job.cbm_gen_topo()
#job.cbm_load_topo()
#job.cbm_display_config()
job.cbm_gen_mesh()
#job.cbm_save()
job.cbm_review_mesh()
#job.cbm_post_2dmesh()
job.cbm_run_vabs()


