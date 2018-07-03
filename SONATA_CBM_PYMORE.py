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

from SONATA.Pymore.marc.marc import MARC

import numpy as np
import SONATA.Pymore.utl.coef as coef


plt.close('all')    
#TODO: Comment the CBM Class and memeber functions properly!
#TODO: include optionflags and Vabs_setup in Configuration

#filename = 'examples/sec_config.yml'
filename = 'jobs/VariSpeed/advanced/sec_config.yml'
config = Configuration(filename)
config.setup['BalanceWeight'] = False

cs = CBM(config)
#job.cbm_save()
cs.cbm_gen_topo()
#job.cbm_load_topo()
#job.cbm_display_config()
cs.cbm_gen_mesh()
#job.cbm_review_mesh()
#job.cbm_post_3dtopo()
#job.cbm_save()
cs.cbm_post_2dmesh()
cs.cbm_run_vabs()
beamProp = np.repeat([cs.cbm_set_DymoreMK()], 2, axis=0)
beamProp[0,-1] = +0.000e+00
beamProp[1,-1] = +7.361e+00
print(beamProp)
#print(beamProp.shape)

mdl_root = 'SONATA/Pymore/dym/mdl/03_rotormodel/05_UH60_rotor_optimization/01_UH60_rotor_snglblade_static/'
mdl = MARC(mdl_root, 'rotor_assembly.dym')
    
nbOfEig = mdl.analysis.sta_get_eigNb()
nbOfNod = mdl.analysis.sta_get_nodNb()
nbOfLoc = 15
RPM_vec = np.linspace(4.3*2*np.pi*0.7, 4.3*2*np.pi*1.1, nbOfLoc)

#beamProp = coef.refBeamProp()
#print(beamProp.shape)
mdl.marc_set_beamProp('BLADE_BP_CG01', beamProp)
result_dir ='SONATA/Pymore/rlt/'
    
#marc.marc_set_beamProp(cs.cbm_set_DymoreMK(x_offset = 0.81786984))
mdl.fanplot(RPM_vec, result_dir)
mdl.fanplot_show(RPM_vec, result_dir)






