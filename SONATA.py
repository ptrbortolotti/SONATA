# -*- coding: utf-8 -*-
"""
SONATA-CBM is a Preprocessor for Parametric Composite Rotor Blade 
Cross-Sections.
The definition of the rotor blade topology is deliberately associated to the 
production of composite rotor blades. Thus, manufacturability is inherent from 
the geometric layup definition. Using orthogonal projection with corner-style 
differentiation the cross-section is discretized and can processed by the 
Variational Asymptotic Beam Sectional Analysis (VABS) afterwards. 

Date: 30/10/2018
@author: T.Pflumm, W.Garre
"""
import numpy as np
from SONATA.cbm.fileIO.configuration import Configuration
from SONATA.cbm.sonata_cbm import CBM

filename = 'jobs/VariSpeed/uh60a_cbm_advanced/sec_config.yml'
# filename = 'examples/sec_config_windturbine.yml'
# filename = 'examples/sec_config_windturbine_root.yml'
config = Configuration(filename)

job = CBM(config)
job.cbm_gen_topo()
job.cbm_gen_mesh(split_quads=True)

job.cbm_review_mesh()

job.cbm_run_vabs()

job.cbm_post_2dmesh(title='NoTitle')
#job.cbm_post_3dtopo()