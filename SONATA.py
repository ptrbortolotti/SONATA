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

Please Use the following Docstring styleguide:
https://numpydoc.readthedocs.io/en/latest/format.html

"""
import numpy as np
from SONATA.cbm.fileIO.configuration import Configuration
from SONATA.cbm.sonata_cbm import CBM

fname = 'jobs/VariSpeed/uh60a_cbm_advanced/sec_config.yml'
#fname = 'jobs/VariSpeed/uh60a_cbm_simple/sec_config.yml'
#fname = 'jobs/AREA/R250/sec_config.yml'
#fname = 'jobs/PBortolotti/sec_config.yml'
config = Configuration(fname)

job = CBM(config)

job.cbm_gen_topo()
job.cbm_gen_mesh(split_quads=True)

job.cbm_review_mesh()
job.cbm_post_2dmesh(title='Hello World!', plotTheta11=True)

job.cbm_run_vabs()
job.cbm_run_anbax()

#job.cbm_post_3dtopo()