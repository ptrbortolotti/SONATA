# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 10:42:32 2018

@author: TPflumm
"""

from concurrent import futures
from time import sleep, time
import copy    
import os
import matplotlib.pyplot as plt

print(os.getcwd())
os.chdir('C://TPflumm_local/work/SONATA')
from SONATA.cbm.fileIO.configuration import Configuration
from SONATA.cbm.sonata_cbm import CBM
from SONATA.cbm.fileIO.hiddenprints import HiddenPrints

def compute_cbm(config):
    job = CBM(config)
    job.cbm_gen_topo()
    job.cbm_gen_mesh()
    job.cbm_run_vabs()
    return job

plt.close('all')    
    
if __name__ == "__main__":
    __spec__ = None
    
    filename = 'examples/sec_config.yml'
    config = Configuration(filename)
    
    configs = []
    for i in range(1,10):
        print(i)
        tmp = copy.deepcopy(config)
        tmp.bw['Material'] = i
        configs.append(tmp)
    
    jobs = []
    with futures.ProcessPoolExecutor(max_workers=8) as e:
        fs = {e.submit(compute_cbm, c): it for it,c in enumerate(configs)}
        print("Alle Aufgaben gestartet.")
        for f in futures.as_completed(fs):
            print("n=%i, MpuS=%f" % (fs[f], f.result().BeamProperties.MpUS))
            jobs.append(f.result())
