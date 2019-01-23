# -*- coding: utf-8 -*-
"""
Created on Tue Apr 04 14:49:03 2017

@author: TPflumm
"""
import os

if __name__ == '__main__':
    os.chdir('../../..')

from SONATA.cbm.fileIO.material import IsotropicMaterial, OrthotropicMaterial, AnisotropicMaterial

def clean_filestring(fname, comments='#', skiprows=0):
    """
    Removes comments, skippes rows and strips '\r\n' newline characters from filestring
    """                     
    b = [x.partition(comments)[0].rstrip()+'\n' for it, x in enumerate(open(fname, 'r').readlines()) if it>=skiprows]
    return ''.join([s for s in b if s.strip("\r\n").strip()])        


#def read_yaml_materialdb(fname):
#    b_string = clean_filestring(fname,comments='#')
#    mdb =  yaml.load(b_string)['Materials']
#    
#    MaterialLst=[]
#    for k,v in mdb.items():
#        v['ID'] = int(k.split()[-1])
#        if v['orth'] == 0:
#            MaterialLst.append(IsotropicMaterial(**v))
#        elif mdb[k]['orth'] == 1:
#            MaterialLst.append(OrthotropicMaterial(**v))
#        elif mdb[k]['orth'] == 2:
#            MaterialLst.append(AnisotropicMaterial(**v))
#        
#    return sorted(MaterialLst, key=lambda x: x.id)

        
if __name__ == '__main__':
    os.chdir('/media/gu32kij/work/TPflumm/SONATA')
    x = clean_filestring('examples/sec_config.yml') 
    #mdb = read_yaml_materialdb('examples/mat_db.yml')