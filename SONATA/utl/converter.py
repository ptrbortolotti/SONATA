#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 09:33:29 2019

@author: gu32kij
"""
import os
import numpy as np
from scipy.interpolate import PchipInterpolator
from collections import OrderedDict
import yaml
from jsonschema import validate
import matplotlib.pyplot as plt


if __name__ == '__main__':
    os.chdir('../..')
    
from SONATA.cbm.classCBMConfig import CBMConfig

def arc_length(x, y):
    """
    Small routine that for given x and y of a profile compute the arc length positions
    
    Parameters
    ----------
    x : np array
        x coordinate of an airfoil, 1-TE 0-LE
    
    y : np array
        y coordinate, positive suction side, negative pressure side
    
    
    Returns
    ----------
    arc : float
        arc position, which can normalized from 0 to 1
        
    """
    
    npts = len(x)
    arc = np.array([0.]*npts)
    for k in range(1, npts):
        arc[k] = arc[k-1] + np.sqrt((x[k] - x[k-1])**2 + (y[k] - y[k-1])**2)

    return arc

def iea37_converter(blade, cs_pos, byml, materials):
    """
    Converts the 
    @author: Pietro Bortolotti
    
    Parameters
    ----------
    blade : Blade
        from classBlade
    
    byml : dict
        yaml data of the blade
    
    materials : 
    
    Returns
    ----------
    cbmconfigs : np.array
        returns a numpy array with grid location and CBMconfig instances
        np.array([[grid, CBMconfig]])
    
    """
    
    # Segments and webs
    tmp0    = byml.get('internal_structure_2d_fem').get('webs')
    #x = blade.x
    x = cs_pos #np.asarray((byml.get('internal_structure_2d_fem').get('positions')))
    
    if tmp0 == None:
        n_webs = 0
    else:
        n_webs  = len(tmp0)
    
    tmp2    = [dict([('position', x[n])]) for n in range(len(x))]
    id_webs = [dict() for n in range(len(x))]

    for i in range(len(x)):
        n_webs_i = 0
        for j in range(n_webs):
            
            if x[i] >=tmp0[j]['start']['grid'][0] and x[i] <=tmp0[j]['start']['grid'][-1]:
                n_webs_i = n_webs_i + 1
                
                id_webs[i][tmp0[j]['name']]       = {}
                id_webs[i][tmp0[j]['name']]['id'] = 2 * j + 2
                
                set_interp = PchipInterpolator(tmp0[j]['start']['grid'], tmp0[j]['start']['values'])
                id_webs[i][tmp0[j]['name']]['start'] = set_interp(x[i])
                
                set_interp = PchipInterpolator(tmp0[j]['end']['grid'], tmp0[j]['end']['values'])
                id_webs[i][tmp0[j]['name']]['end'] = set_interp(x[i])
                
        if n_webs_i == 0:
            tmp2[i]['segments']                 = [{}]
            tmp2[i]['segments'][0]['id']        = 0 
            tmp2[i]['segments'][0]['filler']    = None
            tmp2[i]['segments'][0]['layup']     = [{}]
            
        else:
            tmp2[i]['segments'] = [dict([('id', n),('layup', [{}]),('filler', None)]) for n in range(2 * n_webs_i + 2)]
        
    # Get sections information and init the CBM instances.
    tmp1 = byml.get('internal_structure_2d_fem').get('layers')
    
    # Composite stacking sequences
    thick_web = np.zeros([len(x),n_webs])
    for i in range(len(x)):
        profile         = blade.airfoils[i,1].coordinates
        id_le           = np.argmin(profile[:,0])
        
        if np.mean(profile[0:id_le, 1]) < 0:
            profile = np.flip(profile,0)

        profile_curve   = arc_length(profile[:,0], profile[:,1]) / arc_length(profile[:,0], profile[:,1])[-1]
        
        id_layer = 0
        for idx_sec, sec in enumerate(tmp1):
            if x[i] >= sec['thickness']['grid'][0] and x[i] <= sec['thickness']['grid'][-1]:
                if 'web' not in sec.keys():
                    if idx_sec>0:
                        tmp2[i]['segments'][0]['layup'].append({})
                    tmp2[i]['segments'][0]['layup'][id_layer]['name']          = sec['material'] + '_' + str(x[i])
                    tmp2[i]['segments'][0]['layup'][id_layer]['material_name'] = sec['material']
                    if 'start' in sec.keys():
                        
                        if 'fixed' not in sec['start']:
                            set_interp = PchipInterpolator(sec['start']['grid'], sec['start']['values'])
                            tmp2[i]['segments'][0]['layup'][id_layer]['start']     = set_interp(x[i])
                        else:
                            for j in range(idx_sec):
                                if tmp1[j]['name'] == sec['start']['fixed']:
                                    set_interp = PchipInterpolator(tmp1[j]['end']['grid'], tmp1[j]['end']['values'])
                                    tmp2[i]['segments'][0]['layup'][id_layer]['start']     = set_interp(x[i])
                        
                        if 'fixed' not in sec['end']:
                            set_interp = PchipInterpolator(sec['end']['grid'], sec['end']['values'])
                            tmp2[i]['segments'][0]['layup'][id_layer]['end']       = set_interp(x[i])
                        else:
                            for j in range(idx_sec):
                                if tmp1[j]['name'] == sec['end']['fixed']:
                                    set_interp = PchipInterpolator(tmp1[j]['start']['grid'], tmp1[j]['start']['values'])
                                    tmp2[i]['segments'][0]['layup'][id_layer]['end']     = set_interp(x[i])
                    else:
                        tmp2[i]['segments'][0]['layup'][id_layer]['start']     = 0.
                        tmp2[i]['segments'][0]['layup'][id_layer]['end']       = 1.
                    
                    
                    set_interp = PchipInterpolator(sec['thickness']['grid'], sec['thickness']['values'])
                    tmp2[i]['segments'][0]['layup'][id_layer]['thickness']     = set_interp(x[i])
                    
                    if 'fiber_orientation' in sec.keys():
                        set_interp = PchipInterpolator(sec['fiber_orientation']['grid'], sec['fiber_orientation']['values'])
                        tmp2[i]['segments'][0]['layup'][id_layer]['orientation'] = set_interp(x[i])
                    else:
                        tmp2[i]['segments'][0]['layup'][id_layer]['orientation'] = 0.
                    id_layer = id_layer + 1
                else:
                    id_web          = id_webs[i][(sec['web'])]['id']
                    id_layer_web    = 0
                    for id_mat in range(1,len(materials)):
                        if sec['material'] == materials[id_mat].name:
                        
                            if isinstance(materials[id_mat].E, float): # Isotropic -> fill the web
                                tmp2[i]['segments'][id_web]['filler'] = materials[id_mat].name
                                
                                tmp2[i]['segments'][id_web]['layup'][0]['name']             = 'dummy'
                                tmp2[i]['segments'][id_web]['layup'][0]['material_name']    = sec['material']
                                tmp2[i]['segments'][id_web]['layup'][0]['thickness']        = 1e-3
                                tmp2[i]['segments'][id_web]['layup'][0]['start']            = 0.0
                                tmp2[i]['segments'][id_web]['layup'][0]['end']              = 1.0
                                tmp2[i]['segments'][id_web]['layup'][0]['orientation']      = 0.0
                                
                                set_interp = PchipInterpolator(sec['thickness']['grid'], sec['thickness']['values'])
                                thick_web[i,int(id_web/2-1)] = thick_web[i,int(id_web/2-1)] + set_interp(x[i])
                                                              
                            else:
                                tmp2[i]['segments'][id_web - 1]['layup'][id_layer_web]['name'] = 'web_' + str(int(id_web/2)) + '_' + sec['material']  + '_' + str(x[i])
                                tmp2[i]['segments'][id_web + 1]['layup'][id_layer_web]['name'] = 'web_' + str(int(id_web/2)) + '_' + sec['material']  + '_' + str(x[i])
                                
                                tmp2[i]['segments'][id_web - 1]['layup'][id_layer_web]['material_name'] = sec['material']
                                tmp2[i]['segments'][id_web + 1]['layup'][id_layer_web]['material_name'] = sec['material']
                                
                                set_interp = PchipInterpolator(sec['thickness']['grid'], sec['thickness']['values'])
                                tmp2[i]['segments'][id_web - 1]['layup'][id_layer_web]['thickness'] = set_interp(x[i]) * 0.5
                                tmp2[i]['segments'][id_web + 1]['layup'][id_layer_web]['thickness'] = set_interp(x[i]) * 0.5
                                
                                #flange = 1e-3
                                
                                tmp2[i]['segments'][id_web - 1]['layup'][id_layer_web]['start'] = 0.0#id_webs[i][(sec['web'])]['start'] + flange
                                tmp2[i]['segments'][id_web - 1]['layup'][id_layer_web]['end']   = 1.0 #id_webs[i][(sec['web'])]['end'] - flange
                                
                                tmp2[i]['segments'][id_web + 1]['layup'][id_layer_web]['start'] = 0.0 #id_webs[i][(sec['web'])]['start'] - flange
                                tmp2[i]['segments'][id_web + 1]['layup'][id_layer_web]['end']   = 1.0 #id_webs[i][(sec['web'])]['end'] + flange
                                
                                if 'fiber_orientation' in sec.keys():
                                    set_interp = PchipInterpolator(sec['fiber_orientation']['grid'], sec['fiber_orientation']['values'])
                                    tmp2[i]['segments'][id_web - 1]['layup'][id_layer_web]['orientation'] = set_interp(x[i])
                                    tmp2[i]['segments'][id_web + 1]['layup'][id_layer_web]['orientation'] = set_interp(x[i])
                                else:
                                    tmp2[i]['segments'][id_web - 1]['layup'][id_layer_web]['orientation'] = 0.
                                    tmp2[i]['segments'][id_web + 1]['layup'][id_layer_web]['orientation'] = 0.
                                
                                id_layer_web = id_layer_web + 1
                            

    webs    = [OrderedDict() for n in range(len(x))]
    for i in range(len(x)):    
        for i_web,web in enumerate(tmp0):                
            if x[i] >= web['start']['grid'][0] and x[i] <= web['start']['grid'][-1]:
                
                set_interp      = PchipInterpolator(web['start']['grid'], web['start']['values'])
                start           = set_interp(x[i])
                
                set_interp      = PchipInterpolator(web['end']['grid'], web['end']['values'])
                end             = set_interp(x[i])
                
                profile         = blade.airfoils[i,1].coordinates
                id_le           = np.argmin(profile[:,0])
                
                if np.mean(profile[0:id_le, 1]) < 0:
                    profile     = np.flip(profile,0)

                profile_curve   = arc_length(profile[:,0], profile[:,1]) / arc_length(profile[:,0], profile[:,1])[-1]
                
                set_interp      = PchipInterpolator(profile_curve, profile[:,0])
                x_web_start     = set_interp(start)
                x_web_end       = set_interp(end)
                
                x_web_start_le  = x_web_start - thick_web[i, i_web] / blade.chord[i , 1]
                x_web_start_te  = x_web_start + thick_web[i, i_web] / blade.chord[i , 1]
                
                x_web_end_le    = x_web_end - thick_web[i, i_web] / blade.chord[i , 1]
                x_web_end_te    = x_web_end + thick_web[i, i_web] / blade.chord[i , 1]
                
                # Correction to avoid errors in the interpolation at the edges
                if min(np.diff(np.flip(profile[0:id_le,0]))) < 0:
                    offset_edge = np.argmin(abs(np.diff(profile[0:id_le,0]))) + 1
                else:
                    offset_edge = 1
                    
                set_interp      = PchipInterpolator(np.flip(profile[offset_edge:id_le - offset_edge,0]), np.flip(profile_curve[offset_edge:id_le - offset_edge]))
                web_start_le , web_start_te    = set_interp([x_web_start_le , x_web_start_te])
                
                set_interp      = PchipInterpolator(profile[id_le + offset_edge:-offset_edge,0], profile_curve[id_le + offset_edge:-offset_edge])
                web_end_le , web_end_te        = set_interp([x_web_end_le , x_web_end_te])
                
                w_f = {}
                w_f['id']   = 2 * i_web + 1
                w_f['position'] = [web_start_le, web_end_le]
                w_f['curvature'] = None
                webs[i][2 * i_web + 1] = w_f
                
                w_r = {}
                w_r['id']   = 2 * i_web + 2
                w_r['position'] = [web_start_te, web_end_te]
                w_r['curvature'] = None
                webs[i][2 * i_web + 2] = w_r
                
                
    lst = []
    for (seg,w,x) in zip(tmp2, webs, x):
        tmp = {'webs':list(w.values()), 'segments':seg['segments'], 'position':x, 'mesh_resolution':240}
        lst.append([x, CBMConfig(tmp, materials, iea37=True)])
    
    return np.asarray(lst)
    


#%% MAIN
if __name__ == '__main__':
    from SONATA.classAirfoil import Airfoil
    from SONATA.classBlade import Blade 
    from SONATA.classMaterial import read_IEA37_materials

    with open('./jobs/PBortolotti/IEAonshoreWT.yaml', 'r') as myfile:
        inputs  = myfile.read()
    with open('jobs/PBortolotti/IEAontology_schema.yaml', 'r') as myfile:
        schema  = myfile.read()
    validate(yaml.load(inputs), yaml.load(schema))    
    yml = yaml.load(inputs)
    
    airfoils = [Airfoil(af) for af in yml.get('airfoils')]
    materials = read_IEA37_materials(yml.get('materials'))
    
    job = Blade(name='IEAonshoreWT')
    job.iea37_converter(yml.get('components').get('blade'), airfoils, materials, wt_flag=True)