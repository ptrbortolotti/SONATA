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
# from jsonschema import validate
import matplotlib.pyplot as plt


if __name__ == '__main__':
    os.chdir('../..')

from SONATA.cbm.classCBMConfig import CBMConfig
from SONATA.utl_openmdao.utl_openmdao import calc_axis_intersection

# from jobs.RFeil.utls.utls_openmdao import calc_axis_intersection


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

def converter_WT(blade, cs_pos, byml, materials, mesh_resolution):
    """
    Converts the 
    @author: Pietro Bortolotti, Roland Feil
    
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

    tmp0        = byml.get('internal_structure_2d_fem').get('webs')
    x           = cs_pos # Non dimensional span position of the stations

    if tmp0 == None:
        n_webs = 0
    else:
        n_webs  = len(tmp0) # Max number of webs along span
    web_exist  = np.zeros((len(x), n_webs), dtype=int) # Flag to set whether webs have non zero thickness, initialized at 0

    tmp2    = [dict([('position', x[n])]) for n in range(len(x))]
    id_webs = [dict() for n in range(len(x))]
    
    # Get sections information and init the CBM instances.
    tmp1 = byml.get('internal_structure_2d_fem').get('layers')
    
    # Set the web_exist flag. This checks whether at every station there is at least a non-zero thickness
    # layer defined in the web. If there isn't, webs are not built even if they are defined in terms of start and end positions
    for i in range(len(x)):
        for idx_sec, sec in enumerate(tmp1):
            if x[i] >= sec['thickness']['grid'][0] and x[i] <= sec['thickness']['grid'][-1]:
                set_interp_thick = PchipInterpolator(sec['thickness']['grid'], sec['thickness']['values'])
                thick_i = set_interp_thick(x[i])
                if thick_i > 1.e-6:
                    if 'web' in sec.keys():
                        for j in range(len(tmp0)):
                            if sec['web'] == tmp0[j]['name']:
                                web_exist[i,j]   = 1
                                break


    span_adhesive  = 0.2
    # Determine start and end positions (s-coordinates) of each web
    for i in range(len(x)):
        n_webs_i = 0
        web_config = [dict() for n in range(len(x))]
        for j in range(n_webs):
            if x[i] >=tmp0[j]['start_nd_arc']['grid'][0] and x[i] <=tmp0[j]['start_nd_arc']['grid'][-1] and web_exist[i,j] == 1:
                n_webs_i = n_webs_i + 1

                id_webs[i][tmp0[j]['name']]       = {}
                id_webs[i][tmp0[j]['name']]['id'] = 2 * j + 2

                set_interp = PchipInterpolator(tmp0[j]['start_nd_arc']['grid'], tmp0[j]['start_nd_arc']['values'])
                id_webs[i][tmp0[j]['name']]['start_nd_arc'] = set_interp(x[i])

                set_interp = PchipInterpolator(tmp0[j]['end_nd_arc']['grid'], tmp0[j]['end_nd_arc']['values'])
                id_webs[i][tmp0[j]['name']]['end_nd_arc'] = set_interp(x[i])

        if n_webs_i == 0:
            tmp2[i]['segments']                 = [{}]
            tmp2[i]['segments'][0]['id']        = 0
            tmp2[i]['segments'][0]['filler']    = None
            tmp2[i]['segments'][0]['layup']     = [{}]

        elif x[i] > span_adhesive:
            tmp2[i]['segments'] = [dict([('id', n),('layup', [{}]),('filler', None)]) for n in range(2 * n_webs_i + 3)]  # inits amount of segments in tmp2
        else:
            tmp2[i]['segments'] = [dict([('id', n),('layup', [{}]),('filler', None)]) for n in range(2 * n_webs_i + 2)]  # inits amount of segments in tmp2


    # Composite stacking sequences
    thick_web = np.zeros([len(x),n_webs])
    id_layer_web_le= np.zeros(len(x), dtype=int)
    id_layer_web_te = np.zeros(len(x), dtype=int)
    adhesive_extent = np.zeros(len(x))
    # id_count = np.zeros(len(x), dtype=int)

    for i in range(len(x)):

        id_profile = np.argmin(np.abs(blade.blade_ref_axis[:,0]-x[i]))

        profile         = blade.airfoils[id_profile,1].coordinates
        id_le           = np.argmin(profile[:,0])

        if np.mean(profile[0:id_le, 1]) < 0:
            profile = np.flip(profile,0)

        profile_curve   = arc_length(profile[:,0], profile[:,1]) / arc_length(profile[:,0], profile[:,1])[-1]

        id_layer     = 0
        web_filler_index = False  # introduce web filler index to separate leading and trailing web layups

        for idx_sec, sec in enumerate(tmp1):

            if x[i] >= sec['thickness']['grid'][0] and x[i] <= sec['thickness']['grid'][-1]:
                
                set_interp_thick = PchipInterpolator(sec['thickness']['grid'], sec['thickness']['values'])
                thick_i = float(set_interp_thick(x[i]))  # added float

                if 'start_nd_arc' in sec.keys():
                    set_interp = PchipInterpolator(sec['start_nd_arc']['grid'], sec['start_nd_arc']['values'])
                    start_i     = float(set_interp(x[i]))
                    if start_i>0. and start_i<0.1 and id_layer>0:
                        for kk in range(len(tmp2[i]['segments'][0]['layup'])):
                            if abs(start_i - tmp2[i]['segments'][0]['layup'][kk]['end']) < 1.e-5:
                                start_i += 5.e-3
                else:
                    start_i = 0
                if 'end_nd_arc' in sec.keys():
                    set_interp = PchipInterpolator(sec['end_nd_arc']['grid'], sec['end_nd_arc']['values'])
                    end_i     = float(set_interp(x[i]))
                    # if end_i>0.9 and id_layer>0:
                    #     for kk in range(len(tmp2[i]['segments'][0]['layup'])):
                    #         if end_i == tmp2[i]['segments'][0]['layup'][kk]['start']:
                    #             start_i -= 1.e-2
                else:
                    end_i = 1

                if thick_i > 1.e-6 and abs(start_i - end_i) > 1.e-3:
                    if 'web' not in sec.keys():                        
                        if idx_sec>0:
                            tmp2[i]['segments'][0]['layup'].append({})
                        tmp2[i]['segments'][0]['layup'][id_layer]['thickness']     = thick_i
                        tmp2[i]['segments'][0]['layup'][id_layer]['name']          = sec['material'] + '_' + str(x[i])
                        tmp2[i]['segments'][0]['layup'][id_layer]['material_name'] = sec['material']
                        if 'start_nd_arc' in sec.keys():

                            set_interp = PchipInterpolator(sec['start_nd_arc']['grid'], sec['start_nd_arc']['values'])
                            tmp2[i]['segments'][0]['layup'][id_layer]['start']     = float(set_interp(x[i]))  # added float

                            set_interp = PchipInterpolator(sec['end_nd_arc']['grid'], sec['end_nd_arc']['values'])
                            tmp2[i]['segments'][0]['layup'][id_layer]['end']       = float(set_interp(x[i]))  # added float

                        else:
                            tmp2[i]['segments'][0]['layup'][id_layer]['start']     = 0.
                            tmp2[i]['segments'][0]['layup'][id_layer]['end']       = 1.
                        if 'fiber_orientation' in sec.keys():
                            set_interp = PchipInterpolator(sec['fiber_orientation']['grid'], sec['fiber_orientation']['values'])
                            tmp2[i]['segments'][0]['layup'][id_layer]['orientation'] = float(set_interp(x[i]) * 180 / np.pi)  # added float
                        else:
                            tmp2[i]['segments'][0]['layup'][id_layer]['orientation'] = 0.
                            
                        # # Check consistency
                        # if tmp2[i]['segments'][0]['layup'][id_layer]['end'] < tmp2[i]['segments'][0]['layup'][id_layer]['start']:
                        #     exit('WARNING: Layer ' + tmp2[i]['segments'][0]['layup'][id_layer]['name'] + ' ends before it starts. Check the yaml input file!!')
                        ch = np.interp(x[i], blade.chord[:,0], blade.chord[:,1])
                        adhesive_extent[i] = min([0.04, 0.04 / ch])
                        if x[i] > span_adhesive and tmp2[i]['segments'][0]['layup'][id_layer]['start'] < adhesive_extent[i] and tmp2[i]['segments'][0]['layup'][id_layer]['end'] < 0.5:
                            tmp2[i]['segments'][0]['layup'][id_layer]['start'] = adhesive_extent[i] 
                        elif x[i] > span_adhesive and tmp2[i]['segments'][0]['layup'][id_layer]['end'] > 1. - adhesive_extent[i] and tmp2[i]['segments'][0]['layup'][id_layer]['start'] > 0.5:
                            tmp2[i]['segments'][0]['layup'][id_layer]['end'] = 1. - adhesive_extent[i]

                            # old_start = tmp2[i]['segments'][0]['layup'][id_layer]['start']
                            # old_end   = tmp2[i]['segments'][0]['layup'][id_layer]['end']
                            # if old_start < 1. - adhesive_extent[i]:
                            #     tmp2[i]['segments'][0]['layup'][id_layer]['end']   = old_end
                            #     tmp2[i]['segments'][0]['layup'][id_layer]['start'] = adhesive_extent[i]
                            # else:
                            #     exit('ERROR: The converter is not modeling well the trailing edge')

                            # tmp2[i]['segments'][0]['layup'].append({})
                            # tmp2[i]['segments'][0]['layup'][id_layer + 1]['name']      = tmp2[i]['segments'][0]['layup'][id_layer]['name'] + '_2'
                            # tmp2[i]['segments'][0]['layup'][id_layer + 1]['material_name'] = tmp2[i]['segments'][0]['layup'][id_layer]['material_name']
                            # tmp2[i]['segments'][0]['layup'][id_layer + 1]['thickness']     = tmp2[i]['segments'][0]['layup'][id_layer]['thickness']
                            # tmp2[i]['segments'][0]['layup'][id_layer + 1]['orientation'] = tmp2[i]['segments'][0]['layup'][id_layer]['orientation']
                            # if old_end > adhesive_extent[i]:
                            #     tmp2[i]['segments'][0]['layup'][id_layer + 1]['start'] = old_start
                            #     tmp2[i]['segments'][0]['layup'][id_layer + 1]['end']   = 1. - adhesive_extent[i]
                            # else:
                            #     exit('ERROR: The converter is not modeling well the trailing edge')
                            # id_layer = id_layer + 1
                        
                        id_layer = id_layer + 1
                    else:  # if web in sec.keys():


                        id_seg          = id_webs[i][(sec['web'])]['id']
                        for id_mat in range(1,len(materials)+1):
                            if sec['material'] == materials[id_mat].name:


                                # ------------------------------------------
                                # Use the following code for separated webs
                                # ------------------------------------------
                                # ========= NEW with split yaml files ========= #
                                # Orthotropic at leading edge (_le)
                                # if '_le' in sec['name']:  #  check if on leading edge side of web
                                if (not web_filler_index and not isinstance(materials[id_mat].E, float)):  # begin with leading part of the web BEFORE web has been filled
                                    id_layer_web_le[i] = 0
                                    if tmp2[i]['segments'][id_seg - 1]['layup'] != [{}]:
                                        tmp2[i]['segments'][id_seg - 1]['layup'].append({})
                                        id_layer_web_le[i] = id_layer_web_le[i] + 1

                                    tmp2[i]['segments'][id_seg - 1]['layup'][id_layer_web_le[i]]['name'] = 'web_' + str(int(id_seg / 2)) + '_' + sec['material'] + '_' + str(x[i])
                                    tmp2[i]['segments'][id_seg - 1]['layup'][id_layer_web_le[i]]['material_name'] = sec['material']
                                    # set_interp = PchipInterpolator(sec['thickness']['grid'], sec['thickness']['values'])
                                    tmp2[i]['segments'][id_seg - 1]['layup'][id_layer_web_le[i]]['thickness'] = thick_i  # * 0.5
                                    flange = 0.01
                                    tmp2[i]['segments'][id_seg - 1]['layup'][id_layer_web_le[i]]['start'] = id_webs[i][(sec['web'])]['end_nd_arc'] - flange
                                    tmp2[i]['segments'][id_seg - 1]['layup'][id_layer_web_le[i]]['end']   = id_webs[i][(sec['web'])]['start_nd_arc'] + flange

                                    if 'fiber_orientation' in sec.keys():
                                        set_interp = PchipInterpolator(sec['fiber_orientation']['grid'], sec['fiber_orientation']['values'])
                                        tmp2[i]['segments'][id_seg - 1]['layup'][id_layer_web_le[i]]['orientation'] = float(set_interp(x[i]))
                                    else:
                                        tmp2[i]['segments'][id_seg - 1]['layup'][id_layer_web_le[i]]['orientation'] = 0.

                                    # web_filler_index = False  # <<<<<<<<<<<<<<

                                # Isotropic -> fill the web
                                elif (not web_filler_index and isinstance(materials[id_mat].E, float)):
                                    tmp2[i]['segments'][id_seg]['filler'] = materials[id_mat].name

                                    tmp2[i]['segments'][id_seg]['layup'][0]['name'] = 'dummy'
                                    tmp2[i]['segments'][id_seg]['layup'][0]['material_name'] = sec['material']
                                    tmp2[i]['segments'][id_seg]['layup'][0]['thickness'] = 1e-3
                                    tmp2[i]['segments'][id_seg]['layup'][0]['start'] = 0.0
                                    tmp2[i]['segments'][id_seg]['layup'][0]['end'] = 1.0
                                    tmp2[i]['segments'][id_seg]['layup'][0]['orientation'] = 0.0

                                    # set_interp = PchipInterpolator(sec['thickness']['grid'],sec['thickness']['values'])
                                    thick_web[i, int(id_seg / 2 - 1)] = thick_web[i, int(id_seg / 2 - 1)] + thick_i
                                    if thick_web[i, int(id_seg / 2 - 1)] < 0.01:
                                        thick_web[i, int(id_seg / 2 - 1)] = 0.01
                                        print('WARNING: web filler cannot be thinner than 10mm. This is adjusted here, but please check the input yaml.')
                                    web_filler_index = True  #  changes to true as web has now already been filled with material

                                # Orthotropic at trailing edge (_te)
                                # elif '_te' in sec['name']:  # check if on trailing edge side of web
                                elif (web_filler_index and not isinstance(materials[id_mat].E, float)):  # continue with trailing part of the web AFTER web has been filled

                                    id_layer_web_te[i] = 0
                                    if tmp2[i]['segments'][id_seg + 1]['layup'] != [{}]:
                                        tmp2[i]['segments'][id_seg + 1]['layup'].append({})
                                        id_layer_web_te[i] = id_layer_web_te[i] + 1

                                    tmp2[i]['segments'][id_seg + 1]['layup'][id_layer_web_te[i]]['name'] = 'web_' + str(int(id_seg / 2)) + '_' + sec['material'] + '_' + str(x[i])
                                    tmp2[i]['segments'][id_seg + 1]['layup'][id_layer_web_te[i]]['material_name'] = sec['material']
                                    # set_interp = PchipInterpolator(sec['thickness']['grid'], sec['thickness']['values'])
                                    tmp2[i]['segments'][id_seg + 1]['layup'][id_layer_web_te[i]]['thickness'] = thick_i  # * 0.5
                                    flange = 0.01
                                    tmp2[i]['segments'][id_seg + 1]['layup'][id_layer_web_te[i]]['start'] = id_webs[i][(sec['web'])]['start_nd_arc'] - flange
                                    tmp2[i]['segments'][id_seg + 1]['layup'][id_layer_web_te[i]]['end'] = id_webs[i][(sec['web'])]['end_nd_arc'] + flange

                                    if 'fiber_orientation' in sec.keys():
                                        set_interp = PchipInterpolator(sec['fiber_orientation']['grid'],sec['fiber_orientation']['values'])
                                        tmp2[i]['segments'][id_seg + 1]['layup'][id_layer_web_te[i]]['orientation'] = set_interp(x[i])
                                    else:
                                        tmp2[i]['segments'][id_seg + 1]['layup'][id_layer_web_te[i]]['orientation'] = 0.

                                    web_filler_index = False  #  after completing the te part (this web is finished now!), prepare for next web









                            # ========= OLD without split yaml files ========= #
                            #  NEW YAML FILES INCLUDE ALREADY SEPARATED WEBS  #
                                #
                                # if tmp2[i]['segments'][id_seg - 1]['layup'] != [{}]:
                                #     tmp2[i]['segments'][id_seg - 1]['layup'].append({})
                                #     id_layer_web_le[i] = id_layer_web_le[i] + 1
                                #
                                # if tmp2[i]['segments'][id_seg + 1]['layup'] != [{}]:
                                #     tmp2[i]['segments'][id_seg + 1]['layup'].append({})
                                #     id_layer_web_te[i] = id_layer_web_te[i] + 1
                                #
                                #
                                #
                                # tmp2[i]['segments'][id_seg - 1]['layup'][id_layer_web_le[i]]['name'] = 'web_' + str(int(id_seg/2)) + '_' + sec['material']  + '_' + str(x[i])
                                # tmp2[i]['segments'][id_seg + 1]['layup'][id_layer_web_te[i]]['name'] = 'web_' + str(int(id_seg/2)) + '_' + sec['material']  + '_' + str(x[i])
                                #
                                # tmp2[i]['segments'][id_seg - 1]['layup'][id_layer_web_le[i]]['material_name'] = sec['material']
                                # tmp2[i]['segments'][id_seg + 1]['layup'][id_layer_web_te[i]]['material_name'] = sec['material']
                                #
                                # set_interp = PchipInterpolator(sec['thickness']['grid'], sec['thickness']['values'])
                                # tmp2[i]['segments'][id_seg - 1]['layup'][id_layer_web_le[i]]['thickness'] = set_interp(x[i]) * 0.5
                                # tmp2[i]['segments'][id_seg + 1]['layup'][id_layer_web_te[i]]['thickness'] = set_interp(x[i]) * 0.5
                                #
                                # flange = 0.01
                                #
                                # tmp2[i]['segments'][id_seg - 1]['layup'][id_layer_web_le[i]]['start'] = id_webs[i][(sec['web'])]['end_nd_arc'] - flange
                                # tmp2[i]['segments'][id_seg - 1]['layup'][id_layer_web_le[i]]['end']   = id_webs[i][(sec['web'])]['start_nd_arc'] + flange
                                #
                                #
                                # tmp2[i]['segments'][id_seg + 1]['layup'][id_layer_web_te[i]]['start'] = id_webs[i][(sec['web'])]['start_nd_arc'] - flange
                                # tmp2[i]['segments'][id_seg + 1]['layup'][id_layer_web_te[i]]['end']   = id_webs[i][(sec['web'])]['end_nd_arc'] + flange
                                #
                                # if 'fiber_orientation' in sec.keys():
                                #     set_interp = PchipInterpolator(sec['fiber_orientation']['grid'], sec['fiber_orientation']['values'])
                                #     tmp2[i]['segments'][id_seg - 1]['layup'][id_layer_web_le[i]]['orientation'] = set_interp(x[i])
                                #     tmp2[i]['segments'][id_seg + 1]['layup'][id_layer_web_te[i]]['orientation'] = set_interp(x[i])
                                # else:
                                #     tmp2[i]['segments'][id_seg - 1]['layup'][id_layer_web_le[i]]['orientation'] = 0.
                                #     tmp2[i]['segments'][id_seg + 1]['layup'][id_layer_web_te[i]]['orientation'] = 0.
        # if x[i] > span_adhesive and n_webs > 0:
        if x[i] > span_adhesive and len(tmp2[i]['segments']) > 1:
            # id_seg = n_webs*2 + 2
            tmp2[i]['segments'][-1]['filler'] = 'Adhesive'
            tmp2[i]['segments'][-1]['layup'][0]['name'] = 'dummy'
            tmp2[i]['segments'][-1]['layup'][0]['material_name'] = 'Adhesive'
            tmp2[i]['segments'][-1]['layup'][0]['thickness'] = 5.e-4
            tmp2[i]['segments'][-1]['layup'][0]['start'] = 0.0
            tmp2[i]['segments'][-1]['layup'][0]['end'] = 1.0
            tmp2[i]['segments'][-1]['layup'][0]['orientation'] = 0.0

    # Split webs for determining the
    # print(tmp2)
    # exit()
    webs    = [OrderedDict() for n in range(len(x))]
    for i in range(len(x)):
        for i_web,web in enumerate(tmp0):
            if x[i] >= web['start_nd_arc']['grid'][0] and x[i] <= web['start_nd_arc']['grid'][-1] and web_exist[i,i_web] == 1:

                set_interp      = PchipInterpolator(web['start_nd_arc']['grid'], web['start_nd_arc']['values'])
                start           = set_interp(x[i])

                set_interp      = PchipInterpolator(web['end_nd_arc']['grid'], web['end_nd_arc']['values'])
                end             = set_interp(x[i])

                id_profile      = np.argmin(np.abs(blade.blade_ref_axis[:,0]-x[i]))
                profile         = blade.airfoils[id_profile,1].coordinates
                id_le           = np.argmin(profile[:,0])

                if np.mean(profile[0:id_le, 1]) < 0:
                    profile     = np.flip(profile,0)

                # checks for flatback airfoils
                # profile[0,0] = 1.
                # profile[0,1] = 1.e-4
                # profile[-1,0] = 1.
                # profile[-1,1] = -1.e-4
                if len(np.where(profile[:,0] == 1.)[0]) > 2:
                    profile = profile[:np.where(profile[:,0]==1.)[0][2] , :]

                # if len(profile[:,0]) != len(np.unique(profile[:,0])):
                #     raise Exception('Airfoil at station {:d} does not have unique x points'.format(id_profile))

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


                if 'curvature' in web.keys():
                    set_interp = PchipInterpolator(web['curvature']['grid'], web['curvature']['values'])
                    curve_val = float(set_interp(x[i]))
                else:
                    curve_val = 0.



                w_f = {}
                w_f['id']   = 2 * i_web + 1
                w_f['position'] = [web_start_le, web_end_le]
                w_f['curvature'] = curve_val
                webs[i][2 * i_web + 1] = w_f

                w_r = {}
                w_r['id']   = 2 * i_web + 2
                w_r['position'] = [web_start_te, web_end_te]
                w_r['curvature'] = curve_val
                webs[i][2 * i_web + 2] = w_r
    
    for i in range(len(x)):
        if x[i] > span_adhesive and len(tmp2[i]['segments']) > 1:
            id_segment = sum(web_exist[i,:])*2 + 1
            webs[i][id_segment] = {'curvature': 0.0, 'id': id_segment, 'position': [adhesive_extent[i], 1.-adhesive_extent[i]]}
    lst = []

    for (seg,w,x) in zip(tmp2, webs, x):
        tmp = {'webs':list(w.values()), 'segments':seg['segments'], 'position':x, 'mesh_resolution':mesh_resolution}
        lst.append([x, CBMConfig(tmp, materials)])

    return np.asarray(lst)




#%% MAIN
if __name__ == '__main__':
    from SONATA.classAirfoil import Airfoil
    from SONATA.classBlade import Blade
    from SONATA.classMaterial import read_materials

    with open('./jobs/PBortolotti/IEAonshoreWT.yaml', 'r') as myfile:
        inputs  = myfile.read()
    with open('jobs/PBortolotti/IEAontology_schema.yaml', 'r') as myfile:
        schema  = myfile.read()
    # validate(yaml.load(inputs), yaml.load(schema))
    yml = yaml.load(inputs)

    airfoils = [Airfoil(af) for af in yml.get('airfoils')]
    materials = read_materials(yml.get('materials'))

    job = Blade(name='IEAonshoreWT')
    job.converter_WT(yml.get('components').get('blade'), airfoils, materials, wt_flag=True)
