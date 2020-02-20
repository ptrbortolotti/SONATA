# -*- coding: utf-8 -*-
"""
Created on Wednesday Dec 18 20:43:28 2019

@author: Roland Feil
"""
# ============================================= #

import numpy as np
from scipy.interpolate import PchipInterpolator, interp1d

def utls_openmdao_apply_gains(blade, cs_pos, byml, cbmconfigs, opt_vars):

    # # Determine Arc locations of web
    # af_grid = blade.airfoils[:,0]
    # for i in range(len(cs_pos)):
    #     grid_ind = int(np.where(af_grid == cs_pos[i])[0])  # determine grid index for current radial position
    #     profile_xy         = blade.airfoils[grid_ind,1].coordinates
    #     id_le           = np.argmin(profile_xy[:,0])
    #
    #     if np.mean(profile_xy[0:id_le, 1]) < 0:  # check that profile starts with suction side (first part) and then continues with pressure side (second part); otherwise flip
    #         profile_xy = np.flip(profile_xy,0)
    #
    #     # get pitch axis at current radial location
    #     set_interp = PchipInterpolator(byml.get('outer_shape_bem').get('pitch_axis')['grid'], byml.get('outer_shape_bem').get('pitch_axis')['values'])
    #     pitch_axis_loc = set_interp(cs_pos[i])  # nondim chordwise location of web 1 [-]
    #
    #     # get chord at current radial position
    #     set_interp = PchipInterpolator(byml.get('outer_shape_bem').get('chord')['grid'], byml.get('outer_shape_bem').get('chord')['values'])
    #     chord = set_interp(cs_pos[i])  # chord length [m]
    #
    #
    #     # get info for fore web
    #     set_interp = PchipInterpolator(byml.get('internal_structure_2d_fem').get('webs')[0]['offset_x_pa']['grid'], byml.get('internal_structure_2d_fem').get('webs')[0]['offset_x_pa']['values'])
    #     web0_offset = set_interp(cs_pos[i])  # chordwise location of web 1
    #
    #     set_interp = PchipInterpolator(byml.get('internal_structure_2d_fem').get('webs')[0]['rotation']['grid'], byml.get('internal_structure_2d_fem').get('webs')[0]['rotation']['values'])
    #     web0_rot = set_interp(cs_pos[i])  # rotation of web1 in rad
    #
    #     web_start_nd, web_end_nd = calc_axis_intersection(profile_xy, web0_rot, pitch_axis_loc+web0_offset/chord, [0., 0.],['suction', 'pressure'])







    mod_cbmconfigs = cbmconfigs
    orig_value = cbmconfigs[0][1].webs[1]['Pos1']
    mod_value = cbmconfigs[0][1].webs[1]['Pos1'] + opt_vars[0]
    # optimize chordwise location of the aftweb
    mod_cbmconfigs[0][1].webs[1]['Pos1'] = cbmconfigs[0][1].webs[1]['Pos1'] + opt_vars[0]
    mod_cbmconfigs[0][1].webs[1]['Pos2'] = cbmconfigs[0][1].webs[1]['Pos2'] - opt_vars[0]

    # optimize chordwise location of the foreweb
    # mod_cbmconfigs[0][1].webs[1]['Pos3'] = cbmconfigs[0][1].webs[1]['Pos3'] - gains[0]
    # mod_cbmconfigs[0][1].webs[1]['Pos4'] = cbmconfigs[0][1].webs[1]['Pos4'] + gains[0]

    # print(str(mod_cbmconfigs[0][1].webs[1]['Pos1'] + 0.5))


    return mod_cbmconfigs





def calc_axis_intersection(xy_coord, rotation, offset, p_le_d, side):
    """
    Determine the airfoil intersection points (in s-coordinates) for web placements from its rotation and offset from the pitch axis

    x axis defined pos. from leading to trailing edge (BeamDyn coordinates)
    y axis defined pos. in thickness direction from pressure to suction side
    (optional) z axis follows the blade span

    Local coordinate sys has (0,0) ref location at leading edge of the airfoil

    Inputs:
    xy_coord    -   airfoil outer shape coordinates
    rotation    -   rotation of web
    offset      -   chordwise offset along the ... axis
    p_le_d      -   [0., 0.]
    side        -   ['suction', 'pressure']

    calc_axis_intersection(inputs['coord_xy_dim'][i,:,:], web_rotation[j,i], outputs['web_offset_y_pa'][j,i], [0.,0.], ['suction', 'pressure'])

    """

    # rotation
    offset_x = offset * np.cos(rotation) + p_le_d[0]
    offset_y = offset * np.sin(rotation) + p_le_d[1]

    # m_rot = np.sin(rotation) / np.cos(rotation)  # slope of rotated axis
    # plane_rot = [m_rot, -1 * m_rot * p_le_d[0] + p_le_d[1]]  # coefficients for rotated axis line: a1*x + a0

    m_intersection = np.sin(rotation + np.pi / 2.) / np.cos(rotation + np.pi / 2.)  # slope perpendicular to rotated axis
    plane_intersection = [m_intersection, -1 * m_intersection * offset_x + offset_y]  # coefficients for line perpendicular to rotated axis line at the offset: a1*x + a0

    # intersection between airfoil surface and the line perpendicular to the rotated/offset axis
    y_intersection = np.polyval(plane_intersection, xy_coord[:, 0])

    idx_le = np.argmin(xy_coord[:, 0])
    xy_coord_arc = arc_length(xy_coord[:, 0], xy_coord[:, 1])
    arc_L = xy_coord_arc[-1]
    xy_coord_arc /= arc_L

    # import matplotlib.pyplot as plt
    # fig, ax = plt.subplots()
    # ax.plot(xy_coord[:,0], xy_coord[:,1])

    idx_inter = np.argwhere(np.diff(np.sign(xy_coord[:, 1] - y_intersection))).flatten()  # find closest airfoil surface points to intersection

    midpoint_arc = []
    for sidei in side:
        if sidei.lower() == 'suction':
            tangent_line = np.polyfit(xy_coord[idx_inter[0]:idx_inter[0] + 2, 0], xy_coord[idx_inter[0]:idx_inter[0] + 2, 1], 1)
        elif sidei.lower() == 'pressure':
            tangent_line = np.polyfit(xy_coord[idx_inter[1]:idx_inter[1] + 2, 0], xy_coord[idx_inter[1]:idx_inter[1] + 2, 1], 1)

        midpoint_x = (tangent_line[1] - plane_intersection[1]) / (plane_intersection[0] - tangent_line[0])
        midpoint_y = plane_intersection[0] * (tangent_line[1] - plane_intersection[1]) / (plane_intersection[0] - tangent_line[0]) + plane_intersection[1]

        # convert to arc position
        if sidei.lower() == 'suction':
            x_half = xy_coord[:idx_le + 1, 0]
            arc_half = xy_coord_arc[:idx_le + 1]

        elif sidei.lower() == 'pressure':
            x_half = xy_coord[idx_le:, 0]
            arc_half = xy_coord_arc[idx_le:]

        midpoint_arc.append(remap2grid(x_half, arc_half, midpoint_x, spline=interp1d))

    return midpoint_arc




def arc_length(x, y, z=[]):
    npts = len(x)
    arc = np.array([0.]*npts)
    if len(z) == len(x):
        for k in range(1, npts):
            arc[k] = arc[k-1] + np.sqrt((x[k] - x[k-1])**2 + (y[k] - y[k-1])**2 + (z[k] - z[k-1])**2)
    else:
        for k in range(1, npts):
            arc[k] = arc[k-1] + np.sqrt((x[k] - x[k-1])**2 + (y[k] - y[k-1])**2)

    return arc




def remap2grid(x_ref, y_ref, x, spline=PchipInterpolator):


    try:
        spline_y = spline(x_ref, y_ref)
    except:
        x_ref = np.flip(x_ref, axis=0)
        y_ref = np.flip(y_ref, axis=0)
        spline_y = spline(x_ref, y_ref)

    # error handling for x[-1] - x_ref[-1] > 0 and x[-1]~x_ref[-1]
    try:
        _ = iter(x)
        if x[-1]>max(x_ref) and np.isclose(x[-1], x_ref[-1]):
            x[-1]=x_ref[-1]
    except:
        if np.isclose(x, 0.):
            x = 0.
        if x>max(x_ref) and np.isclose(x, x_ref[-1]):
            x=x_ref[-1]

    y_out = spline_y(x)

    np.place(y_out, y_out < min(y_ref), min(y_ref))
    np.place(y_out, y_out > max(y_ref), max(y_ref))

    return y_out