# -*- coding: utf-8 -*-
"""
Created on Wednesday Dec 18 20:43:28 2019

@author: Roland Feil
"""
# ============================================= #

import numpy as np
from scipy.interpolate import PchipInterpolator, interp1d

def utl_openmdao_apply_gains_mat_thickness(blade, yml, opt_vars):

    # Optimize the thickness of the outer and the inner shell skins

    # ------------------------------------------- #
    # (1) Shell skin optimization of Kevlar_10 - initial values are:
    # ------------------------------------------- #
    # t_outer = yml.get('internal_structure_2d_fem').get('sections')[0]['segments'][0]['layup'][1][2]  # shell skin outer
    # t_inner = yml.get('internal_structure_2d_fem').get('sections')[0]['segments'][0]['layup'][6][2]  # shell skin inner  [10]
    #
    # # replace with optimization variable
    # yml.get('internal_structure_2d_fem').get('sections')[0]['segments'][0]['layup'][1][2] = float(opt_vars[0])
    # yml.get('internal_structure_2d_fem').get('sections')[0]['segments'][0]['layup'][6][2] = float(opt_vars[0])  # [10]




    # ------------------------------------------- #
    # (2) Filament wound optimization of Kevlar_11 - initial values are:
    # ------------------------------------------- #
    thickness = yml.get('internal_structure_2d_fem').get('sections')[0]['segments'][2]['layup'][1][2]  # filament wound thickness  (split - but only optimize inner ply thickness)
    web1_curvature = yml.get('internal_structure_2d_fem').get('sections')[0]['webs'][0]['curvature']  # web 1 curvature
    web2_curvature = yml.get('internal_structure_2d_fem').get('sections')[0]['webs'][1]['curvature']  # web 2 curvature
    ply_orientation_I = yml.get('internal_structure_2d_fem').get('sections')[0]['segments'][2]['layup'][0][3]
    ply_orientation_II = yml.get('internal_structure_2d_fem').get('sections')[0]['segments'][2]['layup'][1][3]

    web1_start = yml.get('internal_structure_2d_fem').get('sections')[0]['webs'][0]['position'][0]
    web1_start = yml.get('internal_structure_2d_fem').get('sections')[0]['webs'][0]['position'][1]
    web2_start = yml.get('internal_structure_2d_fem').get('sections')[0]['webs'][1]['position'][0]
    web2_start = yml.get('internal_structure_2d_fem').get('sections')[0]['webs'][1]['position'][1]

    # replace with optimization variable
    # yml.get('internal_structure_2d_fem').get('sections')[0]['segments'][2]['layup'][0][2] = float(opt_vars[0])  # split elliplis plies in two for better meshing
    yml.get('internal_structure_2d_fem').get('sections')[0]['segments'][2]['layup'][1][2] = float(opt_vars[0])  # split elliplis plies in two for better meshing
    yml.get('internal_structure_2d_fem').get('sections')[0]['webs'][0]['curvature'] = -float(opt_vars[1])
    yml.get('internal_structure_2d_fem').get('sections')[0]['webs'][1]['curvature'] =  float(opt_vars[1])
    yml.get('internal_structure_2d_fem').get('sections')[0]['segments'][2]['layup'][0][3] = float(opt_vars[2])
    yml.get('internal_structure_2d_fem').get('sections')[0]['segments'][2]['layup'][1][3] = float(opt_vars[2])
    # yml.get('internal_structure_2d_fem').get('sections')[0]['segments'][2]['layup'][2][3] = float(opt_vars[2])

    yml.get('internal_structure_2d_fem').get('sections')[0]['webs'][0]['position'][0] = 0.33763813598041315 + float(opt_vars[3])
    yml.get('internal_structure_2d_fem').get('sections')[0]['webs'][0]['position'][1] = 0.6603507085422332 - float(opt_vars[3])
    yml.get('internal_structure_2d_fem').get('sections')[0]['webs'][1]['position'][0] = 0.33763813598041315 - float(opt_vars[3])
    yml.get('internal_structure_2d_fem').get('sections')[0]['webs'][1]['position'][1] = 0.6603507085422332 + float(opt_vars[3])

    # yml.get('internal_structure_2d_fem').get('sections')[0]['segments'][2]['layup'][0][3] = float(opt_vars[0])
    # yml.get('internal_structure_2d_fem').get('sections')[0]['segments'][2]['layup'][1][3] = float(opt_vars[0])


    # ------------------------------------------- #
    # Spar caps thickness optimization - initial values are:
    # ------------------------------------------- #
    # t_ss = yml.get('internal_structure_2d_fem').get('sections')[0]['segments'][0]['layup'][2][2]  # spar cap suction side
    # t_ps = yml.get('internal_structure_2d_fem').get('sections')[0]['segments'][0]['layup'][3][2]  # spar cap pressure side
    #
    # # replace with optimization variable
    # yml.get('internal_structure_2d_fem').get('sections')[0]['segments'][0]['layup'][2][2] = float(opt_vars[0])
    # yml.get('internal_structure_2d_fem').get('sections')[0]['segments'][0]['layup'][3][2] = float(opt_vars[0])

    return yml

def utl_openmdao_apply_gains_web_placement(blade, yml, opt_vars):


    af_grid = blade.airfoils[:,0]
    web_grid = yml.get('internal_structure_2d_fem').get('webs')[0]['start_nd_arc']['grid']

    # web0_start_nd_arc_opt = np.zeros(len(yml.get('internal_structure_2d_fem').get('webs')[0]['start_nd_arc']['values']))
    # web0_end_nd_arc_opt = np.zeros(len(yml.get('internal_structure_2d_fem').get('webs')[0]['start_nd_arc']['values']))

    start_nd_arc = np.zeros(len(yml.get('internal_structure_2d_fem').get('webs')[0]['start_nd_arc']['values']))
    end_nd_arc = np.zeros(len(yml.get('internal_structure_2d_fem').get('webs')[0]['start_nd_arc']['values']))
    # for i in range(len(cs_pos)):
    for i in range(len(web_grid)):

        grid_ind = int(np.where(af_grid == web_grid[i])[0])  # determine grid index for current radial position
        profile_xy         = blade.airfoils[grid_ind,1].coordinates
        id_le           = np.argmin(profile_xy[:,0])

        if np.mean(profile_xy[0:id_le, 1]) < 0:  # check that profile starts with suction side (first part) and then continues with pressure side (second part); otherwise flip
            profile_xy = np.flip(profile_xy,0)

        # Read data from yaml file (doesn't change during optimization)
        set_interp = PchipInterpolator(yml.get('internal_structure_2d_fem').get('webs')[0]['start_nd_arc']['grid'], yml.get('internal_structure_2d_fem').get('webs')[0]['start_nd_arc']['values'])
        web0_start_nd_arc = set_interp(web_grid[i])  # nondim arc position
        set_interp = PchipInterpolator(yml.get('internal_structure_2d_fem').get('webs')[0]['end_nd_arc']['grid'], yml.get('internal_structure_2d_fem').get('webs')[0]['end_nd_arc']['values'])
        web0_end_nd_arc = set_interp(web_grid[i])    # nondim arc position

        web0_nd_offset, web0_angle, web0_line_coeffs = calc_axis_chorwise_offset_angle(profile_xy, web0_start_nd_arc, web0_end_nd_arc)

        # print('Web0 chordwise offset: ' + str(web0_nd_offset))

        web0_nd_offset_opt = web0_nd_offset + opt_vars

        # print('Web0 chordwise offset (opt): ' + str(web0_nd_offset))

        # start_nd_arc[i], end_nd_arc[i] = calc_arc_location(profile_xy, web0_line_coeffs, web0_nd_offset_opt, ['suction', 'pressure'])
        start_nd_arc[i], end_nd_arc[i] = calc_arc_location(profile_xy, web0_line_coeffs, opt_vars, ['suction', 'pressure'])
        # print('Web0 pos. from (orig): ' + str(web0_start_nd_arc) + ' to: ' + str(float(start_nd_arc[i])))
        # print('Web0 pos. from (orig): ' + str(web0_end_nd_arc) + ' to: ' + str(float(end_nd_arc[i])))

        # alternative method:
        # web0_start_nd_arc_opt[i], web0_end_nd_arc_opt[i] = calc_axis_intersection(profile_xy, web0_angle, web0_nd_offset_opt, [0., 0.], ['suction', 'pressure'])



    byml_start_save = yml['internal_structure_2d_fem']['webs'][0]['start_nd_arc']['values']
    byml_end_save = yml['internal_structure_2d_fem']['webs'][0]['end_nd_arc']['values']
    yml['internal_structure_2d_fem']['webs'][0]['start_nd_arc']['values'] = start_nd_arc
    yml['internal_structure_2d_fem']['webs'][0]['end_nd_arc']['values'] = end_nd_arc
    # yml['internal_structure_2d_fem']['webs'][0]['start_nd_arc']['values'] = web0_start_nd_arc_opt
    # yml['internal_structure_2d_fem']['webs'][0]['end_nd_arc']['values'] = web0_end_nd_arc_opt


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







    # mod_cbmconfigs = cbmconfigs
    # orig_value = cbmconfigs[0][1].webs[1]['Pos1']
    # mod_value = cbmconfigs[0][1].webs[1]['Pos1'] + opt_vars[0]
    # # optimize chordwise location of the aftweb
    # mod_cbmconfigs[0][1].webs[1]['Pos1'] = cbmconfigs[0][1].webs[1]['Pos1'] + opt_vars[0]
    # mod_cbmconfigs[0][1].webs[1]['Pos2'] = cbmconfigs[0][1].webs[1]['Pos2'] - opt_vars[0]
    #
    # # optimize chordwise location of the foreweb
    # # mod_cbmconfigs[0][1].webs[1]['Pos3'] = cbmconfigs[0][1].webs[1]['Pos3'] - gains[0]
    # # mod_cbmconfigs[0][1].webs[1]['Pos4'] = cbmconfigs[0][1].webs[1]['Pos4'] + gains[0]
    #
    # # print(str(mod_cbmconfigs[0][1].webs[1]['Pos1'] + 0.5))


    return yml


def calc_axis_chorwise_offset_angle(profile_xy, pointA, pointB):
    """
    Determines offset along the chordwise axis from the connecting line between two points, e.g. the start and the end
    point at the arc of a profile, that are defined in s-coordinates

    Nondimensional coordinates:
    (0, 0) at leading edge
    (0, 1) at trailing edge

    Input:
    profile_xy      -   airfoil outer shape in nondimensonal x and y coordinates
    pointA          -   location of point A, e.g. web0_start_arc or else, in s-coordinates
    pointB          -   location of point B, e.g. web0_start_arc or else, in s-coordinates

    Output:
    Offset          -   offset of connecting line along the chordwise axis
    """
    # Determine arc length and respective s-coordinates
    xy_coord_arc = arc_length(profile_xy[:, 0], profile_xy[:, 1])
    arc_L = xy_coord_arc[-1]    # arc length
    xy_coord_arc /= arc_L       # s-coordinates

    set_interp = PchipInterpolator(xy_coord_arc, profile_xy)
    web0_start_coord = set_interp(pointA)
    web0_end_coord = set_interp(pointB)

    coeffs = np.polyfit((web0_start_coord[0], web0_end_coord[0]), (web0_start_coord[1], web0_end_coord[1]), 1)  # creates linear coefficients through both points; m*x + a
    # m = coeffs[0]  # m
    # a = coeffs[1]  # a

    # === Testplot ===
    # import matplotlib.pyplot as plt
    # # plot airfoil outer shape
    # plt.plot(profile_xy[:,0], profile_xy[:,1])
    # # plot two points
    # plt.plot(web0_start_coord[0], web0_start_coord[1], 'ok')
    # plt.plot(web0_end_coord[0], web0_end_coord[1], 'ok')
    # # plot linear that connects the two points
    # polynomial = np.poly1d(coeffs)
    # x_test = np.arange(0.0, 0.5, 0.01)
    # y_test = np.polyval(polynomial, x_test)
    # plt.plot(x_test, y_test)
    # =================



    offset = np.roots(coeffs)   # determine chordwise offset; nondim

    angle = np.arctan(coeffs[0])  # determine rotational angle; in rad

    return offset, angle, coeffs


def calc_arc_location(profile_xy, coeffs, delta_offset, side):
    """

    """
    # Determine arc length and respective s-coordinates
    idx_le = np.argmin(profile_xy[:, 0])
    xy_coord_arc = arc_length(profile_xy[:, 0], profile_xy[:, 1])
    arc_L = xy_coord_arc[-1]    # arc length
    xy_coord_arc /= arc_L       # s-coordinates

    coeffs_dev = coeffs
    # coeffs_dev[0] = coeffs[0]
    coeffs_dev[1] = coeffs[1] - delta_offset * coeffs[0]   # deviate the y-axis location at x0 for the chordwise offset

    y_intersection = np.polyval(coeffs_dev, profile_xy[:, 0])  # evaluate y-values of linear that is represented through its coefficients (slope:m; y0:a)
    # import matplotlib.pyplot as plt
    # plt.plot(profile_xy[:, 0], y_intersection)

    # polynomial = np.poly1d(coeffs_dev)

    idx_inter = np.argwhere(np.diff(np.sign(profile_xy[:, 1] - y_intersection))).flatten()

    # start_nd_arc = xy_coord_arc[int(idx_inter[0])]
    # end_nd_arc = xy_coord_arc[int(idx_inter[1])]


    midpoint_arc = []
    for sidei in side:
        if sidei.lower() == 'suction':
            tangent_line = np.polyfit(profile_xy[idx_inter[0]:idx_inter[0] + 2, 0], profile_xy[idx_inter[0]:idx_inter[0] + 2, 1], 1)
        elif sidei.lower() == 'pressure':
            tangent_line = np.polyfit(profile_xy[idx_inter[1]:idx_inter[1] + 2, 0], profile_xy[idx_inter[1]:idx_inter[1] + 2, 1], 1)

        midpoint_x = (tangent_line[1] - coeffs_dev[1]) / (coeffs_dev[0] - tangent_line[0])
        midpoint_y = coeffs_dev[0] * (tangent_line[1] - coeffs_dev[1]) / (coeffs_dev[0] - tangent_line[0]) + coeffs_dev[1]

        # convert to arc position
        if sidei.lower() == 'suction':
            x_half = profile_xy[:idx_le + 1, 0]
            arc_half = xy_coord_arc[:idx_le + 1]

        elif sidei.lower() == 'pressure':
            x_half = profile_xy[idx_le:, 0]
            arc_half = xy_coord_arc[idx_le:]

        midpoint_arc.append(remap2grid(x_half, arc_half, midpoint_x, spline=interp1d))




    # === Testplot ===
    # import matplotlib.pyplot as plt
    # # plot airfoil outer shape
    # plt.plot(profile_xy[:,0], profile_xy[:,1])
    #
    # # plot linear that connects the two points based on deviated coefficients
    # polynomial = np.poly1d(coeffs_dev)
    # x_test = np.arange(0.0, 0.5, 0.01)
    # y_test = np.polyval(polynomial, x_test)
    # plt.plot(x_test, y_test)
    #
    # # plot two points
    # plt.plot(profile_xy[idx_inter[0], 0], profile_xy[idx_inter[0], 1], 'or')
    # plt.plot(profile_xy[idx_inter[1], 0], profile_xy[idx_inter[1], 1], 'or')
    #
    #
    # plt.plot(float(midpoint_arc[0]), float(midpoint_arc[1]), 'ob')
    # =================


    return midpoint_arc


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

    # m_intersection = np.sin(rotation + np.pi / 2.) / np.cos(rotation + np.pi / 2.)  # slope perpendicular to rotated axis
    # m_intersection = np.sin(rotation) / np.cos(rotation)  # slope perpendicular to rotated axis
    m_intersection = np.tan(rotation)  # slope perpendicular to rotated axis

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