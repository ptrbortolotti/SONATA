#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wednesday Nov 20 13:33:33 2019

@author: Roland Feil
"""


from builtins import len, range
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
import csv

from SONATA.cbm.cbm_utl import trsf_sixbysix
from jobs.RFeil.utls.utls_sonata2beamdyn import convert_structdef_SONATA_to_beamdyn, write_beamdyn_axis, write_beamdyn_prop

from jobs.RFeil.utls.utls_analytical_rectangle import utls_analytical_rectangle


def beam_struct_eval(flags_dict, cs_pos, job, folder_str, job_str):

    """
    Analyse, transform, evaluate and plot structural results

    Functions:
    beam_struct_eval    - parent function (retrieve & transform data, structure of data evaluation)
    plot_vabs_anbax     - child function to plot the 6x6 stiffness and mass matrices 
    export_beam_struct_properties - csv export of structural beam properties

    Inputs:
    flags_dict          - dictionary containing relevant flags
    cs_pos              - radial station of blade cross sections
    job                 - contains the whole blade data (yaml file content, wires, mesh, etc.)
    folder_str          - name of operating folder
    job_str             - name of job that is currently under investigation 

    Outputs: 
    - None -

    """

    if flags_dict['flag_run_vabs']:
        # ------------------------ #
        # runs both VABS & ANBAX and conducts a verification/ comparison of the results

        # --- VABS --- #
        job.blade_run_vabs()  # run VABS

        if flags_dict['flag_DeamDyn_def_transform']:
            beam_prop = convert_structdef_SONATA_to_beamdyn(cs_pos, job.beam_properties)  # convert to BeamDyn coord sys def
            coordsys = 'BeamDyn'
            str_ext = '_BeamDyn_def'
        else: 
            beam_prop = {}
            beam_prop['beam_section_mass'] = np.zeros([len(cs_pos), 1])
            beam_prop['beam_stiff'] = np.zeros([len(cs_pos), 6, 6])
            beam_prop['beam_inertia'] = np.zeros([len(cs_pos), 6, 6])
            beam_prop['beam_mass_center'] = np.zeros([len(cs_pos), 2])
            beam_prop['beam_neutral_axes'] = np.zeros([len(cs_pos), 2])
            beam_prop['beam_geometric_center'] = np.zeros([len(cs_pos), 2])
            beam_prop['beam_shear_center'] = np.zeros([len(cs_pos), 2])
            for i in range(len(cs_pos)):
                beam_prop['beam_section_mass'][i]       = job.beam_properties[i, 1].m00
                beam_prop['beam_mass_center'][i]        = np.array(job.beam_properties[i, 1].Xm[:])
                beam_prop['beam_neutral_axes'][i]       = np.array(job.beam_properties[i, 1].Xt[:])
                beam_prop['beam_geometric_center'][i]   = np.array(job.beam_properties[i, 1].Xg[:])
                beam_prop['beam_shear_center'][i]       = np.array(job.beam_properties[i, 1].Xs[:])
                for j in range(6):
                    beam_prop['beam_stiff'][i, j, :]    = np.array(job.beam_properties[i, 1].TS[j, :])
                    beam_prop['beam_inertia'][i, j, :]  = np.array(job.beam_properties[i, 1].MM[j, :])


            coordsys = 'VABS/SONATA'
            str_ext = ''


        # --------------------------------------- #
        # Export beam structural properties to csv file
        if flags_dict['flag_csv_export']:
            print('STATUS:\t Export csv files with structural blade characeristics from VABS to: ' + folder_str + 'csv_export/')
            export_beam_struct_properties(folder_str, job_str, cs_pos, coordsys=coordsys, solver='vabs', beam_stiff=beam_prop['beam_stiff'], 
                                          beam_inertia=beam_prop['beam_inertia'], beam_mass_per_length=beam_prop['beam_section_mass'],
                                          beam_mass_center=beam_prop['beam_mass_center'], beam_neutral_axes=beam_prop['beam_neutral_axes'],
                                          beam_geometric_center=beam_prop['beam_geometric_center'], beam_shear_center=beam_prop['beam_shear_center'])

        # --------------------------------------- #
        # write BeamDyn input files
        if flags_dict['flag_write_BeamDyn'] & flags_dict['flag_DeamDyn_def_transform']:
            print('STATUS:\t Write BeamDyn input files')
            refine = int(30/len(cs_pos))  # initiate node refinement parameter
            write_beamdyn_axis(folder_str, job.yml.get('name'), job.yml.get('components').get('blade'), refine)
            write_beamdyn_prop(folder_str, job.yml.get('name'), cs_pos, beam_prop['beam_stiff'], beam_prop['beam_inertia'])

        # --- Plot VABS 6x6 stiffness and mass matrices --- #
        if flags_dict['flag_plot_vabs_struct_characteristics']:
            # plot inertia matrix
            save_path = [folder_str + job_str[0:-5] + '_VABS_mass_matrix' + str_ext + '.png']
            fig_title = 'Mass matrix'
            plot_beam_props_6by6(cs_pos, beam_prop['beam_inertia'], fig_title, save_path)
            # plot stiffness matrix
            save_path = [folder_str + job_str[0:-5] + '_VABS_stiffness_matrix' + str_ext + '.png']
            fig_title = 'Stiffness matrix'
            plot_beam_props_6by6(cs_pos, beam_prop['beam_stiff'], fig_title, save_path)
            # plot axes locations
            save_path = [folder_str + job_str[0:-5] + '_VABS_axes_locations' + str_ext + '.png']
            plot_beam_axes(cs_pos, beam_prop['beam_mass_center'], beam_prop['beam_neutral_axes'], 
                           beam_prop['beam_geometric_center'], beam_prop['beam_shear_center'], save_path)





    # --- ANBAX --- #
    if flags_dict['flag_run_anbax']:
        # --------------------------------------- #
        # clear var before initializing anbax
        job.beam_properties = None


        # --------------------------------------- #
        # --- ANBAX --- #
        job.blade_run_anbax()  # run anbax

        # init used matrices and arrays
        anbax_beam_stiff_init = np.zeros([len(cs_pos), 6, 6])
        anbax_beam_inertia_init = np.zeros([len(cs_pos), 6, 6])
        anbax_beam_stiff = np.zeros([len(cs_pos), 6, 6])
        anbax_beam_inertia = np.zeros([len(cs_pos), 6, 6])
        anbax_beam_section_mass = np.zeros([len(cs_pos), 1])

        # --------------------------------------- #
        # retrieve & allocate ANBAX results
        for i in range(len(job.beam_properties)):
            anbax_beam_section_mass[i] = job.beam_properties[i, 1].m00  # receive mass per unit span
            for j in range(6):
                anbax_beam_stiff_init[i, j, :] = np.array(job.beam_properties[i, 1].TS[j, :])  # receive 6x6 timoshenko stiffness matrix
                anbax_beam_inertia_init[i, j, :] = np.array(job.beam_properties[i, 1].MM[j, :])  # receive 6x6 mass matrix

        # --------------------------------------- #
        #  rotate anbax results from SONATA/VABS def to BeamDyn def coordinate system (for flag_DeamDyn_def_transform = True)
        if flags_dict['flag_DeamDyn_def_transform']:
            print('STATUS:\t Transform to BeamDyn coordinates')
            # B = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]])  # transformation matrix
            B = np.array([[0, 0, 1], [0, -1, 0], [1, 0, 0]])  # transformation matrix
            T = np.dot(np.identity(3), np.linalg.inv(B))
            for n_sec in range(len(cs_pos)):
                anbax_beam_stiff[n_sec, :, :] = trsf_sixbysix(anbax_beam_stiff_init[n_sec, :, :], T)
                anbax_beam_inertia[n_sec, :, :] = trsf_sixbysix(anbax_beam_inertia_init[n_sec, :, :], T)
            str_ext = '_BeamDyn_def'

            print('STATUS:\t Structural characteristics of ANBAX converted from SONATA/VABS to BeamDyn coordinate system definition!')
        else:
            anbax_beam_stiff = anbax_beam_stiff_init
            anbax_beam_inertia = anbax_beam_inertia_init
            str_ext = ''


        # --------------------------------------- #
        # Export beam structural properties to csv file
        # if flags_dict['flag_csv_export']:
        #     print('STATUS:\t Export csv files with structural blade characeristics from ANBAX to: ' + folder_str + 'csv_export/')
        #     export_beam_struct_properties(folder_str, job_str, cs_pos, solver='anbax', beam_stiff=anbax_beam_stiff, beam_inertia=anbax_beam_inertia, beam_mass_per_length=anbax_beam_section_mass)

        # ToDo: also export BeamDyn files for results from anbax as soon as the verification is completed


    # --------------------------------------- #
    # (Optional) - Analytical Rectangle   #
    # --------------------------------------- #
    # rect_data_analytical = utls_analytical_rectangle()

            

    # --------------------------------------- #
    # (Optional) - Code-to-code comparison)   #
    # --------------------------------------- #
    if flags_dict['flag_verify_vabs_anbax'] & flags_dict['flag_run_vabs'] & flags_dict['flag_run_anbax']:
        # ------------------------------------------------------------------------------------------- #
        # --- Plot VABS and anbax 6x6 stiffness and mass matrices for evaluation and verification --- #
        # plot inertia matrix
        save_path = [folder_str + job_str[0:-5] + '_mass_matrix' + str_ext + '.png']
        fig_title = 'Mass matrix'
        print('STATUS:	 Plot Mass Matrices - Comparison between VABS (black) and ANBAX (red)')
        plot_vabs_anbax(cs_pos, beam_prop['beam_inertia'], anbax_beam_inertia, fig_title, save_path)
        # plot stiffness matrix
        save_path = [folder_str + job_str[0:-5] + '_stiffness_matrix' + str_ext + '.png']
        fig_title = 'Stiffness matrix'
        print('STATUS:	 Plot Stiffness Matrices - Comparison between VABS (black) and ANBAX (red)')
        plot_vabs_anbax(cs_pos, beam_prop['beam_stiff'], anbax_beam_stiff, fig_title, save_path)

    else:
        # optionally account for only analyzing VABS or anbax individually?
        # this would be the frame work to put that capability up on
        Warning('Set flag_verify_vabs_anbax = True for plot output !')



# ============================================= #
def plot_beam_props_6by6(cs_pos, data, fig_title, save_path):
    # plots 6x6 matrix
    k = 1
    fig = plt.figure(tight_layout=True, figsize=(14, 10), dpi=80, facecolor='w', edgecolor='k')
    # fig = plt.figure(tight_layout=True, figsize=(9, 6), dpi=80, facecolor='w', edgecolor='k')
    # fig.suptitle(fig_title)
    for i in range(len(data[0, :, 0])):
        for j in range(len(data[0, 0, :])):
            if j >= i:
                ax = fig.add_subplot(len(data[0, :, 0]), len(data[0, 0, :]), k)
                ax.plot(cs_pos, data[:, i, j], '-k')
                plt.ylim(1.1 * min(-1, min(data[:, i, j])), 1.1 * max(1, max(data[:, i, j])))
                ax.set_xlabel('r/R')
                if fig_title == 'Mass matrix':
                    ax.set_title('$m_{%i %i}$' % ((i + 1), (j + 1)))
                elif fig_title == 'Stiffness matrix':
                    ax.set_title('$k_{%i %i}$' % ((i + 1), (j + 1)))
                ax.grid(True)
            k = k + 1
    plt.show()
    fig.savefig(''.join(save_path), dpi=300)

    return None


def plot_beam_axes(cs_pos, vabs_beam_mass_center, vabs_beam_neutral_axes, vabs_beam_geometric_center,
               vabs_beam_shear_center, save_path):

    fig = plt.figure(tight_layout=True, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    ax = plt.axes(projection='3d')
    ax.plot3D(cs_pos, vabs_beam_mass_center[:,0], vabs_beam_mass_center[:,1], label='Mass center')
    ax.plot3D(cs_pos, vabs_beam_neutral_axes[:,0], vabs_beam_neutral_axes[:,1], label='Neutral axes')
    ax.plot3D(cs_pos, vabs_beam_geometric_center[:,0], vabs_beam_geometric_center[:,1], label='Geometric center')
    ax.plot3D(cs_pos, vabs_beam_shear_center[:,0], vabs_beam_shear_center[:,1], label='Shear Center')
    ax.set_xlabel('r/R')
    ax.set_ylabel('chordwise location, m')
    ax.set_zlabel('thickness location, m')
    ax.legend()

    plt.show()
    fig.savefig(''.join(save_path), dpi=300)

    return None
# ============================================= #
def plot_vabs_anbax(cs_pos, vabs_data, anbax_data, fig_title, save_path):
    # plots 6x6 matrix
    k = 1
    fig = plt.figure(tight_layout=True, figsize=(14, 10), dpi=80, facecolor='w', edgecolor='k')
    # fig.suptitle(fig_title)
    for i in range(len(vabs_data[0, :, 0])):
        for j in range(len(vabs_data[0, 0, :])):
            if j >= i:
                ax = fig.add_subplot(len(vabs_data[0, :, 0]), len(vabs_data[0, 0, :]), k)
                ax.plot(cs_pos, vabs_data[:, i, j], '--k')
                ax.plot(cs_pos, anbax_data[:, i, j], ':r')
                plt.ylim(1.1 * min(-1, min(vabs_data[:, i, j]), min(anbax_data[:, i, j])), 1.1 * max(1, max(vabs_data[:, i, j]), max(anbax_data[:, i, j])))
                ax.set_xlabel('r/R')
                if fig_title == 'Mass matrix':
                    ax.set_title('$m_{%i %i}$' % ((i + 1), (j + 1)))
                elif fig_title == 'Stiffness matrix':
                    ax.set_title('$k_{%i %i}$' % ((i + 1), (j + 1)))
                ax.grid(True)
            k = k + 1
    plt.show()
    fig.savefig(''.join(save_path), dpi=300)

    return None



# ============================================= #
def export_beam_struct_properties(folder_str, job_str, radial_stations, coordsys, solver, beam_stiff, beam_inertia, beam_mass_per_length,
                                  beam_mass_center, beam_neutral_axes, beam_geometric_center, beam_shear_center):

    if solver=='vabs':
        export_name_general = 'vabs_beam_properties_general.csv'
        export_name_stiff = 'vabs_beam_properties_stiff_matrices.csv'
        export_name_mass = 'vabs_beam_properties_mass_matrices.csv'
    elif solver=='anbax':
        export_name_general = 'anbax_beam_properties_general.csv'
        export_name_stiff = 'anbax_beam_properties_stiff_matrices.csv'
        export_name_mass = 'anbax_beam_properties_mass_matrices.csv'
    else:
        print('Define correct solver name (vabs or anbax) when calling export_beam_struct_properties')

    # -------------------------------------------------- #
    # Export mass per unit length for the defined radial stations
    with open(''.join([folder_str + 'csv_export/' + job_str[0:-5] + '_' + export_name_general]), mode='w') as csv_file:
        beam_prop_writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        if coordsys == 'BeamDyn':
            beam_prop_writer.writerow(['Coordinate system:', 'Beamdyn coordinates'])
        elif coordsys == 'VABS/SONATA':
            beam_prop_writer.writerow(['Coordinate system:', 'VABS/SONATA coordinates'])
        else:
            beam_prop_writer.writerow(['Coordinate system:', 'to be verified'])

        beam_prop_writer.writerow(['section in r/R', 'Mass per unit length [kg/m]',
                                      'Mass center (chordwise), m', 'Mass center (thickness), m',
                                      'Neutral axes (chordwise), m', 'Neutral axes (thickness), m',
                                      'Geometric center (chordwise), m', 'Geometric center (thickness), m',
                                      'Shear center (chordwise), m', 'Shear center (thickness), m'])
        for i in range(len(beam_mass_per_length)):  # recieve number of radial sections
            beam_prop_writer.writerow([str(radial_stations[i]), str(beam_mass_per_length[i,0]),
                                       str(beam_mass_center[i,0]), str(beam_mass_center[i,1]),
                                       str(beam_neutral_axes[i,0]), str(beam_neutral_axes[i,1]),
                                       str(beam_geometric_center[i,0]), str(beam_geometric_center[i,1]),
                                       str(beam_shear_center[i,0]), str(beam_shear_center[i,1])])

    csv_file.close()
    # -------------------------------------------------- #

    # Export stiffness matrices for the defined radial stations
    with open(''.join([folder_str + 'csv_export/' + job_str[0:-5] + '_' + export_name_stiff]), mode='w') as csv_file:
        beam_prop_writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        for i in range(len(beam_stiff)):  # receive number of radial sections
            beam_prop_writer.writerow(' ')
            beam_prop_writer.writerow(['section in r/R', str(radial_stations[i])])

            for j in range(6):  # number of rows for each matrix
                beam_prop_writer.writerow(beam_stiff[i, j, :])
                # beam_prop_writer.writerow(job.beam_properties[i, 1].TS[j, :])  # can eventually be called as a standalone via the job.beam_properties object
    csv_file.close()

    # -------------------------------------------------- #
    # Export mass matrices for the defined radial stations
    with open(''.join([folder_str + 'csv_export/' + job_str[0:-5] + '_' + export_name_mass]), mode='w') as csv_file:
        beam_prop_writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        for i in range(len(beam_inertia)):  # receive number of radial sections
            beam_prop_writer.writerow(' ')
            beam_prop_writer.writerow(['section in r/R', str(radial_stations[i])])

            for j in range(6):  # number of rows for each matrix
                beam_prop_writer.writerow(beam_inertia[i, j, :])
                # beam_prop_writer.writerow(job.beam_properties[i, 1].TS[j, :])  # can eventually be called as a standalone via the job.beam_properties object
    csv_file.close()
    # -------------------------------------------------- #
    return None


