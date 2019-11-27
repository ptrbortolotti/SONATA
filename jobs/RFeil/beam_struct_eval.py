#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wednesday Nov 20 13:33:33 2019

@author: Roland Feil
"""
# ============================================= #
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
# ============================================= #

from builtins import len, range
import matplotlib.pyplot as plt
import numpy as np
import csv

from SONATA.cbm.cbm_utl import trsf_sixbysix
from jobs.RFeil.write_sonata2beamdyn import write_beamdyn_axis, write_beamdyn_prop



def beam_struct_eval(flags_dict, cs_pos, job, folder_str, job_str):

    if flags_dict['flag_verify_vabs_anbax']:
        # ------------------------ #
        # runs both VABS & ANBAX and conducts a verification/ comparison of the results


        # --- VABS --- #
        job.blade_run_vabs()  # run VABS

        # init used matrices and arrays
        vabs_beam_stiff_init = np.zeros([len(cs_pos), 6, 6])
        vabs_beam_inertia_init = np.zeros([len(cs_pos), 6, 6])
        vabs_beam_stiff = np.zeros([len(cs_pos), 6, 6])
        vabs_beam_inertia = np.zeros([len(cs_pos), 6, 6])
        vabs_beam_section_mass = np.zeros([len(cs_pos), 1])
        vabs_beam_mass_center = np.zeros([len(cs_pos), 2])
        vabs_beam_neutral_axes = np.zeros([len(cs_pos), 2])
        vabs_beam_geometric_center = np.zeros([len(cs_pos), 2])
        vabs_beam_shear_center = np.zeros([len(cs_pos), 2])

        # --------------------------------------- #
        # retrieve & allocate VABS results
        for i in range(len(job.beam_properties)):
            vabs_beam_section_mass[i] = job.beam_properties[i, 1].m00  # mass per unit span
            vabs_beam_mass_center[i] = np.array(job.beam_properties[i, 1].Xm[:])  # Center of gravity (mass center)
            vabs_beam_neutral_axes[i] = np.array(job.beam_properties[i, 1].Xt[:])  # Neutral axes (tension center)
            vabs_beam_geometric_center[i] = np.array(job.beam_properties[i, 1].Xg[:])  # Geometric center
            vabs_beam_shear_center[i] = np.array(job.beam_properties[i, 1].Xs[:])  # Generalized Shear Center of the Cross Section
            for j in range(6):
                vabs_beam_stiff_init[i, j, :] = np.array(job.beam_properties[i, 1].TS[j, :])  # receive 6x6 timoshenko stiffness matrix
                vabs_beam_inertia_init[i, j, :] = np.array(job.beam_properties[i, 1].MM[j, :])  # receive 6x6 mass matrix
                # vabs_beam_inertia_init[i, j, :] = np.array(job.beam_properties[i, 1].MMatMC[j, :])  # receive 6x6 mass matrix at mass center

        # --------------------------------------- #
        #  rotate VABS results from SONATA/VABS def to BeamDyn def coordinate system (for flag_DeamDyn_def_transform = True)
        if flags_dict['flag_DeamDyn_def_transform']:
            print('STATUS:\t Transform to BeamDyn coordinates')
            B = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]])  # transformation matrix
            T = np.dot(np.identity(3), np.linalg.inv(B))
            for n_sec in range(len(cs_pos)):
                vabs_beam_stiff[n_sec, :, :] = trsf_sixbysix(vabs_beam_stiff_init[n_sec, :, :], T)
                vabs_beam_inertia[n_sec, :, :] = trsf_sixbysix(vabs_beam_inertia_init[n_sec, :, :], T)
            str_ext = '_BeamDyn_def'
        else:
            vabs_beam_stiff = vabs_beam_stiff_init
            vabs_beam_inertia = vabs_beam_inertia_init
            str_ext = ''

        # --------------------------------------- #
        # Export beam structural properties to csv file
        if flags_dict['flag_csv_export']:
            print('STATUS:\t Export csv files with structural blade characeristics from VABS to: ' + folder_str + 'csv_export/')
            export_beam_struct_properties(folder_str, job_str, cs_pos, solver='vabs', beam_stiff=vabs_beam_stiff, beam_inertia=vabs_beam_inertia, beam_mass_per_length=vabs_beam_section_mass)

        # --------------------------------------- #
        # write BeamDyn input files
        if flags_dict['flag_write_BeamDyn'] & flags_dict['flag_DeamDyn_def_transform']:

            print('STATUS:\t Write BeamDyn input files')
            write_beamdyn_axis(folder_str, job.yml.get('name'), job.yml.get('components').get('blade'))
            write_beamdyn_prop(folder_str, job.yml.get('name'), cs_pos, vabs_beam_stiff, vabs_beam_inertia)

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
        # retrieve & allocate VABS results
        for i in range(len(job.beam_properties)):
            anbax_beam_section_mass[i] = job.beam_properties[i, 1].m00  # receive mass per unit span
            for j in range(6):
                anbax_beam_stiff_init[i, j, :] = np.array(job.beam_properties[i, 1].TS[j, :])  # receive 6x6 timoshenko stiffness matrix
                anbax_beam_inertia_init[i, j, :] = np.array(job.beam_properties[i, 1].MM[j, :])  # receive 6x6 mass matrix

        # --------------------------------------- #
        #  rotate anbax results from SONATA/VABS def to BeamDyn def coordinate system (for flag_DeamDyn_def_transform = True)
        if flags_dict['flag_DeamDyn_def_transform']:
            print('STATUS:\t Transform to BeamDyn coordinates')
            B = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]])  # transformation matrix
            T = np.dot(np.identity(3), np.linalg.inv(B))
            for n_sec in range(len(cs_pos)):
                anbax_beam_stiff[n_sec, :, :] = trsf_sixbysix(anbax_beam_stiff_init[n_sec, :, :], T)
                anbax_beam_inertia[n_sec, :, :] = trsf_sixbysix(anbax_beam_inertia_init[n_sec, :, :], T)
            str_ext = '_BeamDyn_def'
        else:
            anbax_beam_stiff = anbax_beam_stiff_init
            anbax_beam_inertia = anbax_beam_inertia_init
            str_ext = ''


        # --------------------------------------- #
        # Export beam structural properties to csv file
        if flags_dict['flag_csv_export']:
            print('STATUS:\t Export csv files with structural blade characeristics from ANBAX to: ' + folder_str + 'csv_export/')
            export_beam_struct_properties(folder_str, job_str, cs_pos, solver='anbax', beam_stiff=anbax_beam_stiff, beam_inertia=anbax_beam_inertia, beam_mass_per_length=anbax_beam_section_mass)

        # ToDo: also export BeamDyn files for results from anbax as soon as the verification is completed


        # ------------------------------------------------------------------------------------------- #
        # --- Plot VABS and anbax 6x6 stiffness and mass matrices for evaluation and verification --- #
        # plot inertia matrix
        save_path = [folder_str + job_str[0:-5] + '_mass_matrix' + str_ext + '.png']
        fig_title = 'Mass matrix'
        plot_vabs_anbax(cs_pos, vabs_beam_inertia, anbax_beam_inertia, fig_title, save_path)
        # plot stiffness matrix
        save_path = [folder_str + job_str[0:-5] + '_stiffness_matrix' + str_ext + '.png']
        fig_title = 'Stiffness matrix'
        plot_vabs_anbax(cs_pos, vabs_beam_stiff, anbax_beam_stiff, fig_title, save_path)




    else:
        # optionally account for only analyzing VABS or anbax individually?
        # this would be the frame work to put that capability up on
        Warning('Set flag_verify_vabs_anbax = True for plot output !')



def plot_vabs_anbax(cs_pos, vabs_data, anbax_data, fig_title, save_path):
    # plots 6x6 matrix
    k = 1
    fig = plt.figure(tight_layout=True, figsize=(14, 10), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle(fig_title)
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
    fig.savefig(''.join(save_path), dpi=150)

    return None




def export_beam_struct_properties(folder_str, job_str, radial_stations, solver, beam_stiff, beam_inertia, beam_mass_per_length):

    if solver=='vabs':
        export_name_mass_length = 'vabs_beam_properties_mass_per_unit_length.csv'
        export_name_stiff = 'vabs_beam_properties_stiff_matrices.csv'
        export_name_mass = 'vabs_beam_properties_mass_matrices.csv'
    elif solver=='anbax':
        export_name_mass_length = 'anbax_beam_properties_mass_per_unit_length.csv'
        export_name_stiff = 'anbax_beam_properties_stiff_matrices.csv'
        export_name_mass = 'anbax_beam_properties_mass_matrices.csv'
    else:
        print('Define correct solver name (vabs or anbax) when calling export_beam_struct_properties')

    # -------------------------------------------------- #
    # Export mass per unit length for the defined radial stations
    with open(''.join([folder_str + 'csv_export/' + job_str[0:-5] + '_' + export_name_mass_length]), mode='w') as csv_file:
        beam_prop_writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        beam_prop_writer.writerow(['section in r/R', 'Mass per unit length [kg/m]'])
        for i in range(len(beam_mass_per_length)):  # receive number of radial sections
            beam_prop_writer.writerow([str(radial_stations[i]), str(beam_mass_per_length[i,0])])


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
