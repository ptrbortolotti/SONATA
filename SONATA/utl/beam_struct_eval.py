#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wednesday Nov 20 13:33:33 2019

@author: Roland Feil
"""


from builtins import len, range
import matplotlib.pyplot as plt
import numpy as np
import csv
import os


# from SONATA.utl.analytical_rectangle.utls_analytical_rectangle import utls_analytical_rectangle



def beam_struct_eval(flags_dict, loads_dict, cs_pos, job, folder_str, job_str):

    """
    Analyse, transform, evaluate and plot structural results from ANBAX

    Functions:
    beam_struct_eval                    - parent function (retrieve & transform data, structure of data evaluation)
    plot_beam_props_6by6                - function to plot the 6x6 stiffness and mass matrices
    plot_beam_mass_distribution         - function to plot the beam mass per unit length distribution
    plot_vabs_anbax                     - function to plot the 6x6 stiffness and mass matrices from ANBAX for code-to-code verification
    anbax_export_beam_struct_properties - csv export of structural beam properties

    Inputs:
    flags_dict          - dictionary containing relevant flags
    loads_dict          - dictionary containing the applied loads for recovery analysis
    cs_pos              - radial station of blade cross sections
    job                 - contains the whole blade data (yaml file content, wires, mesh, etc.)
    folder_str          - name of operating folder
    job_str             - name of job that is currently under investigation 

    Outputs: 
    - None -


    Dictionary Examples:
    flag_recovery           = True
    flags_dict = {"flag_recovery": flag_recovery}

    Forces =  np.array([0, 0, 0])  # forces, N (F1: axial force; F2,F3: sectional transverse shear forces)
    Moments =  np.array([0, 5000000, 0])  # moments, Nm
    loads_dict = {"Forces": Forces, "Moments": Moments}

    """
    # --- ANBAX --- #
    # --------------------------------------- #
    # clear var before initializing anbax
    job.beam_properties = None

    if flags_dict['flag_recovery'] == True:
        loads = {  # sonata coord system input converted to anbax coordinates
            "F":    np.array([[0, float(loads_dict["Forces"][1]), float(loads_dict["Forces"][2]), float(loads_dict["Forces"][0])],    [1, float(loads_dict["Forces"][1]), float(loads_dict["Forces"][2]), float(loads_dict["Forces"][0])]]),  # forces, N (F1: shear force in x-direction; F2: shear force in y -direction; F3: axial force)
            "M":    np.array([[0, float(loads_dict["Moments"][1]), float(loads_dict["Moments"][2]), float(loads_dict["Moments"][0])],    [1, float(loads_dict["Moments"][1]), float(loads_dict["Moments"][2]), float(loads_dict["Moments"][0])]])}  # moments, Nm (M1: bending moment around x; M2: bending moment around y; M3: torsional moment)

        job.blade_run_anbax(loads)  # run anbax
    else:
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
    # Export beam structural properties to csv file
    if flags_dict['flag_csv_export']:
        print('STATUS:\t Export csv files with structural blade characeristics from ANBAX to: ' + folder_str + 'csv_export/')
        anbax_export_beam_struct_properties(folder_str, job_str, cs_pos, solver='anbax', beam_stiff=anbax_beam_stiff,
                                        beam_inertia=anbax_beam_inertia, beam_mass_per_length=anbax_beam_section_mass)



    # --------------------------------------- #
    # Export beam structural properties to csv file
    # if flags_dict['flag_csv_export']:
    #     print('STATUS:\t Export csv files with structural blade characeristics from ANBAX to: ' + folder_str + 'csv_export/')
    #     export_beam_struct_properties(folder_str, job_str, cs_pos, solver='anbax', beam_stiff=anbax_beam_stiff, beam_inertia=anbax_beam_inertia, beam_mass_per_length=anbax_beam_section_mass)

    # ToDo: also export BeamDyn files for results from anbax as soon as the verification is completed
    # --------------------------------------- #
    # write BeamDyn input files
    np.savetxt('anbax_BAR00.txt', np.array([cs_pos, anbax_beam_stiff[:, 3, 3], anbax_beam_stiff[:, 4, 4], anbax_beam_stiff[:, 5, 5], anbax_beam_stiff[:, 2, 2], anbax_beam_inertia[:, 0, 0]]).T)

    # --------------------------------------- #
    # (Optional) - Analytical Rectangle   #
    # --------------------------------------- #
    # rect_data_analytical = utls_analytical_rectangle()



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





# ============================================= #
def plot_vabs_anbax(cs_pos, anbax_data, fig_title, save_path):
    # plots 6x6 matrix
    k = 1
    fig = plt.figure(tight_layout=True, figsize=(14, 10), dpi=80, facecolor='w', edgecolor='k')
    # fig.suptitle(fig_title)
    for i in range(len(anbax_data[0, :, 0])):
        for j in range(len(anbax_data[0, 0, :])):
            if j >= i:
                ax = fig.add_subplot(len(anbax_data[0, :, 0]), len(anbax_data[0, 0, :]), k)
                ax.plot(cs_pos, anbax_data[:, i, j], ':r')
                plt.ylim(1.1 * min(-1, min(anbax_data[:, i, j]), min(anbax_data[:, i, j])), 1.1 * max(1, max(anbax_data[:, i, j]), max(anbax_data[:, i, j])))
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

def anbax_export_beam_struct_properties(folder_str, job_str, radial_stations, solver, beam_stiff, beam_inertia, beam_mass_per_length):

    export_name_general = 'anbax_beam_properties_general.csv'
    export_name_stiff = 'anbax_beam_properties_stiff_matrices.csv'
    export_name_mass = 'anbax_beam_properties_mass_matrices.csv'

    # -------------------------------------------------- #
    # Export mass per unit length for the defined radial stations
    if os.path.isdir(folder_str + 'csv_export/') == False:
        os.mkdir(folder_str + 'csv_export/')
    with open(''.join([folder_str + 'csv_export/' + job_str[0:-5] + '_' + export_name_general]), mode='w') as csv_file:
        beam_prop_writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        beam_prop_writer.writerow(['Coordinate system:', 'VABS/SONATA coordinates'])

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
    return None
