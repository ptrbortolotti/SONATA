#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from builtins import len, range
import matplotlib.pyplot as plt
import numpy as np
import csv

def plot_vabs_anbax(cs_pos, vabs_beam_stiff, vabs_beam_inertia, anbax_beam_stiff, anbax_beam_inertia, job_str):


    # ixj (i.e. 6x6) stiffness matrix
    k = 1
    fig = plt.figure(tight_layout=True, figsize=(14,10), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('Stiffness matrix')
    for i in range(len(vabs_beam_stiff[0,:,0])):
        for j in range(len(vabs_beam_stiff[0,0,:])):
            if j >= i:
                ax = fig.add_subplot(len(vabs_beam_stiff[0,:,0]),len(vabs_beam_stiff[0,0,:]),k)
                ax.plot(cs_pos, vabs_beam_stiff[:,i,j])
                ax.plot(cs_pos, anbax_beam_stiff[:,i,j])

                ax.set_xlabel('r/R')
                ax.set_title('$k_{%i %i}$' %((i+1), (j+1)))
                ax.grid(True)

            k = k+1
    plt.show()
    # job_str
    # plt.savefig([folder_str + '/yaml_examples/' + job_str[14:-5] + 'stiffness.png'], dpi=150)
    fig.savefig('yaml_examples/' + job_str[14:-5] + '_stiffness_matrix.png', dpi=150)

    # ixj (i.e. 6x6) mass matrix
    m = 1
    fig = plt.figure(tight_layout=True, figsize=(14,10), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('Mass matrix')
    for i in range(len(vabs_beam_inertia[0,:,0])):
        for j in range(len(vabs_beam_inertia[0,0,:])):
            if j >= i:
                ax = fig.add_subplot(len(vabs_beam_inertia[0,:,0]),len(vabs_beam_inertia[0,0,:]),m)
                ax.plot(cs_pos, vabs_beam_inertia[:,i,j])
                ax.plot(cs_pos, anbax_beam_inertia[:,i,j])

                ax.set_xlabel('r/R')
                ax.set_title('$m_{%i %i}$' %((i+1), (j+1)))
                ax.grid(True)

            m = m+1
    plt.show()
    fig.savefig('yaml_examples/' + job_str[14:-5] + '_mass_matrix.png', dpi=150)
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


# ==============
# Main
# ==============
if __name__ == '__main__':
    import matplotlib.pyplot as plt