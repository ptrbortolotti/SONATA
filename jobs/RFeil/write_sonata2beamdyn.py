# -*- coding: utf-8 -*-
"""
Created on Thursday Oct 10 10:52:28 2019

@author: Roland Feil
"""
# ============================================= #
"""
Convert SONATA structural results to OpenFAST input files

Functions:
write_beamdyn_axis  - writes '*_BeamDyn.dat' file
write_beamdyn_prop  - writes '*_BeamDyn_Blade.dat' file

Inputs:
wt_name         - name of currently investigated concept; i.e. wind turbine 
byml            - yml.get('components').get('blade')
radial_stations - radial stations along the blade span for inertia and stiffness characteristics
beam_stiff      - matrix containing beam stiffness values
beam_inertia    - matrix containing beam inertia values

Outputs: 
- None -

"""
# ============================================= #


import os
import numpy as np
from scipy.interpolate import interp1d


if __name__ == '__main__':
    os.chdir('../..')

# --- Write BeamDyn file with blade reference line locations ---#
def write_beamdyn_axis(folder, wt_name, byml):
    yaml_ref_axis = byml.get('internal_structure_2d_fem').get('reference_axis')
    yaml_twist = byml.get('outer_shape_bem').get('twist')

    n_pts = 50
    grid = np.linspace(0, 1, n_pts)


    # >>> The following is defined correctly when WT definition is used within yaml file <<<
    # The x of yaml corresponds to x in BeamDyn
    f_interp = interp1d(yaml_ref_axis.get('x').get('grid'), yaml_ref_axis.get('x').get('values'))
    kp_xr = f_interp(grid)
    # The y of yaml corresponds to y in BeamDyn
    f_interp = interp1d(yaml_ref_axis.get('y').get('grid'), yaml_ref_axis.get('y').get('values'))
    kp_yr = f_interp(grid)
    # The z of yaml corresponds to z in BeamDyn
    f_interp = interp1d(yaml_ref_axis.get('z').get('grid'), yaml_ref_axis.get('z').get('values'))
    kp_zr = f_interp(grid)
    # Twist
    f_interp = interp1d(yaml_twist.get('grid'), yaml_twist.get('values'))
    twist    = f_interp(grid)*180./np.pi
    # twist = np.zeros_like(kp_yr)  # twist defined to zero

    data = np.vstack((kp_xr, kp_yr, kp_zr, twist)).T

    file = open(folder + wt_name + '_BeamDyn.dat', 'w')
    file.write('--------- BEAMDYN with OpenFAST INPUT FILE -------------------------------------------\n')
    file.write('%s blade\n' % (wt_name))
    file.write('---------------------- SIMULATION CONTROL --------------------------------------\n')
    file.write('False         Echo            - Echo input data to "<RootName>.ech" (flag)\n')
    file.write('True          QuasiStaticInit - Use quasistatic pre-conditioning with centripetal accelerations in initialization (flag) [dynamic solve only]\n')
    file.write(' 0            rhoinf          - Numerical damping parameter for generalized-alpha integrator\n')
    file.write(' 2            quadrature      - Quadrature method: 1=Gaussian; 2=Trapezoidal (switch)\n')
    file.write('DEFAULT       refine          - Refinement factor for trapezoidal quadrature (-). DEFAULT = 1 [used only when quadrature=2]\n')
    file.write('DEFAULT       n_fact          - Factorization frequency (-). DEFAULT = 5\n')
    file.write('DEFAULT       DTBeam          - Time step size (s).\n')
    file.write('DEFAULT       load_retries    - Number of factored load retries before quitting the aimulation\n')
    file.write('DEFAULT       NRMax           - Max number of iterations in Newton-Ralphson algorithm (-). DEFAULT = 10\n')
    file.write('DEFAULT       stop_tol        - Tolerance for stopping criterion (-)\n')
    file.write('DEFAULT       tngt_stf_fd     - Flag to use finite differenced tangent stiffness matrix (-)\n')
    file.write('DEFAULT       tngt_stf_comp   - Flag to compare analytical finite differenced tangent stiffness matrix  (-)\n')
    file.write('DEFAULT       tngt_stf_pert   - perturbation size for finite differencing (-)\n')
    file.write('DEFAULT       tngt_stf_difftol- Maximum allowable relative difference between analytical and fd tangent stiffness (-)\n')
    file.write('True          RotStates       - Orient states in the rotating frame during linearization? (flag) [used only when linearizing]\n')
    file.write('---------------------- GEOMETRY PARAMETER --------------------------------------\n')
    file.write('          1   member_total    - Total number of members (-)\n')
    file.write('         %u   kp_total        - Total number of key points (-) [must be at least 3]\n' % (n_pts))
    file.write('     1     %u                 - Member number; Number of key points in this member\n' % (n_pts))
    file.write('\t\t kp_xr \t\t\t kp_yr \t\t\t kp_zr \t\t initial_twist\n')
    file.write('\t\t  (m)  \t\t\t  (m)  \t\t\t  (m)  \t\t   (deg)\n')
    for i in range(n_pts):
        file.write('\t %.5e \t %.5e \t %.5e \t %.5e \n' % (data[i, 0], data[i, 1], data[i, 2], data[i, 3]))
    file.write('---------------------- MESH PARAMETER ------------------------------------------\n')
    file.write('          5   order_elem     - Order of interpolation (basis) function (-)\n')
    file.write('---------------------- MATERIAL PARAMETER --------------------------------------\n')
    file.write('"%s"    BldFile - Name of file containing properties for blade (quoted string)\n' % (wt_name + '_BeamDyn_Blade.dat'))
    file.write('---------------------- PITCH ACTUATOR PARAMETERS -------------------------------\n')
    file.write('False         UsePitchAct - Whether a pitch actuator should be used (flag)\n')
    file.write('        200   PitchJ      - Pitch actuator inertia (kg-m^2) [used only when UsePitchAct is true]\n')
    file.write('      2E+07   PitchK      - Pitch actuator stiffness (kg-m^2/s^2) [used only when UsePitchAct is true]\n')
    file.write('     500000   PitchC      - Pitch actuator damping (kg-m^2/s) [used only when UsePitchAct is true]\n')
    file.write('---------------------- OUTPUTS -------------------------------------------------\n')
    file.write('True          SumPrint       - Print summary data to "<RootName>.sum" (flag)\n')
    file.write('"ES10.3E2"    OutFmt         - Format used for text tabular output, excluding the time channel.\n')
    file.write('          0   NNodeOuts      - Number of nodes to output to file [0 - 9] (-)\n')
    file.write('          1,          2,          3,          4,          5,          6    OutNd          - Nodes whose values will be output  (-)\n')
    file.write('          OutList            - The next line(s) contains a list of output parameters. See OutListParameters.xlsx for a listing of available output channels, (-)\n')
    file.write('"RootFxr, RootFyr, RootFzr"\n')
    file.write('"RootMxr, RootMyr, RootMzr"\n')
    file.write('"TipTDxr, TipTDyr, TipTDzr"\n')
    file.write('"TipRDxr, TipRDyr, TipRDzr"\n')
    file.write('END of input file (the word "END" must appear in the first 3 columns of this last OutList line)\n')
    file.write('---------------------------------------------------------------------------------------\n')

    file.close()

    print('Finished writing BeamDyn File')

    return None

# --- Write BeamDyn_Blade file with blade properties ---#
def write_beamdyn_prop(folder, wt_name, radial_stations, beam_stiff, beam_inertia):
    n_pts = len(radial_stations)

    file = open(folder + wt_name + '_BeamDyn_Blade.dat', 'w')
    file.write(' ------- BEAMDYN V1.00.* INDIVIDUAL BLADE INPUT FILE --------------------------\n')
    file.write(' Test Format 1\n')
    file.write(' ---------------------- BLADE PARAMETERS --------------------------------------\n')
    file.write('%u   station_total    - Number of blade input stations (-)\n' % (n_pts))
    file.write(' 1   damp_type        - Damping type: 0: no damping; 1: damped\n')
    file.write('  ---------------------- DAMPING COEFFICIENT------------------------------------\n')
    file.write('   mu1        mu2        mu3        mu4        mu5        mu6\n')
    file.write('   (-)        (-)        (-)        (-)        (-)        (-)\n')
    file.write('1.0E-02    1.0E-02    1.0E-02    1.0E-02    1.0E-02    1.0E-02\n')
    file.write(' ---------------------- DISTRIBUTED PROPERTIES---------------------------------\n')
    for i in range(n_pts):
        file.write('\t %.6f \n' % (radial_stations[i]))
        for j in range(6):
            file.write('\t %.16e \t %.16e \t %.16e \t %.16e \t %.16e \t %.16e\n' % (
            beam_stiff[i, j, 0], beam_stiff[i, j, 1], beam_stiff[i, j, 2], beam_stiff[i, j, 3], beam_stiff[i, j, 4],
            beam_stiff[i, j, 5]))
        file.write('\n')
        for j in range(6):
            file.write('\t %.16e \t %.16e \t %.16e \t %.16e \t %.16e \t %.16e\n' % (
            beam_inertia[i, j, 0], beam_inertia[i, j, 1], beam_inertia[i, j, 2], beam_inertia[i, j, 3],
            beam_inertia[i, j, 4], beam_inertia[i, j, 5]))
        file.write('\n')
        # ToDO: check correct translation of stiffness and mass matrices from VABS and anbax !!!
    file.close()

    print('Finished writing BeamDyn_Blade File')

    return None


# ==============
# Main
# ==============

if __name__ == '__main__':
    import yaml
    from SONATA.classAirfoil import Airfoil
    #from SONATA.classBlade import Blade

    from SONATA.classAirfoil import Airfoil
    from SONATA.classMaterial import read_IEA37_materials

    # provide primary path; used for providing yaml input file as well as output directory
    folder = '/Users/rfeil/work/6_SONATA/SONATA/jobs/RFeil/'
    filename = (folder + 'IEAonshoreWT_BAR_005a.yaml')
    with open(filename, 'r') as myfile:
        inputs = myfile.read()
        yml = yaml.load(inputs, Loader=yaml.FullLoader)

    airfoils = [Airfoil(af) for af in yml.get('airfoils')]
    materials = read_IEA37_materials(yml.get('materials'))
    # Blade.read_IEA37(yml.get('components').get('blade'), airfoils, **kwargs)
    wt_name = yml.get('name')

    # radial_stations = [0.0, 0.5, 1.0] # must include 0 and 1!
    radial_stations = np.linspace(0.0, 1.0, 11)

    beam_stiff = np.zeros([len('radial_stations'), 6, 6])
    beam_inertia = np.zeros([len('radial_stations'), 6, 6])

    write_beamdyn_axis(folder, wt_name, yml.get('components').get('blade'))
    write_beamdyn_prop(folder, wt_name, radial_stations, beam_stiff, beam_inertia)

