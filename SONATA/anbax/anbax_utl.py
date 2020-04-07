
# Third party modules
import numpy as np

try:
    import dolfin
    from anba4 import material
except ImportError as error:
    print(error.__class__.__name__ + ": " + error.message)
except Exception as exception:
    print(exception, False)

def build_mat_library(cbm_materials):

    matLibrary = []

    maxE = 0.
    matid = 0
    matdict = {}
    for m in cbm_materials.values():
        if m.orth == 0:
            maxE = max(m.E, maxE)
            matMechanicProp = [m.E, m.nu]
            mat = material.IsotropicMaterial(matMechanicProp, m.rho)
        elif m.orth == 1:
            matMechanicProp = np.zeros((3,3))
            maxE = max(m.E[0], maxE)
            maxE = max(m.E[1], maxE)
            maxE = max(m.E[2], maxE)
            matMechanicProp[0,0] = m.E[1] #Exx
            matMechanicProp[0,1] = m.E[2] #Eyy
            matMechanicProp[0,2] = m.E[0] #Ezz
            matMechanicProp[1,0] = m.G[1] #g_yz
            matMechanicProp[1,1] = m.G[0] #g_xz
            matMechanicProp[1,2] = m.G[2] #g_xy
            matMechanicProp[2,0] = m.nu[1] #nu_zy
            matMechanicProp[2,1] = m.nu[0] #nu_zx
            matMechanicProp[2,2] = m.nu[2] #nu_xy
            mat = material.OrthotropicMaterial(matMechanicProp, m.rho)
        elif m.orth == 2:
            raise ValueError('material type 2 (anysotropic) not supported by Anba')
            
            #mat = material.
        matLibrary.append(mat)
        matdict[m.id] = matid
        matid = matid + 1
    return matLibrary, matdict, maxE

def build_dolfin_mesh(cbm_mesh, cbm_nodes, cbm_materials):
    """function to generate the dolfin.Mesh from a SONATA-CBM definition to run
    with anbax

    Parameters
    ----------
    cbm_mesh : list of cell instances
        from the SONATA-CBM preprocessor

    cbm_nodes : list of nodes
        from the SONATA-CBM preprocessor


    Returns
    ----------
    mesh : dolfin.Mesh
    matLibrary : vector of anbax materials
    materials : dolfin.MeshFunction definign cell materials
    plane_orientations : dolfin.MeshFunction defining cell plane orientations
    fiber_orientations : dolfin.MeshFunction defining cell material fiber orientation
    maxE : reference elastic modulus for scaling rigid mode constraints


    Notes
    ----------
    the cells of cbm_mesh already contain the nodes. So the information is 
    currently passed twice. But consistent with the export_cells_for_vabs.

    """
    (matLibrary, matdict, maxE) = build_mat_library(cbm_materials)

    mesh = dolfin.Mesh()
    me = dolfin.MeshEditor()
    me.open(mesh, "triangle", 2, 2)

    me.init_vertices(len(cbm_nodes))
    for n in cbm_nodes:
        me.add_vertex(n.id-1, n.coordinates)

    me.init_cells(len(cbm_mesh))
    for c in cbm_mesh:
        me.add_cell(c.id-1, [n.id-1 for n in c.nodes])

    me.close()

    materials = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim())
    fiber_orientations = dolfin.MeshFunction("double", mesh, mesh.topology().dim())
    plane_orientations = dolfin.MeshFunction("double", mesh, mesh.topology().dim())
    
    materials.set_all(0)
    plane_orientations.set_all(0.0)
    fiber_orientations.set_all(0.0)
    
    for c in cbm_mesh:
        materials[c.id-1] = matdict[c.MatID]
        plane_orientations[c.id-1] = c.theta_1[0]  # rotation around x1-axis (equiv. to beam axis) in SONATA/VABS coordinates; Theta_11
        fiber_orientations[c.id-1] = c.theta_3     # rotation around x3-axis in SONATA/VABS coordinates; Theta_3


    return mesh, matLibrary, materials, plane_orientations, fiber_orientations, maxE


def anbax_recovery(anba, n_el, force, moment, voigt_convention, T):
    """
    Function to recover stresses and strains from an applied loading
    Results are generated in the global and local ('M' for material coordinate system) coordinates

    INPUTS:
    anba    -   dolfin construct from anbax
    n_el    -   number of mesh elements
    force   -   Forces in anbax coordinates, [F1, F2, F3], e.g. force = [2.2, 3.4, 1.1]
    moment  -   Moments in anbax coordinates, [M1, M2, M3], e.g. moment = [4.2, 5.7, 6.2]
    voigt_convention    -   "anba" with [s_xx, s_yy, s_zz, s_yz, s_xz, s_xy] or "paraview" with [s_xx, s_yy, s_zz, s_xy, s_yz, s_xz]
    T       -   Transformation matrix to convert results from ANBA to SONATA/VABS coordinates

    OUTPUTS:
    *_tran issues that outputs were converted to the SONATA/VABS coordinates
    remove '_tran' from the following output names when needed in ANBA coordinates

    tmp_StressF_tran    -   global stress field
    tmp_StressF_M_tran  -   local stress field
    tmp_StrainF_tran    -   global strain field
    tmp_StrainF_M_tran  -   local strain field



    """


    # Stress field
    anba.stress_field(force, moment, reference="global", voigt_convention=voigt_convention)  # get stress field in global sys
    tmp_StressF_vec = np.array(anba.STRESS.vector().vec())  # global stress field
    anba.stress_field(force, moment, reference="local", voigt_convention=voigt_convention)  # get stress field in local sys (material coordinates)
    tmp_StressF_M_vec = np.array(anba.STRESS.vector().vec())  # local stress field

    # Strain field
    anba.strain_field(force, moment, reference="global", voigt_convention=voigt_convention)  # get strain field in global sys
    tmp_StrainF_vec = np.array(anba.STRAIN.vector().vec())  # global strain field
    anba.strain_field(force, moment, reference="local", voigt_convention=voigt_convention)  # get strain field in local sys (material coordinates)
    tmp_StrainF_M_vec = np.array(anba.STRAIN.vector().vec())  # local strain field

    # cd = anba.STRESS.function_space().dofmap().cell_dofs  # index numbers of cells from dolfin mesh that was used for stress recovery (each cell has 6 dofs)

    s_11 = np.zeros(n_el)
    s_22 = np.zeros(n_el)
    s_33 = np.zeros(n_el)
    s_23 = np.zeros(n_el)
    s_13 = np.zeros(n_el)
    s_12 = np.zeros(n_el)
    tmp_StressF = np.zeros((n_el, 3, 3))
    tmp_StressF_tran = np.zeros((n_el, 3, 3))

    s_11_M = np.zeros(n_el)
    s_22_M = np.zeros(n_el)
    s_33_M = np.zeros(n_el)
    s_23_M = np.zeros(n_el)
    s_13_M = np.zeros(n_el)
    s_12_M = np.zeros(n_el)
    tmp_StressF_M = np.zeros((n_el, 3, 3))
    tmp_StressF_M_tran = np.zeros((n_el, 3, 3))

    e_11 = np.zeros(n_el)
    e_22 = np.zeros(n_el)
    e_33 = np.zeros(n_el)
    e_23 = np.zeros(n_el)
    e_13 = np.zeros(n_el)
    e_12 = np.zeros(n_el)
    tmp_StrainF = np.zeros((n_el, 3, 3))
    tmp_StrainF_tran = np.zeros((n_el, 3, 3))

    e_11_M = np.zeros(n_el)
    e_22_M = np.zeros(n_el)
    e_33_M = np.zeros(n_el)
    e_23_M = np.zeros(n_el)
    e_13_M = np.zeros(n_el)
    e_12_M = np.zeros(n_el)
    tmp_StrainF_M = np.zeros((n_el, 3, 3))
    tmp_StrainF_M_tran = np.zeros((n_el, 3, 3))

    # cell_id = np.zeros((n_el, 6))

    if voigt_convention == "anba":  # [s_xx, s_yy, s_zz, s_yz, s_xz, s_xy]
        for i in range(n_el):
            # stresses in "global" system
            s_11[i] = tmp_StressF_vec[i * 6]
            s_22[i] = tmp_StressF_vec[i * 6 + 1]
            s_33[i] = tmp_StressF_vec[i * 6 + 2]
            s_23[i] = tmp_StressF_vec[i * 6 + 3]  # equiv to s_23
            s_13[i] = tmp_StressF_vec[i * 6 + 4]  # equiv to s_31
            s_12[i] = tmp_StressF_vec[i * 6 + 5]  # equiv to s_21
            tmp_StressF[i, :, :] = np.array([[s_11[i], s_12[i], s_13[i]], [s_12[i], s_22[i], s_23[i]], [s_13[i], s_23[i], s_33[i]]])
            tmp_StressF_tran[i, :, :] = np.dot(np.dot(T.T, tmp_StressF[i]), T)  # transform to sonata coordinate system
            # stresses in "local" system
            s_11_M[i] = tmp_StressF_M_vec[i * 6]
            s_22_M[i] = tmp_StressF_M_vec[i * 6 + 1]
            s_33_M[i] = tmp_StressF_M_vec[i * 6 + 2]
            s_23_M[i] = tmp_StressF_M_vec[i * 6 + 3]  # equiv to s_23_M
            s_13_M[i] = tmp_StressF_M_vec[i * 6 + 4]  # equiv to s_31_M
            s_12_M[i] = tmp_StressF_M_vec[i * 6 + 5]  # equiv to s_21_M
            tmp_StressF_M[i, :, :] = np.array([[s_11_M[i], s_12_M[i], s_13_M[i]], [s_12_M[i], s_22_M[i], s_23_M[i]], [s_13_M[i], s_23_M[i], s_33_M[i]]])
            tmp_StressF_M_tran[i, :, :] = np.dot(np.dot(T.T, tmp_StressF_M[i]), T)  # transform to sonata coordinate system

            # strains in "global" system
            e_11[i] = tmp_StrainF_vec[i * 6]
            e_22[i] = tmp_StrainF_vec[i * 6 + 1]
            e_33[i] = tmp_StrainF_vec[i * 6 + 2]
            e_23[i] = tmp_StrainF_vec[i * 6 + 3]  # equiv to e_23
            e_13[i] = tmp_StrainF_vec[i * 6 + 4]  # equiv to e_31
            e_12[i] = tmp_StrainF_vec[i * 6 + 5]  # equiv to e_21
            tmp_StrainF[i, :, :] = np.array([[e_11[i], e_12[i], e_13[i]], [e_12[i], e_22[i], e_23[i]], [e_13[i], e_23[i], e_33[i]]])
            tmp_StrainF_tran[i, :, :] = np.dot(np.dot(T.T, tmp_StrainF[i]), T)  # transform to sonata coordinate system
            # strains in "local" system
            e_11_M[i] = tmp_StrainF_M_vec[i * 6]
            e_22_M[i] = tmp_StrainF_M_vec[i * 6 + 1]
            e_33_M[i] = tmp_StrainF_M_vec[i * 6 + 2]
            e_23_M[i] = tmp_StrainF_M_vec[i * 6 + 3]  # equiv to e_23_M
            e_13_M[i] = tmp_StrainF_M_vec[i * 6 + 4]  # equiv to e_31_M
            e_12_M[i] = tmp_StrainF_M_vec[i * 6 + 5]  # equiv to e_21_M
            tmp_StrainF_M[i, :, :] = np.array([[e_11_M[i], e_12_M[i], e_13_M[i]], [e_12_M[i], e_22_M[i], e_23_M[i]], [e_13_M[i], e_23_M[i], e_33_M[i]]])
            tmp_StrainF_M_tran[i, :, :] = np.dot(np.dot(T.T, tmp_StrainF_M[i]), T)  # transform to sonata coordinate system

            # cell_id[i, :] = cd(i)
    elif voigt_convention == "paraview":  # different ordering compared to "anba"; "paraview" ordering: [s_xx, s_yy, s_zz, s_xy, s_yz, s_xz]
        print("ToDo - Process to paraview output")

    # Export to Paraview format (to be tested!)
    # file_res = do.XDMFFile('output_filename.xdmf')
    # file_res.parameters['functions_share_mesh'] = True
    # file_res.parameters['rewrite_function_mesh'] = False
    # file_res.parameters["flush_output"] = True
    # file_res.write(anba.STRESS, t=2)  # t=unique_number


    return tmp_StressF_tran, tmp_StressF_M_tran, tmp_StrainF_tran, tmp_StrainF_M_tran

