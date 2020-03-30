
# Third party modules
import numpy as np

try
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
    fiber_orientations.set_all(0.0)
    plane_orientations.set_all(0.0)
    
    for c in cbm_mesh:
        materials[c.id-1] = matdict[c.MatID]
        plane_orientations[c.id-1] = c.theta_1[0]  # rotation around x1-axis (equiv. to beam axis) in SONATA/VABS coordinates
        fiber_orientations[c.id-1] = c.theta_3     # rotation around x3-axis in SONATA/VABS coordinates


    return mesh, matLibrary, materials, plane_orientations, fiber_orientations, maxE
