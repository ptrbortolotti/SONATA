import dolfin

def build_dolfin_mesh(cbm_mesh, cbm_nodes):
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
        
    
    Notes
    ----------
    the cells of cbm_mesh already contain the nodes. So the information is 
    currently passed twice. But consistent with the export_cells_for_vabs.
    
    """
    
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
    
    return mesh

