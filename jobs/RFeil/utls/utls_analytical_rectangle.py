import sys
sys.path.insert(0,'/Users/rfeil/work/6_SONATA/SONATA')  # import sys path to import 'SONATA' & 'job' modules

import math
import numpy as np
# import MassStiffTransformation as mst

import jobs.RFeil.utls.utls_analytical_IsotropicRectangle as ir 
import jobs.RFeil.utls.utls_analytical_MassStiffTransformation as mst 
import quadpy

from SONATA.cbm.cbm_utl import trsf_sixbysix

import matplotlib.pyplot as plt


from dolfin import *

# # orig anbax
# # sys.path.append('/Users/rfeil/work/7_anbax/anba_v4')  # orig anbax v4
# # from anba4.anbax import anbax, material

# # revised anbax from MM
# sys.path.append('/Users/rfeil/work/7_anbax/anba_v4_revised')  # revised anbax v4 version from Marco Morandini (2019-11-27)
# from anba4.anbax_MMtest import anbax, material




# Stand-alone script for analytically determining the structural properties (stiffness & mass matrices) of a rectangular beam
def utls_analytical_rectangle():    
    """
    Stand-alone script for analytically determining the structural properties (stiffness & mass matrices) of a rectangular beam

    inputs: -

    Outputs:
        data    -   analytically determined structural properties of a rectangular beam (transferred to SONATA/VABS coord sys)


    Notes:
    - If desired to determine a different structure, e.g. a circular beam, revise fcts in utls_analytical_IsotropicRectangle
    """

    # ---------------- #
    # User Input (hard coded)
    E = 73800.E6
    nu = 0.33
    a = 0.48
    b = 0.24
    rho = 2700
    dx = np.array([-0.14, -0.06])
    # dx = np.array([0., 0.])  # displacement of center axis in the two cross sectional directions
    theta_in = 0   # input - angle of rotation in deg

    theta = theta_in * math.pi / 180.  # conversion of angle of rotation from deg to rad

    # ---------------- #
    # Analytic Approach
    # ---------------- #
    (M, S, J) = ir.MassMatrix(a, b, rho)
    #print(M)
    #print('----------')
    #print(J)
    #print('----------')
    (Mt, St, Jt) = mst.TransformMass(dx, theta, M, S, J)
    K = ir.StiffnessMatrix(E, nu, a, b)
    Kt = mst.TransformStiffness(dx, theta, K)
    #print(M)
    #print('----------')
    #print(S)
    #print('----------')
    #print(J)
    #print('----------')
    #print(K)
    MM = np.block([[M, S],[S.T, J]])
    MMt = np.block([[Mt, St],[St.T, Jt]])

    # ---------------- #
    # (Optionally) Plot mesh of shape for verification
    # 
    # ---------------- #
    mesh = RectangleMesh(Point(-a/2., -b/2.), Point(a/2., b/2.), 20, 40)
    mesh.rotate(theta * 180. / math.pi)
    mesh.translate(Point(dx))
    plot(mesh)
    plt.show()


    # The following is only for calling anbax
    # matMechanicProp = [E, nu]
    # materials = MeshFunction("size_t", mesh, mesh.topology().dim())
    # fiber_orientations = MeshFunction("double", mesh, mesh.topology().dim())
    # plane_orientations = MeshFunction("double", mesh, mesh.topology().dim())
    # materials.set_all(0)
    # fiber_orientations.set_all(0.0)
    # plane_orientations.set_all(90.0)

    # # Build material property library.
    # mat1 = material.IsotropicMaterial(matMechanicProp, rho)
    # matLibrary = []
    # matLibrary.append(mat1)

    # anba = anbax(mesh, 2, matLibrary, materials, plane_orientations, fiber_orientations)
    # AnbaK = np.mat(anba.compute().getValues(range(6), range(6)))
    # AnbaM = np.mat(anba.inertia().getValues(range(6), range(6)))




    #Define transformation T (from ANBA to SONATA/VABS coordinates)
    # B = np.array([[0,0,1],[1,0,0],[0,1,0]])
    B = np.array([[0,0,-1],[-1,0,0],[0,1,0]])  # new
    T = np.dot(np.identity(3),np.linalg.inv(B))

    K_VABS_coords = trsf_sixbysix(K,T)
    Kt_VABS_coords = trsf_sixbysix(Kt,T)
    MM_VABS_coords = trsf_sixbysix(MM,T)
    MMt_VABS_coords = trsf_sixbysix(MMt,T)

    print('Analytical Rectangular Properties')
    print('Width: ' + str(a))
    print('Height: ' + str(b))
    print('Trans in width dir: ' + str(dx[0]) + ' and height dir: ' + str(dx[1]))
    print('Width: ' + str(theta))
    np.set_printoptions(precision=4)
    print('***************************')
    print('********* RESULTS *********')
    print('***************************')
    print('Stiffness (K) Matrix - Analytical')
    print(Kt_VABS_coords)
    print('***************************')
    print('Mass (M) Matrix - Analytical')
    print(MMt_VABS_coords)
    

    

    data = {}
    data['K'] = K_VABS_coords       # original stiffness matrix
    data['Kt'] = Kt_VABS_coords     # transformed (translational & rotational) stiffness matrix
    data['K'] = MM_VABS_coords      # original mass matrix
    data['K'] = MMt_VABS_coords     # transformed (translational & rotational) mass matrix


    return data


if __name__ == "__main__":

    aRect = utls_analytical_rectangle()


# EOF

