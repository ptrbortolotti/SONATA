
import sys
sys.path.insert(0,'/Users/rfeil/work/6_SONATA/SONATA')  # import sys path to import 'SONATA' & 'job' modules


# import IsotropicRectangle as ir
# import MassStiffTransformation as mst
import jobs.RFeil.anbax_studies.IsotropicRectangle as ir 
import jobs.RFeil.anbax_studies.MassStiffTransformation as mst 
import numpy as np
import math
import matplotlib.pyplot as plt

from dolfin import *
# from anba4 import *

# orig anbax
# sys.path.append('/Users/rfeil/work/7_anbax/anba_v4')  # orig anbax v4
# from anba4.anbax import anbax, material

# revised anbax from MM
sys.path.append('/Users/rfeil/work/7_anbax/anba_v4_revised')  # revised anbax v4 version from Marco Morandini (2019-11-27)
from anba4.anbax_MMtest import anbax, material

if __name__ == "__main__":
    E = 73800.E6
    nu = 0.33
    a = 0.12
    b = 0.24
    rho = 2700
    dx = np.array([0.3, 0.1])
    # dx = np.array([0.15, 0.3])
    # dx = np.array([0., 0.])  # displacement of center axis in the two cross sectional directions
    # theta = 10 * math.pi / 180.  # angle of rotation in rad
    # theta= 0

    theta = math.pi * 0.27  

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
    # Numeric (anbax) Approach
    # ---------------- #
    mesh = RectangleMesh(Point(-a/2., -b/2.), Point(a/2., b/2.), 20, 40)
    mesh.rotate(theta * 180. / math.pi)
    mesh.translate(Point(dx))
    # plot(mesh)
    # plt.show()
    matMechanicProp = [E, nu]
    materials = MeshFunction("size_t", mesh, mesh.topology().dim())
    fiber_orientations = MeshFunction("double", mesh, mesh.topology().dim())
    plane_orientations = MeshFunction("double", mesh, mesh.topology().dim())
    materials.set_all(0)
    fiber_orientations.set_all(0.0)
    plane_orientations.set_all(90.0)

    # Build material property library.
    mat1 = material.IsotropicMaterial(matMechanicProp, rho)
    matLibrary = []
    matLibrary.append(mat1)

    anba = anbax(mesh, 2, matLibrary, materials, plane_orientations, fiber_orientations)
    AnbaK = np.mat(anba.compute().getValues(range(6), range(6)))
    AnbaM = np.mat(anba.inertia().getValues(range(6), range(6)))


    np.set_printoptions(precision=4)
    print('***************************')
    print('********* RESULTS *********')
    print('***************************')
    print('Stiffness (K) Matrix - AnbaX')
    print(AnbaK)
    print('---------------------------')
    print('Stiffness (K) Matrix - Analytical')
    print(Kt)
    print('---------------------------')
    print('Difference between K Matrices')
    # print((Kt-AnbaK)@np.linalg.inv(K))
    print((Kt-AnbaK))
    # print('Matrix multiply')
    # print(np.linalg.inv(K))

    print('***************************')
    print('Mass (M) Matrix - AnbaX')
    print(AnbaM)
    print('---------------------------')
    print('Mass (M) Matrix - Analytical')
    print(MMt)
    print('---------------------------')
    print('Difference between M Matrices')
    # print((MMt-AnbaM)@np.linalg.inv(MM))
    print(MMt-AnbaM)
    
    pass
