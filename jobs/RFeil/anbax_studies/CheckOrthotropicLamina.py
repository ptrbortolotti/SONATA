from dolfin import *
import numpy as np
import math
# import MassStiffTransformation as mst
import jobs.RFeil.utls.utls_analytical_MassStiffTransformation as mst
# from anba4 import *

import sys
sys.path.append('/Users/rfeil/work/7_anbax/anba_v4_revised')  # revised anbax v4 version from Marco Morandini (2019-11-27)
from anba4.anbax_MMtest import anbax, material
#
# import matplotlib.pyplot as plt


e_xx = 9.8e9
e_yy = 9.8e9
e_zz = 1.42e11
g_xy = 4.8e9
g_xz = 6.0e9
g_yz = 6.0e9
nu_xy = 0.34
nu_zx = 0.3
nu_zy = 0.3
#Assmble into material mechanical property Matrix.
matMechanicProp = np.zeros((3,3))
matMechanicProp[0,0] = e_xx
matMechanicProp[0,1] = e_yy
matMechanicProp[0,2] = e_zz
matMechanicProp[1,0] = g_yz
matMechanicProp[1,1] = g_xz
matMechanicProp[1,2] = g_xy
matMechanicProp[2,0] = nu_zy
matMechanicProp[2,1] = nu_zx
matMechanicProp[2,2] = nu_xy

# Meshing domain.
sectionWidth = 3.0023e-2
sectionHeight = 1.9215e-3
nPly = 1 # t = 0.2452mm per ply.
#mesh = RectangleMesh.create([Point(0., 0.), Point(sectionWidth, sectionHeight)], [30, 32], CellType.Type.quadrilateral)
mesh = RectangleMesh(Point(0., 0.), Point(sectionWidth, sectionHeight), 30, 32, 'crossed')
# plot(mesh)
# plt.show()
ALE.move(mesh, Constant([-sectionWidth/2.0, -sectionHeight/2.0]))
mesh_points=mesh.coordinates()
#print(mesh_points)

# CompiledSubDomain
materials = MeshFunction("size_t", mesh, mesh.topology().dim())
fiber_orientations = MeshFunction("double", mesh, mesh.topology().dim())
plane_orientations = MeshFunction("double", mesh, mesh.topology().dim())
#isActive = MeshFunction("bool", mesh, mesh.topology().dim())
tol = 1e-14

# Rotate mesh.
rotation_angle = 35.#23.0       # equiv. theta_1
theta = 23.                     # equiv. theta_3
dx = np.array([0., 0.])
materials.set_all(0)
mesh.rotate(rotation_angle)
fiber_orientations.set_all(theta)
plane_orientations.set_all(rotation_angle)

#subdomain_0_p20.mark(materials, 1)

# Build material property library.
mat1 = material.OrthotropicMaterial(matMechanicProp)

matLibrary = []
matLibrary.append(mat1)


anba = anbax(mesh, 1, matLibrary, materials, plane_orientations, fiber_orientations)
Kbase = np.mat(anba.compute().getValues(range(6), range(6)))
np.set_printoptions(precision=4)
print(Kbase)
print()


print('CHECK THE FOLLOWING:')

# rotate mesh.
rotation_angle = 35.0
theta = 23.

mesh.rotate(rotation_angle)

fiber_orientations.set_all(theta)
plane_orientations.set_all(rotation_angle)

anba = anbax(mesh, 1, matLibrary, materials, plane_orientations, fiber_orientations)
K = np.mat(anba.compute().getValues(range(6), range(6)))
Kr = mst.TransformStiffness(dx, rotation_angle / 180. * math.pi, Kbase)
print(K)
print()

print((Kr-K)@np.linalg.inv(Kbase))

