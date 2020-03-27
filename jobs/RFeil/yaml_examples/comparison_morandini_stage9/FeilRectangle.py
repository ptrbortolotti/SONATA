#
# Copyright (C) 2018 Marco Morandini
#
#----------------------------------------------------------------------
#
#    This file is part of Anba.
#
#    Anba is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Hanba is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Anba.  If not, see <https://www.gnu.org/licenses/>.
#
#----------------------------------------------------------------------
#

from dolfin import *
# from dolfin import compile_extension_module
import time
import math
import numpy as np
from petsc4py import PETSc
import os
import matplotlib.pyplot as plt

# from anba4 import *
import mshr

# from voight_notation import stressVectorToStressTensor, stressTensorToStressVector, strainVectorToStrainTensor, strainTensorToStrainVector
# from material import material
# from anbax import anbax

import sys
sys.path.append('/Users/rfeil/work/7_anbax/anba_v4_revised')  # revised anbax v4 version from Marco Morandini (2019-11-27)
from anba4.anbax_MMtest import anbax

parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize_flags"] = "-O2"
parameters["form_compiler"]["quadrature_degree"] = 2

# Basic material parameters. 9 is needed for orthotropic materials.

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

thickness = 0.05
width = 0.48 
height = 0.24
Square1 = mshr.Rectangle(Point(-width/2., -height/2., 0.), Point(width/2., height/2., 0.))
Square2 = mshr.Rectangle(Point(-width/2.+thickness, -height/2.+thickness, 0.), Point(width/2.-thickness, height/2.-thickness, 0.))
C_shape = Square1 - Square2
mesh = mshr.generate_mesh(C_shape, 64)
mesh_points=mesh.coordinates()
#plot(mesh)
#plt.show()

# CompiledSubDomain
materials = MeshFunction("size_t", mesh, mesh.topology().dim())
fiber_orientations = MeshFunction("double", mesh, mesh.topology().dim())
plane_orientations = MeshFunction("double", mesh, mesh.topology().dim())
tol = 1e-14


# Rotate mesh.
theta = 23.#0.0
rotation_angle = 0.
materials.set_all(0)
fiber_orientations.set_all(theta)
plane_orientations.set_all(rotation_angle)


# rotate mesh.
rotate = Expression(("x[0] * (cos(rotation_angle)-1.0) - x[1] * sin(rotation_angle)",
    "x[0] * sin(rotation_angle) + x[1] * (cos(rotation_angle)-1.0)"), rotation_angle = rotation_angle * np.pi / 180.0,
    degree = 1)

ALE.move(mesh, rotate)

# Build material property library.
mat1 = material.OrthotropicMaterial(matMechanicProp)

matLibrary = []
matLibrary.append(mat1)


anba = anbax(mesh, 1, matLibrary, materials, plane_orientations, fiber_orientations, 1.E9)
stiff = anba.compute()
stiff.view()
