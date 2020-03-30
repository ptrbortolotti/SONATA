#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 13:02:57 2018

@author: gu32kij
"""

# Core Library modules
# def read_VABS_Results(filename):
import os
import platform
import shutil
import subprocess
from datetime import datetime

# Third party modules
import numpy as np
from OCC.Core.gp import gp_Pnt, gp_Pnt2d, gp_Vec2d

# First party modules
from SONATA.cbm.display.display_mesh import plot_cells
from SONATA.cbm.fileIO.read_yaml_input import read_yaml_materialdb
from SONATA.cbm.mesh.cell import Cell
from SONATA.cbm.mesh.consolidate_mesh import consolidate_mesh_on_web
from SONATA.cbm.mesh.mesh_core import gen_core_cells
from SONATA.cbm.mesh.mesh_intersect import map_mesh_by_intersect_curve2d
from SONATA.cbm.mesh.mesh_utils import (find_node_by_ID,
                                        grab_nodes_of_cells_on_BSplineLst,
                                        merge_nodes_if_too_close,
                                        sort_and_reassignID,)
from SONATA.cbm.mesh.node import Node
from SONATA.vabs.strain import Strain
from SONATA.vabs.stress import Stress
from SONATA.vabs.VABS_interface import (VABS_config, XSectionalProperties,
                                        export_cells_for_VABS,)

if __name__ == "__main__":
    os.chdir("../../..")  # print(os.getcwd())





def read_nastran_bulkdata(filename):
    with open(filename) as f:
        nodes = []
        mesh = []
        for line in f:
            line = line.partition("$")[0]
            line = line.rstrip()

            if "GRID" in line:
                tmp = list(map("".join, zip(*[iter(line)] * 8)))[1:]
                tmp = [t.strip() for t in tmp if not t.isalpha()]
                tmp = list(filter(None, tmp))
                tmp = [float(t) for t in tmp]
                node = Node(gp_Pnt2d(tmp[1], tmp[2]))
                node.id = int(tmp[0])
                nodes.append(node)

            elif "CTRIA3" in line:
                tmp = list(map("".join, zip(*[iter(line)] * 8)))[1:]
                tmp = [t.strip() for t in tmp if not t.isalpha()]
                tmp = list(filter(None, tmp))

                nlst = [tmp[2], tmp[4], tmp[3]]  # counterclockwise
                c_nodeLst = [find_node_by_ID(nodes, int(n)) for n in nlst]
                cell = Cell(c_nodeLst)
                cell.id = int(tmp[0])
                cell.MatID = int(tmp[1])
                # Important See CTRIA3 Element Geometry and Coordinate System in NASTRAN Reference Guide
                nas_theta = float(tmp[5])
                v0 = gp_Vec2d(gp_Pnt2d(0, 0), gp_Pnt2d(1, 0))
                v1 = gp_Vec2d(cell.nodes[0].Pnt2d, cell.nodes[2].Pnt2d)
                theta_11 = (v0.Angle(v1)) * 180 / np.pi + nas_theta
                if theta_11 < 0:
                    theta_11 = 360 + theta_11
                cell.theta_1[0] = theta_11
                cell.theta_1[1] = 540
                cell.theta_3 = 0
                mesh.append(cell)

            elif "CQUAD4" in line:
                tmp = list(map("".join, zip(*[iter(line)] * 8)))[1:]
                tmp = [t.strip() for t in tmp if not t.isalpha()]
                tmp = list(filter(None, tmp))

                nlst = [tmp[2], tmp[5], tmp[4], tmp[3]]  # counterclockwise
                c_nodeLst = [find_node_by_ID(nodes, int(n)) for n in nlst]
                cell = Cell(c_nodeLst)
                cell.id = int(tmp[0])
                cell.MatID = int(tmp[1])

                # Important See CQUAD4 Element Geometry and Coordinate System in NASTRAN Reference Guide
                nas_theta = float(tmp[6])
                v0 = gp_Vec2d(gp_Pnt2d(0, 0), gp_Pnt2d(1, 0))
                v1 = gp_Vec2d(cell.nodes[0].Pnt2d, cell.nodes[3].Pnt2d)
                theta_11 = (v0.Angle(v1)) * 180 / np.pi + nas_theta
                if theta_11 < 0:
                    theta_11 = 360 + theta_11
                cell.theta_1[0] = theta_11
                cell.theta_1[1] = 540
                cell.theta_3 = 0
                mesh.append(cell)
    return nodes, mesh


if __name__ == "__main__":
    filename = "jobs/JBlaut/Querschnitt.nas"
    nodes, mesh = read_nastran_bulkdata(filename)
    plot_cells(mesh, nodes, "MatID")
    plot_cells(mesh, nodes, "theta_3", plotTheta11=True)
