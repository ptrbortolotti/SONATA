# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 09:31:30 2018

@author: TPflumm
"""

#hello.py
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
print "hello world from process ", rank