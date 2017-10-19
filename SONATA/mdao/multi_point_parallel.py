# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 16:29:42 2017

@author: TPflumm
"""
if __name__ == "__main__":
    from openmdao.api import Problem
    from openmdao.core.mpi_wrap import MPI

    if MPI: # pragma: no cover
        # if you called this script with 'mpirun', then use the
        # petsc data passing
        from openmdao.core.petsc_impl import PetscImpl as impl
    else:
        # if you didn't use `mpirun`, then use the numpy data passing
        from openmdao.api import BasicImpl as impl

    def mpi_print(prob, *args):
        """ helper function to only print on rank 0"""
        if prob.root.comm.rank == 0:
            print(*args)

    prob = Problem(impl=impl) #set the implementation

    size = 10 #number of points

    adders = np.arange(size)/10.
    scalars = np.arange(size, 2*size)/10.

    prob.root = ParallelMultiPoint(adders, scalars)

    #turning off setup checking to avoid getting 10 sets of printouts to the screen
    prob.setup(check=False)

    st = time.time()

    prob['x'] = np.ones(size) * 0.7
    st = time.time()
    mpi_print(prob, "run started")
    prob.run()
    mpi_print(prob, "run finished", time.time() - st)

    mpi_print(prob, prob['aggregate.total'])