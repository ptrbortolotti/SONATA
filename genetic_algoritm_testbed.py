# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 14:09:03 2018

@author: TPflumm
"""
import numpy as np
import copy
from concurrent import futures

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from time import sleep, time
from datetime import datetime

from six.moves import range

from pyDOE import lhs
from openmdao.utils.concurrent import concurrent_eval


class GeneticAlgorithm():
    """
    Simple Genetic Algorithm.

    This is the Simple Genetic Algorithm implementation based on 2009 AAE550: MDO Lecture notes of
    Prof. William A. Crossley. It can be used standalone or as part of the OpenMDAO Driver.


    Attributes
    ----------
    comm : MPI communicator or None
        The MPI communicator that will be used objective evaluation for each generation.
    elite : bool
        Elitism flag.
    lchrom : int
        Chromosome length.
    npop : int
        Population size.
    objfun : function
        Objective function callback.
    """

    def __init__(self, objfun, run_parallel=False, max_workers=8):
        """
        Initialize genetic algorithm object.

        Parameters
        ----------
        objfun : function
            Objective callback function.
        comm : MPI communicator or None
            The MPI communicator that will be used objective evaluation for each generation.
        max_workers: 8
        """
        self.objfun = objfun
        self.run_parallel = run_parallel
        self.max_workers = max_workers
        self.startTime = datetime.now()
        self.lchrom = 0
        self.npop = 0
        self.elite = True

    def execute_ga(self, vlb, vub, bits, pop_size, max_gen):
        """
        Perform the genetic algorithm.

        Parameters
        ----------
        vlb : ndarray
            Lower bounds array.
        vub : ndarray
            Upper bounds array.
        bits : ndarray
            Number of bits to encode the design space for each element of the design vector.
        pop_size : int
            Number of points in the population.
        max_gen : int
            Number of generations to run the GA.

        Returns
        -------
        ndarray
            Best design point
        float
            Objective value at best design point.
        int
            Number of successful function evaluations.
        """
        xopt = copy.deepcopy(vlb)
        fopt = np.inf
        self.lchrom = int(np.sum(bits))

        if np.mod(pop_size, 2) == 1:
            pop_size += 1
        self.npop = int(pop_size)
        fitness = np.zeros((self.npop, ))

        Pc = 0.5
        Pm = (self.lchrom + 1.0) / (2.0 * pop_size * np.sum(bits))
        elite = self.elite

        # TODO: from an user-supplied intial population
        # new_gen, lchrom = encode(x0, vlb, vub, bits)
        new_gen = np.round(lhs(self.lchrom, self.npop, criterion='center'))

        # Main Loop
        nfit = 0
        
        if self.run_parallel:
            with futures.ProcessPoolExecutor(max_workers=self.max_workers) as e:
            
            
                for generation in range(max_gen + 1):
                    print(generation, (datetime.now() - self.startTime))
                    old_gen = copy.deepcopy(new_gen)
                    x_pop = self.decode(old_gen, vlb, vub, bits)
        
        
                    # Evaluate points in this generation.
                    if self.run_parallel:
        
                        res = list(e.map(self.objfun,x_pop))
                        #print(results)
                        #
                        #cases = [((item, ii), None) for ii, item in enumerate(x_pop)]
                        #results = concurrent_eval(self.objfun, cases, self.comm, allgather=True)
        
                        fitness[:] = np.inf
                        for ii,r in enumerate(res):
                            if r[1]==True:
                                fitness[ii] = r[0] 
                                nfit += 1
                            else:
                                 print('A case failed:')
        
                    else:
                        # Serial
                        for ii in range(self.npop):
                            x = x_pop[ii]
                            
                            fitness[ii], success, _ = self.objfun(x, 0)
        
                            if success:
                                nfit += 1
                            else:
                                fitness[ii] = np.inf
        
                    # Elitism means replace worst performing point with best from previous generation.
                    if elite and generation > 0:
                        max_index = np.argmax(fitness)
                        old_gen[max_index] = min_gen
                        x_pop[max_index] =   min_x
                        fitness[max_index] = min_fit
        
                    # Find best performing point in this generation.
                    min_fit = np.min(fitness)
                    min_index = np.argmin(fitness)
                    min_gen = old_gen[min_index]
                    min_x = x_pop[min_index]
        
                    if min_fit < fopt:
                        fopt = min_fit
                        xopt = min_x
        
                    # Evolve new generation.
                    new_gen = self.tournament(old_gen, fitness)
                    new_gen = self.crossover(new_gen, Pc)
                    new_gen = self.mutate(new_gen, Pm)

        return xopt, fopt, nfit

    def tournament(self, old_gen, fitness):
        """
        Apply tournament selection and keep the best points.

        Parameters
        ----------
        old_gen : ndarray
            Points in current generation

        fitness : ndarray
            Objective value of each point.

        Returns
        -------
        ndarray
            New generation with best points.
        """
        new_gen = []
        idx = np.array(range(0, self.npop - 1, 2))
        for j in range(2):
            old_gen, i_shuffled = self.shuffle(old_gen)
            fitness = fitness[i_shuffled]

            # Each point competes with its neighbor; save the best.
            i_min = np.argmin(np.array([[fitness[idx]], [fitness[idx + 1]]]), axis=0)
            selected = i_min + idx
            new_gen.append(old_gen[selected])

        return np.concatenate(np.array(new_gen), axis=1).reshape(old_gen.shape)

    def crossover(self, old_gen, Pc):
        """
        Apply crossover to the current generation.

        Crossover flips two adjacent genes.

        Parameters
        ----------
        old_gen : ndarray
            Points in current generation

        Pc : float
            Probability of crossover.

        Returns
        -------
        ndarray
            Current generation with crossovers applied.
        """
        new_gen = copy.deepcopy(old_gen)
        num_sites = self.npop // 2
        sites = np.random.rand(num_sites, self.lchrom)
        idx, idy = np.where(sites < Pc)
        for ii, jj in zip(idx, idy):
            i = 2 * ii
            j = i + 1
            new_gen[i][jj] = old_gen[j][jj]
            new_gen[j][jj] = old_gen[i][jj]
        return new_gen

    def mutate(self, current_gen, Pm):
        """
        Apply mutations to the current generation.

        A mutation flips the state of the gene from 0 to 1 or 1 to 0.

        Parameters
        ----------
        current_gen : ndarray
            Points in current generation

        Pm : float
            Probability of mutation.

        Returns
        -------
        ndarray
            Current generation with mutations applied.
        """
        temp = np.random.rand(self.npop, self.lchrom)
        idx, idy = np.where(temp < Pm)
        current_gen[idx, idy] = 1 - current_gen[idx, idy]
        return current_gen

    def shuffle(self, old_gen):
        """
        Shuffle (reorder) the points in the population.

        Used in tournament selection.

        Parameters
        ----------
        old_gen : ndarray
            Old population.

        Returns
        -------
        ndarray
            New shuffled population.
        ndarray(dtype=np.int)
            Index array that maps the shuffle from old to new.
        """
        temp = np.random.rand(self.npop)
        index = np.argsort(temp)
        return old_gen[index], index

    def decode(self, gen, vlb, vub, bits):
        """
        Decode from binary array to real value array.

        Parameters
        ----------
        gen : ndarray
            Population of points, encoded.
        vlb : ndarray
            Lower bound array.
        vub : ndarray
            Upper bound array.
        bits : ndarray
            Number of bits for decoding.

        Returns
        -------
        ndarray
            Decoded design variable values.
        """
        num_desvar = len(bits)
        interval = (vub - vlb) / (2**bits - 1)
        x = np.empty((self.npop, num_desvar))
        sbit = 0
        ebit = 0
        for jj in range(num_desvar):
            exponents = 2**np.array(range(bits[jj] - 1, -1, -1))
            ebit += bits[jj]
            fact = exponents * (gen[:, sbit:ebit])
            x[:, jj] = np.einsum('ij->i', fact) * interval[jj] + vlb[jj]
            sbit = ebit
        return x

    def encode(self, x, vlb, vub, bits):
        """
        Encode array of real values to array of binary arrays.

        Parameters
        ----------
        x : ndarray
            Design variable values.
        vlb : ndarray
            Lower bound array.
        vub : ndarray
            Upper bound array.
        bits : int
            Number of bits for decoding.

        Returns
        -------
        ndarray
            Population of points, encoded.
        """
        # TODO : We need this method if we ever start with user defined initial sampling points.
        pass


def func(v, tmp=0):
    X = v[0]
    Y = v[1]
    Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
    Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
    Z = (10.0 * (Z2 - Z1))    # difference of Gaussians
    #sleep(1.5)
    return (Z, 1, 1)


if __name__ == "__main__":
    __spec__ = None
    ga = GeneticAlgorithm(func, run_parallel=True, max_workers=8)
    vlb = np.array([-3,-2])
    vub = np.array([3,2])
    bits = np.array([10,10])
    pop = 20
    max_gen = 15
    
    xopt, fopt, nfit = ga.execute_ga(vlb, vub, bits, pop, max_gen)
    
    # Plot the surface.
    X = np.linspace(-3,3,200)
    Y = np.linspace(-2,2,200)
    X, Y = np.meshgrid(X, Y)
    Z,sucess,_ = func([X,Y])
    plt.figure()
    CS = plt.contour(X, Y, Z, 60, cmap='plasma')
    plt.clabel(CS, inline=1, fontsize=9)
    plt.plot(xopt[0],xopt[1],'o')
    
    plt.show()