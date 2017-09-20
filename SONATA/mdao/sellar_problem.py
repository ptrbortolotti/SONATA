# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 15:12:45 2017

@author: TPflumm
"""

from __future__ import print_function
import numpy as np
from openmdao.api import Component
from openmdao.api import ExecComp, IndepVarComp, Group, NLGaussSeidel, \
                         ScipyGMRES
from openmdao.api import Problem, ScipyOptimizer

class SellarDis1(Component):
    """Component containing Discipline 1."""
    
    def __init__(self):
        super(SellarDis1, self).__init__()
        self.add_param('z', val=np.zeros(2))    # Global Design Variable
        self.add_param('x', val=0.)             # Local Design Variable
        self.add_param('y2', val=1.0)           # Coupling parameter
        self.add_output('y1', val=1.0)          # Coupling output

    def solve_nonlinear(self, params, unknowns, resids):
        """Evaluates the equation
        y1 = z1**2 + z2 + x1 - 0.2*y2"""
        
        z1 = params['z'][0]
        z2 = params['z'][1]
        x1 = params['x']
        y2 = params['y2']
        unknowns['y1'] = z1**2 + z2 + x1 - 0.2*y2

    def linearize(self, params, unknowns, resids):
        """ Jacobian for Sellar discipline 1."""
        J = {}
        J['y1','y2'] = -0.2
        J['y1','z'] = np.array([[2*params['z'][0], 1.0]])
        J['y1','x'] = 1.0
        return J


class SellarDis2(Component):
    """Component containing Discipline 2."""

    def __init__(self):
        super(SellarDis2, self).__init__()

        self.add_param('z', val=np.zeros(2))    # Global Design Variable
        self.add_param('y1', val=1.0)           # Coupling parameter
        self.add_output('y2', val=1.0)          # Coupling output
        
    def solve_nonlinear(self, params, unknowns, resids):
        """Evaluates the equation
        y2 = y1**(.5) + z1 + z2"""

        z1 = params['z'][0]
        z2 = params['z'][1]
        y1 = params['y1']

        # Note: this may cause some issues. However, y1 is constrained to be
        # above 3.16, so lets just let it converge, and the optimizer will
        # throw it out
        y1 = abs(y1)

        unknowns['y2'] = y1**.5 + z1 + z2

    def linearize(self, params, unknowns, resids):
        """ Jacobian for Sellar discipline 2."""
        J = {}
        J['y2', 'y1'] = .5*params['y1']**-.5
        J['y2', 'z'] = np.array([[1.0, 1.0]])       #Extra set of brackets below ensure we have a 2D array instead of a 1D array
        # for the Jacobian;  Note that Jacobian is 2D (num outputs x num inputs).    
        return J
    


class SellarDerivatives(Group):
    """ Group containing the Sellar MDA. This version uses the disciplines
    with derivatives."""

    def __init__(self):
        super(SellarDerivatives, self).__init__()

        self.add('px', IndepVarComp('x', 1.0), promotes=['x'])
        self.add('pz', IndepVarComp('z', np.array([5.0, 2.0])), promotes=['z'])

        self.add('d1', SellarDis1(), promotes=['z', 'x', 'y1', 'y2'])
        self.add('d2', SellarDis2(), promotes=['z', 'y1', 'y2'])

        self.add('obj_cmp', ExecComp('obj = x**2 + z[1] + y1 + exp(-y2)', z=np.array([0.0, 0.0]), x=0.0, y1=0.0, y2=0.0), promotes=['obj', 'z', 'x', 'y1', 'y2'])

        self.add('con_cmp1', ExecComp('con1 = 3.16 - y1'), promotes=['y1', 'con1'])
        self.add('con_cmp2', ExecComp('con2 = y2 - 24.0'), promotes=['con2', 'y2'])

        self.nl_solver = NLGaussSeidel()
        self.nl_solver.options['atol'] = 1.0e-12

        self.ln_solver = ScipyGMRES()
        
        
        
if __name__ == "__main__": 
    top = Problem()
    top.root = SellarDerivatives()
    
    top.driver = ScipyOptimizer()
    top.driver.options['optimizer'] = 'SLSQP'
    top.driver.options['tol'] = 1.0e-8
    
    top.driver.add_desvar('z', lower=np.array([-10.0, 0.0]),
                         upper=np.array([10.0, 10.0]))
    top.driver.add_desvar('x', lower=0.0, upper=10.0)
    
    top.driver.add_objective('obj')
    top.driver.add_constraint('con1', upper=0.0)
    top.driver.add_constraint('con2', upper=0.0)
    
    top.setup()
    
    # Setting initial values for design variables
    top['x'] = 1.0
    top['z'] = np.array([5.0, 2.0])
    
    top.run()
    
    print("\n")
    print( "Minimum found at (%f, %f, %f)" % (top['z'][0], \
                                             top['z'][1], \
                                             top['x']))
    print("Coupling vars: %f, %f" % (top['y1'], top['y2']))
    print("Minimum objective: ", top['obj'])