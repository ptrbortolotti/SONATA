# -*- coding: utf-8 -*-
from __future__ import print_function

from openmdao.api import IndepVarComp, Component, Problem, Group
from openmdao.api import ScipyOptimizer
from openmdao.api import ExecComp
from openmdao.api import view_tree, view_connections

class Paraboloid(Component):
    """ Evaluates the equation f(x,y) = (x-3)^2 + xy + (y+4)^2 - 3 """

    def __init__(self):
        super(Paraboloid, self).__init__()

        self.add_param('x', val=0.0)
        self.add_param('y', val=0.0)

        self.add_output('f_xy', shape=1)

    def solve_nonlinear(self, params, unknowns, resids):
        """f(x,y) = (x-3)^2 + xy + (y+4)^2 - 3
        """

        x = params['x'] #dictionary
        y = params['y']

        unknowns['f_xy'] = (x-3.0)**2 + x*y + (y+4.0)**2 - 3.0  #dictionary

    def linearize(self, params, unknowns, resids):
        """ Jacobian for our paraboloid."""

        x = params['x']
        y = params['y']
        J = {} #J, should be a dictionary whose keys are tuples of the form (‘unknown’, ‘param’) and whose values are n-d arrays or scalars. 

        J['f_xy', 'x'] = 2.0*x - 6.0 + y
        J['f_xy', 'y'] = 2.0*y + 8.0 + x
        return J
    
    
if __name__ == "__main__":

    top = Problem() #An instance of an OpenMDAO Problem is always the top object for running a model.

    root = top.root = Group() #Each Problem in OpenMDAO must contain a root Group. A Group is a System that contains other Components or Groups.

    root.add('p1', IndepVarComp('x', 3.0)) # IndepVarComp is a Component that provides the source for a variable which we can later give to a Driver as a design variable to control.
    root.add('p2', IndepVarComp('y', -4.0)) # The numbers 3.0 and -4.0 are values chosen for each as starting points for the optimizer.
    root.add('p', Paraboloid()) #Then we add the paraboloid using the same syntax as before, giving it the name ‘p’.)

    # Constraint Equation
    root.add('con', ExecComp('c = x-y')) #An ExecComp is a shortcut that lets us easily create a component that defines a simple expression for us.

    root.connect('p1.x', 'p.x') #Then we connect up the outputs of the IndepVarComps to the parameters of the Paraboloid. Notice the dotted naming convention used to refer to variables. 
    root.connect('p2.y', 'p.y')
    root.connect('p.x', 'con.x')
    root.connect('p.y', 'con.y')

    top.driver = ScipyOptimizer()
    top.driver.options['optimizer'] = 'SLSQP'

    top.driver.add_desvar('p1.x', lower=-50, upper=50)
    top.driver.add_desvar('p2.y', lower=-50, upper=50)
    top.driver.add_objective('p.f_xy')  #Finally, we add the objective. You can use any unknown in your model as the objective.
    top.driver.add_constraint('con.c', lower=15.0, upper=16.0)

    top.setup() #Before we can run our model we need to do some setup. This is done using the setup method on the Problem. This method performs all the setup of vector storage, data transfer, etc., necessary to perform calculations. Calling setup is required before running the model.
            
    top.run() #Now we can run the model using the run method of Problem.
    
    top.setup(check=False)
    view_tree(top, show_browser=True)
    view_connections(top.root, show_browser=True)
 
    print('Minimum of %f found at (%f, %f)' % (top['p.f_xy'], top['p.x'], top['p.y']))