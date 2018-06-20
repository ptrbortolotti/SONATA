from __future__ import print_function
from six.moves import range
import time
import numpy as np

from openmdao.api import Component, Group, ParallelGroup, IndepVarComp, ExecComp


class Plus(Component):
    """
    adder: float
        value that is added to every element of the x array parameter
    """

    def __init__(self, adder):
        super(Plus, self).__init__()
        self.add_param('x', 0.0)
        self.add_output('f1', shape=1)
        self.adder = float(adder)

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['f1'] = params['x'] + self.adder

        #sleep to slow things down a bit, to show the parallelism better
        time.sleep(.1)

class Times(Component):
    """
    scalar: float
        every element of the x array parameter is multiplied by this value
    """

    def __init__(self, scalar):
        super(Times, self).__init__()
        self.add_param('f1', 0.0)
        self.add_output('f2', shape=1)
        self.scalar = float(scalar)

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['f2'] = params['f1'] * self.scalar

class Point(Group):
    """
    Single point combining Plus and Times. Multiple copies will be made for
    a multi-point problem
    """

    def __init__(self, adder, scalar):
        super(Point, self).__init__()

        self.add('plus', Plus(adder), promotes=['x', 'f1'])
        self.add('times', Times(scalar), promotes=['f1', 'f2'])


class Summer(Component):
    """
    Aggregating component that takes all the multi-point values
    and adds them together
    """

    def __init__(self, size):
        super(Summer, self).__init__()
        self.size = size
        self.vars = []
        for i in range(size):
            v_name = 'f2_%d'%i
            self.add_param(v_name, 0.)
            self.vars.append(v_name)

        self.add_output('total', shape=1)

    def solve_nonlinear(self, params, unknowns, resids):
        tot = 0
        for v_name in self.vars:
            tot += params[v_name]
        unknowns['total'] = tot

class ParallelMultiPoint(Group):

    def __init__(self, adders, scalars):
        super(ParallelMultiPoint, self).__init__()

        size = len(adders)
        self.add('desvar', IndepVarComp('x', val=np.zeros(size)),
                                        promotes=['x'])

        self.add('aggregate', Summer(size))

        # a ParallelGroup works just like a Group if it's run in serial,
        # so using a ParallelGroup here will make our ParallelMultipoint
        # class work well in serial and under MPI.
        pg = self.add('multi_point', ParallelGroup())

        #This is where you stamp out all the points you need
        for i,(a,s) in enumerate(zip(adders, scalars)):
            c_name = 'p%d'%i
            pg.add(c_name, Point(a,s))
            self.connect('x', 'multi_point.%s.x'%c_name, src_indices=[i])
            self.connect('multi_point.%s.f2'%c_name,'aggregate.f2_%d'%i)


from openmdao.api import Problem


prob = Problem()

size = 10 #number of points

adders = np.arange(size)/10.
scalars = np.arange(size, 2*size)/10.

prob.root = ParallelMultiPoint(adders, scalars)

prob.setup()

st = time.time()

prob['x'] = np.ones(size) * 0.7
st = time.time()
print("run started")
prob.run()
print("run finished", time.time() - st)

print(prob['aggregate.total'])