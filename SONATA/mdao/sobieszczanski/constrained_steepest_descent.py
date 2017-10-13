# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 11:14:04 2017

@author: TPflumm


Sobieszczanski-Sobieski, J. - Multidiciplinary Design Optimization Supported by Knowledge Based Engineering 1st Edition (2015, Wiley)
4.5.3 Interior Panelty Method


"""
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import inv

plt.close('all')
def f(x):
    return x[0]**2+2*x[1]**2

def g(x):
    return 10-3*x[0]-4*x[1]

def B(x):
    return -np.log(-g(x))

def phi(x,r):
    #return f(x)+r*B(x)
    return x[0]**2+2*x[1]**2-r*np.log(3*x[0]+4*x[1]-10)

def nabla_phi(x,r):
    return np.array([2*x[0]-3*r/(3*x[0]+3*x[1]-10), 4*x[1]-4*r/(3*x[0]+3*x[1]-10)]) 

def Hessian_phi(x,r):
    return np.array([[2+9*r/(3*x[0]+3*x[1]-10)**2, 12*r/(3*x[0]+3*x[1]-10)**2],
             [12*r/(3*x[0]+3*x[1]-10)**2,  4+16*r/(3*x[0]+3*x[1]-10)**2]])
    
def nabla_f(x):
    return np.array([2*x[0],4*x[1]])

def N(x):
    return np.array([-3,-4])
    
#Plot f(x) and g(x):
x1 = np.linspace(-4, 4, 100)
x2 = np.linspace(-4, 4, 100)
x1,x2 = np.meshgrid(x1, x2)
CS = plt.contour(x1, x2, f([x1,x2]), colors='k')
plt.clabel(CS, fontsize=9, inline=1)
CS = plt.contour(x1, x2, g([x1,x2]),20, colors='r')
plt.clabel(CS, fontsize=9, inline=1)


x0 = np.array([5,0])
r = 0.1
plt.plot(x0[0],x0[1],'ko')


x = x0
alpha = 0.2
k = 0 #iteration counter.
crit = abs(np.linalg.norm(nabla_phi(x,r)))
#while crit>=1e-6:
#   k=k+1
#   print k,x,crit
#   plt.plot(x[0],x[1],'b.')
#   dx = -(I-N(x)*)
#   x_new = x+alpha*dx
#   #x_new = x-alpha*nabla_phi(x,r)
#   x=x_new
#   crit = abs(np.linalg.norm(nabla_phi(x,r)))
#plt.show()