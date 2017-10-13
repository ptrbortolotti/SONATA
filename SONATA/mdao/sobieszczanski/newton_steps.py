# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 11:14:04 2017

@author: TPflumm


Sobieszczanski-Sobieski, J. - Multidiciplinary Design Optimization Supported by Knowledge Based Engineering 1st Edition (2015, Wiley)
4.4.1 Unconstrained Quadratic Search Method Employing Newton Steps


"""
import matplotlib.pyplot as plt
import numpy as np


#3D:
# Make data.
x1 = np.arange(-2, 4, 0.1)
x2 = np.arange(-3, 3, 0.1)
x1, x2 = np.meshgrid(x1, x2)
f = (x1-1)**2 + x2**2
plt.subplot(212)
CS = plt.contour(x1, x2, f, colors='k')
plt.clabel(CS, fontsize=9, inline=1)

def f2(x):
    return (x[0]-1)**2 + x[1]**2

def gradient_f2(x):
    return np.array([2*(x[0]-1),2*x[1]])
    
alpha = 0.1
epsilon = 1e-8

x = np.array([-1,2])
B = 1
k = 0 #iteration counter.
while B>=epsilon:
   k=k+1
   print k,x,B
   plt.plot(x[0],x[1],'b.')
   x_new = x-alpha*gradient_f2(x)
   B=abs(np.linalg.norm(gradient_f2(x)))
   x=x_new
   
plt.show()