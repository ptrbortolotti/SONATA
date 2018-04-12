# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 16:37:54 2018

@author: TPflumm
"""
from concurrent import futures
import matplotlib.mlab as mlab
from time import sleep, time

def func(v, tmp=0):
    X = v[0]
    Y = v[1]
    Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
    Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
    Z = (-10.0 * (Z2 - Z1))    # difference of Gaussians
    sleep(abs(v[0])*2)
    return (Z, 1, 1)

if __name__ == "__main__":
    __spec__ = None
    test = [[1,1],[0,0],[0.5,0.2],[2,1.3]]
    with futures.ProcessPoolExecutor(max_workers=4) as e:
#        fs = {e.submit(func, n): n for n in test}
#        print('Alle Aufgaben gestartet')
#        for f in futures.as_completed(fs):
#            print(f.result())
         #for i in zip(test, e.map(func,test)):
         res = list(e.map(func,test))
        