import math
import numpy as np

def cross2d(pos2d):
# m = [  0   -v(3)   v(2)
#       v(3)   0    -v(1)
#      -v(2)  v(1)    0  ];
    return np.array([[0, 0, pos2d[1]],[0., 0., -pos2d[0]],[-pos2d[1], pos2d[0], 0.]])

def rot2d(theta):
# m = [  0   -v(3)   v(2)
#       v(3)   0    -v(1)
#      -v(2)  v(1)    0  ];
    return np.array([[math.cos(theta), -math.sin(theta), 0.],[math.sin(theta), math.cos(theta), 0.],[0., 0., 1.]])

def TransformMass(dx, theta, M, S, J):
    R = rot2d(theta)
    dxcross = cross2d(dx)
    Sr = -dxcross@M + R@S@R.T
    Jr = R@J@R.T + R@S@R.T@dxcross + dxcross@R@S@R.T - dxcross@M@dxcross
    return (M, Sr, Jr)

def TransformStiffness(dx, theta, K):
    R = rot2d(theta)
    dxcross = cross2d(dx)
    H1 = np.block([[R, np.zeros((3,3))],[dxcross@R, R]])
    #H2 = np.block([[R, dxcross@R],[np.zeros((3,3)), R]])
    return H1@K@H1.T


