
import sys
sys.path.insert(0,'/Users/rfeil/work/6_SONATA/SONATA')  # import sys path to import 'SONATA' & 'job' modules

import math
import numpy as np
# import MassStiffTransformation as mst
import jobs.RFeil.anbax_studies.MassStiffTransformation as mst 

#import matplotlib.pyplot as plt
#from scipy.integrate import dblquad

import quadpy


def MassMatrix(a, b, density):
    Area = a * b
    Mass = Area * density
    M = np.diag([Mass, Mass, Mass])
    S = np.zeros((3,3))
    jxx = 1./12. * a * b**3 * density
    jyy = 1./12. * b * a**3 * density
    jzz = 1./12. * a * b * (a**2 + b**2) * density
    J = np.diag([jxx, jyy, jzz])
    return (M, S, J)

def k_shear_hutch(a, b, nu):
    # J. R. Hutchinson, Shear Coefficients for Timoshenko Beam Theory, J. Appl. Mechanics, 
    # Jan 2001, Vol 68, pp. 87--92, https://doi.org/10.1115/1.1349417
    a = a / 2.
    b = b / 2.
    C4 = 4./45. * a**3 * b *( -12 * a**2 - 15 * nu * a**2 + 5 * nu * b**2)
    for n in range(1,1001):
        C4 += 16. * nu**2 *b**5 * (n * math.pi * a - b * math.tanh( n * math.pi * a / b )) / (( n * math.pi )**5 * ( 1. + nu ))
    k = - 2. * ( 1 + nu ) / (9. / ( 4. * a**5 *b ) * C4 + nu * (1. - b**2 / a**2))
    return k

def w_shear_Tim1(xp, yp, a, b, nu):
    #expression for tau_xz and tau_yz taken from
    #Timoshenko and Goodier, pp. 323-324
    #print('chiamata ',x, y)
    
    J = 1. / 12. * b * a**3

    a = a / 2.
    b = b / 2.

    w = np.zeros(np.size(xp))
    n = np.array(range(1,101))
    for i in range(np.size(xp)):
        x = xp[i]
        y = yp[i]
        tauxz = np.zeros(100)#1 / 2 / J * (a**2 - x**2);
        tauyz = np.zeros(100)
        for m in range(0,101):
            #for n in range (1,101):
                # a = 0.12 b = 0.24
                #tauxz = tauxz -19119271143784725/2251799813685248*(-1)**(m+n-1)*np.cos((50/3*m+25/3)*math.pi*x)*np.cos(25/3*n*math.pi*y)*math.pi/(2*m+1)/((2*m+1)**2+n**2);
                #tauyz = tauyz -2294312537254167/2251799813685248*(-1)**(m+n-1)*np.sin((50/3*m+25/3)*math.pi*x)*(50/3*m+25/3)*math.pi*np.sin(25/3*n*math.pi*y)/(2*m+1)/n/((2*m+1)**2+n**2);
                # a = 0.24 b = 0.12
                tauxz += -(8*(-1)**(m + n - 1)*b**2*nu*np.cos((x*math.pi*(2*m + 1))/(2*a))*np.cos((math.pi*n*y)/b))/(J*math.pi**3*(2*m + 1)*(n**2 + (b**2*(2*m + 1)**2)/(4*a**2))*(nu + 1))
                tauyz += -(4*(-1)**(m + n - 1)*b**3*nu*np.sin((x*math.pi*(2*m + 1))/(2*a))*np.sin((math.pi*n*y)/b))/(J*a*n*math.pi**3*(n**2 + (b**2*(2*m + 1)**2)/(4*a**2))*(nu + 1))
        tauxz = np.sum(tauxz) + 1 / 2 / J * (a**2 - x**2)
        tauyz = np.sum(tauyz)
        w[i] = (tauxz**2 + tauyz**2) / (2.)
    return w

# def w_shear_Tim2(xp, yp, a, b, nu):
#     #expression for tau_xz and tau_yz taken from
#     #Timoshenko and Goodier, pp. 323-324
#     print(a,b)
#     J = 1 / 12 * b * a**3
# 
#     a = a / 2;
#     b = b / 2;
#     tauxz = np.zeros(100)#1 / 2 / J * (a**2 - x**2)
#     tauyz = np.zeros(100)#np.zeros(np.size(x))
# 
#     w = np.zeros(np.size(xp))
#     n = np.array(range(1,101))
#     count = -1
#     for i in range(np.size(xp)):
#         count = count + 1
#         x = xp[i]
#         y = yp[i]
#         tauxz = np.zeros(100)#1 / 2 / J * (a**2 - x**2);
#         tauyz = np.zeros(100)
#         for m in range(0,101):
#             #for n in range (1,101):
#                 #a = 0.12 b = 0.24
#                 tauxz += -(140737488355328*(-1)**(m + n - 1)*b**2*nu*math.pi*np.cos((x*math.pi*(2*m + 1))/(2*a))*np.cos((math.pi*n*y)/b))/(1713638851887625*J*(2*m + 1)*(n**2 + (b**2*(2*m + 1)**2)/(4*a**2))*(nu + 1))
#                 tauyz += -(70368744177664*(-1)**(m + n - 1)*b**3*nu*math.pi*np.sin((x*math.pi*(2*m + 1))/(2*a))*np.sin((math.pi*n*y)/b))/(1713638851887625*J*a*n*(n**2 + (b**2*(2*m + 1)**2)/(4*a**2))*(nu + 1))
#                 # a = 0.24 b = 0.12
#                 #tauxz = tauxz -19119271143784725/36028797018963968*(-1)**(m+n-1)*np.cos((25/3*m+25/6)*math.pi*x)*np.cos(50/3*n*math.pi*y)*math.pi/(2*m+1)/(1/16*(2*m+1)**2+n**2);
#                 #tauyz = tauyz -2294312537254167/72057594037927936*(-1)**(m+n-1)*np.sin((25/3*m+25/6)*math.pi*x)*(25/3*m+25/6)*math.pi*np.sin(50/3*n*math.pi*y)/(2*m+1)/n/(1/16*(2*m+1)**2+n**2);
#         tauxz = np.sum(tauxz) + 1 / 2 / J * (a**2 - x**2)
#         tauyz = np.sum(tauyz)
#         w[count] = (tauxz**2 + tauyz**2) / (2.);
#     return w


def Kshear_Tim(a, b, G, nu):
    scheme = quadpy.quadrilateral.sommariva_20()
    scheme1d = quadpy.line_segment.gauss_legendre(20)
    #scheme = quadpy.quadrilateral.product(scheme1d)
    #scheme.plot(quadpy.quadrilateral.rectangle_points([-a/2., a/2.], [-b/2., b/2.]), True)
    #print(quadpy.quadrilateral.rectangle_points([-a/2., a/2.], [-b/2., b/2.]))
    #plt.show()
#    w1 = dblquad(lambda x, y: w_shear_Tim1(x, y, a, b, nu), -b/2., b/2., lambda x: -a/.2, lambda x: a/2.) / G
#    w2 = dblquad(lambda x, y: w_shear_Tim2(x, y, a, b, nu),  -a/2., a/2., lambda x: -b/.2, lambda x: b/2.) / G
    w1 = scheme.integrate(lambda x: w_shear_Tim1(x[0], x[1], a, b, nu), quadpy.quadrilateral.rectangle_points([-a/2., a/2.], [-b/2., b/2.])) / G
    w2 = scheme.integrate(lambda x: w_shear_Tim1(x[0], x[1], b, a, nu), quadpy.quadrilateral.rectangle_points([-b/2., b/2.], [-a/2., a/2.])) / G
    return (1./2./w1, 1./2./w2)

def k_Mt(a, b):
    a = a / 2.
    b = b / 2.
    #Fraeijes pp 194
    #k1 = 8/3*b*a**3;
    #Timoshenko Eq. 154 pp. 300
    k2 = 16./3.*b*a**3
    for i in range(101):
        #k1 = k1 - 1024/math.pi**5/(2*i+1)**5*tanh(((2*i+1)*math.pi*b)/(2*a));
        k2 += 16./3.*b*a**3 * (- 192./math.pi**5*a/b/(2.*i+1.)**5*math.tanh(((2.*i+1.)*math.pi*b)/(2.*a)))
    return k2

def StiffnessMatrix(E, nu, a, b):
    G = E / (2 * ( 1 + nu ) )
    Area = a * b
    EA = E * Area
    #k1 = k_shear_hutch(a, b, nu)
    #k2 = k_shear_hutch(b, a, nu)
    #GA1 = G * Area * k1
    #GA2 = G * Area * k2
    (GA1, GA2) = Kshear_Tim(a, b, G, nu)
    Jxx = 1./12. * b**3 * a
    Jyy = 1./12. * a**3 * b
    EJxx = E * Jxx
    EJyy = E * Jyy
    GJ = G * k_Mt(a, b)
    Stiff = np.diag([GA1, GA2, EA, EJxx, EJyy, GJ])
    return Stiff


if __name__ == "__main__":
    E = 73800.E6
    nu = 0.33
    a = 0.12
    b = 0.24
    rho = 2700

    (M,S, J) = MassMatrix(a, b, rho)
    # (M, S, J) = mst.TransformMass((0,0), 0., M, J)
    (M, S, J) = mst.TransformMass(dx=(0.2,0.15), theta=0.17, M=M, S=S, J=J)
    K = StiffnessMatrix(E, nu, a, b)
    K = mst.TransformStiffness(dx=(0.2,0.15), theta=0.17, K=K)

    print('Mass (M) Matrix:')
    print(M)
    print('----------')
    print('S Matrix:')
    print(S)
    print('----------')
    print('Inertia (J) Matrix:')
    print(J)
    print('----------')
    print('Stiffness (K) Matrix:')
    print(K)
    pass
