#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 13:59:15 2019

@author: gu32kij
"""
import numpy as np
from io import StringIO 
import os

if __name__ == '__main__':
    os.chdir('../..')

from SONATA.cbm.fileIO.readinput import read_FLOATrowSTR
from SONATA.cbm.cbm_utl import trsf_sixbysix, trsf_coords
from SONATA.vabs.vabs_utl import transfer_matrix_unitconvertion, read_VABS_Results, \
                                grab_str_segment, read_PIA

class BeamSectionalProps(object):
    """
    this class stores the beam cross-sectional data and has methods to 
    read the result files from vabs
    
    The coordinate system follows the definition of VABS and DYMORE. 
    For more information see the the vabs and dymore user manual.
    the x/x1 axis is in the direction of the beam. the y/x2 axis points towards
    the leading edge, and z/x3 points accordingly upwards.
    
    The 6x6 sectional stiffness matrix, TS (Timoshenko Stiffness Matrix)
    (1-extension; 2,3-shear, 4-twist; 5,6-bending) relates the sectional axial 
    strain, epsilon1, transverse shearing strains, epsilon2 and epsilon3, 
    twisting curvatures, kappa1 and two bending curvatures, kappa2 and kappa3, 
    to the axial force, F1, transverse shear forces, F2 and F3, twisting 
    moment, M1, and two bending moments, M2 and M3. The relationship between 
    these sectional strains and sectional stress resultants takes the form of a 
    symmetric, 6x6 matrix
    | F1 |   | k11 k12 k13 k14 k15 k16 |   | epsilon1 | 
    | F2 |   | k12 k22 k23 k24 k25 k26 |   | epsilon2 |
    | F3 |   | k13 k23 k33 k34 k35 k36 |   | epsilon3 |
    | M1 | = | k14 k24 k34 k44 k45 k46 | * | kappa1   |
    | M2 |   | k15 k25 k35 k45 k55 k56 |   | kappa2   |
    | M3 |   | k16 k26 k36 k46 k56 k66 |   | kappa3   |

    The 6x6 sectional mass matrix, MM, relates the sectional linear velocities, 
    denoted v1, v2 and v3, and angular velocities, denoted omega1, omega2 and 
    omgea3, to the sectional linear momenta, denoted p1, p2 and p3, and angular
    momenta, denoted h1, h2 and h3. The relationship between these sectional 
    velocities and sectional momenta takes the form of a symmetric, 6x6 matrix:
    | p1 |   | M11 M12 M13 M14 M15 M16 |   | v1     | 
    | p2 |   | M12 M22 M23 M24 M25 M26 |   | v2     |
    | p3 |   | M13 M23 M33 M34 M35 M36 |   | v3     |
    | h1 | = | M14 M24 M34 M44 M45 M46 | * | omega1 |
    | h2 |   | M15 M25 M35 M45 M55 M56 |   | omega2 |
    | h3 |   | M16 M26 M36 M46 M56 M66 |   | omega3 |
    
    Due to the nature of the problem, many of these coefficients vanish, and 
    the remaining entries are written as:
    | p1 |   | m00        0        0         0  m00*Xm3  −m00*Xm2 |   | v1    | 
    | p2 |   | 0        m00        0  −m00*Xm3        0         0 |   | v2    |
    | p3 |   | 0          0      m00   m00*Xm2        0         0 |   | v3    |
    | h1 | = | 0   −m00*Xm3  m00*Xm2       m11        0         0 | * | omega1|
    | h2 |   | m00*Xm3    0        0        0       m22      −m23 |   | omega2|
    | h3 |   | −m00*Xm2   0        0        0      −m23       m33 |   | omega3|
    
    where m00 is the sectional mass per unit span and Xm2 and Xm3 the 
    coordinates of the sectional center of mass. m22 and m33 are the section 
    mass moments of inertia per unit span with respect to unit vectors i2 and 
    i3, respectively, and m23 the sectional cross-product of inertia per unit 
    span. The polar moment of inertia per unit span, denoted m11, is given by
    m11 = m22 + m33.


    Attributes
    ----------
    TS : ndarray
        Timoshenko Stiffness Matrix (1-extension; 2,3-shear, 4-twist; 5,6-bending)
    MM : ndarray
        6x6 Mass Matrix    
    Xg : ndarray
        Geometric Center of the Cross Section
    Xt : ndarray
        Neutral Axes (or Tension Center) of the Cross Section
    Xs : ndarray
        The Generalized Shear Center of the Cross Section
    VS: ndarray
        Vlasov Stiffness Matrix (1-extension; 2-twist; 3,4-bending; 5-twist rate)
    Ag: ndarray
        Trapeze Effects Related Matrices
    Bk: ndarray
        Trapeze Effects Related Matrices
    Ck: ndarray
        Trapeze Effects Related Matrices
    Dk: ndarray
        Trapeze Effects Related Matrices   


    Properties
    ----------
    CS : ndarray
        Classical Stiffness Matrix (1-extension; 2-twist; 3,4-bending)
    CF : ndarray
        Classical Flexibility Matrix (1-extension; 2-twist; 3,4-bending)
    TF : ndarray
        Timoshenko Flexibility Matrix (1-extension; 2,3-shear, 4-twist; 5,6-bending)
    MMatMC : ndarray
        6x6 Mass Matrix at Mass Center
    Xm : ndarray
        Center of Gravity (Mass Center)
    m00 : float
        mass per unit span
    m11 : float
        polar moment of inertia per unit span
    m22 : float
        mass moments of inertia per unit span with respect to unit vectors i2 
    m33 : float
        mass moments of inertia per unit span with respect to unit vectors i3 
    VF: ndarray
        Vlasov Flexibility Matrix (1-extension; 2-twist; 3,4-bending; 5-twist rate)

    ToDo
    ----------
    - Implement function to rotate the Vlasov Stiffness matrix and the Trapez 
        Effect related terms
    - Make one unit convertion method that converts the complete instance
    """    
    
    __slots__ = ('TS', 'MM', 'Xg', 'Xt', 'Xs', 'PIA', 'VS', 'Ag', \
                 'Bk', 'Ck', 'Dk', 'ELE', 'U')
    

    def __init__(self, fname=None):
        self.TS = None
        self.MM = None
        self.Xg = None
        self.Xt = None
        self.Xs = None
        self.PIA = None
        self.VS = None
        self.Ag = None
        self.Bk = None
        self.Ck = None
        self.Dk = None

        self.ELE = None
        self.U = None

        if fname:
            self.read_vabs_K(fname)

    @property
    def m00(self):
        return self.MM[0,0]
    
    @property
    def m11(self):
        return self.MM[3,3]
    
    @property
    def m22(self):
        return self.MM[4,4]
    
    @property
    def m33(self):
        return self.MM[5,5]
    
    @property
    def m23(self):
        return -self.MM[4,5]
    
    @property
    def Xm(self):
        return np.array([-self.MM[0,5]/self.m00, self.MM[0,4]/self.m00])
    
    @property
    def rg(self):
        return np.sqrt(self.m11/self.m00)
    
    @property
    def MMatMC(self):
        mask = np.zeros(self.MM.shape, dtype=bool)
        np.fill_diagonal(mask, True)
        mask[4,5] = True
        mask[5,4] = True
        return self.MM*mask
    
#    @property
#    def CS(self):
#        tmp = np.delete(self.TS,[1,2],axis=0)
#        return np.delete(tmp,[1,2],axis=1)

    @property
    def TF(self):
        return np.linalg.inv(self.TS)

#    @property
#    def CF(self):
#        return np.linalg.inv(self.CS)
    
    @property
    def VF(self):
        return np.linalg.inv(self.VS)
    

    
    def read_vabs_K(self, fname):    
        """
        reads the vabs .K output file
        
        Parameters
        ----------
        fname : str
            filename (incl. filepath) to the vabs .K output file
        """
    
        a = ''
        with open(fname) as f:
            for line in f:
                line = line.partition('#')[0]
                line = line.rstrip()
                a += line 
                a += '\n'
             
        STR = ''.join([s for s in a.strip().splitlines(True) if s.strip("\r\n").strip()])
    
        i_MM = STR.find('The 6X6 Mass Matrix')
        i_MC = STR.find('The Mass Center of the Cross Section')
        i_MMatMC = STR.find('The 6X6 Mass Matrix at the Mass Center')   
        i_MP = STR.find('The Mass Properties with respect to Principal Inertial Axes')
        i_GC = STR.find('The Geometric Center of the Cross Section')
        i_CS = STR.find('Classical Stiffness Matrix (1-extension; 2-twist; 3,4-bending)')
        i_CF = STR.find('Classical Flexibility Matrix (1-extension; 2-twist; 3,4-bending)')
        i_NA = STR.find('The Neutral Axes (or Tension Center) of the Cross Section')
        i_TS = STR.find('Timoshenko Stiffness Matrix (1-extension; 2,3-shear, 4-twist; 5,6-bending)')
        i_TF =STR.find('Timoshenko Flexibility Matrix (1-extension; 2,3-shear, 4-twist; 5,6-bending)')
        i_GSC = STR.find('The Generalized Shear Center of the Cross Section in the User Coordinate System')        
        i_VS = STR.find('Vlasov Stiffness Matrix (1-extension; 2-twist; 3,4-bending; 5-twist rate)')
        i_VF = STR.find('Vlasov Flexibility Matrix (1-extension; 2-twist; 3,4-bending; 5-twist rate)'
                        )
        # Trapeze Effects Related Matrices
        i_Ag = STR.find('Ag1--Ag1--Ag1--Ag1')
        i_Bk = STR.find('Bk1--Bk1--Bk1--Bk1')
        i_Ck = STR.find('Ck2--Ck2--Ck2--Ck2')
        i_Dk = STR.find('Dk3--Dk3--Dk3--Dk3')
        
        splitpoints = [i_MM,i_MC,i_MMatMC,i_MP,i_GC,i_CS,i_CF,i_NA,i_TS,i_TF,i_GSC,i_VS,i_VF,i_Ag,i_Bk,i_Ck,i_Dk]
        splitpoints.sort()
                
        MM = grab_str_segment(STR,i_MM,splitpoints)
        #MC = grab_str_segment(STR,i_MC,splitpoints)
        #MMatMC = grab_str_segment(STR,i_MMatMC,splitpoints)
        MP = grab_str_segment(STR,i_MP,splitpoints)
        GC = grab_str_segment(STR,i_GC,splitpoints)
        #CS = grab_str_segment(STR,i_CS,splitpoints)
        #CF = grab_str_segment(STR,i_CF,splitpoints)
        NA = grab_str_segment(STR,i_NA,splitpoints)
        TS = grab_str_segment(STR,i_TS,splitpoints)
        #TF = grab_str_segment(STR,i_TF,splitpoints)
        GSC = grab_str_segment(STR,i_GSC,splitpoints)
        VS = grab_str_segment(STR,i_VS,splitpoints)
        #VF = grab_str_segment(STR,i_VF,splitpoints)
        Ag = grab_str_segment(STR,i_Ag,splitpoints)
        Bk = grab_str_segment(STR,i_Bk,splitpoints)
        Ck = grab_str_segment(STR,i_Ck,splitpoints)
        Dk = grab_str_segment(STR,i_Dk,splitpoints)

        # MM_new = np.zeros((6,6))
        # #The 6X6 Mass Matrix:
        # MM_temp = np.loadtxt(StringIO(MM),skiprows=2)
        # for i in range(6):
        #     for j in range(6):
        #         MM_new[i,j] = float(MM_temp[i,j])
        # ToDO !!!! change the way of loading the textfile content

        self.MM = np.loadtxt(StringIO(MM),skiprows=2)
               
        #The 6X6 Mass Matrix at the Mass Center:
#        try:
#            self.MMatMC = np.loadtxt(StringIO(MMatMC),skiprows=2)
#        except:
#            self.MMatMC = None
            
        try:
            self.PIA = read_PIA(MP, 'The Principal Inertial Axes Rotated from User Coordinate System by')
        except:
            if  'The user coordinate axes are the principal inertial axes.' in MP:
                self.PIA = float(0)
        
        self.Xg = np.array([read_FLOATrowSTR(GC,'Xg2'),read_FLOATrowSTR(GC,'Xg3')])
        self.Xt = np.array([read_FLOATrowSTR(NA,'Xt2'),read_FLOATrowSTR(NA,'Xt3')])

        #Timoshenko Stiffness and Flexibility Matrix (1-extension; 2,3-shear, 4-twist; 5,6-bending)
        try:
            self.TS = np.loadtxt(StringIO(TS),skiprows=2)
        except:
            self.TS =None
                        
        #The Generalized Shear Center of the Cross Section in the User Coordinate System    
        try:
            self.Xs = np.array([read_FLOATrowSTR(GSC,'Xs2'),read_FLOATrowSTR(GSC,'Xs3')])
        except:
            self.Xs = None

        #Vlasov Stiffness and Flexibility Matrix (1-extension; 2-twist; 3,4-bending; 5-twist rate)
        try:
            self.VS = np.loadtxt(StringIO(VS),skiprows=2)
        except:
            self.VS = None
        
        #Trapeze Effects Related Matrices
        try:    
            self.Ag = np.loadtxt(StringIO(Ag),skiprows=2)
            self.Bk = np.loadtxt(StringIO(Bk),skiprows=2)
            self.Ck = np.loadtxt(StringIO(Ck),skiprows=2)
            self.Dk = np.loadtxt(StringIO(Dk),skiprows=2)
        except:
            self.Ag=None
            self.Bk=None
            self.Ck=None
            self.Dk=None
            

    def read_all_VABS_Results(self, filename='vabs_filename'):
        """
        reads the strains, stresses and displacements from VABS result files
        """
        # self.ELE = read_VABS_Results(self.filename.replace('.K','.ELE'))
        # self.U = read_VABS_Results(self.filename.replace('.K','.U'))
        self.ELE = read_VABS_Results(filename+'.ELE')
        self.U = read_VABS_Results(filename+'.U')

    
    def MM_convert_units(self,in_dct = {'mass': 'kg', 'length': 'm'},out_dct = {'mass':'kg', 'length':'m'}):
        """
        method to convert units of the 6x6 Mass Matrix
        
        Parameters
        ----------
        in_dct : dict, 
            input units, default = {'mass': 'g', 'length': 'mm'}
            
        out_dct : dict,
            output units,  default = {'mass':'kg', 'length':'m'}
            
        Returns
        ----------
        MM : np.array
            the converted 6x6 Mass Matrix
        """
        
        MMunits = np.array([['mass/length','','','','mass','mass'],
                           ['','mass/length','', 'mass','',''],
                           ['','','mass/length', 'mass','',''],
                           ['','mass','mass','mass*length','',''],
                           ['mass','','','','mass*length','mass*length'],
                           ['mass','','','','mass*length','mass*length']])
        MM_TM = transfer_matrix_unitconvertion(MMunits,in_dct,out_dct)
        return np.multiply(self.MM, MM_TM)
        
    
    def TS_convert_units(self,in_dct = {'force': 'N', 'length': 'm'}, out_dct = {'force': 'N', 'length': 'm'}):
        """
        method to convert units of the 6x6 Timoshenko Stiffness Matrix
        
        Parameters
        ----------
        in_dct : dict, 
            input units, default = {'force': 'N', 'length': 'mm'}
            
        out_dct : dict,
            output units,  default = {'force': 'N', 'length': 'm'}
            
        Returns
        ----------
        TS : np.array
            the converted 6x6 Timoshenko Stiffness Matrix
        
        """
        
        TSunits = np.array([['force','force','force','force*length','force*length','force*length'],
                           ['force','force','force','force*length','force*length','force*length'],
                           ['force','force','force','force*length','force*length','force*length'],
                           ['force*length','force*length','force*length','force*length*length','force*length*length','force*length*length'],
                           ['force*length','force*length','force*length','force*length*length','force*length*length','force*length*length'],
                           ['force*length','force*length','force*length','force*length*length','force*length*length','force*length*length']])
    
        TS_TM = transfer_matrix_unitconvertion(TSunits,in_dct,out_dct)
        return np.multiply(self.TS, TS_TM)


    def rotate(self, Theta=0, copy=True):
        """
        rotates the VabsSectionalProps by angle theta. 
        if copy
        
        Parameters
        ----------
        Theta : float
            angle of rotation about x1 in radians 
        copy : bool
            turn true if you want a new copied instance returned,
            default is False

        Returns
        ----------
        tmp : BeamSectionalProps
            copy and rotated VABSSectionalProps of self,
        """
        
        #defuine rotation matrix
        T = np.array([[1, 0 ,0],
                      [0, np.cos(Theta), -np.sin(Theta)],
                      [0, np.sin(Theta), np.cos(Theta)]])
        
        #transform Stiffness and Mass Matrix
        if copy:
            tmp = BeamSectionalProps()
        else:
            tmp = self 

        tmp.TS = trsf_sixbysix(self.TS,T)
        tmp.MM = trsf_sixbysix(self.MM,T)
        
        #transform coordinates
        try:
            tmp.Xg = trsf_coords(np.hstack((np.zeros((1,)),self.Xg)),T)
            tmp.Xt = trsf_coords(np.hstack((np.zeros((1,)),self.Xt)),T)
            tmp.Xs = trsf_coords(np.hstack((np.zeros((1,)),self.Xs)),T)
        except:
            pass
        
        #TODO: transform PIA
        tmp.PIA = None
        
        #TODO: transform Vlasov Stiffness Matrix and Trapez Effect related Matrices
        tmp.VS = None
        tmp.Ag = None
        tmp.Bk = None
        tmp.Ck = None
        tmp.Dk = None
        
        if copy:
            return tmp
        else:
            return None
    
#======================================================
#       MAIN
#======================================================       
if __name__ == '__main__':
    filename = 'cbm_noname_2019-04-02_145252576501.vab.K'
    bp1 = BeamSectionalProps(filename)
    bp2 = bp1.rotate(np.pi/8, copy=True)