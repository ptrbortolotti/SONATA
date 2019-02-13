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
from SONATA.vabs.vabs_utl import transfer_matrix_unitconvertion, read_VABS_Results, \
                                grab_str_segment, read_PIA

class VABSSectionalProps(object):
    """
    this class stores the vabs cross-sectional result data and has methods to 
    read the result files
    """    
    def __init__(self, filename=None):
        if filename is None:
            pass
        elif filename is not None:
            self.read_vabs_K(filename)
            self.filename = filename
 
    def read_vabs_K(self, filename):    

        a = ''
        with open(filename) as f:
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
        MC = grab_str_segment(STR,i_MC,splitpoints)
        MMatMC = grab_str_segment(STR,i_MMatMC,splitpoints)
        MP = grab_str_segment(STR,i_MP,splitpoints)
        GC = grab_str_segment(STR,i_GC,splitpoints)
        CS = grab_str_segment(STR,i_CS,splitpoints)
        CF = grab_str_segment(STR,i_CF,splitpoints)
        NA = grab_str_segment(STR,i_NA,splitpoints)
        TS = grab_str_segment(STR,i_TS,splitpoints)
        TF = grab_str_segment(STR,i_TF,splitpoints)
        GSC = grab_str_segment(STR,i_GSC,splitpoints)
        VS = grab_str_segment(STR,i_VS,splitpoints)
        VF = grab_str_segment(STR,i_VF,splitpoints)
        Ag = grab_str_segment(STR,i_Ag,splitpoints)
        Bk = grab_str_segment(STR,i_Bk,splitpoints)
        Ck = grab_str_segment(STR,i_Ck,splitpoints)
        Dk = grab_str_segment(STR,i_Dk,splitpoints)
        
        #The 6X6 Mass Matrix:
        self.MM = np.loadtxt(StringIO(MM),skiprows=2)
        #The Mass Center of the Cross Section:
        self.Xm2 = read_FLOATrowSTR(MC,'Xm2')
        self.Xm3 = read_FLOATrowSTR(MC,'Xm3')
        #The 6X6 Mass Matrix at the Mass Center:
        self.MMatMC = np.loadtxt(StringIO(MMatMC),skiprows=2)
        #The Mass Properties with respect to Principal Inertial Axes:
        self.MpUS = read_FLOATrowSTR(MP,'Mass Per Unit Span')
        self.I1 = read_FLOATrowSTR(MP,'Mass Moments of Intertia about x1 axis')
        self.I2 = read_FLOATrowSTR(MP,'Mass Moments of Intertia about x2 axis')
        self.I3 = read_FLOATrowSTR(MP,'Mass Moments of Intertia about x3 axis')
        
        try:
            self.PIA = read_PIA(MP, 'The Principal Inertial Axes Rotated from User Coordinate System by')
        except:
            if  'The user coordinate axes are the principal inertial axes.' in MP:
                self.PIA = float(0)
                
        self.Rg = read_FLOATrowSTR(MP,'The mass-weighted radius of gyration')
        #The Geometric Center of the Cross Section:
        self.Xg2 = read_FLOATrowSTR(GC,'Xg2')
        self.Xg3 = read_FLOATrowSTR(GC,'Xg3')
        # The Neutral Axes (or Tension Center) of the Cross Section:
        self.Xt2 = read_FLOATrowSTR(NA,'Xt2')
        self.Xt3 = read_FLOATrowSTR(NA,'Xt3')
        
        
        # Classical Stiffness and Flexibility Matrix (1-extension; 2-twist; 3,4-bending)
        try:
            self.CS = np.loadtxt(StringIO(CS),skiprows=2)
            self.CF = np.loadtxt(StringIO(CF),skiprows=2)
        except:
            self.CS =None
            self.CF =None
    
        # Timoshenko Stiffness and Flexibility Matrix (1-extension; 2,3-shear, 4-twist; 5,6-bending)
        try:
            self.TS = np.loadtxt(StringIO(TS),skiprows=2)
            self.TF = np.loadtxt(StringIO(TF),skiprows=2)
        except:
            self.TS =None
            self.TF =None
                        
        #The Generalized Shear Center of the Cross Section in the User Coordinate System    
        try:
            self.Xs2 = read_FLOATrowSTR(GSC,'Xs2')
            self.Xs3 = read_FLOATrowSTR(GSC,'Xs3')
        except:
            self.Xs2 =None
            self.Xs3 =None
        
        #Vlasov Stiffness and Flexibility Matrix (1-extension; 2-twist; 3,4-bending; 5-twist rate)
        try:
            self.VS = np.loadtxt(StringIO(VS),skiprows=2)
            self.VF = np.loadtxt(StringIO(VF),skiprows=2)
        except:
            self.VS =None
            self.VF =None
        
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
            

    def read_all_VABS_Results(self):
        self.ELE = read_VABS_Results(self.filename.replace('.K','.ELE'))
        self.U = read_VABS_Results(self.filename.replace('.K','.U'))

    
    def MM_convert_units(self,in_dct = {'mass': 'g', 'length': 'mm'},out_dct = {'mass':'kg', 'length':'m'}):
        MMunits = np.array([['mass/length','','','','mass','mass'],
                           ['','mass/length','', 'mass','',''],
                           ['','','mass/length', 'mass','',''],
                           ['','mass','mass','mass*length','',''],
                           ['mass','','','','mass*length','mass*length'],
                           ['mass','','','','mass*length','mass*length']])
        MM_TM = transfer_matrix_unitconvertion(MMunits,in_dct,out_dct)
        return np.multiply(self.MM, MM_TM)
        
    
    def TS_convert_units(self,in_dct = {'force': 'N', 'length': 'mm'}, out_dct = {'force': 'N', 'length': 'm'}):
        TSunits = np.array([['force','force','force','force*length','force*length','force*length'],
                           ['force','force','force','force*length','force*length','force*length'],
                           ['force','force','force','force*length','force*length','force*length'],
                           ['force*length','force*length','force*length','force*length*length','force*length*length','force*length*length'],
                           ['force*length','force*length','force*length','force*length*length','force*length*length','force*length*length'],
                           ['force*length','force*length','force*length','force*length*length','force*length*length','force*length*length']])
    
        TS_TM = transfer_matrix_unitconvertion(TSunits,in_dct,out_dct)
        return np.multiply(self.TS, TS_TM)


#======================================================
#       MAIN
#======================================================       
if __name__ == '__main__':
    filename = 'naca0012_cspar_mesh.vab.K'
    BeamProperties = VABSSectionalProps(filename)
    filename.replace('.K','.ELE')
    BeamProperties.read_VABS_ELE(filename)