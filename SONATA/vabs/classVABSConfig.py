"""
Created on Mon Mar 20 14:15:44 2017

@author: TPflumm
"""
import numpy as np

class VABSConfig(object):
    
    """
    this class contains the Configuration for a VABS (Variational Asymptotic 
    Beam Sectional Analysis)
    """
    
    def __init__(self, **kw):
        self.format_flag = 0     #0:old format, 1:new format
        self.nlayer = 0          #nlayer is not used in the old format. nlayer shoud be always given a value greater than one if format_flag=1
        self.Timoshenko_flag = 1 #VABS will construct both the classical model and the generalized Timoshenko model. If it is 0, it will only construct the classical model. 
        self.recover_flag = 0    #0: VABS will carry out the consitutive modeling. 1: non-linear beam theory, 3D stress,displacement,strain recovery, 2:linear beam theory 3D recovery
        self.thermal_flag = 0    #Either 0 or 3: 0:pure mechanical analysis. 3:one-way coupled thermoelastic analyis.
        self.curve_flag = 0      #To model initially curved or twisted beams. if 1: three real numbers for twist(k1) and curvatures (k2 and k3) shoud be provided in the very next line    
        self.trapeze_flag = 0    #to obtain the trapeze effect trapeze_flag is 1
        self.Vlasov_flag = 0     #To obtain a generalized Vlasov model, can be 1 only if Timoshenko_flag is 1. For more information see VABS USER MANUAL
        self.oblique_flag = 0    #to model oblique cross sections. See VABS USER MANUAL
        
        if 'format_flag' in kw:     self.format_flag = kw['format_flag']  #Currently no need to use this!
        if 'nlayer' in kw:          self.nlayer = kw['nlayer']  
        if 'Timoshenko_flag' in kw: self.Timoshenko_flag = kw['Timoshenko_flag']
        if 'recover_flag' in kw:    self.recover_flag = kw['recover_flag']
        if 'thermal_flag' in kw:    self.thermal_flag = kw['thermal_flag']
        if 'trapeze_flag' in kw:    self.trapeze_flag = kw.get('trapeze_flag')
        if 'Vlasov_flag'  in kw:    self.Vlasov_flag = kw['Vlasov_flag']
        
        if 'curve_flag' in kw:      
            self.curve_flag = kw['curve_flag']
            self.k1 = kw['k1']
            self.k2 = kw['k2']
            self.k3 = kw['k3']
        if 'oblique_flags' in kw:
            self.oblique_flag = kw['oblique_flag']     #to model oblique cross sections. See VABS USER MANUAL
            self.oblique_cosine1 = kw['oblique_cosine1']             #Angle between beam axis x1 and and oblique axis y1. See VABS USER MANUAL Figure 6
            self.oblique_cosine2 = kw['oblique_cosine2']             #Angle between beam axis x1 and and oblique axis y2. See VABS USER MANUAL Figure 6
            self.Timoshenko_flag = 0                             #THis feature is only enabled in the classical beam model
            
        self.u = [0,0,0]                                    
        self.Cij = np.array([[0,0,0],[0,0,0],[0,0,0]])
        self.F = [0,0,0]                #Timoshenko_flag == 0 -> F = [0]    
        self.M = [0,0,0]
        self.f = [0,0,0]                #Timoshenko_flag == 1:
        self.df = [0,0,0]               #Timoshenko_flag == 1:
        self.ddf =[0,0,0]               #Timoshenko_flag == 1:
        self.dddf =[0,0,0]              #Timoshenko_flag == 1:
        self.m = [0,0,0]                #Timoshenko_flag == 1:
        self.dm = [0,0,0]               #Timoshenko_flag == 1:
        self.ddm =[0,0,0]               #Timoshenko_flag == 1:
        self.dddm =[0,0,0]              #Timoshenko_flag == 1:
        self.gamma11 = [0]              #Vlasov_flag == 1:
        self.kappa = [0,0,0]            #Vlasov_flag == 1:
        self.dkappa1 = 0                #Vlasov_flag == 1:
        self.ddkapp1 = 0                #Vlasov_flag == 1:
        

if __name__ == '__main__':
    test = VABSConfig(recover_flag=1)