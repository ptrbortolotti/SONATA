"""
Created on Mon Mar 20 14:15:44 2017

@author: TPflumm
"""
import numpy as np

class VABSConfig(object):
    
    """
    this class contains the Configuration for a VABS (Variational Asymptotic 
    Beam Sectional Analysis) for detaile information please read the VABS user 
    manual
    
    Attributes:
    ----------
    format_flag : int, optional
        0 = old format, 1 = new format. (default = 0)
    
    nlayer : int, optional
        nlayer is not used in the old format. nlayer shoud be always given a 
        value greater than one if format_flag=1.
    
    Timoshenko_flag : int, optional
        VABS will construct both the classical model and the generalized 
        Timoshenko model. If it is 0, it will only construct the classical model. 
        (default = 1)
    
    recover_flag : int, optional
        0: VABS will carry out the consitutive modeling. 1: non-linear beam 
        theory, 3D stress, displacement, strain recovery, 2:linear beam theory 
        3D recovery. (default = 0)
        
    thermal_flag : int, optional
        Either 0 or 3: 0:pure mechanical analysis. 3:one-way coupled 
        thermoelastic analyis. (default = 0)
        
    curve_flag :  int 
        To model initially curved or twisted beams. if 1: three real numbers 
        for twist(k1) and curvatures (k2 and k3) shoud be provided. (default = 0)
        
    k1 : float 
        twist about x1 axis, see curve_flag
        
    k2 : float
        curvature about x2 axis, see curve_flag
    
    k3 : float
        curvature about x3 axis, see curve_flag
    
    trapeze_flag : int 
         To obtain the trapeze effect, trapeze flag is 1. (default = 0)    
    
    Vlasov_flag : int  
        To obtain a generalized Vlasov model, can be 1 only if Timoshenko_flag 
        is 1. For more information see VABS user manual. (default = 0)
        
    oblique_flag : int
        to model oblique cross sections. If oblique_flag is 1, two real numbers 
        are needed in the following line to specify the orientation of an 
        oblique reference cross section, see Figure 6 in the VABS user manual 
        for a sketch of such a cross section. (default = 0)
        This feature is only enabled in the classical beam model 
        (Timoshenko_flag = 0)
        
    oblique_cosine1 : float
        is the cosine of the angle between normal of the oblique section y1
        and beam axis x1.
    
    oblique_cosine2 : float
        is the cosine of the angle between y2 of the oblique section and beam 
        axis x1.
    
    u : nparray([u1, u2, u3])
        u1, u2 , and u3 are the 1D beam displacements along x1 , x2 , x3 , 
        respectively.ui and Cij are needed only for recovering 3D 
        displacements. If the user is not interested in 3D displacements,
        these values can be arbitrary real numbers. 
        
    Cij : nparray([C11, C12, C13],[C21,C22,C23],...] 
        The matrix Cij, with i = 1, 2, 3 and j = 1, 2, 3, is the direction 
        cosine matrix defined as Bi = Ci1b1 + Ci2 b2 + Ci3b3 with i = 1, 2, 3 
        where B1, B2, and B3 are the base vectors of the deformed triad and b1,
        b2, and b3 are the base vectors of the undeformed triad. Details of 
        this definition can be found in Ref. [10]. ui and Cij are needed only 
        for recovering 3D displacements. If the user is not interested in 3D 
        displacements, these values can be arbitrary real numbers. 
    
    F : nparray([F1, F2, F3]]) 
        F1 is the sectional aial force, F2 and F3 are the sectional transverse 
        shear forces along x2 and x3 respectively
        
    M : nparray([M1, M2, M3]]) 
        M1 is the sectional torque, M2 is the setional bending moment around x2 
        and M3 is the sectional bending moment around x3.
        
    f : nparray([f1, f2, f2]])
        f1, f2 , f3 are distributed forces (including both applied forces and 
        inertial forces) per unit span along x1, x2, x3 respectively.
         
    m : nparray([m1, m2, m2]])   
        m1 , m2 , m3 are distributed moment (including both applied and 
        inertial moments) per unit span along x1, x2, x3 respectively. 
    
    df : nparray([[f1', f2', f3']])
        first derivative of distributed loads with respect to beam axis
        The prime denotes derivative with respect to beam axis, 
        
    dm : nparray([[f1', f2', f3']])
        first derivative of distributed loads with respect to beam axis
        The prime denotes derivative with respect to beam axis, 
        
    ddf : nparray([[ f1'', f2'', f3'']])
        second derivative of distributed loads with respect to beam axis
        
    ddm : nparray([[m1'', m2'', m3'']])
        second derivative of distributed loads with respect to beam axis
        
    dddf : nparray([[ f1''', f2''', f3''']])
        third derivative of distributed loads with respect to beam axis
        
    dddf : nparray([[ m1''', m2''', m3''']])
        third derivative of distributed loads with respect to beam axis
    
    gamma11 : float
        is the beam axial stretching strain measure.
        
    kappa : nparray([[κ̄1 κ̄2 κ̄3]])
        κ̄1 is the twist measure, κ̄2 and κ̄3 are then curvature measures around 
        x2 and x3 respectively.
        
    dkappa1 : float
        The prime denotes derivative with respect to beam axis x1 of the twist 
        measure 
        
    ddkappa1 : float
        The primes denotes the second derivative with respect to beam axis x1 
        of the twist measure     
    """
    
    def __init__(self, **kw):
        self.format_flag = 0     
        self.nlayer = 0          
        self.Timoshenko_flag = 1 
        self.recover_flag = 0    
        self.thermal_flag = 0    
        self.curve_flag = 0     
        self.trapeze_flag = 0    
        self.Vlasov_flag = 0     
        self.oblique_flag = 0
        
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
            self.Timoshenko_flag = 0
            
        self.u = [0,0,0]                                    
        self.Cij = np.array([[0,0,0],[0,0,0],[0,0,0]])
        self.F = np.array([0,0,0])                #Timoshenko_flag == 0 -> F = [0]    
        self.M = np.array([0,0,0])
        self.f = np.array([0,0,0])                #Timoshenko_flag == 1:
        self.df = np.array([0,0,0])               #Timoshenko_flag == 1:
        self.ddf = np.array([0,0,0])               #Timoshenko_flag == 1:
        self.dddf = np.array([0,0,0])               #Timoshenko_flag == 1:
        self.m = np.array([0,0,0])                #Timoshenko_flag == 1:
        self.dm = np.array([0,0,0])               #Timoshenko_flag == 1:
        self.ddm = np.array([0,0,0])               #Timoshenko_flag == 1:
        self.dddm = np.array([0,0,0])               #Timoshenko_flag == 1:
        self.gamma11 = 0              #Vlasov_flag == 1:
        self.kappa = [0,0,0]            #Vlasov_flag == 1:
        self.dkappa1 = 0                #Vlasov_flag == 1:
        self.ddkappa1 = 0                #Vlasov_flag == 1:
        
if __name__ == '__main__':
    test = VABSConfig(recover_flag=1)