#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 13:04:03 2019

@author: gu32kij

Utility Module to calculate the approximate elastic properties for 
unidirectional laminates according to the semi-empiric functions of Schuermann

"""


import yaml

#TODO: Implement more models for elastic properties of Composites \cite{Lamers -Review on Micromechanical Modelling (1999, Twente)}
#TODO: Halpin-Tsai and CCA model!
#TODO: Implement model for linear thermal expasion properties
#TODO: Make Puck, Halpin-Tsai and CCA Model a function that needs Fiber and Matrix instances + FVC they return a Orthotropic Material !!!

#Fiber:
class Fiber(object):
    def __init__(self,rho = 1.74, E_1 = 230, E_2 = 28, G_12 = 50, nu_21 = 0.23, name = 'cfht', description = 'ht carbon fiber', source = 'Schuermann p.42'):
       self.rho = rho
       self.E_1 = E_1 #GPa
       self.E_2 = E_2 #GPa
       self.G_12 = G_12 #-
       self.nu_21 = nu_21
       self.name = name
       self.description = description
       self.source = source
    
    def __repr__(self):
        """__repr__ is the built-in function used to compute the "official" 
        string reputation of an object, """
        return 'Fiber: '+ self.name
            
#Matrix:
class Resin(object):
    def __init__(self, rho = 1.23, E = 3.4, nu=0.35, name = 'ep', description = 'epoxy resin', source = 'nosource'):
       self.rho = rho
       self.E = E #GPa
       self.nu = nu #-
       self.name = name
       self.description = description
       self.source = source
       
    def __repr__(self):
        """__repr__ is the built-in function used to compute the "official" 
        string reputation of an object, """
        return 'Resin: '+ self.name
       
    @property
    def G(self): 
       return self.E/(2*(1+self.nu))

#Composite: 
class Puck_UDLaminate(object):
    def __init__(self, fiber, resin, phi = 0.6, description = 'nodescription', source = 'nosource'):

        self.fiber = fiber
        self.resin = resin
        self.phi = phi
        self.description = description
        self.source = source
      
    def __repr__(self):
        """__repr__ is the built-in function used to compute the "official" 
        string reputation of an object, """
        return 'UDLaminate: '+ self.name
        
        
    @property
    def name(self):
        str_k = '%.2f' % self.phi
        return 'ud_'+self.fiber.name+'_'+self.resin.name+'_'+str_k.replace('.','')
    
    @property
    def rho(self):
        m = self.resin 
        f = self.fiber
        return self.phi*f.rho + (1-self.phi)*m.rho
        
    #EModuli
    @property
    def E_1(self):
        return self.fiber.E_1*self.phi + (self.resin.E*(1-self.phi))
    
    @property
    def E_2(self):
       m = self.resin 
       f = self.fiber
       E_2 = (m.E / (1-m.nu**2)) * (1 + 0.85*self.phi**2)/( (1-self.phi)**1.25 + (m.E / ((1-m.nu**2)*f.E_2))*self.phi)
       return E_2
   
    @property
    def E_3(self):
        return self.E_2
   
    #Shear-Modul
    @property
    def G_12(self):
        m = self.resin 
        f = self.fiber
        G_12 = m.G * (1+0.4*self.phi**0.5)/((1-self.phi)**1.45 + (m.G/f.G_12)*self.phi)
        return G_12
    
    @property
    def G_13(self):
        return self.G_12
    
    @property
    def G_23(self):
        return self.E_2/(1*(1+self.nu_23))
    
    #Poisson-Ratios
    @property
    def nu_12(self):
        nu_12 = self.E_2 * self.nu_21 / self.E_1
        return nu_12
    
    @property
    def nu_21(self):
        m = self.resin 
        f = self.fiber
        nu_21 = self.phi * f.nu_21 + (1-self.phi)*m.nu
        return nu_21
    
    @property
    def nu_13(self):
        return self.nu_12
    
    @property
    def nu_23(self):
        m = self.resin 
        f = self.fiber
        nu_meff = m.nu*((1+m.nu - self.nu_21*m.E/self.E_1)/(1 - m.nu**2 + m.nu*self.nu_21*m.E/self.E_1))
        nu_23 = self.phi*f.nu_21 + (1-self.phi)*nu_meff
        return nu_23
    
    def write_yaml(self, fname=None):
        yml = {}
        yml['name'] = self.name
        yml['description'] = self.description
        yml['source'] = self.source
        yml['rho'] = self.rho
        yml['E_1'] = self.E_1
        yml['E_2'] = self.E_2
        yml['E_3'] = self.E_3
        
        yml['G_12'] = self.G_12
        yml['G_13'] = self.G_13
        yml['G_23'] = self.G_23
        
        yml['nu_12'] = self.nu_12
        yml['nu_13'] = self.nu_13
        yml['nu_23'] = self.nu_23

        if fname:            
            with open(fname, 'w') as f:
                yaml.dump(yml, f, default_flow_style=False)
            return None
        else:
            return yml
    
    
if __name__ == '__main__':
    def float_representer(dumper, value):
        text = '{0:.3f}'.format(value)
        return dumper.represent_scalar(u'tag:yaml.org,2002:float', text)
    yaml.add_representer(float, float_representer)
    
    
    
    ep = Resin(rho = 1.23, E = 3.4, nu=0.35, name = 'ep', description = 'epoxy resin', source = 'Schuermann')
    cf_ht = Fiber(rho = 1.74, E_1 = 230, E_2 = 28, G_12 = 50, nu_21 = 0.23, name = 'cfht', description = 'ht carbon fiber', source = 'Schuermann p.42')
    cf_hm = Fiber(rho=1.81, E_1=392, E_2 = 15.2, G_12 = 28.6, nu_21 = 0.2, name='cf_hm',  source = 'Schuermann p.42')
    cf_im = Fiber(rho=1.8, E_1=294, name='cf_im',  source = 'Schuermann p.42')
    gf_e =  Fiber(rho=2.54, E_1=73, E_2 = 73, G_12 = 29.92, nu_21= 0.22, name='gf_e')
    gf_s =  Fiber(rho=2.49, E_1=86.81, E_2 = 86.81, G_12 = 35.578, nu_21= 0.22, name='gf_s')
    af_hm = Fiber(rho=1.45, E_1=130, E_2 = 5.4, G_12=1.45, nu_21=0.32,  name='af_hm', source = 'Schuermann p.45')
    
    cf_ht_ep = Puck_UDLaminate(cf_ht, ep).write_yaml('cf_ht_ep.yml')
    cf_hm_ep = Puck_UDLaminate(cf_hm, ep).write_yaml('cf_hm_ep.yml')
    cf_im_ep = Puck_UDLaminate(cf_im, ep).write_yaml('cf_im_ep.yml')
    gf_e_ep = Puck_UDLaminate(gf_e, ep).write_yaml('gf_e_ep.yml')
    gf_s_ep = Puck_UDLaminate(gf_s, ep).write_yaml('gf_s_ep.yml')
    af_hm_ep = Puck_UDLaminate(af_hm, ep).write_yaml('af_hm_ep.yml')