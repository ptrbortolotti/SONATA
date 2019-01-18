#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 09:38:19 2018

@author: Tobias Pflumm
"""
import os
from OCC.gp import gp_Ax2, gp_Pnt, gp_Dir, gp_Ax1
import numpy as np
import math
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

if __name__ == '__main__':
    os.chdir('..')
    
from SONATA.classComponent import Component
from SONATA.classAirfoil import Airfoil
from SONATA.cbm.topo.wire_utils import rotate_wire, mirror_wire_pnt_dir, translate_wire, scale_wire
from SONATA.cbm.bladegen.bladegen_utils import make_loft

class Blade(Component):

    __slots__ = ('coordinates', 'chord', 'twist', 'airfoil_position', 'blade_matrix', 'airfoilLst')
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args,**kwargs)

    def __repr__(self):
        """__repr__ is the built-in function used to compute the "official" 
        string reputation of an object, """
        return 'Blade: '+ self.name
    
    def read_IAE37(self, yml):
        """
        reads the IAE Wind Task 37 style Blade dictionary and assigsn them to
        the class attribute      
        """
        self.coordinates = {}
        tmp = byml.get('bem_aero')
        self.coordinates['x'] = np.asarray((tmp.get('coordinates').get('x').get('grid'),tmp.get('coordinates').get('x').get('values'))).T
        self.coordinates['y'] = np.asarray((tmp.get('coordinates').get('y').get('grid'),tmp.get('coordinates').get('y').get('values'))).T
        self.coordinates['z'] = np.asarray((tmp.get('coordinates').get('z').get('grid'),tmp.get('coordinates').get('z').get('values'))).T
        self.twist = np.asarray((tmp.get('twist').get('grid'),tmp.get('twist').get('values'))).T
        self.chord = np.asarray((tmp.get('chord').get('grid'),tmp.get('chord').get('values'))).T
        self.airfoil_position = (tmp.get('airfoil_position').get('grid'),tmp.get('airfoil_position').get('labels'))
        
    def gen_blade_matrix(self, airfoils):
        tmp = []
        for an in self.airfoil_position[1]: 
            tmp.append(next((x for x in airfoils if x.name == an), None).id)
        arr = np.asarray([self.airfoil_position[0],tmp]).T
            
        f_chord = interp1d(self.chord[:,0], self.chord[:,1],bounds_error=False,fill_value='extrapolate')
        f_twist = interp1d(self.twist[:,0], self.twist[:,1],bounds_error=False,fill_value='extrapolate')
        f_coordinates_x = interp1d(self.coordinates['x'][:,0], self.coordinates['x'][:,1],bounds_error=False,fill_value='extrapolate')
        f_coordinates_y = interp1d(self.coordinates['y'][:,0], self.coordinates['y'][:,1],bounds_error=False,fill_value='extrapolate')
        f_coordinates_z = interp1d(self.coordinates['z'][:,0], self.coordinates['z'][:,1],bounds_error=False,fill_value='extrapolate')
        #f_airfoil = interp1d(arr[:,0], arr[:,1],bounds_error=False,fill_value='extrapolate')
        
        x = np.unique(np.sort(np.hstack((self.chord[:,0],self.twist[:,0],self.coordinates['x'][:,0],self.coordinates['y'][:,0],self.coordinates['z'][:,0], arr[:,0]))))
        self.blade_matrix = np.transpose(np.unique(np.array([x,f_chord(x),f_twist(x),f_coordinates_x(x),f_coordinates_y(x),f_coordinates_z(x)]), axis=1))
        self.airfoilLst = [interp_airfoil_position(self.airfoil_position, airfoils, x) for x in x]
        return self.blade_matrix, self.airfoilLst
        
    def plot_blade_matrix(self, airfoils):
        plt.close("all")   
        plt.figure()
        plt.subplot(321)
        plt.plot(self.blade_matrix[:,0], self.blade_matrix[:,1], 'r.-')
        plt.ylabel('chord [mm]')
        
        plt.subplot(323)
        plt.plot(self.blade_matrix[:,0], self.blade_matrix[:,2], 'g.-')
        plt.ylabel('twist [deg]')
        
        plt.subplot(325)
        plt.plot(self.blade_matrix[:,0], self.blade_matrix[:,6], 'b.-')
        plt.ylabel('airfoil')
        strLst = []
        idLst = []
        for af in airfoils:
            idLst.append(af.id)
            strLst.append(af.name)
        plt.yticks(idLst, strLst)
        
        plt.subplot(322)
        plt.plot(self.blade_matrix[:,0], self.blade_matrix[:,3], 'k.-')
        plt.ylabel('x of 1/4-line [mm]')
        
        plt.subplot(324)
        plt.plot(self.blade_matrix[:,0], self.blade_matrix[:,4], 'k.-')
        plt.ylabel('y of 1/4-line [mm]')
        
        plt.subplot(326)
        plt.plot(self.blade_matrix[:,0], self.blade_matrix[:,5], 'k.-')
        plt.ylabel('z of 1/4-line [mm]')
        
        plt.show()
    
    def gen_surface(self):
        wireLst = []
        for bm,afl in zip(self.blade_matrix, self.airfoilLst):
            if afl.wire == None:
                afl.gen_OCCtopo()
            
            wire = afl.wire
            wire = rotate_wire(wire, gp_Ax1(gp_Pnt(0,0,0), gp_Dir(1,0,0)), np.pi/2)
            wire = rotate_wire(wire, gp_Ax1(gp_Pnt(0,0,0), gp_Dir(0,0,1)), -np.pi/2)
            wire = translate_wire(wire, gp_Pnt(0,0,0), gp_Pnt(0,0.25,0))
            wire = scale_wire(wire, gp_Pnt(0,0,0), bm[1])

            wire = rotate_wire(wire, gp_Ax1(gp_Pnt(0,0,0), gp_Dir(1,0,0)), bm[2])  
            wire = translate_wire(wire, gp_Pnt(0,0,0), gp_Pnt(bm[3],bm[4],bm[5]))
            display.DisplayShape(wire, color='BLACK')
            wireLst.append(wire)
            
        return make_loft(wireLst[:12], ruled=False, tolerance=1e-5, continuity=4, check_compatibility=True)
    
def interp_airfoil_position(airfoil_position, airfoils, grid_loc):
    #TBD: Extrapolation
    if grid_loc in airfoil_position[0]:
        afname = airfoil_position[1][airfoil_position[0].index(grid_loc)]
        return next((x for x in airfoils if x.name == afname), None)

    #find closest value:
    min_idx = np.argmin([abs(x-grid_loc) for x in airfoil_position[0]])
    min_val = airfoil_position[0][min_idx]
    
    if grid_loc > min_val:
        iv_idx = (min_idx,min_idx+1)
    else:
        iv_idx = (min_idx-1, min_idx)
    
    iv_val = tuple(airfoil_position[0][iv_idx[0]:iv_idx[1]+1])
    iv_af = tuple(airfoil_position[1][iv_idx[0]:iv_idx[1]+1])
    k = (grid_loc-iv_val[0]) / (iv_val[1]-iv_val[0])
    
    #select af from airfoils
    af1 = next((x for x in airfoils if x.name == iv_af[0]), None)
    af2 = next((x for x in airfoils if x.name == iv_af[1]), None)
    
    if af1 == af2:
        return af1
    
    #return transformed airfoil
    return af1.transformed(af2,k,200)
    
if __name__ == '__main__':

    import yaml
    from OCC.Display.SimpleGui import init_display
    from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeEdge, BRepBuilderAPI_Transform
    from OCC.gp import gp_Vec, gp_Trsf2d, gp_Pnt2d
    from OCC.Geom import Geom_Plane
    from SONATA.cbm.display.display_utils import show_coordinate_system
    
    display, start_display, add_menu, add_function_to_menu = init_display('qt-pyqt5')
    display.Context.SetDeviationAngle(1e-5)       # 0.001 default. Be careful to scale it to the problem.
    display.Context.SetDeviationCoefficient(1e-5) #
        
    with open('jobs/PBortolotti/IEAonshoreWT.yaml', 'r') as myfile:
        inputs  = myfile.read()
    with open('jobs/PBortolotti/IEAturbine_schema.yaml', 'r') as myfile:
        schema  = myfile.read()
    
    yml = yaml.load(inputs)
    yml_schema = yaml.load(schema)
    
    airfoils = [Airfoil(af) for af in yml.get('airfoils')]
    
    byml = yml.get('components').get('blade')
    B = Blade(gp_Pnt(0,10,0), gp_Dir(0,0,1), name='TestBlade')
    B.read_IAE37(byml)
    bm,afl = B.gen_blade_matrix(airfoils)
    loft = B.gen_surface()

    display.DisplayShape(loft)
    show_coordinate_system(display,1)

    display.FitAll()
    start_display()   
    