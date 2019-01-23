# -*- coding: utf-8 -*-
"""
Created on Mon Sep 04 10:45:43 2017

@author: TPflumm
testbed for the bladegen tool!
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import math

#PythonOCC Libraries
from OCC.Display.SimpleGui import init_display
from OCC.gp import gp_Pnt, gp_Dir, gp_Ax1

#SONATA Modules
if __name__ == '__main__':
    os.chdir('../../..')


from SONATA.cbm.bladegen.airfoil import Airfoil
from SONATA.blade_utl import make_loft
from SONATA.cbm.topo.wire_utils import rotate_wire, translate_wire,scale_wire, mirror_wire_pnt_dir
from SONATA.cbm.display.display_utils import show_coordinate_system 
from SONATA.cbm.fileIO.CADinput import BSplineLst_from_intersect_shape


class Blade(object):
    def __init__(self, folder, name=None, FLAG_SHOW_3D_TOPO=False, FLAG_SHOW_2D_DATA=False):
        self.folder = folder
        if name:
            self.name = name
        else:
            self.name = folder
        self.FLAG_SHOW_3D_TOPO = FLAG_SHOW_3D_TOPO
        self.FLAG_SHOW_2D_DATA = FLAG_SHOW_2D_DATA
        self.airfoilLst = []
        self.blade_matrix = self.read_and_interpolate_blade_matrix()
        self.surface = self.design_surface()
        self.BSplineLst = []
        
    def read_and_interpolate_blade_matrix(self):
        #READ DATA FROM FILES and interpolat        
        twist = np.loadtxt(self.folder + '/twist.dat')
        chord = np.loadtxt(self.folder + '/chord.dat')
        quarterline = np.loadtxt(self.folder + '/quarter_line.dat')
        foildata = np.loadtxt(self.folder + '/airfoil.dat',dtype='str')
        self.airfoilLst = []
        airfoil_tmp = Airfoil()
        airfoil_tmp.restore_counter()
        for i,af in enumerate(np.unique(foildata[:,1])):
            self.airfoilLst.append(Airfoil(af))        
       
        for i,af in enumerate(foildata):
            for x in self.airfoilLst:
                if x.name == af[1]:
                    foildata[i,1] = x.id
    
        airfoil = np.asarray(foildata,dtype=float)
        f_chord = interp1d(chord[:,0], chord[:,1],bounds_error=False,fill_value='extrapolate')
        f_twist = interp1d(twist[:,0], twist[:,1],bounds_error=False,fill_value='extrapolate')
        f_quarterline_y = interp1d(quarterline[:,0], quarterline[:,1],bounds_error=False,fill_value='extrapolate')
        f_quarterline_z = interp1d(quarterline[:,0], quarterline[:,2],bounds_error=False,fill_value='extrapolate')
        f_airfoil = interp1d(airfoil[:,0], airfoil[:,1],bounds_error=False,fill_value='extrapolate')
        
        x = np.sort(np.hstack((chord[:,0],twist[:,0],airfoil[:,0])))
        return np.transpose(np.unique(np.array([x,f_chord(x),f_twist(x),f_quarterline_y(x),f_quarterline_z(x),f_airfoil(x)]), axis=1))

    
    def design_surface(self):
    #DESIGN HOMOGENEOUS SECTION   
        wireLst = []
        for x in self.blade_matrix:
            wire = self.airfoilLst[int(x[5])].wire
            wire = mirror_wire_pnt_dir(wire,gp_Pnt(0,0,0),gp_Dir(0,0,1))
            wire = translate_wire(wire,gp_Pnt(0,0,0),gp_Pnt(0,0.25,0))
            wire = scale_wire(wire,gp_Pnt(0,0,0),x[1])
            wire = rotate_wire(wire,gp_Ax1(gp_Pnt(0,0,0),gp_Dir(1,0,0)),math.radians(x[2]))  
            wire = translate_wire(wire,gp_Pnt(0,0,0),gp_Pnt(x[0],x[3],x[4]))
            wireLst.append(wire)
    
        return make_loft(wireLst,ruled=True, tolerance=1e-6, continuity=1, check_compatibility=True)
    
    def get_crosssection(self,R,scale_factor=1):
        aResShape = self.surface
        self.BSplineLst =  BSplineLst_from_intersect_shape(aResShape,R,scale_factor,self.get_Theta(R))
        return self.BSplineLst
            
    
    def get_Theta(self,R):
        f_twist = interp1d(self.blade_matrix[:,0], self.blade_matrix[:,2],bounds_error=False,fill_value='extrapolate')
        return f_twist(R)
        
        
    def display(self,display):
        #==========PLOT=====================
        if self.FLAG_SHOW_2D_DATA:
            plt.close("all")   
            plt.figure()
            plt.subplot(321)
            plt.plot(self.blade_matrix[:,0], self.blade_matrix[:,1], 'r-o')
            plt.ylabel('chord [mm]')
            
            plt.subplot(323)
            plt.plot(self.blade_matrix[:,0], self.blade_matrix[:,2], 'g-o')
            plt.ylabel('twist [deg]')
            
            plt.subplot(325)
            plt.plot(self.blade_matrix[:,0], self.blade_matrix[:,5], 'g-o')
            plt.ylabel('airfoil')
            strLst = []
            idLst = []
            for af in self.airfoilLst:
                idLst.append(af.id)
                strLst.append(af.name)
            plt.yticks(idLst, strLst)
            
            plt.subplot(322)
            plt.plot(self.blade_matrix[:,0], self.blade_matrix[:,3], 'b-o')
            plt.ylabel('y of 1/4-line [mm]')
            
            plt.subplot(324)
            plt.plot(self.blade_matrix[:,0], self.blade_matrix[:,4], 'b-o')
            plt.ylabel('z of 1/4-line [mm]')
            
            plt.show()
        
        #OCC DISPLAY
        if self.FLAG_SHOW_3D_TOPO:
            display.Context.SetDeviationAngle(1e-5)       # 0.001 default. Be careful to scale it to the problem.
            display.Context.SetDeviationCoefficient(1e-5) # 0.001 default. Be careful to scale it to the problem. 
                        
            display.DisplayShape(self.surface)
            display.set_bg_gradient_color(20,6,111,200,200,200)
            show_coordinate_system(display,50)
            
            for spline in self.BSplineLst:
                display.DisplayShape(spline)

        return None
    
#===========================MAIN===============================================
if __name__ == '__main__':
    UH60A_blade = Blade('examples/UH-60A','UH-60A',True,True)
    display, start_display, add_menu, add_function_to_menu = init_display('qt-pyqt5')
    BSplineLst = UH60A_blade.get_crosssection(7600)
    UH60A_blade.display(display)
    TopoDS_Shape = UH60A_blade.surface
    
    display.DisplayMessage(gp_Pnt(100,0,0),'UH-60A Main Rotor Blade')
            
    display.FitAll()
    start_display()   