import numpy as np

from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge,BRepBuilderAPI_MakeWire
from OCC.Geom import Geom_Plane
from OCC.gp import gp_Pnt, gp_Pnt2d, gp_Pln, gp_Dir, gp_Vec

from utils import  calc_DCT_angles, TColgp_HArray1OfPnt2d_from_nparray
from BSplineLst_utils import get_BSpline_length, get_BSplineLst_length, find_BSplineLst_coordinate, get_BSplineLst_Pnt2d, trim_BSplineLst, seg_boundary_from_dct  
from wire_utils import trim_wire, build_wire_from_BSplineLst
from readinput import UIUCAirfoil2d, AirfoilDat2d
class Segment(object):
    ''' 
    The Segment object is constructed from multiple Layers obejcts. 
    Each Segment has one Segment boundary.
    ''' 
   
    def __init__(self, ID = 0, **kwargs): #gets called whenever we create a new instance of a class
        ''' Initialize with BSplineLst:             Segment(item, Layup = Configuration.Layup[1], CoreMaterial = Configuration.SEG_CoreMaterial[i], OCC=True,  Boundary = BSplineLst)
            Initialize with airfoil database:       Segment(item, Layup = Configuration.Layup[1], CoreMaterial = Configuration.SEG_CoreMaterial[i], OCC=False, airfoil = 'naca23012')
            Initialize from file:                   Segment(item, Layup = Configuration.Layup[1], CoreMaterial = Configuration.SEG_CoreMaterial[i], OCC=False, filename = 'naca23012.dat')   
        #empty initialization with no Boundary: Segment(item, Layup = Configuration.Layup[1], CoreMaterial = Configuration.SEG_CoreMaterial[i])'''
        self.BSplineLst = None
        self.ID = ID
        self.Layup = kwargs.get('Layup') 
        self.CoreMaterial = kwargs.get('CoreMaterial') 
        self.OCC = kwargs.get('OCC')

        if self.OCC == True:
            self.BSplineLst = kwargs.get('Boundary')

        elif self.OCC == False:
            if kwargs.get('airfoil') != None:
                self.BSplineLst = self.BSplineLst_from_airfoil_database(kwargs.get('airfoil'),140) 
                
            elif kwargs.get('filename') != None:
                self.BSplineLst = self.BSplineLst_from_file(self,kwargs.get('filename'),140)  
                
        #self.wire  = []
        self.LayerLst = []

    def __repr__(self):
        return '{}: {}'.format(self.__class__.__name__,self.ID) 
       
                  
    def get_length(self): #Determine and return Legth of Layer self
         self.length = get_BSplineLst_length(self.BSplineLst)
         return self.length
        
    def pnt2d(self,S): #Return, gp_Pnt2d of argument S of Segment self  
        return get_BSplineLst_Pnt2d(self.BSplineLst,S)
    
    def build_wire(self): #Builds TopoDS_Wire from connecting BSplineSegments and returns it                  
        self.wire = build_wire_from_BSplineLst(self.BSplineLst)           
        
    def trim(self,S1,S2, start, end): #Trims layer between S1 and S2
        return trim_BSplineLst(self.BSplineLst, S1, S2, start, end)
    
    def trim_SEGwire(self, S1, S2):
        return trim_wire(self.wire, S1, S2)
        

    def BSplineLst_from_airfoil_database(self,string,min_degree=140):
        '''
        string: 'naca23012'
        min_degree: allowed angle in discrete representation before starting to split
        '''
        DCT_data = UIUCAirfoil2d(string).T
        self.BSplineLst = seg_boundary_from_dct(DCT_data,min_degree)
        return seg_boundary_from_dct(DCT_data,min_degree)
 
            
    def BSplineLst_from_file(self,filename,min_degree=140):
        '''
        filename: 'naca23012.dat'
        min_degree: allowed angle in discrete representation before starting to split default
        '''
        DCT_data = AirfoilDat2d(filename).T
        self.BSplineLst = seg_boundary_from_dct(DCT_data)
        return seg_boundary_from_dct(DCT_data)
        
    
    def set_to_origin(self):
        pass

    
    def show(self): #display the layer with pythonocc viewer module
        """
        TBD: renders the topological entity in the viewer: 
        """
#        if not self.display_set:
#            display.Display()(self, *args, **kwargs)
#        else:
#            self.disp.DisplayShape(*args, **kwargs)

            
