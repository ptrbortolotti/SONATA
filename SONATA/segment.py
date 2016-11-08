import numpy as np

from OCC.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_CompCurve
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge,BRepBuilderAPI_MakeWire
from OCC.BRepLib import BRepLib_MakeFace 
from OCC.GCPnts import GCPnts_AbscissaPoint
from OCC.Geom import Geom_Plane
from OCC.Geom2dAPI import Geom2dAPI_Interpolate
from OCC.gp import gp_Pnt, gp_Pnt2d, gp_Pln, gp_Dir, gp_Vec
from OCC.IntCurvesFace import IntCurvesFace_Intersector

from utils import (UIUCAirfoil2d, calc_DCT_angles, TColgp_HArray1OfPnt2d_from_nparray,
                   AirfoilDat2d)
from explorer import WireExplorer
from segment_utils import *
  

class Segment(object):
    ''' 
    The Segment object is constructed from multiple Layers obejcts. 
    Each Segment has one Segment boundary.
    ''' 
   
    def __init__(self): #gets called whenever we create a new instance of a class
        self.dct   = []
        self.BSplineLst  = []
        self.wire  = []


    #def __repr__(self): #we can tell Python how to prepresent an object of our calss (when using a print statement)
        #pass
            
    def length(self): #Determine and return Legth of Segment self
        pass
 
    def pnt(self,S): #Return, gp_Pnt of argument S of Segment self  
        pass
    
    def pnt2d(self,S): #Return, gp_Pnt2d of argument S of Segment self  
        pass
    
    def build_wire(self): #Builds TopoDS_Wire from connecting BSplineSegments and returns it        
        tmp_wire = 	BRepBuilderAPI_MakeWire()
        for i,item in enumerate(self.BSplineLst):
            #Convert to 3D SPACE onto x,y plane
            P = gp_Pnt(0,0,0)
            V = gp_Dir(gp_Vec(0,0,1))
            Plane = Geom_Plane(P, V)
            tmp_edge = BRepBuilderAPI_MakeEdge(item.GetHandle(),Plane.GetHandle())
            tmp_wire.Add(tmp_edge.Edge())
            
        self.wire = tmp_wire.Wire()
            
        
    def trim(self,S1,S2): #Trims layer between S1 and S2
        pass
    
    def first():
        pass
    
    def last():
        pass
    
 
    def check():
        pass
              
    
    def BSplineLst_from_airfoil_database(self,string,min_degree=140):
        '''
        string: 'naca23012'
        min_degree: allowed angle in discrete representation before starting to split
        '''
        DCT_data = UIUCAirfoil2d(string).T
        self.BSplineLst = seg_boundary_from_dct(DCT_data,min_degree)
       
            
    def BSplineLst_from_file(self,filename,min_degree=140):
        '''
        filename: 'naca23012.dat'
        min_degree: allowed angle in discrete representation before starting to split default
        '''
        DCT_data = AirfoilDat2d(filename).T
        self.BSplineLst = seg_boundary_from_dct(DCT_data)
        
    
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

            
