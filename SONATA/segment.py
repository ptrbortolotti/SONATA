import numpy as np

from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge,BRepBuilderAPI_MakeWire
from OCC.Geom import Geom_Plane
from OCC.gp import gp_Pnt, gp_Pnt2d, gp_Pln, gp_Dir, gp_Vec

from utils import  calc_DCT_angles, TColgp_HArray1OfPnt2d_from_nparray
from BSplineLst_utils import get_BSpline_length, get_BSplineLst_length, \
                            find_BSplineLst_coordinate, get_BSplineLst_Pnt2d, \
                            trim_BSplineLst, seg_boundary_from_dct, set_BSplineLst_to_Origin, \
                            copy_BSplineLst
from wire_utils import trim_wire, build_wire_from_BSplineLst
from readinput import UIUCAirfoil2d, AirfoilDat2d
from projection import layup_to_projection
from layer import Layer
class Segment(object):
    ''' 
    The Segment object is constructed from multiple Layers obejcts. 
    Each Segment has one Segment boundary.
    ''' 
   
    def __init__(self, ID = 0, **kwargs): #gets called whenever we create a new instance of a class
        ''' Initialize with BSplineLst:             Segment(ID, Layup = Configuration.Layup[1], CoreMaterial = Configuration.SEG_CoreMaterial[i], OCC=True, Boundary = BSplineLst)
            Initialize with airfoil database:       Segment(ID, Layup = Configuration.Layup[1], CoreMaterial = Configuration.SEG_CoreMaterial[i], OCC=False, airfoil = 'naca23012')
            Initialize from file:                   Segment(ID, Layup = Configuration.Layup[1], CoreMaterial = Configuration.SEG_CoreMaterial[i], OCC=False, filename = 'naca23012.dat')   
        #empty initialization with no Boundary: Segment(item, Layup = Configuration.Layup[1], CoreMaterial = Configuration.SEG_CoreMaterial[i])'''
        self.BSplineLst = None
        self.ID = ID
        self.Layup = kwargs.get('Layup') 
        self.CoreMaterial = kwargs.get('CoreMaterial') 
        self.OCC = kwargs.get('OCC')
        #self.wire  = []
        self.LayerLst = []

        if self.OCC == True:
            self.BSplineLst = kwargs.get('Boundary')
            #self.BSplineLst = set_BSplineLst_to_Origin(BSplineLst_tmp)
            
        elif self.OCC == False:
            if kwargs.get('airfoil') != None:
               BSplineLst_tmp = self.BSplineLst_from_airfoil_database(kwargs.get('airfoil'),30) 
            elif kwargs.get('filename') != None:
                BSplineLst_tmp = self.BSplineLst_from_file(kwargs.get('filename'),30)  
            self.BSplineLst = set_BSplineLst_to_Origin(BSplineLst_tmp)

        self.Projection = layup_to_projection(self.Layup)
    
    def __repr__(self):
        return '{}: {}'.format(self.__class__.__name__,self.ID) 
        
    def build_layers(self):
      
        for i in range(1,len(self.Layup)+1):
            print "STATUS:\t Building Segment %d, Layer: %d" % (self.ID,i)
            new_Boundary_BSplineLst = []
            #get_boundary_layer
            if i == 1:
                new_Boundary_BSplineLst = self.BSplineLst
            
            else:
                new_Boundary_BSplineLst += trim_BSplineLst(self.LayerLst[-1].Boundary_BSplineLst, 0, self.LayerLst[-1].S1, 0, 1)  #start und ende der lage
                new_Boundary_BSplineLst += copy_BSplineLst(self.LayerLst[-1].BSplineLst)
                new_Boundary_BSplineLst += trim_BSplineLst(self.LayerLst[-1].Boundary_BSplineLst, self.LayerLst[-1].S2, 1, 0, 1)  #start und ende der lage
                   
                new_Boundary_BSplineLst = set_BSplineLst_to_Origin(new_Boundary_BSplineLst)
        
            #CREATE LAYER Object
            tmp_Layer = Layer(i,new_Boundary_BSplineLst, self.Layup[i-1][0], self.Layup[i-1][1],self.Layup[i-1][2],self.Layup[i-1][3],self.Layup[i-1][4],cutoff_style= 1, join_style=1, name = 'test')   
            tmp_Layer.build_layer() 
            tmp_Layer.build_wire()
            self.LayerLst.append(tmp_Layer)     
    
    def copy(self):
        BSplineLstCopy =  copy_BSplineLst(self.BSplineLst)
        SegmentCopy = Segment(self.ID, Layup = self.Layup, CoreMaterial = self.CoreMaterial, OCC = True, Boundary = BSplineLstCopy)
        return SegmentCopy
                  
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
        

    def BSplineLst_from_airfoil_database(self,string,angular_deflection=30):
        '''
        string: 'naca23012'
        angular deflection: allowed angluar deflection in discrete representation before starting to split
        '''
        DCT_data = UIUCAirfoil2d(string).T
        self.BSplineLst = seg_boundary_from_dct(DCT_data,angular_deflection)
        return seg_boundary_from_dct(DCT_data,angular_deflection)
                                     
    def BSplineLst_from_file(self,filename,angular_deflection=30):
        '''
        filename: 'naca23012.dat'
        angular_deflection allowed angular deflection in discrete representation before starting to split default
        '''
        DCT_data = AirfoilDat2d(filename).T
        self.BSplineLst = seg_boundary_from_dct(DCT_data,angular_deflection)
        return seg_boundary_from_dct(DCT_data,angular_deflection)
        
    
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

            