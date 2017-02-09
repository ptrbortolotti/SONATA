import numpy as np

from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge,BRepBuilderAPI_MakeWire
from OCC.Geom import Geom_Plane
from OCC.Geom2dAPI import Geom2dAPI_PointsToBSpline
from OCC.gp import gp_Pnt, gp_Pnt2d, gp_Pln, gp_Dir, gp_Vec

from utils import  calc_DCT_angles, TColgp_HArray1OfPnt2d_from_nparray, point2d_list_to_TColgp_Array1OfPnt2d
from BSplineLst_utils import get_BSpline_length, get_BSplineLst_length, \
                            find_BSplineLst_coordinate, get_BSplineLst_Pnt2d, \
                            trim_BSplineLst, seg_boundary_from_dct, set_BSplineLst_to_Origin, \
                            copy_BSplineLst, trim_BSplineLst_by_Pnt2d
from wire_utils import trim_wire, build_wire_from_BSplineLst
from readinput import UIUCAirfoil2d, AirfoilDat2d
from projection import layup_to_projection
from layer import Layer
from para_Geom2d_BsplineCurve import ParaLst_from_BSplineLst, BSplineLst_from_ParaLst

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
        self.Theta = kwargs.get('Theta') 
        self.scale_factor  = kwargs.get('scale_factor') 
        self.LayerLst = []

        if self.OCC == True:
            self.BSplineLst = kwargs.get('Boundary')
            #self.BSplineLst = set_BSplineLst_to_Origin(BSplineLst_tmp)
            
        elif self.OCC == False:
            if kwargs.get('airfoil') != None:
               BSplineLst_tmp = self.BSplineLst_from_airfoil_database(kwargs.get('airfoil'),30,self.scale_factor) 
            elif kwargs.get('filename') != None:
                BSplineLst_tmp = self.BSplineLst_from_file(kwargs.get('filename'),30,self.scale_factor)  
            self.BSplineLst = set_BSplineLst_to_Origin(BSplineLst_tmp)

        self.Projection = layup_to_projection(self.Layup)
    
    def __repr__(self):
        return '{}: {}'.format(self.__class__.__name__,self.ID) 
    

    def __getstate__(self):
        """Return state values to be pickled."""
        self.Para_BSplineLst = ParaLst_from_BSplineLst(self.BSplineLst)
        return (self.ID, self.Layup, self.CoreMaterial, self.OCC, self.Theta, self.scale_factor, self.Projection, self.LayerLst, self.Para_BSplineLst)   
    
    
    def __setstate__(self, state):
        """Restore state from the unpickled state values."""
        self.ID, self.Layup, self.CoreMaterial, self.OCC, self.Theta, self.scale_factor, self.Projection, self.LayerLst,  self.Para_BSplineLst = state
        self.BSplineLst = BSplineLst_from_ParaLst(self.Para_BSplineLst)



    
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
                new_Boundary_BSplineLst = set_BSplineLst_to_Origin(new_Boundary_BSplineLst,self.Theta)
        
            #CREATE LAYER Object
            tmp_Layer = Layer(i,new_Boundary_BSplineLst, self.Layup[i-1][0], self.Layup[i-1][1],self.Layup[i-1][2],self.Layup[i-1][3],self.Layup[i-1][4],cutoff_style= 2, join_style=1, name = 'test')   
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
        

    def BSplineLst_from_airfoil_database(self,string,angular_deflection=30,scale_factor=1.0):
        '''
        string: 'naca23012'
        angular deflection: allowed angluar deflection in discrete representation before starting to split
        '''
        DCT_data = UIUCAirfoil2d(string).T
        DCT_data = np.multiply(DCT_data,scale_factor)
        self.BSplineLst = seg_boundary_from_dct(DCT_data,angular_deflection)
        return self.BSplineLst 
                                     
    def BSplineLst_from_file(self,filename,angular_deflection=30,scale_factor=1.0):
        '''
        filename: 'naca23012.dat'
        angular_deflection allowed angular deflection in discrete representation before starting to split default
        '''
        DCT_data = AirfoilDat2d(filename).T
        DCT_data = np.multiply(DCT_data,scale_factor)
        self.BSplineLst = seg_boundary_from_dct(DCT_data,angular_deflection)
        return seg_boundary_from_dct(DCT_data,angular_deflection)
            
    def determine_final_boundary(self):
        new_Boundary_BSplineLst = []
        new_Boundary_BSplineLst += trim_BSplineLst(self.LayerLst[-1].Boundary_BSplineLst, 0, self.LayerLst[-1].S1, 0, 1)  #start und ende der lage
        new_Boundary_BSplineLst += copy_BSplineLst(self.LayerLst[-1].BSplineLst)
        new_Boundary_BSplineLst += trim_BSplineLst(self.LayerLst[-1].Boundary_BSplineLst, self.LayerLst[-1].S2, 1, 0, 1)  #start und ende der lage
        new_Boundary_BSplineLst = set_BSplineLst_to_Origin(new_Boundary_BSplineLst)
        self.final_Boundary_BSplineLst = new_Boundary_BSplineLst
        
    def build_segment_boundary_from_WebLst(self,WebLst,Segment0_final_Boundary_BSplineLst):
        """Input has to be a complete WebLst"""
        print 'STATUS: \t Building Segment Boundaries %s' %(self.ID)
        NbOfWebs = len(WebLst)
        i = self.ID - 1    
            
        if self.ID == 0:
            None
        
        if self.ID == 1:
            #CREATE SEGMENT BOUNDARY 1
            P1 = WebLst[i].IntPnts_Pnt2d[0]
            P2 = WebLst[i].IntPnts_Pnt2d[1]
            trimmed_Boundary = trim_BSplineLst_by_Pnt2d(Segment0_final_Boundary_BSplineLst,P1,P2)
            Boundary_WEB_BSplineLst = [Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([P2,P1])).Curve().GetObject()]
            Boundary_BSplineLst = trimmed_Boundary + Boundary_WEB_BSplineLst                                  
            Boundary_BSplineLst = set_BSplineLst_to_Origin(Boundary_BSplineLst)
            self.BSplineLst = Boundary_BSplineLst
            self.build_wire()   
            
        elif self.ID == NbOfWebs+1:
            #CREATE LAST BOUNDARY
            P1 = WebLst[i-1].IntPnts_Pnt2d[0]
            P2 = WebLst[i-1].IntPnts_Pnt2d[1]
            trimmed_Boundary = trim_BSplineLst_by_Pnt2d(Segment0_final_Boundary_BSplineLst,P2,P1)
            Boundary_WEB_BSplineLst = [Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([P1,P2])).Curve().GetObject()]
            Boundary_BSplineLst = trimmed_Boundary + Boundary_WEB_BSplineLst                                  
            Boundary_BSplineLst = set_BSplineLst_to_Origin(Boundary_BSplineLst)
            self.BSplineLst = Boundary_BSplineLst
            self.build_wire()   
            
        else:
            #CREATE INTERMEDIATE BOUNDARIES
            P1 = WebLst[i-1].IntPnts_Pnt2d[0]
            P2 = WebLst[i-1].IntPnts_Pnt2d[1]
            P3 = WebLst[i].IntPnts_Pnt2d[0]
            P4 = WebLst[i].IntPnts_Pnt2d[1]
            Boundary_WEB_BSplineLst_1 = [Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([P1,P2])).Curve().GetObject()]
            Boundary_WEB_BSplineLst_2 = [Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([P4,P3])).Curve().GetObject()]
            trimmed_Boundary1 = trim_BSplineLst_by_Pnt2d(Segment0_final_Boundary_BSplineLst,P3,P1)
            trimmed_Boundary2 = trim_BSplineLst_by_Pnt2d(Segment0_final_Boundary_BSplineLst,P2,P4)
            Boundary_BSplineLst = trimmed_Boundary1 + Boundary_WEB_BSplineLst_1 + trimmed_Boundary2 + Boundary_WEB_BSplineLst_2
            Boundary_BSplineLst = set_BSplineLst_to_Origin(Boundary_BSplineLst)
            self.BSplineLst = Boundary_BSplineLst
            self.build_wire()
        