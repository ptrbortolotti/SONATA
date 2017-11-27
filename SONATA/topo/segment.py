import numpy as np

from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge,BRepBuilderAPI_MakeWire
from OCC.Geom import Geom_Plane
from OCC.Geom2dAPI import Geom2dAPI_PointsToBSpline
from OCC.gp import gp_Pnt, gp_Pnt2d, gp_Pln, gp_Dir, gp_Vec


from SONATA.fileIO.readinput import UIUCAirfoil2d, AirfoilDat2d
from SONATA.bladegen.blade import Blade
from SONATA.fileIO.CADinput import import_2d_stp, import_3d_stp

from SONATA.topo.utils import  getID  
from SONATA.topo.utils import  calc_DCT_angles, TColgp_HArray1OfPnt2d_from_nparray, point2d_list_to_TColgp_Array1OfPnt2d
from SONATA.topo.BSplineLst_utils import get_BSpline_length, get_BSplineLst_length, \
                            find_BSplineLst_coordinate, get_BSplineLst_Pnt2d, \
                            trim_BSplineLst, seg_boundary_from_dct, set_BSplineLst_to_Origin, \
                            copy_BSplineLst, trim_BSplineLst_by_Pnt2d
from SONATA.topo.wire_utils import trim_wire, build_wire_from_BSplineLst
from SONATA.topo.projection import cummulated_layup_boundaries, relevant_cummulated_layup_boundaries,\
                                    plot_layup_projection, inverse_relevant_cummulated_layup_boundaries
from SONATA.topo.layer import Layer
from SONATA.topo.para_Geom2d_BsplineCurve import ParaLst_from_BSplineLst, BSplineLst_from_ParaLst

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
        self.cells = []

        if self.OCC == True:
            self.BSplineLst = kwargs.get('Boundary')
            #self.BSplineLst = set_BSplineLst_to_Origin(BSplineLst_tmp)
            
        elif self.OCC == False:
            if kwargs.get('airfoil') != None:
               BSplineLst_tmp = self.BSplineLst_from_airfoil_database(kwargs.get('airfoil'),30,self.scale_factor) 
            elif kwargs.get('filename') != None:
                BSplineLst_tmp = self.BSplineLst_from_file(kwargs.get('filename'),30,self.scale_factor)  
            self.BSplineLst = set_BSplineLst_to_Origin(BSplineLst_tmp)

        self.Projection = relevant_cummulated_layup_boundaries(self.Layup)
    
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

    
    def ivLst_to_BSplineLst(self,ivLst):
        '''The member function ivLst_to_BSplineLst generates the 
        BSplineLst from the InvervalLst definitions. It loops through all 
        intervals,trims them accordingly and assembles them into the 
        iv_BSplineLst, which is returned.
        
        Args:   self: the attributes of the Segment class
                iVLst: the iVLst is a entry of the Layup Projection. 
                It has the following form: array([[ 0.2  ,  0.3  ,  6.   ],
                                                  [ 0.3  ,  0.532,  8.   ]])
        returns: iv_BSplineLst: (list of BSplines)
        '''  
        
        iv_BSplineLst = []
        for iv in ivLst:
            #print iv[0],iv[1],iv[2]
            if int(iv[2]) == 0:
                BSplineLst = self.BSplineLst
                start = 0.0
                end = 1.0
            else:
                layer = self.LayerLst[int(iv[2])-1]
                BSplineLst = layer.BSplineLst
                start = layer.S1
                end = layer.S2

            
            S1 = iv[0]
            S2 = iv[1]
            iv_BSplineLst.extend(trim_BSplineLst(BSplineLst,S1,S2,start,end))
        return iv_BSplineLst
    
    
    def build_layers(self):
        '''The build_layers member function of the class Segment generates all Layer objects and it's associated wires
        and return the relevant_boundary_BSplineLst'''
        #plot_layup_projection(self.Layup)
        for i in range(1,len(self.Layup)+1):
            print "STATUS:\t Building Segment %d, Layer: %d" % (self.ID,i)
            relevant_boundary_BSplineLst = self.ivLst_to_BSplineLst(self.Projection[i-1])
          
            #CREATE LAYER Object
            tmp_Layer = Layer(i,relevant_boundary_BSplineLst, self.Layup[i-1][0], 
                              self.Layup[i-1][1],self.Layup[i-1][2],self.Layup[i-1][3],
                              self.Layup[i-1][4],cutoff_style= 2, join_style=1, name = 'test')   
            tmp_Layer.build_layer() 
            tmp_Layer.ivLst = self.Projection[i-1]
            tmp_Layer.inverse_ivLst = inverse_relevant_cummulated_layup_boundaries(self.Layup)[i-1]
            tmp_Layer.build_wire()
            self.LayerLst.append(tmp_Layer)     
    
        return relevant_boundary_BSplineLst
                
    def determine_final_boundary(self):
        '''The member function determin_final_boundary2 generates the 
        BSplineLst that encloses all Layers of the Segement. This final 
        boundary is needed for the generation of the subordinate segments, 
        where the final boundary BSplineLst is intersected with the Webs.
        
        Args:   self: only the attributes of the Segment class itself.
        
        returns: None, but assignes the final_Boundary_BSplineLst class 
            attribute
        '''    
        projectionlist=cummulated_layup_boundaries(self.Layup)
        self.final_Boundary_ivLst = projectionlist[-1]
        self.final_Boundary_BSplineLst = self.ivLst_to_BSplineLst(projectionlist[-1])
        return None


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
            

        
    def build_segment_boundary_from_WebLst(self,WebLst,Segment0_final_Boundary_BSplineLst):
        """Input has to be a complete WebLst"""
        print 'STATUS:\t Building Segment Boundaries %s' %(self.ID)
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
       

def generate_SegmentLst(Configuration):
    
    #generate Segment Lst
    SegmentLst = []   #List of Segment Objects
    for i,item in enumerate(Configuration.SEG_ID):
        if item == 0:        
            if Configuration.SETUP_input_type == 0:   #0) Airfoil from UIUC Database  --- naca23012
                SegmentLst.append(Segment(item, Layup = Configuration.SEG_Layup[i], CoreMaterial = Configuration.SEG_CoreMaterial[i], scale_factor = Configuration.SETUP_scale_factor, Theta = Configuration.SETUP_Theta, OCC=False, airfoil = Configuration.SETUP_datasource))
            
            elif Configuration.SETUP_input_type == 1: #1) Geometry from .dat file --- AREA_R250.dat
                SegmentLst.append(Segment(item, Layup = Configuration.SEG_Layup[i], CoreMaterial = Configuration.SEG_CoreMaterial[i], scale_factor = Configuration.SETUP_scale_factor,  Theta = Configuration.SETUP_Theta, OCC=False, filename = Configuration.SETUP_datasource))
            
            elif Configuration.SETUP_input_type == 2: #2)2d .step or .iges  --- AREA_R230.stp
                BSplineLst = import_2d_stp(Configuration.SETUP_datasource, Configuration.SETUP_scale_factor,Configuration.SETUP_Theta)
                SegmentLst.append(Segment(item, Layup = Configuration.SEG_Layup[i], CoreMaterial = Configuration.SEG_CoreMaterial[i], Theta = Configuration.SETUP_Theta, OCC=True, Boundary = BSplineLst))
            
            elif Configuration.SETUP_input_type == 3: #3)3D .step or .iges and radial station of crosssection --- AREA_Blade.stp, R=250
                BSplineLst = import_3d_stp(Configuration.SETUP_datasource,Configuration.SETUP_radial_station,Configuration.SETUP_scale_factor,Configuration.SETUP_Theta)
                SegmentLst.append(Segment(item, Layup = Configuration.SEG_Layup[i], CoreMaterial = Configuration.SEG_CoreMaterial[i], Theta = Configuration.SETUP_Theta, OCC=True, Boundary = BSplineLst))  
    
            elif Configuration.SETUP_input_type == 4: #4)generate 3D-Shape from twist,taper,1/4-line and airfoils, --- examples/UH-60A, R=4089, theta is given from twist distribution
                genblade = Blade(Configuration.SETUP_datasource,Configuration.SETUP_datasource,False,False)
                BSplineLst = genblade.get_crosssection(Configuration.SETUP_radial_station,Configuration.SETUP_scale_factor)
                SegmentLst.append(Segment(item, Layup = Configuration.SEG_Layup[i], CoreMaterial = Configuration.SEG_CoreMaterial[i], Theta = genblade.get_Theta(Configuration.SETUP_radial_station), OCC=True, Boundary = BSplineLst))  
                
            else:
                print 'ERROR:\t WRONG input_type'
     
        else:
            SegmentLst.append(Segment(item, Layup = Configuration.SEG_Layup[i], CoreMaterial = Configuration.SEG_CoreMaterial[i],Theta = Configuration.SETUP_Theta))
    sorted(SegmentLst, key=getID)  
    return SegmentLst
