import matplotlib.pyplot as plt

from OCC.gp import gp_Pnt2d

from SONATA.fileIO.CADinput import order_BSplineLst_Head2Tail

from SONATA.topo.BSplineLst_utils import get_BSplineLst_length, get_BSplineLst_Pnt2d, trim_BSplineLst, copy_BSplineLst, BSplineLst_from_dct, discretize_BSplineLst 
from SONATA.topo.wire_utils import build_wire_from_BSplineLst
from SONATA.topo.cutoff import cutoff_layer
from SONATA.topo.offset import shp_parallel_offset
from SONATA.topo.para_Geom2d_BsplineCurve import ParaLst_from_BSplineLst, BSplineLst_from_ParaLst

from SONATA.mesh.mesh_utils import equidistant_nodes_on_BSplineLst, grab_nodes_on_BSplineLst
class Layer(object):
    ''' 
    The layer object is constructed from multiple BSplineCurveSegments. It is the basis for all future operations. 
    The object can be constructed from either a discrete formulation of point tables or from an existing TopoDS_Wire.
    ''' 
   
    def __init__(self, ID, Boundary_BSplineLst, globalStart, globalEnd, thickness, Orientation = 0, MatID = 1, **kwargs):
        self.ID = ["%04d" % ID] 	                          #First single Digit: SegmentNb, Last 3 Digits: LayerNb; e.g.: 1029, Segment1, Layer29
        self.Boundary_BSplineLst = Boundary_BSplineLst        #List of Geom2d_BSplineCurve, Geom_BSplineCurve 
        self.S1 = globalStart	                              #Starting Point in S coordinates
        self.S2 = globalEnd		                              #End Point in S coordinates
        self.thickness = thickness   	                      #in mm
        self.Orientation = Orientation                        #in deg
        self.MatID = MatID
        self.cells = []
        self.ivLst = []
        self.inverse_ivLst = []
        self.a_nodes = []
        self.b_nodes = []
        self.cells = [] #Container to collet all cells that are composing this layere

        #KWARGS:
        if kwargs.get('name') == None:
             self.name = 'DEFAULT'                             
        else:
            self.name = kwargs.get('name')
        
        if (kwargs.get('cutoff_style') == None) or (type(kwargs.get('cutoff_style')) is not int):   #cutoff_style (step, linear, smooth_bezier)
             self.cutoff_style = 2                             
        else:
            self.cutoff_style = kwargs.get('cutoff_style')      
        
        
        if (kwargs.get('join_style') == None) or (type(kwargs.get('join_style')) is not int):      #offset algorithm join_style = 1#( 1:round,2:mitre,3:bevels)
             self.join_style = 1                             
        else:
            self.join_style = kwargs.get('join_style')          
            

        #self.wire = []                                     #Make Wire from BSplineSegments
        #self.display_set = False  
    

    
    @property
    def StartPoint(self): #gp_Pnt2d
        return self.BSplineLst[0].StartPoint()
    
    @property
    def EndPoint(self): #gp_Pnt2d
        return self.BSplineLst[-1].EndPoint()
    
    @property
    def a_BSplineLst(self): #gp_Pnt2d
        return self.BSplineLst
    
    @property
    def b_BSplineLst(self): #gp_Pnt2d
        return self.Boundary_BSplineLst

    @property
    def IsClosed(self):
        return self.a_BSplineLst[0].StartPoint().IsEqual(self.a_BSplineLst[-1].EndPoint(),1e-5)
    
    def __str__(self): 
        #we can tell Python how to prepresent an object of our class (when using a print statement) for general purposes use  __repr__(self): 
        return  str('LayerID: \tStart[-]: \tEnd[-]: \tthickness[-]: \tOrientation[deg]: \tMatID \tName:\t\n' \
                    '%s, \t%s, \t%s, \t%s, \t\t%s, \t\t%s, \t%s, ' % (self.ID, self.S1, self.S2, self.thickness, self.Orientation, self.MatID, self.name))
     
        
    def __getstate__(self):
        """Return state values to be pickled."""

        self.Para_BSplineLst = ParaLst_from_BSplineLst(self.BSplineLst)
        self.Para_Boundary_BSplineLst = ParaLst_from_BSplineLst(self.Boundary_BSplineLst)
        
        return (self.ID, self.S1, self.S2, self.thickness, self.Orientation, self.MatID, self.Para_Boundary_BSplineLst, self.Para_BSplineLst)   
    
    
    def __setstate__(self, state):
        """Restore state from the unpickled state values."""
        self.ID, self.S1, self.S2, self.thickness, self.Orientation, self.MatID, self.Para_Boundary_BSplineLst, self.Para_BSplineLst = state
        self.Boundary_BSplineLst = BSplineLst_from_ParaLst(self.Para_Boundary_BSplineLst)
        self.BSplineLst = BSplineLst_from_ParaLst(self.Para_BSplineLst)
        
        
    def copy(self):
        BSplineLstCopy =  copy_BSplineLst(self.BSplineLst)
        namecopy = self.name + '_Copy'
        LayerCopy = Layer(self.ID,BSplineLstCopy,self.globalStart,self.globalEnd,self.thickness, self.Orientation, self.MatID, namecopy)
        return LayerCopy
  
    
    def get_length(self): #Determine and return Legth of Layer self
         self.length = get_BSplineLst_length(self.BSplineLst)
         return self.length
         
     
    def get_pnt2d(self,S,start,end): #Return, gp_Pnt2d of argument S of layer self  
        return get_BSplineLst_Pnt2d(self.BSplineLst,S,start,end)
            
    
    def build_wire(self): #Builds TopoDS_Wire from connecting BSplineSegments and returns it  
        self.wire = build_wire_from_BSplineLst(self.BSplineLst)   
        
    def trim(self,S1,S2,start, end): #Trims layer between S1 and S2
        return trim_BSplineLst(self.BSplineLst, S1, S2,  start, end)
        
    def trim_to_coords(self, start, end):
        self.BSplineLst = trim_BSplineLst(self.BSplineLst, self.globalStart, self.globalEnd,  start, end)
        return self.BSplineLst
            
#    def build_layer(self):
#        Obsolete
#        Trimmed_BSplineLst = trim_BSplineLst(self.Boundary_BSplineLst, self.S1, self.S2, 0, 1)
#        npArray = discretize_BSplineLst(Trimmed_BSplineLst, 1e-3) 
#        offlinepts = shp_parallel_offset(npArray,self.thickness,self.join_style)
#        OffsetBSplineLst = BSplineLst_from_dct(offlinepts)
#        OffsetBSplineLst = cutoff_layer(Trimmed_BSplineLst,OffsetBSplineLst,self.S1,self.S2,self.cutoff_style)
#        self.BSplineLst = OffsetBSplineLst
    
    def build_layer(self):
        npArray = discretize_BSplineLst(self.Boundary_BSplineLst, 1e-3) 
        offlinepts = shp_parallel_offset(npArray,self.thickness,self.join_style)
        OffsetBSplineLst = BSplineLst_from_dct(offlinepts)
        OffsetBSplineLst = cutoff_layer(self.Boundary_BSplineLst,OffsetBSplineLst,self.S1,self.S2,self.cutoff_style)
        self.BSplineLst = OffsetBSplineLst
         
        
    def determine_a_nodes(self,LayerLst,global_minLen):
        nLayers = len(LayerLst)
        new_a_nodes=[]
        for iv_counter,iv in enumerate(self.inverse_ivLst):
            if int(iv[2])==nLayers: 
                #print iv, "equidistand nodes on BsplineLst of LayerLst entry",
                eq_nodes = []
                BSplineLst = self.a_BSplineLst             
                iv_BSplineLst = trim_BSplineLst(BSplineLst,iv[0],iv[1],self.S1,self.S2  )
                if iv_counter==0 and len(self.inverse_ivLst)>1: #first but not last
                    IncStart=True
                    IncEnd=False
                
                elif iv_counter==0 and len(self.inverse_ivLst)==1: #first and last
                    IncStart=True
                    IncEnd=True

                elif iv_counter==len(self.inverse_ivLst)-1 and len(self.inverse_ivLst)>1: #last but not first
                    if  iv_counter==len(self.inverse_ivLst)-1 and iv[1]==1 and self.inverse_ivLst[0][0]==0:
                        IncStart=False
                        IncEnd=False
                    else:
                        IncStart=False
                        IncEnd=True
                                        
                else:
                    IncStart=False
                    IncEnd=False
                
                eq_nodes = equidistant_nodes_on_BSplineLst(iv_BSplineLst, True, IncStart, IncEnd, minLen = global_minLen, LayerID = self.ID[0])
                new_a_nodes.extend(eq_nodes)

                
            else:
                #print iv, "use nodes b_nodes of layer", int(iv[2])
                tmp_layer = LayerLst[int(iv[2])]
                iv_BSplineLst = trim_BSplineLst(tmp_layer.b_BSplineLst,iv[0],iv[1],tmp_layer.S1,tmp_layer.S2 )
                isClosed = iv_BSplineLst[0].StartPoint().IsEqual(iv_BSplineLst[-1].EndPoint(),1e-5)
               
                if isClosed:
                    tmp_nodes = tmp_layer.b_nodes
                else: 
                    tmp_nodes = [tmp_layer.a_nodes[0]]+tmp_layer.b_nodes+[tmp_layer.a_nodes[-1]]
                
                disco_nodes = grab_nodes_on_BSplineLst(tmp_nodes,iv_BSplineLst)
                new_a_nodes.extend(disco_nodes)
        
        self.a_nodes = new_a_nodes

    
    def show(self): #display the layer with pythonocc viewer module
        """
        TBD: renders the topological entity in the viewer: 
        """
        if not self.display_set:
            display.Display()(self, *args, **kwargs)
        else:
            self.disp.DisplayShape(*args, **kwargs)


    
   
#execute the following code if this file is executed as __main__   
if __name__ == '__main__':   
     L1 = Layer()
     pass