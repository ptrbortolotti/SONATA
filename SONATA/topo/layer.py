import matplotlib.pyplot as plt

from SONATA.topo.BSplineLst_utils import get_BSplineLst_length, get_BSplineLst_Pnt2d, trim_BSplineLst, copy_BSplineLst, BSplineLst_from_dct, discretize_BSplineLst 
from SONATA.topo.wire_utils import build_wire_from_BSplineLst
from SONATA.topo.cutoff import cutoff_layer
from SONATA.topo.offset import shp_parallel_offset
from SONATA.topo.para_Geom2d_BsplineCurve import ParaLst_from_BSplineLst, BSplineLst_from_ParaLst


class Layer(object):
    ''' 
    The layer object is constructed from multiple BSplineCurveSegments. It is the basis for all future operations. 
    The object can be constructed from either a discrete formulation of point tables or from an existing TopoDS_Wire.
    ''' 
   
    def __init__(self, ID, Boundary_BSplineLst, globalStart, globalEnd, thickness, Orientation = 0, MatID = 1, **kwargs):
        self.ID = ["%04d" % ID] 	                            #First single Digit: SegmentNb, Last 3 Digits: LayerNb; e.g.: 1029, Segment1, Layer29
        self.Boundary_BSplineLst = Boundary_BSplineLst        #List of Geom2d_BSplineCurve, Geom_BSplineCurve
        self.S1 = globalStart	                      #Starting Point in S coordinates
        self.S2 = globalEnd		                      #End Point in S coordinates
        self.thickness = thickness   	                      #in
        self.Orientation = Orientation
        self.MatID = MatID
        self.cells = []
        
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
        #self.DCTArrayIN = []		#Discrete point distribution over the layer from polygon offset
        #self.DCTArrayOUT = []		#Discrete point distribution for the next layer from polygon offset	
    
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
         
    def get_pnt2d(self,S,start,end): #Return, gp_Pnt of argument S of layer self  
        return get_BSplineLst_Pnt2d(self.BSplineLst,S,start,end)
            
    def build_wire(self): #Builds TopoDS_Wire from connecting BSplineSegments and returns it  
        self.wire = build_wire_from_BSplineLst(self.BSplineLst)   
        
    def trim(self,S1,S2,start, end): #Trims layer between S1 and S2
        return trim_BSplineLst(self.BSplineLst, S1, S2,  start, end)
        
    def trim_to_coords(self):
        self.BSplineLst = trim_BSplineLst(self.BSplineLst, self.globalStart, self.globalEnd,  start, end)
        return self.BSplineLst
            
    def build_layer(self):
        Trimmed_BSplineLst = trim_BSplineLst(self.Boundary_BSplineLst, self.S1, self.S2, 0, 1)
        npArray = discretize_BSplineLst(Trimmed_BSplineLst, 1e-3) 
        offlinepts = shp_parallel_offset(npArray,self.thickness,self.join_style)
        OffsetBSplineLst = BSplineLst_from_dct(offlinepts)
        OffsetBSplineLst = cutoff_layer(Trimmed_BSplineLst,OffsetBSplineLst,self.S1,self.S2,self.cutoff_style)
        self.BSplineLst = OffsetBSplineLst
        
    
    def get_last():
        pass
    
 
    def check_layer():
        pass
    
    def layer_from_wire(self,TopoDS_wire):
        pass
    
    def layer_from_dctData(self,dctData):
        pass    

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