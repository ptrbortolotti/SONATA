from BSplineLst_utils import get_BSplineLst_length, get_BSplineLst_Pnt2d, trim_BSplineLst
from wire_utils import build_wire_from_BSplineLst

class Layer(object):
    ''' 
    The layer object is constructed from multiple BSplineCurveSegments. It is the basis for all future operations. 
    The object can be constructed from either a discrete formulation of point tables or from an existing TopoDS_Wire.
    ''' 
   
    def __init__(self, ID, BSplineLst, globalStart, globalEnd, thickness, Orientation = 0, MatID = 1, name='DEFAULT'): #gets called whenever we create a new instance of a class
        self.BSplineLst = BSplineLst        #List of Geom2d_BSplineCurve, Geom_BSplineCurve
        self.ID = ID                        # Segment 1, Layer 29,
        self.name = name
        #self.wire = []                      #Make Wire from BSplineSegments
        #self.First = []			    #Starting Point in S coordinates
        #self.Last = []				    #End Point in S Coordinates
        self.display_set = False
        self.globalStart = globalStart	    #Starting Point in S coordinates
        self.globalEnd = globalEnd		    #End Point in S coordinates
        self.thickness = thickness   	    #in
        self.Orientation = Orientation
        self.MatID = MatID
        #self.DCTArrayIN = []		#Discrete point distribution over the layer from polygon offset
        #self.DCTArrayOUT = []		#Discrete point distribution for the next layer from polygon offset	
    
    def __str__(self): 
        #we can tell Python how to prepresent an object of our class (when using a print statement) for general purposes use  __repr__(self): 
        return  str('LayerID: \tStart[-]: \tEnd[-]: \tthickness[mm]: \tOrientation[deg]: \tMatID \tName:\t\n' \
                    '%s, \t%s, \t%s, \t%s, \t\t%s, \t\t%s, \t%s, ' % (self.ID, self.globalStart, self.globalEnd, self.thickness, self.Orientation, self.MatID, self.name))
                   
    def get_length(self): #Determine and return Legth of Layer self
         self.length = get_BSplineLst_length(self.BSplineLst)
         return self.length
         
    def get_pnt2d(self,S): #Return, gp_Pnt of argument S of layer self  
        return get_BSplineLst_Pnt2d(self.BSplineLst,S)
            
    def build_wire(self): #Builds TopoDS_Wire from connecting BSplineSegments and returns it  
        self.wire = build_wire_from_BSplineLst(self.BSplineLst)   
        
    def trim(self,S1,S2): #Trims layer between S1 and S2
        return trim_BSplineLst(self.BSplineLst, S1, S2)
    
    def get_first():
        pass
    
    def get_last():
        pass
    
 
    def check_layer():
        pass
    
    def layer_from_wire(self,TopoDS_wire):
        pass
    
    def layer_from_dctData(self,dctData):
        pass    
    
    
    def offset(self):
        #return the new layer
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