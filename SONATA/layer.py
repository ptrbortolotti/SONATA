from BSplineLst_utils import get_BSplineLst_length, get_BSplineLst_Pnt2d


class Layer(object):
    ''' 
    The layer object is constructed from multiple BSplineCurveSegments. It is the basis for all future operations. 
    The object can be constructed from either a discrete formulation of point tables or from an existing TopoDS_Wire.
    ''' 
   
    def __init__(self, BSplineLst, globalStart, globalEnd, thickness, Orientation = 0, MatID = 1): #gets called whenever we create a new instance of a class
        self.BSplineLst = BSplineLst        #List of Geom2d_BSplineCurve, Geom_BSplineCurve
        #self.Wire = []                      #Make Wire from BSplineSegments
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
		
    #def __repr__(self): #we can tell Python how to prepresent an object of our calss (when using a print statement)
        #pass
            
    def get_length(self): #Determine and return Legth of Layer self
         self.length = get_BSplineLst_length(self.BSplineLst)
         return self.length
         
    def get_pnt2d(self,S): #Return, gp_Pnt of argument S of layer self  
        return get_BSplineLst_Pnt2d(self.BSplineLst,S)
            
    def build_wire(self): #Builds TopoDS_Wire from connecting BSplineSegments and returns it  
        pass   
        
    def trim(self,S1,S2): #Trims layer between S1 and S2
        pass
    
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
    
    def set_to_origin(self):
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