class Layer(object):
    ''' 
    The layer object is constructed from multiple BSplineCurveSegments. It is the basis for all future operations. 
    The object can be constructed from either a discrete formulation of point tables or from an existing TopoDS_Wire.
    ''' 
   
    def __init__(self): #gets called whenever we create a new instance of a class
        self.BSplineSegments = []       #List of Geom2d_BSplineCurve, Geom_BSplineCurve
        self.Wire = []                  #Make Wire from BSplineSegments
        self.Length = []
        self.First = []
        self.Last = []
        self.display_set = False
        self.globalStart = []
        self.globalEnd = []
        self.thickness = []    #in
        self.Orientation =[]
        self.MaterialID = []
 
    #def __repr__(self): #we can tell Python how to prepresent an object of our calss (when using a print statement)
        #pass
            
    def length(self): #Determine and return Legth of Layer self
        pass
 
    def pnt(self,S): #Return, gp_Pnt of argument S of layer self  
        pass
    
    def pnt2d(self,S): #Return, gp_Pnt2d of argument S of layer self  
        pass
    
    def build_wire(self): #Builds TopoDS_Wire from connecting BSplineSegments and returns it  
        pass   
        
    def trim(self,S1,S2): #Trims layer between S1 and S2
        pass
    
    def first():
        pass
    
    def last():
        pass
    
 
    def check():
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