import os
import numpy as np
import pickle
import matplotlib as plt

from OCC.AIS import AIS_Shape
from OCC.BRep import BRep_Builder, BRep_Tool
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeFace
from OCC.BRepAdaptor import BRepAdaptor_CompCurve
from OCC.gp import gp_Pnt2d, gp_Pnt, gp_Vec2d,gp_Lin2d, gp_Dir2d
from OCC.GCPnts import GCPnts_QuasiUniformAbscissa
from OCC.Geom2d import Geom2d_Line
from OCC.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.Geom2dAPI import Geom2dAPI_ProjectPointOnCurve
from OCC.Display.SimpleGui import init_display
from OCC.Quantity import Quantity_Color
from OCC.Graphic3d import Graphic3d_EF_PDF, Graphic3d_EF_SVG, Graphic3d_EF_TEX, Graphic3d_EF_PostScript, Graphic3d_EF_EnhPostScript
from OCC.BRepLib import breplib_BuildCurves3d

from OCC.TopoDS import TopoDS_Compound, topods_Face, topods_Edge


from BSplineLst_utils import get_BSplineLst_length, get_BSpline_length, trim_BSplineLst
from wire_utils import build_wire_from_BSplineLst,get_wire_length
from utils import Pnt2dLst_to_npArray, unique_rows
from display_mesh import plot_mesh
#=======================DISPLAY FUCTIONS....===================================
def display_SONATA_SegmentLst(SegmentLst):
    # transfer shapes and display them in the viewer
    display.DisplayShape(SegmentLst[0].wire, color="BLACK")
    display.DisplayShape(SegmentLst[0].BSplineLst[0].StartPoint())
    
    for i,seg in enumerate(SegmentLst):
        display.DisplayShape(seg.wire,color="BLACK")
        k = 0
        for j,layer in enumerate(seg.LayerLst):
            [R,G,B,T] =  plt.cm.jet(k*50)
            
            if i==0:
                display.DisplayColoredShape(layer.wire, Quantity_Color(R, G, B, 0),update=True)
                
            elif i==1:
                display.DisplayColoredShape(layer.wire, Quantity_Color(R, G, B, 0),update=True)
    
            else:
                display.DisplayColoredShape(layer.wire, Quantity_Color(R, G, B, 0),update=True)
    
            k = k+1;
            if k>5:
                k = 0
    return None


def display_custome_shape(shape,linewidth,transparency,RGB):
    s = shape
    ais_shp = AIS_Shape(s)
    ais_shp.SetWidth(linewidth)
    ais_shp.SetTransparency(transparency)
    ais_shp.SetColor(Quantity_Color(RGB[0], RGB[1], RGB[2], 0))
    ais_context = display.GetContext().GetObject()
    ais_context.Display(ais_shp.GetHandle())
    return None



def export_to_PDF(event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('capture_pdf%s.pdf' % i):
        i += 1
    f.Export('capture_pdf%s.pdf' % i, Graphic3d_EF_PDF)
    print "EXPORT: \t Screencapture exported to capture_pdf%s.pdf" % i
    display.set_bg_gradient_color(20,6,111,200,200,200)
    
def export_to_SVG(event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('capture_svg%s.svg' % i):
        i += 1
    f.Export('capture_svg_%s.svg' % i, Graphic3d_EF_SVG)
    print "EXPORT: \t Screencapture exported to capture_svg%s.svg" % i
    display.set_bg_gradient_color(20,6,111,200,200,200)
    
def export_to_PS(event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('capture_ps%s.ps' % i):
        i += 1
    f.Export('capture_ps%s.ps' % i, Graphic3d_EF_PostScript)
    print "EXPORT: \t Screencapture exported to capture_ps%s.ps" % i
    display.set_bg_gradient_color(20,6,111,200,200,200)

def export_to_EnhPS(event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('capture_Enh_ps%s.ps' % i):
        i += 1
    f.Export('capture_Enh_ps%s.ps' % i, Graphic3d_EF_EnhPostScript)
    print "EXPORT: \t Screencapture exported to capture_Enh_ps%s.ps" % i
    display.set_bg_gradient_color(20,6,111,200,200,200)
    
def export_to_TEX(event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('capture_tex%s.tex' % i):
        i += 1
    f.Export('capture_tex%s.tex' % i, Graphic3d_EF_TEX)
    print "EXPORT: \t Screencapture exported to capture_tex%s.tex" % i
    display.set_bg_gradient_color(20,6,111,200,200,200)
    
def export_to_BMP(event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('capture_bmp%s.bmp' % i):
        i += 1
    display.View.Dump('capture_bmp%s.bmp' % i)
    print "EXPORT: \t Screencapture exported to capture_bmp%s.bmp" % i
    display.set_bg_gradient_color(20,6,111,200,200,200)

def export_to_PNG(event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('capture_png%s.png' % i):
        i += 1
    display.View.Dump('capture_png%s.png' % i)
    print "EXPORT: \t Screencapture exported to capture_png%s.bmp" % i
    display.set_bg_gradient_color(20,6,111,200,200,200)

def export_to_JPEG(event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('capture_jpeg%s.jpeg' % i):
        i += 1
    display.View.Dump('capture_jpeg%s.jpeg' % i)
    print "EXPORT: \t Screencapture exported to capture_jpeg%s.jpeg" % i
    display.set_bg_gradient_color(20,6,111,200,200,200)

def export_to_TIFF(event=None):
    display.set_bg_gradient_color(255,255,255,255,255,255)
    i = 0
    while os.path.exists('capture_tiff%s.tiff' % i):
        i += 1
    display.View.Dump('capture_tiff%s.tiff' % i)
    print "EXPORT: \t Screencapture exported to capture_tiff%s.tiff" % i
    display.set_bg_gradient_color(20,6,111,200,200,200)



#====================INIT DISPLAY:==========================================================
display, start_display, add_menu, add_function_to_menu = init_display('wx')
display.Context.SetDeviationAngle(0.000001)       # 0.001 default. Be careful to scale it to the problem.
display.Context.SetDeviationCoefficient(0.000001) # 0.001 default. Be careful to scale it to the problem. 
display.set_bg_gradient_color(20,6,111,200,200,200) 
    

#====================NOTES:==========================================================
'''
1) Distribute Equidistant Points on Layer_i
2) Find Cell where the Max.Error between Cell and the Layer is larger than allowed -> Split Cell at max deviation
3) Determine Area of Cells and Cell Angles -> Improve Cells 





'''
#===================LOAD CROSSSECTION==========================================
#LOAD .pkl data with SegmentLst
filename = 'sec_config.pkl'
with open(filename, 'rb') as handle:
    SegmentLst = pickle.load(handle)
    

#Build wires for each layer and segment
for seg in SegmentLst:
    seg.build_wire()
    for layer in seg.LayerLst:
        layer.build_wire()



class Node(object):
    class_counter= 1
    def __init__(self, Pnt2d, parameters=['0',0,0]):
        self.id= self.__class__.class_counter
        self.__class__.class_counter += 1
        
        self.Pnt2d = Pnt2d  #gp_Pnt2d
        self.Pnt = gp_Pnt(self.Pnt2d.X(),self.Pnt2d.Y(),0)  #gp_Pnt
        self.coordinates = [self.Pnt2d.X(),self.Pnt2d.Y()]  #[x,y]
        self.parameters = parameters    #[LayerID, idx, U]
        self.corner = False
        self.face_pointer = []
        
    def __repr__(self): 
        #we can tell Python how to prepresent an object of our class (when using a print statement) for general purposes use  __repr__(self): 
        return  str('Node: %s @ [%.3f,%.3f]' % (self.id, self.coordinates[0],self.coordinates[1]))
            
  

def Pnt2dLst_to_NodeLst(Pnt2dLst):
    NodeLst = []
    for Pnt2d in Pnt2dLst:
        NodeLst.append(Node(Pnt2d))
    return NodeLst

 
            
class Cell(object):
    class_counter= 1
    def __init__(self,nodeLst,Orientation,MatID):                  #int
        self.id = self.__class__.class_counter
        self.__class__.class_counter += 1
        
        self.nodes = nodeLst              #[node,node,node,nodes]      !!!counterclockwise direction!!!
        self.face  = []                     #[rear,top,front,bottom]        !!!counterclockwise direction!!!       
        self.wire  = self.build_wire()  #TopoDs_wire
        self.neighbours = []                #-[Cell_ID,CELL_ID... ]
        self.theta_1 = self.calc_theta_1()  #Ply coordinate system is formed by rotating the global coordinate system in the right-hand sense about the amount 0<Theta_1<260.
        self.theta_3 = Orientation          #The Ply coordiate system is rotated about y3 in the right hand sense by the amount -90<Theta_3<90 to for the material system. 
        self.MatID  = MatID                 #material id, int
        
    def __repr__(self): 
        #we can tell Python how to prepresent an object of our class (when using a print statement) for general purposes use  __repr__(self): 
        STR = '\n'
        STR += str('Cell %s w. nodes:\t' % (self.id))
        for n in self.nodes:
            STR += str('%i, ' % (n.id))
            
        return  STR
        
    def calc_theta_1(self):  
        if len(self.nodes) == 3:
            theta_1 = 69
        elif len(self.nodes) == 4:
            #points[0].coordinates+points[1].coordinates
            theta_1 = 69
        
        return theta_1
            

    def build_wire(self):
        WireBuilder = BRepBuilderAPI_MakeWire()
        for i in range(0,len(self.nodes)-1):
            me = BRepBuilderAPI_MakeEdge(self.nodes[i].Pnt, self.nodes[i+1].Pnt)
            if me.IsDone():
                WireBuilder.Add(me.Edge())
        
        me = BRepBuilderAPI_MakeEdge(self.nodes[-1].Pnt, self.nodes[0].Pnt)
        if me.IsDone():
            WireBuilder.Add(me.Edge())         
        
        return WireBuilder.Wire()

def equidistant_nodes_on_BSplineLst(BSplineLst, IC=False, **kwargs): #include corners
    ''' minLen = the minimum distance between points. Should be determined by the segment0.BSplineLstLenght/Resolution
        IC (Include Corners): True or False, Note: IC=True the the Points are not equidistant placed on BSplineLst!!!
        
        (BSplineLst, IC=False, either NbPoints or MinLen)
        (BSplineLst, IC=True,  MinLen)
        
    ''' 
    if IC==True:
        
        #KWARGS:
        if kwargs.get('minLen') !=  None:
            minLen = kwargs.get('minLen')

        if kwargs.get('LayerID') !=  None:
            LayerID = kwargs.get('LayerID')
        
        nodes = []
        for idx,item in enumerate(BSplineLst):
            Adaptor = Geom2dAdaptor_Curve(item.GetHandle())
            length = get_BSpline_length(item)
            NbPoints = int(length//minLen)+2  
            discretization = GCPnts_QuasiUniformAbscissa(Adaptor,NbPoints)
            
            for j in range(1, NbPoints):
                    para = discretization.Parameter(j)
                    Pnt = gp_Pnt2d()
                    item.D0(para,Pnt)
                    node = Node(Pnt,[LayerID,idx,para])
                    nodes.append(node)
            
        #add last point 
        para = discretization.Parameter(j+1)
        Pnt = gp_Pnt2d()
        item.D0(para,Pnt)
        node = Node(Pnt, [LayerID,idx,para])
        nodes.append(node)       
                    

    else:
        wire = build_wire_from_BSplineLst(BSplineLst)
        wire_length = get_wire_length(wire)
        
        #KWARGS:
        if  kwargs.get('NbPoints') !=  None:
            NbPoints = kwargs.get('NbPoints')
            
        elif kwargs.get('minLen') !=  None and kwargs.get('NbPoints') == None:
            NbPoints = int(wire_length//kwargs.get('minLen'))+2  
                          
        AdaptorComp = BRepAdaptor_CompCurve(wire, True)
        discretization = GCPnts_QuasiUniformAbscissa(AdaptorComp,NbPoints)
        
        nodes = []
        for j in range(1, NbPoints+1):
            para = discretization.Parameter(j)
            P = AdaptorComp.Value(para)
            node = Node(gp_Pnt2d(P.X(),P.Y()))
            nodes.append(nodes)
       
    return nodes


def corners_of_BSplineLst(BSplineLst):
    corners = [] 
    for item in BSplineLst:
        corners.append(item.EndPoint())
        
    corners.pop(-1)
    return corners #gp_Pnt2d Lst
    

def project_Pnt2d_on_BSplineLst(a_BSplineLst,a_nodes,b_BSplineLst,layer_thickness,minLen, tol=2e-3):
    LayerID = 'T_' + a_nodes[0].parameters[0]
    b_nodes = []
    cellLst = []
    distance = (1+tol)*layer_thickness
               
    #==================PROJECT POINTS ON LOWER BOUNDARY =======================            
    for i,node in enumerate(a_nodes[1:-1], start=1):
        Pnt2d = node.Pnt2d
        pPnts = []
        pPara = []
        pIdx = []
        
        for idx,item in enumerate(b_BSplineLst):
            projection = Geom2dAPI_ProjectPointOnCurve(Pnt2d,item.GetHandle())
            
            for j in range(1,projection.NbPoints()+1):
                if projection.Distance(j)<=distance:
                    pPnts.append(projection.Point(j))
                    pPara.append(projection.Parameter(j))
                    pIdx.append(idx)
                else: None   
            
        #==================DETECT CORNERS====================================== 
        if len(pPnts) == 1: 
            b_nodes.append(Node(pPnts[0],[LayerID,pIdx[0],pPara[0]]))    


        elif len(pPnts) == 2:
            #display.DisplayShape(Pnt2d,color='GREEN')
            v1 = gp_Vec2d(Pnt2d,pPnts[0])
            v2 = gp_Vec2d(Pnt2d,pPnts[1])
            angle = (180-v1.Angle(v2)*180/np.pi)
            crit_angle = 100

            if angle < crit_angle:
                node.corner = True
                b_nodes.append(Node(pPnts[0],[LayerID,pIdx[0],pPara[0]]))
                b_nodes.append(Node(b_BSplineLst[pIdx[0]].EndPoint(),[LayerID,pIdx[0],b_BSplineLst[pIdx[0]].LastParameter()]))
                b_nodes.append(Node(pPnts[1],[LayerID,pIdx[1],pPara[1]]))    
                #display.DisplayShape(pPnts[0],color='YELLOW')        
                #display.DisplayShape(b_BSplineLst[pIdx[0]].EndPoint(),color='YELLOW')
                #display.DisplayShape(pPnts[1],color='YELLOW')
             
            else:
                b_nodes.append(Node(b_BSplineLst[pIdx[0]].EndPoint(),[LayerID,pIdx[0],b_BSplineLst[pIdx[0]].LastParameter()]))
                #display.DisplayShape(b_BSplineLst[pIdx[0]].EndPoint(),color='YELLOW')
                
        else: 
            print 'Projection Error, number of projection points: ', len(pPnts)
            print projection.NbPoints()
            #display.DisplayShape(Pnt2d,color='RED')  
    
    
    #==============REVERSED PROJECTION=========================================
    leftover_exterior = [] 
    leftover_exterior_para = []
    for i,item in enumerate(b_BSplineLst[:-1]):
        spline1 = item
        spline2 = b_BSplineLst[i+1]
        u1,p1,v1 = spline1.LastParameter(),gp_Pnt2d(),gp_Vec2d()
        u2,p2,v2  = spline2.FirstParameter(),gp_Pnt2d(),gp_Vec2d()
        spline1.D1(u1,p1,v1)
        spline2.D1(u2,p2,v2)
        
        Angle = abs(v1.Angle(v2))*180/np.pi       
        if Angle>0.5:
            leftover_exterior.append(item.EndPoint())
            leftover_exterior_para.append([LayerID,i,u1])  
    
    #find exterior corner Points that are not part of b_nodes
    to_delete = []
    LinearTolerance = 1e-3
    for idx,corn in enumerate(leftover_exterior):
        for node in b_nodes:
            if node.Pnt2d.IsEqual(corn, LinearTolerance):
                to_delete.append(idx)
                break                
    
    for offset,idx in enumerate(to_delete):
        idx -= offset
        del leftover_exterior[idx]
        del leftover_exterior_para[idx]
    
    print len(leftover_exterior)
    #do the reversed projection! -> the original Pnt2dLst must be modified and be returned as well!
    leftover_interior = []
    leftover_interior_para = []
    for Pnt2d in leftover_exterior:
        pPnts = []
        pIdx = []
        pPara = []
        for idx,item in enumerate(a_BSplineLst):
            projection = Geom2dAPI_ProjectPointOnCurve(Pnt2d,item.GetHandle())
            for i in range(1,projection.NbPoints()+1):
                if projection.Distance(i)<=distance:
                    pPnts.append(projection.Point(i))
                    pPara.append(projection.Parameter(i))
                    pIdx.append(idx)
              
        if len(pPnts) == 1:
            leftover_interior.append(pPnts[0])
            leftover_interior_para.append([a_nodes[0].parameters[0],pIdx[0],pPara[0]]) 


    
    leftover_exterior_nodes = []
    leftover_interior_nodes = []
    for i,p in enumerate(leftover_exterior):
        leftover_exterior_nodes.append(Node(leftover_exterior[i],leftover_exterior_para[i]))
        leftover_interior_nodes.append(Node(leftover_interior[i],leftover_interior_para[i]))
    

    #======INSERT LEFTOVER NODES ==============================================
    newlist = a_nodes + leftover_interior_nodes       
    a_nodes =  sorted(newlist, key=lambda Node: (Node.parameters[1],Node.parameters[2]))
       
    newlist = b_nodes + leftover_exterior_nodes       
    b_nodes =  sorted(newlist, key=lambda Node: (Node.parameters[1],Node.parameters[2]))
    
#    for n in b_nodes:
#        print n,n.parameters

    
    #==============CREATE CELLS PROJECTION=========================================
    
    #Last Cell as Triangle:
    b = 0   #b_nodes idx
    for a,node in enumerate(a_nodes[1:-1], start=1):
        if a == 1: #Start Triangle
            cellLst.append(Cell([a_nodes[a],a_nodes[a-1],b_nodes[b]],69,69))
        
        elif a==len(a_nodes)-2: #End Triangle
            cellLst.append(Cell([a_nodes[a-1],b_nodes[b-1],b_nodes[b],a_nodes[a]],69,69))
            cellLst.append(Cell([a_nodes[a],b_nodes[b],a_nodes[a+1]],69,69))
        
        else: #Regular Cell Creation
            if a_nodes[a].corner == True:
                cellLst.append(Cell([a_nodes[a-1],b_nodes[b-1],b_nodes[b],a_nodes[a]],69,69))
                b += 2
                cellLst.append(Cell([a_nodes[a],b_nodes[b-2],b_nodes[b-1],b_nodes[b]],69,69))

            else:   
                cellLst.append(Cell([a_nodes[a-1],b_nodes[b-1],b_nodes[b],a_nodes[a]],69,69))
        
        b += 1
        
        
        
    #==============DISPLAY POINTS AND CELLS=========================================    

    #print cellLst
    for cell in cellLst:
        display.DisplayShape(cell.wire,color='BLACK')
        
    for n in b_nodes:
        display.DisplayShape(n.Pnt2d,color='GREEN')
    
    for n in a_nodes:
        display.DisplayShape(n.Pnt2d,color='ORANGE')  
        
        
    return a_nodes, b_nodes, cellLst




#==============================================================================
Projection = SegmentLst[-1].Projection
Resolution = 3000 # Nb of Points on Segment0
length = get_BSplineLst_length(SegmentLst[0].BSplineLst)
global_minLen = round(length/Resolution,5)
layer = SegmentLst[-1].LayerLst[-1]
BSplineLst = layer.BSplineLst
Trimmed_BSplineLst = trim_BSplineLst(layer.Boundary_BSplineLst, layer.S1, layer.S2, 0, 1)


#Point Projection on BSPLINELST 
a_BSplineLst = BSplineLst
a_nodes = equidistant_nodes_on_BSplineLst(BSplineLst, True, minLen = global_minLen, LayerID = layer.ID[0])

b_BSplineLst = Trimmed_BSplineLst
a_nodes, b_nodes, cells = project_Pnt2d_on_BSplineLst(a_BSplineLst,a_nodes,b_BSplineLst,layer.thickness,global_minLen)



















 
##CREATE mesh
#nodes = np.vstack([a,b])
#s_idx = len(a)
#len_b = len(b)
#elements = []
#for i in range(0, len(b)-1):
#    elements.append([i+3,i+2,s_idx+(i+1),s_idx+(i+2)])
#data = np.zeros(len(elements))
#
##plot_mesh(nodes,elements,data,False,False)
#
#
#
##DisplayMesh(nodex,elements)
#builder = BRep_Builder()
#comp = TopoDS_Compound()
#WireBuilder = BRepBuilderAPI_MakeWire()
#builder.MakeCompound(comp)
#WireLst = []
#
#for ele in elements:
#    WireBuilder = BRepBuilderAPI_MakeWire()
#    for i in range(0,len(ele)-1):
#        node1 = ele[i]-1
#        node2 = ele[i+1]-1
#        Pnt1 = gp_Pnt(nodes[node1][0],nodes[node1][1],0)
#        Pnt2 = gp_Pnt(nodes[node2][0],nodes[node2][1],0)
#        me = BRepBuilderAPI_MakeEdge(Pnt1, Pnt2)
#        if me.IsDone():
#            WireBuilder.Add(me.Edge())
#    
#    node1 = ele[-1]-1
#    node2 = ele[0]-1         
#    Pnt1 = gp_Pnt(nodes[node1][0],nodes[node1][1],0)
#    Pnt2 = gp_Pnt(nodes[node2][0],nodes[node2][1],0)
#    me = BRepBuilderAPI_MakeEdge(Pnt1, Pnt2)
#    if me.IsDone():
#        WireBuilder.Add(me.Edge())         
#        
#    WireLst.append(WireBuilder.Wire())
#
#FaceLst = [] 
#for wire in WireLst:
#    Face = BRepBuilderAPI_MakeFace(wire)
#    FaceLst.append(Face.Face())


#OPTIMIZE MESH!!!!!


#NEXT LAYER
# - if: is Pnt2d on BSplineLst???
#       then....




#====================DISPLAY===================================================

#for p in a_Pnt2dLst:
#    display.DisplayShape(p,color='BLACK')
#
#for p in b_Pnt2dLst:
#    display.DisplayShape(p,color='GREEN')


    
display_SONATA_SegmentLst(SegmentLst)


#for p in c_Pnt2dLst:
#    display.DisplayShape(p,color='CYAN')
    
#for p in d_Pnt2dLst:
#    display.DisplayShape(p,color='ORANGE')

#for w in WireLst:
#    display.DisplayShape(w)


#==============================================================================
'''CREATE AXIS SYSTEM for Visualization'''
O  = gp_Pnt(0., 0., 0.)
p1 = gp_Pnt(10.0,0.,0.)
p2 = gp_Pnt(0.,10.0,0.)
p3 = gp_Pnt(0.,0.,10.0)

h1 = BRepBuilderAPI_MakeEdge(O,p1).Shape()
h2 = BRepBuilderAPI_MakeEdge(O,p2).Shape()
h3 = BRepBuilderAPI_MakeEdge(O,p3).Shape()

display.DisplayShape(O,color='BLACK')
display.DisplayShape(h1,color='RED')
display.DisplayShape(h2,color='GREEN')
display.DisplayShape(h3,color='BLUE')
    
    
f = display.View.View().GetObject()

display.set_bg_gradient_color(20,6,111,200,200,200)
add_menu('screencapture')
add_function_to_menu('screencapture', export_to_PDF)
add_function_to_menu('screencapture', export_to_SVG)
add_function_to_menu('screencapture', export_to_PS)
add_function_to_menu('screencapture', export_to_EnhPS)
add_function_to_menu('screencapture', export_to_TEX)
add_function_to_menu('screencapture', export_to_BMP)
add_function_to_menu('screencapture', export_to_PNG)
add_function_to_menu('screencapture', export_to_JPEG)
add_function_to_menu('screencapture', export_to_TIFF)

display.View_Top()
display.FitAll()
start_display()

