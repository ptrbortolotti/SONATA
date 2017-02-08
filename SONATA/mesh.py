import pickle
import matplotlib as plt

from OCC.Display.SimpleGui import init_display
from OCC.Quantity import Quantity_Color


#LOAD .pkl data with SegmentLst
filename = 'sec_config.pkl'
with open(filename, 'rb') as handle:
    SegmentLst = pickle.load(handle)
    

#Build wires for each layer and segment
for seg in SegmentLst:
    seg.build_wire()
    for layer in seg.LayerLst:
        layer.build_wire()

#==============================================================================


# HERE WILL BE THE MESHING ALGORITHM!






#====================DISPLAY===================================================
display, start_display, add_menu, add_function_to_menu = init_display('wx')
display.Context.SetDeviationAngle(0.000001)       # 0.001 default. Be careful to scale it to the problem.
display.Context.SetDeviationCoefficient(0.000001) # 0.001 default. Be careful to scale it to the problem. 
display.set_bg_gradient_color(20,6,111,200,200,200) 
    
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

display.View_Top()
display.FitAll()
start_display()
