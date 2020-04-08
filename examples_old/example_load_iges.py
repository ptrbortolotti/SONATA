
from OCC.Display.SimpleGui import init_display
from OCC.Extend.DataExchange import read_iges_file

shapes = read_iges_file('wt.iges')

display, start_display, add_menu, add_function_to_menu = init_display()
display.DisplayShape(shapes, update=True)
