#!/usr/bin/env python
"""
Launch HelloWorld wxFormBuilder App
"""
 
import numpy as np
import wx
import HelloWorld
 
import matplotlib
matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
        FigureCanvasWxAgg as FigCanvas, \
        NavigationToolbar2WxAgg as NavigationToolbar
 
class MyFrameSub(HelloWorld.MyFrame1):
    """
   Subclass the form created by wxFormBuilder.  We use the generated output
   from wxFormBuilder verbatim and add all of the functionality here.
   """
    # Overload methods in GUI
    # Perhaps event handlers.
    def __init__(self, parent):
        # If we overload __init__ we still need to make sure the parent
        # is initialized.
        # This also sets up GUI elements with in self namespace
        HelloWorld.MyFrame1.__init__(self, parent)
 
        # Events can also be bound in the GUI builder arguably is easier
        # because you don't have to remember the WX events.
        # The gui builder will bind the events and generate a virtual method
        # that we can override in the derived class.
 
        #self.m_Calculate.Bind(wx.EVT_BUTTON, self.OnCalculate)
        self.CreatePlot()
 
        return
 
    def CreatePlot(self):
        self.figure = Figure()#figsize=(6, 4), dpi=100)
        self.axes = self.figure.add_subplot(111)
        x = np.arange(0, 6, .01)
        y = np.sin(x**2)*np.exp(-x)
        self.axes.plot(x, y)
        # Add it to the panel created in wxFormBuilder
        self.canvas = FigCanvas(self.m_panel1, wx.ID_ANY, self.figure)
 
        return
 
    def OnCalculate(self, event):
        # Overload calculate button
        # Do something interesting here
        print "%s" % self.m_value1.GetValue()
        return
 
app = wx.App()
window = MyFrameSub(None)
window.Show(True)
app.MainLoop()
 