#!/usr/bin/env python

# Copyright 2009-2016 Thomas Paviot (tpaviot@gmail.com)
##
# This file is part of pythonOCC.
##
# pythonOCC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
##
# pythonOCC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
##
# You should have received a copy of the GNU Lesser General Public License
# along with pythonOCC.  If not, see <http://www.gnu.org/licenses/>.
import sys
import wx
from OCC.Display.wxDisplay import wxViewer3d
from OCC import VERSION
from OCC.Display.backend import load_backend, get_qt_modules

def check_callable(_callable):
    assert callable(_callable), 'the function supplied is not callable'

class AppFrame(wx.Frame):

        def __init__(self, parent):
            wx.Frame.__init__(self, parent, -1, "SONATA GUI - pythonOCC-%s 3d viewer" % VERSION, style=wx.DEFAULT_FRAME_STYLE, size=(1024, 768))
            self.canva = wxViewer3d(self)
            self.menuBar = wx.MenuBar()
            self._menus = {}
            self._menu_methods = {}
            self.sb = self.CreateStatusBar( 1, wx.ST_SIZEGRIP, wx.ID_ANY )
            
        def add_menu(self, menu_name):
            _menu = wx.Menu()
            self.menuBar.Append(_menu, "&" + menu_name)
            self.SetMenuBar(self.menuBar)
            self._menus[menu_name] = _menu

        def add_function_to_menu(self, menu_name, _callable):
            # point on curve
            _id = wx.NewId()
            check_callable(_callable)
            try:
                self._menus[menu_name].Append(_id, _callable.__name__.replace('_', ' ').lower())
            except KeyError:
                raise ValueError('the menu item %s does not exist' % menu_name)
            self.Bind(wx.EVT_MENU, _callable, id=_id)

            
            
def init_display():
    """ This function loads and initialize a GUI using either wx, pyq4, pyqt5 or pyside.
    If ever the environment variable PYTHONOCC_SHUNT_GUI, then the GUI is simply ignored.
    It can be useful to test some algorithms without being polluted by GUI statements.
    This feature is used for running the examples suite as tests for
    pythonocc-core development.
    """    
    # wxPython based simple GUI
    
    app = wx.App(False)
    win = AppFrame(None)
    win.Show(True)
    wx.SafeYield()
    win.canva.InitDriver()
    app.SetTopWindow(win)
    display = win.canva._display
    display.SetBackgroundImage('default_background.png')    

    def add_menu(*args, **kwargs):
        win.add_menu(*args, **kwargs)

    def add_function_to_menu(*args, **kwargs):
        win.add_function_to_menu(*args, **kwargs)

    def start_display():
        app.MainLoop()

    return display, start_display, add_menu, add_function_to_menu


#=================================================================================================
if __name__ == '__main__':
    display, start_display, add_menu, add_function_to_menu = init_display()
    from OCC.BRepPrimAPI import BRepPrimAPI_MakeSphere, BRepPrimAPI_MakeBox

    def sphere(event=None):
        display.DisplayShape(BRepPrimAPI_MakeSphere(100).Shape(), update=True)

    def cube(event=None):
        display.DisplayShape(BRepPrimAPI_MakeBox(1, 1, 1).Shape(), update=True)

    def exit(event=None):
        sys.exit()

    add_menu('File')
    add_menu('Edit')
    add_menu('Tools')
    add_menu('primitives')
    add_menu('Help')
    add_function_to_menu('primitives', sphere)
    add_function_to_menu('primitives', cube)
    add_function_to_menu('primitives', exit)
    start_display()