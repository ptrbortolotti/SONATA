
from OCC.Core.gp import gp_OX2d
from OCC.Core.GCE2d import GCE2d_MakeEllipse
from OCC.Core.Geom2d import Geom2d_TrimmedCurve
from OCC.Core.Geom2dConvert import geom2dconvert_CurveToBSplineCurve
from OCC.Core.Convert import Convert_TgtThetaOver2

from OCC.Display.SimpleGui import init_display
display, start_display, add_menu, add_function_to_menu = init_display()


def curves2d_from_curves():
    major, minor = 12, 4
    axis = gp_OX2d()
    ellipse = GCE2d_MakeEllipse(axis, major, minor).Value()
    trimmed_curve = Geom2d_TrimmedCurve(ellipse, -1, 2, True)
    bspline = geom2dconvert_CurveToBSplineCurve(trimmed_curve, Convert_TgtThetaOver2)
    display.DisplayShape(bspline, update=True)


if __name__ == '__main__':
    curves2d_from_curves()
    start_display()