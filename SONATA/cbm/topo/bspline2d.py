from OCC.gp import gp_Pnt2d
from OCC.GCPnts import GCPnts_AbscissaPoint
from OCC.Geom2dAPI import Geom2dAPI_PointsToBSpline
from OCC.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.Geom2d import Geom2d_BSplineCurve, Handle_Geom2d_BSplineCurve_DownCast

import os
os.chdir('C:\\TPflumm_local\\work\\SONATA\\')
from SONATA.cbm.topo.utils import point2d_list_to_TColgp_Array1OfPnt2d
from SONATA.cbm.fileIO.readinput import UIUCAirfoil2d, AirfoilDat2d


class BSpline2d(Geom2d_BSplineCurve):
    '''
    Describes a list of OCC.Geom2d.Geom2d_BSplineCurve curve.
    This object is inherited from the OCC.Geom2d.Geom2d_BSplineCurve Class
    '''
    
    def __init__(self, instance):
        assert isinstance(instance, Geom2d_BSplineCurve), 'need a Geom2d_BSplineCurve, got a %s' % instance.__class__
        super(BSpline2d, self).__init__()

    def getLength(self,tolerance = 1e-7):
        first = self.FirstParameter()
        last = self.LastParameter()
        Adaptor = Geom2dAdaptor_Curve(self)
        length = GCPnts_AbscissaPoint().Length(Adaptor, first, last, tolerance)
        return length
  

    
if __name__ == '__main__':
   
    p1 = gp_Pnt2d(0,0)
    p2 = gp_Pnt2d(2,0)
    tmp = Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([p1,p2])).Curve()
    test = BSpline2d(tmp)
    test.getLength()