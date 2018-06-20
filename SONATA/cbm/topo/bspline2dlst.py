from OCC.gp import gp_Pnt2d
from OCC.GCPnts import GCPnts_AbscissaPoint
from OCC.Geom2dAPI import Geom2dAPI_PointsToBSpline, Geom2dAPI_ProjectPointOnCurve
from OCC.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.Geom2d import Geom2d_BSplineCurve, Handle_Geom2d_BSplineCurve_DownCast

import os
import numpy as np
os.chdir('C:\\TPflumm_local\\work\\SONATA\\')
from SONATA.cbm.topo.utils import point2d_list_to_TColgp_Array1OfPnt2d
from SONATA.cbm.fileIO.readinput import UIUCAirfoil2d, AirfoilDat2d

from SONATA.cbm.topo.BSplineLst_utils import get_BSpline_length, get_BSplineLst_length, \
                            find_BSplineLst_coordinate, get_BSplineLst_Pnt2d, \
                            trim_BSplineLst, seg_boundary_from_dct, set_BSplineLst_to_Origin, \
                            copy_BSplineLst, trim_BSplineLst_by_Pnt2d, findPnt_on_2dcurve, \
                            isPnt_on_2dcurve

from SONATA.cbm.topo.utils import calc_DCT_angles, TColgp_HArray1OfPnt2d_from_nparray, Pnt2dLst_to_npArray, \
                    discrete_stepsize, curvature_of_curve, isclose, unique_rows, \
                    P2Pdistance, PolygonArea, TColgp_Array1OfPnt2d_from_nparray


class BSpline2dLst(list):
    '''
    Describes a list of OCC.Geom2d.Geom2d_BSplineCurve and is a subclass of list
    '''
    
    def __init__(self, lst):
        #TODO: make sure that lst is of type list with BSplines!
        list.__init__(self,lst)
        self.a = 1
        
        
    def corners(self):
        corners = [] 
        for item in self:
            corners.append(item.EndPoint())        
        corners.pop(-1)
        return corners #gp_Pnt2d Lst
    
    
    def contains_Pnt2d(self,Pnt2d,tolerance=1e-6):
        Bool = False
        for item in self:
            if isPnt_on_2dcurve(Pnt2d,item.GetHandle(),tolerance):
                Bool = True
                break
        return Bool
    
    
    def find_Pnt2d(self,Pnt2d):
        for i,item in enumerate(self):
            if isPnt_on_2dcurve(Pnt2d,item.GetHandle()):
                u = findPnt_on_2dcurve(Pnt2d,item.GetHandle())
                coordinates = (i,u)
                break
            else:
                coordinates = None
        return coordinates
        
    
    def project_Pnt2d(self, Pnt2d, tolerance_distance = 100):
        ''' returns tuple of the closest orthogonal projects of Pnt2d on BSplineLst that are 
            within the tolerance_distance as list
        '''
        poclst = []
        for idx,item in enumerate(self,):
            poc = Geom2dAPI_ProjectPointOnCurve(Pnt2d,item.GetHandle())
            for j in range(1,poc.NbPoints()+1):
                poclst.append((poc.NearestPoint(), idx, poc.LowerDistanceParameter(), poc.LowerDistance()))
        
        poclst = np.asarray(poclst)
        min_index = poclst[:,3].argmin() 
        return tuple(poclst[min_index,:])
        
    
    def distance_btw_coords(self,para1,para2):
        
        
    #    para1 = findPnt_on_BSplineLst(P1,BSplineLst)
    #    para2 = findPnt_on_BSplineLst(P2,BSplineLst)
        
    #    if closed:
    #        print para1,para2
    #        para2 = [len(BSplineLst)-1, BSplineLst[-1].LastParameter()]
    #        print para1,para2
        tol=1e-7
        #print para1,para2
        
        Distance = 0
        for i,item in enumerate(self):
            Adaptor = Geom2dAdaptor_Curve(item.GetHandle())
            First = item.FirstParameter() 
            Last =  item.LastParameter() 
            if para1[0] == i and para2[0] != i:
                Length = GCPnts_AbscissaPoint().Length(Adaptor, para1[1], Last, tol)
                Distance += Length
                 
            elif (para1[0] != i and para2[0] != i) and (para1[0] < i and para2[0] > i):
                Length = GCPnts_AbscissaPoint().Length(Adaptor, First, Last, tol)
                Distance += Length
                 
            elif para1[0] == i and para2[0] == i:
                Length = GCPnts_AbscissaPoint().Length(Adaptor, para1[1], para2[1], tol)
                Distance += Length
                break
        
        #NOTE: Somehow the execption if para2 is close to Last doesn't work properly!!!!!!! Similar to the trimm function     
            elif para1[0] != i and para2[0] == i:
                if isclose(para2[1],Last):
                     Length = GCPnts_AbscissaPoint().Length(Adaptor, First, Last, tol)
                     Distance += Length
                else:
                     Length = GCPnts_AbscissaPoint().Length(Adaptor, First, para2[1], tol)
                     Distance += Length
                break
    
        return Distance	


    

if __name__ == '__main__': 
    DCT_data = UIUCAirfoil2d('naca23012').T
    DCT_data = np.multiply(DCT_data,1.0)
    testlst = seg_boundary_from_dct(DCT_data)
    a_Bspline2dLst = BSpline2dLst(testlst)
    print('a_Bspline2dLst:', a_Bspline2dLst)
    print('.corners:', a_Bspline2dLst.corners())
    p = gp_Pnt2d(0,0)
    print('.contains_Pnt2d(gp_Pnt2d(0,0)):', a_Bspline2dLst.contains_Pnt2d(p))
    print('.find_Pnt2d(gp_Pnt2d(0,0)):', a_Bspline2dLst.find_Pnt2d(p))
    print('.project_Pnt2d(gp_Pnt2d(0,0.001)),', a_Bspline2dLst.project_Pnt2d(gp_Pnt2d(0,0o01)))
    p = gp_Pnt2d(0,1)
    print('.contains_Pnt2d(gp_Pnt2d(0,0)):', a_Bspline2dLst.contains_Pnt2d(p))
    print('.find_Pnt2d(gp_Pnt2d(0,0)):', a_Bspline2dLst.find_Pnt2d(p))
    p = gp_Pnt2d(0,1)
