# -*- coding: utf-8 -*-
"""
Created on Wed Jan 04 13:52:57 2017

@author: TPflumm

"""


from OCC.Geom2d import Geom2d_BezierCurve
from OCC.Geom2dAPI import Geom2dAPI_PointsToBSpline
from OCC.Geom2dConvert import geom2dconvert_CurveToBSplineCurve

from SONATA.cbm.topo.BSplineLst_utils import get_BSplineLst_Pnt2d, get_BSplineLst_D2, trim_BSplineLst, get_BSplineLst_length
from SONATA.cbm.topo.utils import point2d_list_to_TColgp_Array1OfPnt2d


def cutoff_layer(Trimmed_BSplineLst,OffsetBSplineLst,S1,S2,cutoff_style=2):
    ''' cutoff_layer creates a cutoff to connect the offsetBSpline List with the original BSplineList '''
    #cutoff_style: 0 step, 1 linear, 2 round, 3 Bezier
    
    Offset_StartPnt = get_BSplineLst_Pnt2d(OffsetBSplineLst,S1,S1,S2)
    Offset_EndPnt = get_BSplineLst_Pnt2d(OffsetBSplineLst,S2,S1,S2)
    
    Org_StartPnt = get_BSplineLst_Pnt2d(Trimmed_BSplineLst,S1,S1,S2)   
    Org_EndPnt = get_BSplineLst_Pnt2d(Trimmed_BSplineLst,S2,S1,S2)
    
    dist = Org_StartPnt.Distance(Offset_StartPnt) #layer thickness
    
    if S1>S2:
        S1, S2 = S2, S1
            

    length = get_BSplineLst_length(OffsetBSplineLst)   
    paralength = abs(S2-S1)
    cutoff_depth = 1.2*dist/length*paralength
    
    #check if OffsetBSplineLst is closed:
    if not Offset_StartPnt.IsEqual(Offset_EndPnt, 1e-4) and not(S1 == 0.0 and S2==1.0):
        if cutoff_style == 0: # STEP-CUTOFF
            Offset_StartPnt = get_BSplineLst_Pnt2d(OffsetBSplineLst,S1,S1,S2)
            Offset_EndPnt = get_BSplineLst_Pnt2d(OffsetBSplineLst,S2,S1,S2)
            
            Start_bspline = Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d[Org_StartPnt,Offset_StartPnt]).Curve()
            End_bspline = Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([Offset_EndPnt,Org_EndPnt])).Curve()
                       
        elif cutoff_style == 1: #LINEAR-CUTOFF
            #cutoff_depth = 1.05*dist
            OffsetBSplineLst = trim_BSplineLst(OffsetBSplineLst, S1+cutoff_depth, S2-cutoff_depth, S1, S2)
            Offset_StartPnt = get_BSplineLst_Pnt2d(OffsetBSplineLst,S1,S1,S2)
            Offset_EndPnt = get_BSplineLst_Pnt2d(OffsetBSplineLst,S2,S1,S2)
            
            Start_bspline = Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([Org_StartPnt,Offset_StartPnt])).Curve()
            End_bspline = Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([Offset_EndPnt,Org_EndPnt])).Curve()
        
        elif cutoff_style == 2: #ROUND-CUTOFF
            #cutoff_depth = 1.05*dist
            OffsetBSplineLst = trim_BSplineLst(OffsetBSplineLst, S1+cutoff_depth, S2-cutoff_depth, S1, S2)
            Org_StartPnt, Org_V1Start, Org_V2Start = get_BSplineLst_D2(Trimmed_BSplineLst,S1,S1,S2)
            Org_EndPnt,  Org_V1End, Org_V2End = get_BSplineLst_D2(Trimmed_BSplineLst,S2,S1,S2)
                    
            Offset_StartPnt, Offset_V1Start, Offset_V2Start = get_BSplineLst_D2(OffsetBSplineLst,S1,S1,S2)
            Offset_EndPnt, Offset_V1End, Offset_V2End = get_BSplineLst_D2(OffsetBSplineLst,S2,S1,S2)
                    
            Offset_V1Start.Normalize()
            Offset_V1End.Normalize()
    
            kappa_bezier = 0.5*dist  
            Offset_V1Start.Multiply(-kappa_bezier)
            Offset_V1End.Multiply(kappa_bezier)
                   
            Start_Bezier_PntList = [Org_StartPnt, Offset_StartPnt.Translated(Offset_V1Start), Offset_StartPnt]
            End_Bezier_PntList = [Offset_EndPnt,Offset_EndPnt.Translated(Offset_V1End),Org_EndPnt]                        
            
            Start_Bezier = Geom2d_BezierCurve(point2d_list_to_TColgp_Array1OfPnt2d(Start_Bezier_PntList))
            End_Bezier = Geom2d_BezierCurve(point2d_list_to_TColgp_Array1OfPnt2d(End_Bezier_PntList)) 
                    
            Start_bspline = geom2dconvert_CurveToBSplineCurve(Start_Bezier)  
            End_bspline = geom2dconvert_CurveToBSplineCurve(End_Bezier)
            
        elif cutoff_style == 3: #BEZIER-CUTOFF
            #cutoff_depth = 1.05*dist
            OffsetBSplineLst = trim_BSplineLst(OffsetBSplineLst, S1+cutoff_depth, S2-cutoff_depth, S1, S2)
            Org_StartPnt, Org_V1Start, Org_V2Start = get_BSplineLst_D2(Trimmed_BSplineLst,S1,S1,S2)
            Org_EndPnt,  Org_V1End, Org_V2End = get_BSplineLst_D2(Trimmed_BSplineLst,S2,S1,S2)
                    
            Offset_StartPnt, Offset_V1Start, Offset_V2Start = get_BSplineLst_D2(OffsetBSplineLst,S1,S1,S2)
            Offset_EndPnt, Offset_V1End, Offset_V2End = get_BSplineLst_D2(OffsetBSplineLst,S2,S1,S2)
                    
            Org_V1Start.Normalize()
            Org_V1End.Normalize()
            Offset_V1Start.Normalize()
            Offset_V1End.Normalize()
    
            kappa_bezier = cutoff_depth
            Org_V1Start.Multiply(kappa_bezier) 
            Org_V1End.Multiply(-kappa_bezier)      
            Offset_V1Start.Multiply(-kappa_bezier)
            Offset_V1End.Multiply(kappa_bezier)
                   
            Start_Bezier_PntList = [Org_StartPnt, Org_StartPnt.Translated(Org_V1Start), Offset_StartPnt.Translated(Offset_V1Start), Offset_StartPnt]
            End_Bezier_PntList = [Offset_EndPnt,Offset_EndPnt.Translated(Offset_V1End),Org_EndPnt.Translated(Org_V1End),Org_EndPnt]                        
            
            Start_Bezier = Geom2d_BezierCurve(point2d_list_to_TColgp_Array1OfPnt2d(Start_Bezier_PntList))
            End_Bezier = Geom2d_BezierCurve(point2d_list_to_TColgp_Array1OfPnt2d(End_Bezier_PntList)) 
                    
            Start_bspline = geom2dconvert_CurveToBSplineCurve(Start_Bezier)   
            End_bspline = geom2dconvert_CurveToBSplineCurve(End_Bezier)  
            

        OffsetBSplineLst.insert(0,Start_bspline)
        OffsetBSplineLst.append(End_bspline)
        
    else:
        None

    return OffsetBSplineLst

if __name__ == '__main__':
    exec(compile(open("SONATA.py").read(), "SONATA.py", 'exec'))
