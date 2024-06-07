# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 15:52:58 2015

 LIFTINGSURFACE.PY ============================================================
 This module contains the definition of the class of 3d lifting surfaces.
 This class can be instantiated to generate wings, tailplanes, fins, propeller-
 or rotor blades, etc.

 This is an OCC_AirCONICS file, based on the Rhino 'AirCONICS' plugin
 by A. Sobester: https://github.com/sobester/AirCONICS
 ==============================================================================

@author: pchambers
"""
import numpy as np
# from base import AirconicsShape
# from primitives import Airfoil
# import AirCONICStools as act


import sys
sys.path.insert(0,'/Users/rfeil/work/6_SONATA/SONATA')  # import sys path to import 'SONATA' & 'job' modules
from SONATA.airconics_blade_cad.airconics.base import AirconicsShape
import SONATA.airconics_blade_cad.airconics.AirCONICStools as act


from OCC.Core.gp import gp_XOY, gp_Ax3, gp_Dir
from OCC.Core.GeomAbs import GeomAbs_C2
from OCC.Core.Geom import Geom_Plane

class LiftingSurface(AirconicsShape):
    """Airconics class for defining lifting surface shapes

    Parameters
    ----------
    AirfoilFunct - function
        function defining the sectional Airfoil (see primitives.Airfoil)
        vs epsilon spanwise variable coordinate between 0 and 1
        (curvilinear attached).
        Updating will rebuild the geometry.

    NSegments - int (default = 11)
        Number of segments to sample the wing defined by input functions.
        Updating will rebuild the geometry.

    TipRequired - bool (default = False)
        TODO: Not yet used
        adds the wing tip face to components if true

    max_degree - (default = 8)
        maximum degree of the fitted NURBS surface

    continuity - OCC.GeomAbs.GeomAbs_XX Type
        the order of continuity i.e. C^0, C^1, C^2... would be
        GeomAbs_C0, GeomAbs_C1, GeomAbs_C2 ...

    construct_geometry : bool
        If true, Build method will be called on construction

    Attributes
    ----------
    self['Surface'] : TopoDS_Shape
        The generated lifting surface

    Sections : list of airconics.primitives.Airfoil objects
        The rib curves from which the main surface is lofted. Updating any of
        the spanwise functions (ChordFunct, TwistFunct, etc...), SpanFactor,
        ChordFactor Raises an
        error if attempting to write over it manually.

    Notes
    -----
    * Output surface is stored in self['Surface']
    *
    * See airconics.examples.wing_example_transonic_airliner for
      example input functions

    See also
    --------
    airconics.primitives.Airfoil,
    airconics.examples.wing_example_transonic_airliner
    """

    def __init__(self,
                 AirfoilFunct=False,
                 NSegments=11,
                 TipRequired=False,
                 max_degree=8,
                 continuity=GeomAbs_C2,
                 construct_geometry=True,
                 ):
        # convert ApexPoint from list if necessary

        self.CreateConstructionGeometry()

#        Initialise the components using base class:
        super(LiftingSurface, self).__init__(components={},
                                             _AirfoilFunct=AirfoilFunct,
                                             _Sections=[],
                                             _NSegments=NSegments,
                                             max_degree=max_degree,
                                             Cont=continuity,
                                             construct_geometry=construct_geometry
                                             )

    # Properties
    # ----------------------------------------------------------------------
    # Airfoil
    @property
    def AirfoilFunct(self):
        return self._AirfoilFunct

    @AirfoilFunct.setter
    def AirfoilFunct(self, newAirfoilFunct):
        # Maybe add some tests here
        self._AirfoilFunct = newAirfoilFunct
        if self.construct_geometry:
            self.Build()

    @property
    def NSegments(self):
        return self._NSegments

    @NSegments.setter
    def NSegments(self, newNSegments):
        self._NSegments = newNSegments
        if self.construct_geometry:
            self.Build()

    @property
    def Sections(self):
        # Note: there is no setter for this, as it should not be directly
        # interacted with
        return self._Sections
    # ----------------------------------------------------------------------

    def CreateConstructionGeometry(self):
        """
        Creates the plane and vector used for projecting wetted area
        """
        self.XoY_Plane = Geom_Plane(gp_Ax3(gp_XOY()))
        self.ProjVectorZ = gp_Dir(0, 0, 1)

    def Build(self):
        """Builds the section curves and lifting surface using the current

        Uses the current ChordFactor, ScaleFactor, NSegments, ApexPoint, and
        all spanwise variation functions (e.g. ChordFunct) defined in self to
        produce a surface

        Notes
        -----
        Called on initialisation of a lifting surface class.

        :Example:
            >>> Wing = liftingsurface.LiftingSurface(
                                                myAirfoilFunction)
            >>> Surface = Wing['Surface']

        See Also
        --------
        airconics.examples.wing_example_transonic_airliner
        """
        super(LiftingSurface, self).Build()
        assert(self.AirfoilFunct), 'No Airfoil function defined'

        # The first time all functions are successfully defined, change the
        # build flag to true so that further changes to these functions will
        # rebuild geometry via the setter property functions
        self.construct_geometry = True

        self.GenerateSectionCurves()
        self.GenerateLiftingSurface()

        # Also store the tip leading edge point (for fitting tip devices)?
        # self.TipLE = self.Sections[-1].Chord.

    def GenerateSectionCurves(self):
        """Generates the loft section curves  based on the current
        functional parameters and ChordFactor of the object.

        Uses self._AirfoilFunct, _ChordFunct etc. and other attributes to
        update the content of self._Sections.

        Returns
        -------
        None
        """
        # Empty the current geometry
        self._Sections = []

        Eps = np.linspace(0, 1, self.NSegments + 1)

        for i, eps in enumerate(Eps):
            Af = self.AirfoilFunct(Epsilon=eps)
            self._Sections.append(Af)

    def GenerateLiftingSurface(self):
        """Builds a lifting surface (wing, tailplane, etc.) with the Chord and
        Scale factors, and section list defined in self.

        This function should be called after GenerateSectionCurves. Note that
        both operations are performed with `Build', which should be used
        preferentially

        Returns
        -------
        None

        Notes
        -----
        Adds a ('Surface': Shape) key value pair to self.
        """
        LS = act.AddSurfaceLoft(self._Sections,
                                max_degree=self.max_degree,
                                continuity=self.Cont,
                                solid=False)

        if LS is None:
            raise ValueError("Failed to generate Lofted Surface")
        #######################################################################
            # TODO: backup surface loft is not yet implemented for OCC
            # Version of Airconics (this is legacy from the Rhino plugin)
        #######################################################################

        #  Update instance components:
        self.AddComponent(LS, 'Surface')
        return None
