#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 16:00:30 2019

@author: gu32kij
"""
# Core Library modules
import sys
from datetime import datetime

# Third party modules
import numpy as np
import yaml
from openmdao.api import ExplicitComponent

# First party modules
import Pymore.utl.coef as coef
from Pymore.app.marc_flight_analysis import flight_analysis
from Pymore.app.marc_flight_analysis_utl import (extract_blade_data,
                                                 fade_beam_props, interp_dui,
                                                 load_pmfa_config, plot_dui,
                                                 plot_rotor_polarcontour,
                                                 plot_sensors,)
from Pymore.app.marc_reference_utl import plot_ref_data
from SONATA.cbm.fileIO.hiddenprints import HiddenPrints
from SONATA.openmdao_utl.doe_utl import filename_generator


class ExplComp_Hover_Analysis(ExplicitComponent):
    """
    A simple Frequency Analysis ExplicitComponent that computes the the fanplot 
    for the rotor with the input beam properties of the blade
    """

    def __init__(self, savepath="uh60a_hover.pkl"):
        super().__init__()
        self.savepath = savepath

    def setup(self):
        self.counter = 0
        self.startTime = datetime.now()
        self.set_input()
        self.set_output()
        self.set_partials()
        (self.path, self.dui, self.sensors, beamprops, massprops) = load_pmfa_config("../../jobs/MonteCarlo/uh60a_hover.yml")

    def set_input(self):
        self.add_input("BeamProps", val=np.zeros((13, 29)), desc="Massterms(6), Stiffness(21), damping(1) and coordinate(1)")

    def set_output(self):
        self.add_output("steady_mean_elastic_tip_response", val=np.zeros((3)))
        self.add_output("steady_stdvs_elastic_tip_response", val=np.zeros((3)))
        self.add_output("steady_mean_bearing_response", val=np.zeros((3)))
        self.add_output("steady_stdvs_bearing_response", val=np.zeros((3)))

    def set_partials(self):
        self.declare_partials("*", "*", method="fd", step=0.05)  # finite differences all partials

    def compute(self, inputs, outputs):
        refProps = np.load("../../jobs/MonteCarlo/baseline_beamprops.npy")
        dynBeamProps = fade_beam_props(refProps, inputs["BeamProps"], self.dui, t0=0.0, t1=1.0)
        BeamProps = {"BLADE_BP_AB01": dynBeamProps}
        with HiddenPrints():
            data = flight_analysis(self.path, self.dui, beamprops=BeamProps, sensors=self.sensors, savepath=self.savepath)
        data["Psi"] = np.deg2rad(-data["SHAFT_SS_POS"][:, 3] + 180)

        data = {k: v[-465:] for k, v in data.items()}  # last two rotor rotations

        # Elastic Response!
        mean_elastic_flap_response = np.mean(data["BLADE_SS_POSTIP_BLADEMF_BE01"][:, [2]])
        stdvs_elastic_flap_response = np.std(data["BLADE_SS_POSTIP_BLADEMF_BE01"][:, [2]])

        mean_elastic_lag_response = np.mean(data["BLADE_SS_POSTIP_BLADEMF_BE01"][:, [1]])
        stdvs_elastic_lag_response = np.std(data["BLADE_SS_POSTIP_BLADEMF_BE01"][:, [1]])

        mean_elastic_torsion_response = np.mean(data["BLADE_SS_POSTIP_BLADEMF_BE01"][:, [3]]) * 180 / np.pi
        stdvs_elastic_torsion_response = np.std(data["BLADE_SS_POSTIP_BLADEMF_BE01"][:, [3]]) * 180 / np.pi

        # Bearing Response
        mean_flap_response = np.mean(data["ATTACHMENT_SS_FLAP_BE01"][:, [5]]) * 180 / np.pi
        stdvs_flap_response = np.std(data["ATTACHMENT_SS_FLAP_BE01"][:, [5]]) * 180 / np.pi

        mean_lag_response = np.mean(data["ATTACHMENT_SS_LAG_BE01"][:, [5]]) * 180 / np.pi
        stdvs_lag_response = np.std(data["ATTACHMENT_SS_LAG_BE01"][:, [5]]) * 180 / np.pi

        mean_feathering_response = np.mean(data["ATTACHMENT_SS_FEATHERING_BE01"]) * 180 / np.pi
        stdvs_feathering_response = np.std(data["ATTACHMENT_SS_FEATHERING_BE01"]) * 180 / np.pi

        outputs["steady_mean_elastic_tip_response"] = np.array([mean_elastic_flap_response, mean_elastic_lag_response, mean_elastic_torsion_response])
        outputs["steady_stdvs_elastic_tip_response"] = np.array([stdvs_elastic_flap_response, stdvs_elastic_lag_response, stdvs_elastic_torsion_response])
        outputs["steady_mean_bearing_response"] = np.array([mean_flap_response, mean_lag_response, mean_feathering_response])
        outputs["steady_stdvs_bearing_response"] = np.array([stdvs_flap_response, stdvs_lag_response, stdvs_feathering_response])
