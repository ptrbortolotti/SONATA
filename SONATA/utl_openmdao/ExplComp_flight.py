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
from Pymore.app.fft import extract_vibratory_hubloads
from Pymore.app.marc_flight_analysis import flight_analysis
from Pymore.app.marc_flight_analysis_utl import (extract_blade_data,
                                                 fade_beam_props,
                                                 interp_azimuth, interp_dui,
                                                 load_pmfa_config, plot_dui,
                                                 plot_rotor_polarcontour,
                                                 plot_sensors,)
from Pymore.app.marc_reference_utl import plot_ref_data
from SONATA.cbm.fileIO.hiddenprints import HiddenPrints
from SONATA.openmdao_utl.doe_utl import filename_generator


class ExplComp_Flight_Analysis(ExplicitComponent):
    """

    """

    def __init__(self, savepath="uh60a_c8513.pkl"):
        super().__init__()
        self.savepath = savepath

    def setup(self):
        self.counter = 0
        self.startTime = datetime.now()
        self.set_input()
        self.set_output()
        self.set_partials()
        (self.path, self.dui, self.sensors, beamprops, massprops) = load_pmfa_config("jobs/MonteCarlo/uh60a_c8513.yml")

    def set_input(self):
        self.add_input("BeamProps", val=np.zeros((13, 29)), desc="Massterms(6), Stiffness(21), damping(1) and coordinate(1)")

    def set_output(self):
        self.add_output("vibratory_hubloads", val=np.zeros((6)))
        self.add_output("mean_elastic_tip_response", val=np.zeros((4, 181)))
        self.add_output("stdvs_elastic_tip_response", val=np.zeros((4, 181)))

    def set_partials(self):
        self.declare_partials("*", "*", method="fd", step=0.05)  # finite differences all partials

    def compute(self, inputs, outputs):

        refProps = np.load("jobs/MonteCarlo/baseline_beamprops.npy")
        dynBeamProps = fade_beam_props(refProps, inputs["BeamProps"], self.dui, t0=0.0, t1=1.0)
        BeamProps = {"BLADE_BP_AB01": dynBeamProps}
        with HiddenPrints():
            data = flight_analysis(self.path, self.dui, beamprops=BeamProps, sensors=self.sensors, savepath=self.savepath)
        hubforces = data["SHAFT_SS_HUBFORCES"]
        time = data["time"]
        samples = 465

        # ========================= Normalized Vibratory Hubloads  ============
        peaks = extract_vibratory_hubloads(time, hubforces, samples=samples)
        steady_thrust = -np.mean(data["SHAFT_SS_HUBFORCES"][:, 2][-samples:])
        seady_torque = np.mean(data["SHAFT_SS_HUBFORCES"][:, 5][-samples:])

        # ======================== Blade Tip-Response over azimuth =============
        psi = data["Psi"][-samples:]

        y = data["BLADE_SS_POSTIP_BLADEMF_BE01"][-samples:, 2]
        (xint, mean_elastic_flap_response, stdvs_elastic_flap_response) = interp_azimuth(psi, y)

        y = data["BLADE_SS_POSTIP_BLADEMF_BE01"][-samples:, 1]
        (xint, mean_elastic_lag_response, stdvs_elastic_lag_response) = interp_azimuth(psi, y)

        y = data["BLADE_SS_POSTIP_BLADEMF_BE01"][-samples:, 3] * 180 / np.pi
        (xint, mean_elastic_torsion_response, stdvs_elastic_torsion_response) = interp_azimuth(psi, y)

        outputs["vibratory_hubloads"] = np.hstack((peaks[:3, 1] / steady_thrust, peaks[3:, 1] / seady_torque))
        outputs["mean_elastic_tip_response"] = np.array([xint, mean_elastic_flap_response, mean_elastic_lag_response, mean_elastic_torsion_response])
        outputs["stdvs_elastic_tip_response"] = np.array([xint, stdvs_elastic_flap_response, stdvs_elastic_lag_response, stdvs_elastic_torsion_response])
