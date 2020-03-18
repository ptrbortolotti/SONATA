# -*- coding: utf-8 -*-

"""
CLI driver for generating wind turbine CAD models

"""

import os
import argparse
from .airfoil_cst import TurbineCAD
from .airconics.liftingSurface import LiftingSurface

class TurbineCLI:
    """Entry point for the command-line driver for CAD generation"""

    #: Description of the CLI app
    description = "Generate a CAD model of a wind turbine based on YAML definitions"

    #: Epilog for help messages
    epilog = "Wind Turbine CAD Generator"

    def __init__(self):
        self.parser = argparse.ArgumentParser(
            description=self.description,
            epilog=self.epilog)
        self.cli_options()
        self.args = self.parser.parse_args()

    def cli_options(self):
        """Set up the CLI options and arguments"""
        parser = self.parser
        parser.add_argument(
            '-o', '--output-file', default='turbine_cad.iges',
            help="IGES output file name for the CAD model")
        parser.add_argument(
            '-i', '--input-file', default="turbine_cad.yaml",
            help="Input file containing the turbine description")

    def __call__(self):
        """Run the command"""
        args = self.args

        inpfile = os.path.abspath(args.input_file)
        if not os.path.exists(inpfile):
            print("ERROR: Input file does not exist = ", args.input_file)
            sys.parser.exit(1)

        turbine = TurbineCAD(inpfile)
        print("==> Creating blade CST")
        turbine.create_blade_cst()
        print("==> Lofting blade")
        blade = LiftingSurface(turbine.airfoil_func)
        print("==> Finished lofting blade")

        outfile = args.output_file
        status = blade.Write(outfile)
        print("==> IGES file written to file %s with status %d"%(outfile, status))

def main():
    """Entry point for CLI app"""
    cli = TurbineCLI()
    cli()
