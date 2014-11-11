#!/usr/bin/env python2.6
# Authors: Darren Kavanagh and Brad McRae

"""Climate Linkage Mapper configuration module.

Assigns input parameters from ToolBox to variables, and sets constants

"""

import os
import sys


def nullstring(arg_string):
    """Convert ESRI nullstring to Python null"""
    if arg_string == "#":
        arg_string = None
    return arg_string


class ClimateConfig():
    """Class container to hold Climate Tool global variables"""
    def __init__(self):
        """Init class (empty)"""
        pass

    def configure(self, arg):
        """Setup input parameters"""
        # Input files
        self.proj_dir = arg[1]  # Project directory
        self.core_fc = arg[2]  # Core area feature class
        self.core_fld = arg[3]  # Core area field name
        self.climate_rast = arg[4]  # Climate raster (+ path)
        self.resist_rast = nullstring(arg[5])  # Resistance raster (+ path)

        # Setup GRASS base folder environmental setting
        self.gisbase = os.environ['GISBASE'] = arg[6]  # GRASS path
        sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))

        # Tool settings
        self.min_euc_dist = float(arg[7])  # Min distance between core pairs
        self.max_euc_dist = float(arg[8])  # Max distance between core pairs
        self.climate_threshold = float(arg[9])  # Climate threshold
        self.climate_cost = float(arg[10])  # Cost incurred by climate variable

        # Prune network settings
        # For Linkage Mapper. No type conversions necessary
        self.prune_network = arg[11]
        self.max_nn = arg[12]  # No of connected nearest neighbors
        self.nn_unit = arg[13]  # NN Unit
        self.keep_constelations = arg[14]

        # Setup model global variables
        self.code_dir = os.path.dirname(os.path.abspath(__file__))
        self.out_dir = os.path.join(self.proj_dir, "clm_cor")  # CC directory
        self.tmp_dir = os.path.join(self.proj_dir, "tmp")
        self.inputs_gdb = os.path.join(self.out_dir, "LM_inputs.gdb")

        # Define project area files
        # Using gdb seems to best best avoids dll conflict
        self.prj_core_fc = (os.path.join(self.inputs_gdb, "cc_cores"))
        self.prj_climate_rast = os.path.join(self.inputs_gdb, "climate")
        self.prj_core_rast = os.path.join(self.inputs_gdb, "cores")
        self.prj_resist_rast = os.path.join(self.inputs_gdb, "cc_resist")
        self.simplify_cores = True
        self.core_simp = os.path.join(cc_env.out_dir, "coresim.shp")

cc_env = ClimateConfig()
