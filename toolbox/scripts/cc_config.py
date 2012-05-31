#!/usr/bin/env python2.6
# Author: Darren Kavanagh

"""Climate Corridor configuration module.

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
    """Class container to holdreu Climate Tool global variables"""

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
        self.min_euc_dist = arg[7]  # Min distance between core pairs
        self.max_euc_dist = arg[8]  # Max distance between core pairs
        self.climate_threashold = arg[9]  # Temperature threashold
        self.climate_cost = arg[10]  # Cost incurred by climate variable

        # Setup model global variables
        path = os.path.abspath(__file__)
        self.code_dir = os.path.dirname(os.path.abspath(__file__))
        self.out_dir = os.path.join(self.proj_dir, "cc")  # Output directory
        self.prj_area_rast = os.path.join(self.out_dir, "projarea")
        self.prj_core_fc = (os.path.join(
                            self.out_dir, "cores.shp")) # Proj core area  name
        prj_climate_rast, ext =  os.path.splitext(os.path.basename(
            self.climate_rast))   # Climate raster file name
        self.prj_climate_rast = os.path.join(self.out_dir,
                                             prj_climate_rast[:8])
        self.simplyfy_cores = True

        # Use project area raster if no resistance layer is provided
        if self.resist_rast is not None:
            prj_resist_rast, ext =  os.path.splitext(os.path.basename(
                self.resist_rast))   # Climate raster file name
            self.prj_resist_rast = os.path.join(self.out_dir,
                                                prj_resist_rast[:8])
        else:
            self.prj_resist_rast = self.prj_area_rast


cc_env = ClimateConfig()
