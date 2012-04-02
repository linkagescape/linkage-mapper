#!/usr/bin/env python2.6
# Author: Darren Kavanagh

"""Climate Corridor configuration module.

Assigns input parameters from ToolBox to variables, and sets constants

""" 

# $Revision$

import os.path as path


def nullstring(arg_string):
    """Convert ESRI nullstring to Python null"""
    if arg_string == "#":
        arg_string = None
    return arg_string

class ClimateConfig():
    """Class container to holdreu Climate Tool global variables"""

    def configure(self, arg):
        """Setup input parameters"""
        self.proj_dir = arg[1]  # Project directory
        self.core_fc = arg[2]  # Core area feature class
        self.core_fld = arg[3]  # Core area field name
        self.climate_rast = arg[4]  # Climate raster (+ path)
        self.resist_rast = nullstring(arg[5])  # Resistance raster (+ path)
        self.search_radius = arg[6]  # Temperature threashold
        self.min_euc_dist = arg[7]  # Min distance between core pairs
        self.max_euc_dist = arg[8]  # Max distance between core pairs
        self.climate_threashold = arg[9]  # Temperature threashold
        self.climate_cost = arg[10]  # Cost incurred by climate variable

        self.out_dir = path.join(self.proj_dir, "cc")  # Output directory
        self.prj_area_rast = path.join(self.out_dir, "parea_rast")
        self.prj_core_fc = path.join(self.out_dir, "pcores.shp")  # Proj core area  name
        
        prj_climate_rast, ext =  path.splitext(path.basename(
            self.climate_rast))   # Climate raster file name
        self.prj_climate_rast = path.join(self.out_dir, prj_climate_rast)
        
        # Use project area raster if no resistance layer is provided
        if self.resist_rast is not None:
            prj_resist_rast, ext =  path.splitext(path.basename(
                self.resist_rast))   # Climate raster file name
            self.prj_resist_rast = path.join(self.out_dir, prj_resist_rast)            
        else:
            self.prj_resist_rast = self.prj_area_rast
                

cc_env = ClimateConfig()
