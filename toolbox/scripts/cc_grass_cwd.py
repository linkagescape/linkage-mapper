#!/usr/bin/env python2.6
# Author: Darren Kavanagh

"""Climate Corridor grass module.

"""

import os
import sys
import shutil
import subprocess
#import pdb

import arcpy
import grass.script as grass
import grass.script.setup as gsetup

from cc_config import cc_env
import cc_util

def main(core_list):
    """ """
    try:
        arcpy.AddMessage("Running Grass to create cost-weighted distance  rasters")
        gisbase = os.path.join(cc_env.proj_dir, "gwksp")
        climate_asc = os.path.join(cc_env.out_dir, "cc_climate.asc")
        resist_asc = os.path.join(cc_env.out_dir, "cc_resist.asc")
        location = "gcwd"
        climate_lyr = "climate"
        resist_lyr = "resist"
        core_lyr = "cores"

        arcpy.RasterToASCII_conversion(cc_env.prj_climate_rast, climate_asc)
        arcpy.RasterToASCII_conversion(cc_env.prj_resist_rast, resist_asc)

        cur_path = subprocess.Popen("echo %PATH%", stdout=subprocess.PIPE,
                                    shell=True).stdout.read()

        setup_wrkspace(gisbase, location, climate_asc)

        grass.run_command("r.in.arc", input=climate_asc, output=climate_lyr)
        grass.run_command("r.in.arc", input=resist_asc, output=resist_lyr)
        grass.run_command("v.in.ogr", dsn=cc_env.out_dir, output=core_lyr)

        slope_factor = "1"
        walk_coeff_flat = "1"
        walk_coeff_uphill = str(cc_env.climate_cost)
        walk_coeff_downhill = str(cc_env.climate_cost * -1)
        walk_coeff = (walk_coeff_flat + "," + walk_coeff_uphill + ","
                      + walk_coeff_downhill + "," + walk_coeff_downhill)

        core = "core"
        core_rast = "core_rast"
        gcwd = "gcwd"
        gback = "gback"
        core_points = "corepoints"

        ascii_fld = cc_util.mk_proj_dir("cwdascii")

        for core_no in core_list:
            grass.run_command("v.extract", input=core_lyr,
                output=core, where= cc_env.core_fld +  " = " + core_no,
                overwrite=True)
            grass.run_command("v.to.rast", input=core, output=core_rast,
                              use="val", overwrite=True)
            grass.run_command("r.to.vect", flags="z", input=core_rast,
                output=core_points, feature="point", overwrite=True)
            grass.run_command("r.walk", elevation=climate_lyr,
                friction=resist_lyr, output=gcwd, outdir=gback,
                start_points=core_points, walk_coeff=walk_coeff,
                slope_factor=slope_factor, overwrite=True)

            cwd_ascii = os.path.join(ascii_fld, "cwd_" + core_no + ".asc")
            back_ascii = os.path.join(ascii_fld, "back_" + core_no + ".asc")
            grass.run_command("r.out.arc", input=gcwd, output=cwd_ascii,
                              overwrite=True)
            grass.run_command("r.out.arc", input=gback, output=back_ascii,
                              overwrite=True)
    except Exception:
        raise
    finally:
        os.environ['PATH'] = cur_path  # Revert to original windows path
        shutil.rmtree(gisbase, True)


def setup_wrkspace(gisdbase, location, geo_file):
    """Setup GRASS workspace and modify windows path fGRASS GDAL"""
    gisbase = os.environ['GISBASE']
    gisbase_etc = os.path.join(gisbase, "etc")
    mapset   = "PERMANENT"

    # Update PATH so that GRASS GDAL is called instead of ArcGIS version
    env_list = os.environ['PATH'].split(';')
    env_list.insert(0, gisbase_etc)
    os.environ['PATH'] = ';'.join(env_list)

    grass.create_location(gisdbase, location, filename=geo_file)
    gsetup.init(gisbase, gisdbase, location, mapset)


if __name__ == "__main__":
    # options, flags = grass.parser()
    sys.exit(main())
