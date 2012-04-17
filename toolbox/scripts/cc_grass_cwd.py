#!/usr/bin/env python2.6
# Author: Darren Kavanagh

"""Climate Corridor grass module.

"""

import os
import sys
import shutil
import subprocess
import pdb

import arcpy
import grass.script as grass
import grass.script.setup as gsetup

from cc_config import cc_env
import cc_util

def main(core_list):
    """ """
    cur_path = subprocess.Popen("echo %PATH%", stdout=subprocess.PIPE,
                                shell=True).stdout.read()
    gisdbase = os.path.join(cc_env.proj_dir, "gwksp")

    try:
        arcpy.AddMessage("\nRUNNING GRASS TO CREATE COST-WEIGHTED DISTANCE "
                         "RASTERS")
        gisdbase = os.path.join(cc_env.proj_dir, "gwksp")
        climate_asc = os.path.join(cc_env.out_dir, "cc_climate.asc")
        resist_asc = os.path.join(cc_env.out_dir, "cc_resist.asc")
        location = "gcwd"
        climate_lyr = "climate"
        resist_lyr = "resist"
        core_lyr = "cores"

        arcpy.RasterToASCII_conversion(cc_env.prj_climate_rast, climate_asc)
        arcpy.RasterToASCII_conversion(cc_env.prj_resist_rast, resist_asc)

        setup_wrkspace(gisdbase, location, climate_asc)

        arcpy.AddMessage("Importing raster files into GRASS")
        grass.run_command("r.in.arc", input=climate_asc, output=climate_lyr)
        grass.run_command("r.in.arc", input=resist_asc, output=resist_lyr)
        arcpy.AddMessage("Importing cores feature class into GRASS")
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
        gbackrc = "gbackrc"
        core_points = "corepoints"

        ascii_fld = cc_util.mk_proj_dir("cwdascii")

        arcpy.AddMessage("Generating cost-weighted-distance rasters for "
                         + str(len(core_list)) + " cores")

        for position, core_no in enumerate(core_list):
            arcpy.AddMessage("\nCore ID:" + core_no)
            arcpy.AddMessage("Core No:" + str(position + 1))
            grass.run_command("v.extract", input=core_lyr,
                output=core, where= cc_env.core_fld +  " = " + core_no)
            grass.run_command("v.to.rast", input=core, output=core_rast,
                              use="val")
            grass.run_command("r.to.vect", flags="z", input=core_rast,
                output=core_points, feature="point")
            arcpy.AddMessage("Running r.walk to create CWD raster")
            grass.run_command("r.walk", elevation=climate_lyr,
                friction=resist_lyr, output=gcwd, outdir=gback,
                start_points=core_points, walk_coeff=walk_coeff,
                slope_factor=slope_factor)
            
            # Reclassify from the directional degree output from GRASS to
            # Arc's 1 to 8 directions format            
            path = os.path.abspath(__file__)
            dir_path = os.path.dirname(path)
            rc_rules = os.path.join(dir_path, "cc_bkrast.rcf")            
            grass.run_command("r.reclass", input=gback, output=gbackrc,
                               rules=rc_rules)
            
            arcpy.AddMessage("Exporting CWD and back rasters to ASCII grids")
            cwd_ascii = os.path.join(ascii_fld, "cwd_" + core_no + ".asc")
            back_ascii = os.path.join(ascii_fld, "back_" + core_no + ".asc")
            grass.run_command("r.out.arc", input=gcwd, output=cwd_ascii)
            grass.run_command("r.out.arc", input=gbackrc, output=back_ascii)

    except Exception:
        raise
    finally:
        os.environ['PATH'] = cur_path  # Revert to original windows path
        shutil.rmtree(gisdbase, True)


def setup_wrkspace(gisdbase, location, geo_file):
    """Setup GRASS workspace and modify windows path for GRASS GDAL"""
    arcpy.AddMessage("Creating GRASS workspace")
    gisbase = cc_env.gisbase

    os.environ['GISRC'] = os.path.join(os.environ['USERPROFILE'], ".grassrc6")
    os.environ['LD_LIBRARY_PATH'] = os.path.join(gisbase, "lib")
    os.environ['GRASS_SH'] = os.path.join(gisbase, "msys", "bin", "sh.exe")
    # os.environ['PYTHONLIB']

    env_list = os.environ['PATH'].split(';')
    env_list.insert(0, os.path.join(gisbase, "msys", "bin"))
    env_list.insert(0, os.path.join(gisbase, "extralib"))
    env_list.insert(0, os.path.join(gisbase, "bin"))
    env_list.insert(0, os.path.join(gisbase, "lib"))
    env_list.insert(0, os.path.join(gisbase, "etc", "python"))
    env_list.insert(0, os.path.join(gisbase, "etc"))
    os.environ['PATH'] = ';'.join(env_list)

    mapset   = "PERMANENT"

    gdal = subprocess.Popen("where gdal*", stdout=subprocess.PIPE,
                                shell=True).stdout.read()
    arcpy.AddMessage("GDAL: " + gdal)

    grass.create_location(gisdbase, location, filename=geo_file)
    gsetup.init(gisbase, gisdbase, location, mapset)
    grass.run_command("g.gisenv", set="OVERWRITE=1")
    os.environ['GRASS_VERBOSE'] = "0"


if __name__ == "__main__":
    # options, flags = grass.parser()
    sys.exit(main())
