#!/usr/bin/env python2.6
# Author: Darren Kavanagh

"""Climate Corridor grass module.

"""

import os
import shutil
import subprocess
# import pickle


import arcpy

# sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))
import grass.script as grass
import grass.script.setup as gsetup

from cc_config import cc_env
# from lm_config import tool_env as lm_env
import cc_util
import lm_util


def grass_cwd(core_list):
# def grass_cwd(pfile_nm):
    """Creating CWD and Back rasters using GRASS r.walk function"""
    # with open(pfile_nm, "rb") as pfile:
        # argv, lm_arg, core_list = pickle.load(pfile)
    # cc_env.configure(argv)
    # lm_env.configure(lm_env.TOOL_CC, lm_arg)

    cur_path = subprocess.Popen("echo %PATH%", stdout=subprocess.PIPE,
                                shell=True).stdout.read()
    gisdbase = os.path.join(cc_env.proj_dir, "gwksp")
    ccr_grassrc = os.path.join(cc_env.proj_dir, "ccr_grassrc")
    climate_asc = os.path.join(cc_env.out_dir, "cc_climate.asc")
    resist_asc = os.path.join(cc_env.out_dir, "cc_resist.asc")
    climate_lyr = "climate"
    resist_lyr = "resist"
    core_lyr = "cores"

    try:
        arcpy.AddMessage("\nRUNNING GRASS TO CREATE COST-WEIGHTED DISTANCE "
                         "RASTERS")

        # Convert input GRID rasters to ASCII
        arcpy.AddMessage("Converting ARCINFO GRID rasters to ASCII")
        arcpy.RasterToASCII_conversion(cc_env.prj_climate_rast, climate_asc)
        arcpy.RasterToASCII_conversion(cc_env.prj_resist_rast, resist_asc)

        # Create resource file and setup workspace
        write_grassrc(ccr_grassrc)
        grass_version = setup_wrkspace(gisdbase, ccr_grassrc, climate_asc)

        # Make cwd folder for Linkage Mapper
        lm_util.make_cwd_paths(max(core_list))

        # Import files into GRASS
        arcpy.AddMessage("Importing raster files into GRASS")
        run_grass_cmd("r.in.arc", input=climate_asc, output=climate_lyr)
        run_grass_cmd("r.in.arc", input=resist_asc, output=resist_lyr)
        arcpy.AddMessage("Importing cores feature class into GRASS\n")
        run_grass_cmd("v.in.ogr", dsn=cc_env.out_dir, output=core_lyr)

        # Generate CWD and Back rasters
        gen_cwd_back(grass_version, core_list, climate_lyr, resist_lyr,
                     core_lyr)
    except Exception:
        raise
    finally:
        os.environ['PATH'] = cur_path  # Revert to original windows path
        try:
            shutil.rmtree(gisdbase, True)
        except OSError:
            arcpy.AddWarning("Unable to delete temporary GRASS folder. "
                             "Program will contine")
        cc_util.delete_feature(climate_asc)
        cc_util.delete_feature(resist_asc)
        cc_util.delete_feature(ccr_grassrc)


def write_grassrc(ccr_grassrc):
    """"Write GRASS resource file to project folder"""
    with open(ccr_grassrc, 'w') as f:
        f.write("GISDBASE: <UNKNOWN>\n")
        f.write("LOCATION_NAME: <UNKNOWN>\n")
        f.write("MAPSET: <UNKNOWN>\n")        

        
def setup_wrkspace(gisdbase, ccr_grassrc, geo_file):
    """Setup GRASS workspace and modify windows path for GRASS GDAL"""
    arcpy.AddMessage("Creating GRASS workspace")
    gisbase = cc_env.gisbase
    location = "gcwd"
    mapset = "PERMANENT"

    os.environ['GISRC'] = ccr_grassrc
    os.environ['LD_LIBRARY_PATH'] = os.path.join(gisbase, "lib")
    os.environ['GRASS_SH'] = os.path.join(gisbase, "msys", "bin", "sh.exe")

    env_list = os.environ['PATH'].split(';')
    env_list.insert(0, os.path.join(gisbase, "msys", "bin"))
    env_list.insert(0, os.path.join(gisbase, "extralib"))
    env_list.insert(0, os.path.join(gisbase, "bin"))
    env_list.insert(0, os.path.join(gisbase, "lib"))
    env_list.insert(0, os.path.join(gisbase, "etc", "python"))
    env_list.insert(0, os.path.join(gisbase, "etc"))
    os.environ['PATH'] = ';'.join(env_list)

    # Code to check GDAL dlls and system path
    # gdal = subprocess.Popen("where gdal*", stdout=subprocess.PIPE,
    #                         shell=True).stdout.read()
    # arcpy.AddMessage("GDAL DLL/s: " + gdal)

    # os_path = subprocess.Popen("echo %PATH%", stdout=subprocess.PIPE,
    #                         shell=True).stdout.read()
    # arcpy.AddMessage("Path: " + os_path)

    grass.create_location(gisdbase, location, filename=geo_file)
    gsetup.init(gisbase, gisdbase, location, mapset)
    run_grass_cmd("g.gisenv", set="OVERWRITE=1")
    os.environ['GRASS_VERBOSE'] = "0"  # only errors and warnings are printed
    return grass.version()['version']


def gen_cwd_back(grass_version, core_list, climate_lyr, resist_lyr, core_lyr):
    """"Generate CWD and back rasters using r.walk in GRASS"""
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

    try:
        for position, core_no in enumerate(core_list):
            core_no_txt = str(core_no)
            arcpy.AddMessage("Generating CWD and back rasters for"
                " Core " + core_no_txt + " (" + str(position + 1) + "/" +
                str(len(core_list)) + ")")

            # Extracting current core
            run_grass_cmd("v.extract", flags="t", input=core_lyr,
                          type="area", output=core,
                          where=cc_env.core_fld + " = " + core_no_txt)

            # Converting core vector to raster
            run_grass_cmd("v.to.rast", input=core, output=core_rast,
                          use="val")

            # Converting raster core to point feature
            if grass_version.startswith('7'):  # Command different in GRASS 7
                run_grass_cmd("r.to.vect", flags="z", input=core_rast,
                              output=core_points, type="point")
            else:
                run_grass_cmd("r.to.vect", flags="z", input=core_rast,
                              output=core_points, feature="point")

            # Running r.walk to create CWD and back raster
            run_grass_cmd("r.walk", elevation=climate_lyr,
                          friction=resist_lyr, output=gcwd, outdir=gback,
                          start_points=core_points, walk_coeff=walk_coeff,
                          slope_factor=slope_factor)

            # Reclassify from the directional degree output from GRASS to
            # Arc's 1 to 8 directions format
            rc_rules = os.path.join(cc_env.code_dir, "ccr_rcbkrast")
            run_grass_cmd("r.reclass", input=gback, output=gbackrc,
                          rules=rc_rules)

            # Exporting CWD and back rasters to ASCII grids
            cwd_ascii = os.path.join(cc_env.out_dir,
                                     "cwd_" + core_no_txt + ".asc")
            cwd_grid = lm_util.get_cwd_path(core_no)
            back_ascii = os.path.join(cc_env.out_dir,
                                      "back_" + core_no_txt + ".asc")
            back_grid = cwd_grid.replace("cwd_", "back_")

            run_grass_cmd("r.out.arc", input=gcwd, output=cwd_ascii)
            run_grass_cmd("r.out.arc", input=gbackrc, output=back_ascii)

            # Take grass cwd and back asciis and write them as ARCINFO grids
            arcpy.ASCIIToRaster_conversion(cwd_ascii, cwd_grid, "FLOAT")
            os.remove(cwd_ascii)
            arcpy.ASCIIToRaster_conversion(back_ascii, back_grid, "INTEGER")
            os.remove(back_ascii)
    except Exception:
        raise


def run_grass_cmd(*args, **kwargs):
    """ Run inputed GRASS command.

    Checks stderr for error and raises exception if it finds one.
    """
    kwargs['stdout'] = grass.PIPE
    kwargs['stderr'] = grass.PIPE
    ps = grass.start_command(*args, **kwargs)
    stderr = ps.communicate()[1]
    if 'ERROR:' in stderr:
        raise Exception("GRASSS ERROR: %s" % stderr[8:])
    # elif stderr:
        # arcpy.AddMessage(stderr)
