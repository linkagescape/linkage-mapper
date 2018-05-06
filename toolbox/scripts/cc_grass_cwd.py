#!/usr/bin/env python2.6
# Authors: Darren Kavanagh and Brad McRae

"""Climate Linkage Mapper grass module.

"""

import os

import arcpy

import grass.script as grass
import grass.script.setup as gsetup

from cc_config import cc_env
import cc_util
import lm_util


def grass_cwd(core_list):
    """Creating CWD and Back rasters using GRASS r.walk function"""
    out_fldr = cc_env.scratch_dir

    gisdbase = os.path.join(out_fldr, "cc_grass")

    ccr_grassrc = os.path.join(out_fldr, "cc_grassrc")
    climate_asc = os.path.join(out_fldr, "cc_gclimate.asc")
    resist_asc = os.path.join(out_fldr, "cc_gresist.asc")
    core_asc = os.path.join(out_fldr, "cc_gcores.asc")
    climate_lyr = "climate"
    resist_lyr = "resist"
    core_lyr = "cores"

    try:
        lm_util.gprint("\nRUNNING GRASS TO CREATE COST-WEIGHTED DISTANCE "
                       "RASTERS")

        start_path = os.environ["PATH"]

        # Convert input GRID rasters to ASCII
        lm_util.gprint("Converting ARCINFO GRID rasters to ASCII")
        arcpy.RasterToASCII_conversion(cc_env.prj_climate_rast, climate_asc)
        arcpy.RasterToASCII_conversion(cc_env.prj_resist_rast, resist_asc)
        arcpy.RasterToASCII_conversion(cc_env.prj_core_rast, core_asc)

        # Create resource file and setup workspace
        os.environ["PATH"] = cc_env.gpath
        write_grassrc(ccr_grassrc, gisdbase)
        setup_wrkspace(gisdbase, ccr_grassrc, climate_asc)

        # Make cwd folder for Linkage Mapper
        lm_util.make_cwd_paths(max(core_list))

        # Import files into GRASS
        lm_util.gprint("Importing raster files into GRASS")
        run_grass_cmd("r.in.gdal", input=climate_asc, output=climate_lyr)
        run_grass_cmd("r.in.gdal", input=resist_asc, output=resist_lyr)
        run_grass_cmd("r.in.gdal", input=core_asc, output=core_lyr)

        # Generate CWD and Back rasters
        gen_cwd_back(core_list, climate_lyr, resist_lyr, core_lyr)

    except Exception:
        raise
    finally:
        os.environ["PATH"] = start_path
        cc_util.arc_delete(gisdbase, ccr_grassrc,
                           climate_asc, resist_asc, core_asc)


def write_grassrc(ccr_grassrc, gisdbase):
    """"Write GRASS resource file to project folder"""
    with open(ccr_grassrc, 'w') as rc_file:
        rc_file.write("GISDBASE: %s\n" % gisdbase)
        rc_file.write("LOCATION_NAME: <UNKNOWN>\n")
        rc_file.write("MAPSET: <UNKNOWN>\n")


def setup_wrkspace(gisdbase, ccr_grassrc, geo_file):
    """Setup GRASS workspace"""
    lm_util.gprint("Creating GRASS workspace")
    gisbase = cc_env.gisbase
    location = "gcwd"
    mapset = "PERMANENT"

    os.environ['GISRC'] = ccr_grassrc
    os.environ['LD_LIBRARY_PATH'] = os.path.join(gisbase, "lib")
    os.environ['GRASS_SH'] = os.path.join(gisbase, "msys", "bin", "sh.exe")

    try:
        grass.create_location(gisdbase, location, filename=geo_file)
    except:
        warn_msg = ("Cannot create GRASS workspace.\n"
                    "Try rebooting and restarting ArcGIS. If that doesn't\n"
                    "work you can try using the 'cc_demo.py' python script\n"
                    "in the demo/demo_scripts/ directory where the Linkage\n"
                    "Mapper toolbox is installed instead of ArcGIS to call\n"
                    "the tool (see user guide).")
        raise Exception(warn_msg)

    gsetup.init(gisbase, gisdbase, location, mapset)
    run_grass_cmd("g.gisenv", set="OVERWRITE=1")
    os.environ['GRASS_VERBOSE'] = "0"  # Only errors and warnings are printed


def gen_cwd_back(core_list, climate_lyr, resist_lyr, core_lyr):
    """"Generate CWD and back rasters using r.walk in GRASS"""
    slope_factor = "1"
    walk_coeff_flat = "1"
    walk_coeff_uphill = str(cc_env.climate_cost)
    walk_coeff_downhill = str(cc_env.climate_cost * -1)
    walk_coeff = (walk_coeff_flat + "," + walk_coeff_uphill + "," +
                  walk_coeff_downhill + "," + walk_coeff_downhill)

    focal_core_rast = "focal_core_rast"
    gcwd = "gcwd"
    gback = "gback"
    gbackrc = "gbackrc"
    core_points = "corepoints"
    no_cores = str(len(core_list))

    # Map from directional degree output from GRASS to Arc's 1 to 8 directions
    # format. See r.walk source code and ArcGIS's 'Understanding cost distance
    # analysis' help page.
    rc_rules = "180=5\n225=4\n270=3\n315=2\n360=1\n45=8\n90=7\n135=6"

    try:
        for position, core_no in enumerate(core_list):
            core_no_txt = str(core_no)
            lm_util.gprint("Generating CWD and back rasters for Core " +
                           core_no_txt + " (" + str(position + 1) + "/" +
                           no_cores + ")")

            # Pull out focal core for cwd analysis
            write_grass_cmd("r.reclass", input=core_lyr,
                            output=focal_core_rast, overwrite=True,
                            rules="-", stdin=core_no_txt + '=' + core_no_txt)

            # Converting raster core to point feature
            run_grass_cmd("r.to.vect", flags="z", input=focal_core_rast,
                          output=core_points, type="point")

            # Running r.walk to create CWD and back raster
            run_grass_cmd("r.walk", elevation=climate_lyr,
                          friction=resist_lyr, output=gcwd, outdir=gback,
                          start_points=core_points, walk_coeff=walk_coeff,
                          slope_factor=slope_factor)

            # Reclassify back raster directional degree output to ArcGIS format
            write_grass_cmd("r.reclass", input=gback, output=gbackrc,
                            rules="-", stdin=rc_rules)

            # Get spatial reference for defining ARCINFO raster projections
            desc_data = arcpy.Describe(cc_env.prj_core_rast)
            spatial_ref = desc_data.spatialReference

            # Get cwd path (e.g. ..\datapass\cwd\cw\cwd_3)
            cwd_path = lm_util.get_cwd_path(core_no)

            def create_arcgrid(rtype, grass_grid):
                """Export GRASS raster to ASCII grid and then to ARCINFO grid
                """
                ascii_grid = os.path.join(cc_env.scratch_dir,
                                          rtype + core_no_txt + ".asc")
                arc_grid = cwd_path.replace("cwd_", rtype)
                run_grass_cmd("r.out.gdal", input=grass_grid,
                              output=ascii_grid, format="AAIGrid")
                arcpy.CopyRaster_management(ascii_grid, arc_grid)
                arcpy.DefineProjection_management(arc_grid, spatial_ref)
                cc_util.arc_delete(ascii_grid)

            create_arcgrid("cwd_", gcwd)  # Export CWD raster
            create_arcgrid("back_", gbackrc)  # Export reclassified back raster
    except Exception:
        raise


def start_grass_cmd(*args, **kwargs):
    """Calls the grass module's start_command to run the inputed GRASS command.

    Returns Popen object
    """
    kwargs['stdout'] = grass.PIPE
    kwargs['stderr'] = grass.PIPE
    return grass.start_command(*args, **kwargs)


def chk_stderr(stderr):
    """Check process for error and raises exception if one is found"""
    if 'ERROR:' in stderr:
        raise Exception("GRASSS ERROR: %s" % stderr[7:])


def run_grass_cmd(*args, **kwargs):
    """Run inputed GRASS command"""
    return_ps = start_grass_cmd(*args, **kwargs)
    chk_stderr(return_ps.communicate()[1])


def write_grass_cmd(*args, **kwargs):
    """Feeds stdin string to process stdin and runs inputed GRASS command"""
    stdin = kwargs['stdin']
    kwargs['stdin'] = grass.PIPE
    return_ps = start_grass_cmd(*args, **kwargs)
    chk_stderr(return_ps.communicate(input=stdin)[1])
