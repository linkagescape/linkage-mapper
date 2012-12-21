#!/usr/bin/env python2.6
# Author: Darren Kavanagh

"""Climate Linkage Mapper grass module.

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
        lm_util.gprint("\nRUNNING GRASS TO CREATE COST-WEIGHTED DISTANCE "
                         "RASTERS")

        # Convert input GRID rasters to ASCII
        lm_util.gprint("Converting ARCINFO GRID rasters to ASCII")
        arcpy.RasterToASCII_conversion(cc_env.prj_climate_rast, climate_asc)
        arcpy.RasterToASCII_conversion(cc_env.prj_resist_rast, resist_asc)
        # Create resource file and setup workspace
        write_grassrc(ccr_grassrc, gisdbase)
        
        grass_version = setup_wrkspace(gisdbase, ccr_grassrc, climate_asc)

        # Make cwd folder for Linkage Mapper
        lm_util.make_cwd_paths(max(core_list))

        # Import files into GRASS
        lm_util.gprint("Importing raster files into GRASS")
        run_grass_cmd("r.in.arc", input=climate_asc, output=climate_lyr)
        run_grass_cmd("r.in.arc", input=resist_asc, output=resist_lyr)
        lm_util.gprint("Importing cores feature class into GRASS\n")
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


def write_grassrc(ccr_grassrc, gisdbase):
    """"Write GRASS resource file to project folder"""
    with open(ccr_grassrc, 'w') as f:
        f.write("GISDBASE: %s\n" % gisdbase)
        f.write("LOCATION_NAME: <UNKNOWN>\n")
        f.write("MAPSET: <UNKNOWN>\n")


def setup_wrkspace(gisdbase, ccr_grassrc, geo_file):
    """Setup GRASS workspace and modify windows path for GRASS GDAL"""
    lm_util.gprint("Creating GRASS workspace")
    gisbase = cc_env.gisbase
    location = "gcwd"
    mapset = "PERMANENT"

    os.environ['GISRC'] = ccr_grassrc
    os.environ['LD_LIBRARY_PATH'] = os.path.join(gisbase, "lib")
    os.environ['GRASS_SH'] = os.path.join(gisbase, "msys", "bin", "sh.exe")

    cc_util.add_grass_path(gisbase)
    # os_path = subprocess.Popen("echo %PATH%", stdout=subprocess.PIPE,
    #                         shell=True).stdout.read()
    # lm_util.gprint("Path: " + os_path)
    if os.path.exists(gisdbase):
        try:
            shutil.rmtree(gisdbase, True)
        except:
            raise Exception("Cannot delete grass workspace: " + gisdbase)  
        if os.path.exists(gisdbase):
            raise Exception("Cannot delete grass workspace: " + gisdbase)    
    try:
        grass.create_location(gisdbase, location, filename=geo_file)
    except:
        gdal_fail_check('create location failure')
        arcpy.AddWarning("GRASS ERROR. Try rebooting and restarting ArcGIS.")
        arcpy.AddWarning("If that doesn't work you can try using ")
        arcpy.AddWarning("the 'CC Run Script.py' python script in the ")
        arcpy.AddWarning("scripts directory where the Linkage Mapper toolbox")
        arcpy.AddWarning("is installed instead of ArcGIS to call the tool")
        arcpy.AddWarning("(see user guide).")
        raise Exception("GRASS ERROR: Cannot create workspace.")    
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
            core='core' + core_no_txt
            lm_util.gprint("Generating CWD and back rasters for"
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
                                     "cwd_" + core_no_txt + ".asc") #xxxx
            cwd_grid = lm_util.get_cwd_path(core_no)
            back_ascii = os.path.join(cc_env.out_dir,
                                      "back_" + core_no_txt + ".asc") #xxxx
            back_grid = cwd_grid.replace("cwd_", "back_")                   
  
            run_grass_cmd("r.out.arc", input=gcwd, output=cwd_ascii)
            run_grass_cmd("r.out.arc", input=gbackrc, output=back_ascii)

            # # Take grass cwd and back asciis and write them as ARCINFO grids  
            # ASCIIToRaster INVOKES GDAL- replaced with copyraster
            # arcpy.ASCIIToRaster_conversion(cwd_ascii, 
                                           # cwd_grid, "FLOAT") # 
            arcpy.CopyRaster_management(cwd_ascii,cwd_grid)
            os.remove(cwd_ascii)
            # arcpy.ASCIIToRaster_conversion(back_ascii, back_grid, "INTEGER")
            arcpy.CopyRaster_management(back_ascii,back_grid)
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
        raise Exception("GRASSS ERROR: %s" % stderr[7:])
    # elif stderr:
        # lm_util.gprint(stderr)

def gdal_check(msg):
# Code to check GDAL dlls and system path
    gdal = subprocess.Popen("where gdal*", stdout=subprocess.PIPE,
                            shell=True).stdout.read()
    lm_util.gprint("\nGDAL DLL/s at " + msg + ': ' + gdal)
    
def gdal_fail_check(msg):
# Code to check GDAL dlls and system path
    gdal = subprocess.Popen("where gdal*", stdout=subprocess.PIPE,
                            shell=True).stdout.read()
    gdalList = gdal.split('\n')                    
    if 'arcgis' in gdalList[1].lower():
        arcpy.AddWarning("It looks like there is a conflict between ArcGIS")
        arcpy.AddWarning("and GRASS. \nPlease *RESTART ArcGIS* and try again.")
        arcpy.AddWarning("\nIf that doesn't work you can try using ")
        arcpy.AddWarning("the 'CC Run Script.py' python script in the ")
        arcpy.AddWarning("scripts directory where the Linkage Mapper toolbox")
        arcpy.AddWarning("is installed instead of ArcGIS to call the tool")        
        raise Exception("GDAL DLL conflict")    

