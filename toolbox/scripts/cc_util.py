#!/usr/bin/env python2.6
# Author: Darren Kavanagh

"""Climate Linkage Mapper utility module"""

import os
import subprocess
import shutil

import arcpy

from cc_config import cc_env
import lm_util


def mk_proj_dir(in_dir):
    """Create directory off project folder if it does not already exist"""
    new_dir = os.path.join(cc_env.proj_dir, in_dir)
    if not os.path.exists(new_dir):
        os.mkdir(new_dir)
    return new_dir


def delete_features(in_features):
    """Delete temporary feature/s if it exists"""
    for in_feature in in_features:
        if arcpy.Exists(in_feature):
            try:
                arcpy.Delete_management(in_feature)
            except arcpy.ExecuteError:
                arcpy.AddWarning("Error deleting temporary %s. Program will "
                                 "continue." % in_feature)


def add_grass_path(gisbase):
    """Add GRASS directories to system PATH"""
    env_list = os.environ['PATH'].split(';')
    env_list.insert(0, os.path.join(gisbase, "msys", "bin"))
    env_list.insert(0, os.path.join(gisbase, "extrabin"))
    env_list.insert(0, os.path.join(gisbase, "bin"))
    env_list.insert(0, os.path.join(gisbase, "lib"))
    env_list.insert(0, os.path.join(gisbase, "etc", "python"))
    env_list.insert(0, os.path.join(gisbase, "etc"))
    # Path should now prioritize proper gdal dll location
    os.environ['PATH'] = ';'.join(env_list)


def check_cc_project_dir():
    """Checks to make sure path name is not too long.

    Long path names can cause problems with ESRI grids.
    """
    if len(cc_env.proj_dir) > 140:
        msg = ('ERROR: Project directory "' + cc_env.proj_dir +
               '" is too deep.  Please choose a shallow directory'
               r'(something like "C:\PUMA").')
        raise Exception(msg)

    if ("-" in cc_env.proj_dir or " " in cc_env.proj_dir or
            "." in cc_env.proj_dir):
        msg = ('ERROR: Project directory cannot contain spaces, dashes, or '
               'special characters.')
        raise Exception(msg)
    return


def remove_grass_wkspc(gisdbase):
    """Remove GRASS workspace if it exists"""
    if os.path.exists(gisdbase):
        try:
            shutil.rmtree(gisdbase, True)
        except OSError:
            return False
    return True


def gdal_fail_check():
    """Code to check GDAL dlls and system path"""
    gdal = subprocess.Popen("where gdal*", stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell=True).stdout.read()

    if gdal != '':
        gdal_list = gdal.split('\n')
        if 'arcgis' in gdal_list[1].lower():
            lm_util.gprint("\nGDAL DLL/s found at: " + gdal)
            arcpy.AddWarning("It looks like there is a conflict between "
                             "ArcGIS")
            arcpy.AddWarning("and GRASS. This could be the result of a "
                             "previous ")
            arcpy.AddWarning("analysis (like a Linkage Mapper run) or it might"
                             "be")
            arcpy.AddWarning("caused by conflicts with pre-loaded ArcGIS ")
            arcpy.AddWarning("extensions.")
            arcpy.AddWarning("\nThis error often goes away if you run the tool"
                             "in")
            arcpy.AddWarning("the background (see user guide). ")
            arcpy.AddWarning("\nIf that doesn't work, try restarting ArcMap.")
            arcpy.AddWarning("\nIf you still get an error, then restart again")
            arcpy.AddWarning(" and disable any extensions you are not using")
            arcpy.AddWarning("(Click on Customize >> Extensions).")
            arcpy.AddWarning("\nAnd if that doesn't work try closing Arc and ")
            arcpy.AddWarning("instead run the tool using the "
                             "'CC Run Script.py' ")
            arcpy.AddWarning("python script.  This script can be found in "
                             "the ")
            arcpy.AddWarning("'demo' directory, located where the Linkage")
            arcpy.AddWarning("Mapper toolbox is installed.\n")
            raise Exception("ArcGIS-GRASS GDAL DLL conflict")
    else:
        lm_util.gprint("GDAL DLL/s not found in system PATH")
        arcpy.AddWarning("Check if the appropriate version of GRASS is "
                         "correctly installed.")
        raise Exception("GRASS DLL/s not found")
