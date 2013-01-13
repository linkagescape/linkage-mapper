#!/usr/bin/env python2.6
# Author: Darren Kavanagh

"""Climate Linkage Mapper utility module"""

import os

import arcpy

from cc_config import cc_env


def mk_proj_dir(in_dir):
    """Create directory off project folder if it does not already exist"""
    new_dir = os.path.join(cc_env.proj_dir, in_dir)
    if not os.path.exists(new_dir):
        os.mkdir(new_dir)
    return new_dir


def delete_feature(in_feature):
    """Delete temporary feature if already exists"""
    if arcpy.Exists(in_feature):
        try:
            arcpy.Delete_management(in_feature)
        except arcpy.ExecuteError:
            arcpy.AddWarning("Error deleting temporary %s. Program will "
                             "continue." % in_feature)

def add_grass_path(gisbase):
    env_list = os.environ['PATH'].split(';')
    env_list.insert(0, os.path.join(gisbase, "msys", "bin"))
    env_list.insert(0, os.path.join(gisbase, "extralib"))
    env_list.insert(0, os.path.join(gisbase, "bin"))
    env_list.insert(0, os.path.join(gisbase, "lib"))
    env_list.insert(0, os.path.join(gisbase, "etc", "python"))
    env_list.insert(0, os.path.join(gisbase, "etc"))
    os.environ['PATH'] = ';'.join(env_list)  # Path should now prioritize 
                                             # proper gdal dll location
                             
def check_cc_project_dir():
    """Checks to make sure path name is not too long.

    Long path names can cause problems with ESRI grids.
    """
    if len(cc_env.proj_dir) > 140:
        msg = ('ERROR: Project directory "' + cc_env.proj_dir +
               '" is too deep.  Please choose a shallow directory'
               '(something like "C:\PUMA").')
        raise Exception(msg)

    if "-" in cc_env.proj_dir or " " in cc_env.proj_dir or "." in cc_env.proj_dir:
        msg = ('ERROR: Project directory cannot contain spaces, dashes, or '
                'special characters.')
        raise Exception(msg)
    return
    
    
def remove_grass_wkspc(gisdbase):
        # Remove GRASS workspace from earlier run if any
        if os.path.exists(gisdbase):
            import shutil
            try:
                shutil.rmtree(gisdbase, True)
            except:
                arcpy.AddWarning("\nCannot delete GRASS workspace from earlier"
                                " run."
                                "\nPlease choose a new project directory.")
                raise Exception("Cannot delete GRASS workspace: " + gisdbase)              
            if os.path.exists(gisdbase):
                arcpy.AddWarning("\nCannot delete GRASS workspace from earlier"
                                " run."
                                "\nPlease choose a new project directory.")
                raise Exception("Cannot delete GRASS workspace: " + gisdbase)  
    
                             