#!/usr/bin/env python2.6
# Author: Darren Kavanagh

"""Climate Linkage Mapper utility module"""

import os
import shutil
import arcpy

from cc_config import cc_env


def mk_proj_dir(in_dir):
    """Create directory off project folder if it does not already exist"""
    new_dir = os.path.join(cc_env.proj_dir, in_dir)
    if not os.path.exists(new_dir):
        os.mkdir(new_dir)
    return new_dir


def delete_features(*in_features):
    """Delete temporary feature/s if it exists"""
    for in_feature in in_features:
        if arcpy.Exists(in_feature):
            try:
                arcpy.Delete_management(in_feature)
            except arcpy.ExecuteError:
                arcpy.AddWarning("Error deleting temporary %s. Program will "
                                 "continue." % in_feature)


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
