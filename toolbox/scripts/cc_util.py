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


def arc_delete(*items):
    """Delete data from disk using Arcpy."""
    for item in items:
        if arcpy.Exists(item):
            try:
                arcpy.Delete_management(item)
            except arcpy.ExecuteError:
                arcpy.AddWarning("Error deleting temporary object - %s. "
                                 "Program will continue." % item)


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
