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
