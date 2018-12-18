"""Utility functions for configuration modules."""

import os.path as path
import imp


GP_NULL = '#'


def get_code_path():
    """Get full path to tool scripts folder."""
    return path.dirname(path.realpath(__file__))


def set_custom(settings, config_obj):
    """Add attributes for custom file to configuration object."""
    cust_settings = imp.load_source("set_mod", settings)
    for setting in dir(cust_settings):
        if setting == setting.upper():
            setting_value = getattr(cust_settings, setting)
            setattr(config_obj, setting, setting_value)


def nullstring(arg_string):
    """Convert ESRI nullstring to Python null."""
    if arg_string == GP_NULL:
        arg_string = None
    return arg_string


def str2bool(pstr):
    """Convert ESRI boolean string to Python boolean type."""
    return pstr == "true"
