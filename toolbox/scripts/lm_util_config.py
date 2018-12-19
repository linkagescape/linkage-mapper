"""Utility functions for configuration modules."""

GP_NULL = '#'


def nullstring(arg_string):
    """Convert ESRI nullstring to Python null."""
    if arg_string == GP_NULL:
        arg_string = None
    return arg_string


def str2bool(pstr):
    """Convert ESRI boolean string to Python boolean type."""
    return pstr == "true"
