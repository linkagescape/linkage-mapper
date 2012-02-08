#!/usr/bin/env python2.5

"""Master script for barrier analysis in linkage mapper.

Reguired Software:
ArcGIS 9.3 with Spatial Analyst extension
Python 2.5
Numpy

"""

import os.path as path
import sys

import arcgisscripting
import arcpy

from lm_config import tool_env as cfg
import lm_util as lu
import s6_barriers as s6


_filename = path.basename(__file__)


def bar_master():
    """ Experimental code to detect barriers using cost-weighted distance
    outputs from Linkage Mapper tool.

    """
    cfg.configure("barrier_mapper", sys.argv)
    gprint = lu.gprint
    try:
        lu.createfolder(cfg.LOGDIR)
        lu.createfolder(cfg.MESSAGEDIR)

        cfg.logFile = lu.create_log_file(cfg.MESSAGEDIR, cfg.TOOL, cfg.PARAMS)

        # Move adj and cwd results from earlier versions to datapass directory
        lu.move_old_results()

        # Delete final ouptut geodatabase
        lu.delete_dir(cfg.BARRIERGDB)
        if not arcpy.Exists(cfg.BARRIERGDB):
            # Create output geodatabase
            arcpy.CreateFileGDB_management(cfg.OUTPUTDIR,
                                           path.basename(cfg.BARRIERGDB))

        lu.createfolder(cfg.OUTPUTDIR)
        lu.delete_dir(cfg.SCRATCHDIR)
        lu.createfolder(cfg.SCRATCHDIR)
        lu.createfolder(cfg.ARCSCRATCHDIR)

        arcpy.env.extent = cfg.RESRAST_IN
        arcpy.env.snapRaster = cfg.RESRAST_IN

        gprint('\nMaking local copy of resistance raster.')
        lu.delete_data(cfg.RESRAST)
        arcpy.CopyRaster_management(cfg.RESRAST_IN, cfg.RESRAST)

        s6.STEP6_calc_barriers()

        #clean up
        lu.delete_dir(cfg.SCRATCHDIR)
        if not cfg.SAVEBARRIERDIR:
            lu.delete_dir(cfg.BARRIERBASEDIR)
        gprint('\nDONE!\n')

    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.print_geoproc_error(_filename)

    # Return any PYTHON or system specific errors
    except:
        lu.print_python_error(_filename)

if __name__ == "__main__":
    bar_master()



        # desc = arcpy.Describe(cfg.RESRAST_IN)

        # if hasattr(desc, "name"):
            # gprint ("Name:        " + desc.name)

        # if hasattr(desc, "catalogPath"):
            # gprint ("CatalogPath: " + desc.catalogPath)
