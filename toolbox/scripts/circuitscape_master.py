#!/usr/bin/env python2.5

"""Master script for circuitscape analysis in linkage mapper.

Reguired Software:
ArcGIS 9.3 with Spatial Analyst extension
Python 2.5
Numpy

"""

import os.path as path
import shutil
import sys

import arcgisscripting

from lm_config import tool_env as cfg
import lm_util as lu
import s7_centrality as s7
import s8_pinchpoints as s8


_filename = path.basename(__file__)


def circuitscape_master():
    """

    """
    cfg.configure("circuitscape", sys.argv)
    gp = cfg.gp
    gprint = lu.gprint
    try:
        lu.createfolder(cfg.LOGDIR)
        lu.createfolder(cfg.MESSAGEDIR)
        cfg.logFile = lu.create_log_file(cfg.MESSAGEDIR, cfg.TOOL, cfg.PARAMS)

        # Check core ID field.
        lu.check_cores(cfg.COREFC, cfg.COREFN)

        gp.OutputCoordinateSystem = gp.describe(cfg.COREFC).SpatialReference
        gp.pyramid = "NONE"
        gp.rasterstatistics = "NONE"

        # Move adj and cwd results from earlier versions to datapass directory
        lu.move_old_results()

        lu.delete_dir(cfg.SCRATCHDIR)

        if not cfg.DOPINCH and not cfg.DOCENTRALITY:
            msg = ('ERROR: Please choose at least one option: pinch point or\n'
                    'network centrality analysis.')
            lu.raise_error(msg)

        lu.createfolder(cfg.SCRATCHDIR)
        lu.createfolder(cfg.ARCSCRATCHDIR)

        if cfg.DO_ALLPAIRS:
            #  Fixme: move raster path to config
            S5CORRIDORRAS = path.join(cfg.OUTPUTGDB, cfg.PREFIX +
                                         "_lcc_mosaic_int")
            if not gp.Exists(S5CORRIDORRAS):
                msg = ('ERROR: Corridor raster created in step 5 is required'
                        '\nfor all-pair analyses, but was not found.')
                lu.raise_error(msg)
        if cfg.DOPINCH:
            # Make a local grid copy of resistance raster-
            # will run faster than gdb.
            lu.delete_data(cfg.RESRAST)
            if not gp.Exists(cfg.RESRAST_IN):
                msg = ('ERROR: Resistance raster is required for pinch point'
                        ' analyses, but was not found.')
                lu.raise_error(msg)

            arcpy.env.extent = cfg.RESRAST_IN
            arcpy.env.snapRaster = cfg.RESRAST_IN
            gprint('\nMaking local copy of resistance raster.')
            gp.CopyRaster_management(cfg.RESRAST_IN, cfg.RESRAST)

        if cfg.DOCENTRALITY:
            gprint("Creating output folder: " + cfg.CENTRALITYBASEDIR)
            if path.exists(cfg.CENTRALITYBASEDIR):
                shutil.rmtree(cfg.CENTRALITYBASEDIR)
            lu.createfolder(cfg.CENTRALITYBASEDIR)
            gp.CreateFolder_management(cfg.CENTRALITYBASEDIR,
                                        cfg.CIRCUITOUTPUTDIR_NM)
            gp.CreateFolder_management(cfg.CENTRALITYBASEDIR,
                                        cfg.CIRCUITCONFIGDIR_NM)
            lu.clean_out_workspace(cfg.CORECENTRALITYGDB)

            s7.STEP7_calc_centrality()
            if not cfg.SAVECENTRALITYDIR:
                lu.delete_dir(cfg.CENTRALITYBASEDIR)

        if cfg.DOPINCH:
            gprint("Creating output folder: " + cfg.CIRCUITBASEDIR)
            lu.delete_dir(cfg.CIRCUITBASEDIR)
            lu.createfolder(cfg.CIRCUITBASEDIR)
            gp.CreateFolder_management(cfg.CIRCUITBASEDIR,
                                        cfg.CIRCUITOUTPUTDIR_NM)
            gp.CreateFolder_management(cfg.CIRCUITBASEDIR,
                                        cfg.CIRCUITCONFIGDIR_NM)

            lu.clean_out_workspace(cfg.PINCHGDB)
            lu.delete_data(cfg.PINCHGDB)

            s8.STEP8_calc_pinchpoints()

            lu.delete_dir(cfg.SCRATCHDIR)
            if not cfg.SAVECIRCUITDIR:
                lu.delete_dir(cfg.CIRCUITBASEDIR)

        gprint('\nDONE!\n')

    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.print_geoproc_error(_filename)

    # Return any PYTHON or system specific errors
    except:
        lu.print_python_error(_filename)


if __name__ == "__main__":
    circuitscape_master()
