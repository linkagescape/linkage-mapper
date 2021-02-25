#!/usr/bin/env python2
# Author: Brad McRae

"""Master script for circuitscape analysis in linkage mapper.

Reguired Software:
ArcGIS Desktop 10.3+ or ArcGIS Pro with Spatial Analyst extension
Numpy

"""

import os
from os import path
import shutil
import sys

import arcpy

from lm_config import tool_env as cfg
import lm_util as lu

try:
    import s7_centrality as s7
except Exception:
    pass
import s8_pinchpoints as s8

_SCRIPT_NAME = "circuitscape_master.py"


def circuitscape_master(argv=None):
    """

    """
    gprint = lu.gprint

    if argv is None:
        argv = sys.argv

    argv.append(get_cs_path())  # Add Circuitscape path

    cfg.configure(cfg.TOOL_CS, argv)

    try:
        lu.create_dir(cfg.LOGDIR)
        lu.create_dir(cfg.MESSAGEDIR)
        cfg.logFilePath = lu.create_log_file(cfg.MESSAGEDIR, cfg.TOOL,
                                           cfg.PARAMS)

        if cfg.CSPATH is None:
            lu.raise_error("Cannot find an installation of Circuitscape"
                           "\nin your Program Files directory.")

        lu.print_drive_warning()
        # Check core ID field.
        lu.check_cores(cfg.COREFC, cfg.COREFN)

        arcpy.env.outputCoordinateSystem = arcpy.Describe(cfg.COREFC).SpatialReference
        # Set data frame spatial reference to coordinate system of input data
        lu.set_dataframe_sr()

        arcpy.env.pyramid = "NONE"
        arcpy.env.rasterStatistics = "NONE"

        # Move adj and cwd results from earlier versions to datapass directory
        lu.move_old_results()

        if cfg.CWDCUTOFF > 0:
            lu.delete_dir(cfg.SCRATCHDIR)

        # restart code- in progress
        if cfg.CWDCUTOFF < 0:
            cfg.CWDCUTOFF = cfg.CWDCUTOFF * -1

        if not cfg.DOPINCH and not cfg.DOCENTRALITY:
            msg = ('ERROR: Please choose at least one option: pinch point or\n'
                    'network centrality analysis.')
            lu.raise_error(msg)

        lu.create_dir(cfg.SCRATCHDIR)
        lu.create_dir(cfg.ARCSCRATCHDIR)

        if cfg.DO_ALLPAIRS:
            #  Fixme: move raster path to config
            S5CORRIDORRAS = path.join(cfg.OUTPUTGDB,cfg.PREFIX + "_corridors")
            if not arcpy.Exists(S5CORRIDORRAS):
                S5CORRIDORRAS = path.join(cfg.OUTPUTGDB, cfg.PREFIX +
                                         "_lcc_mosaic_int")
            if not arcpy.Exists(S5CORRIDORRAS):
                msg = ('ERROR: Corridor raster created in step 5 is required'
                        '\nfor all-pair analyses, but was not found.')
                lu.raise_error(msg)
        if cfg.DOPINCH:
            if cfg.CWDCUTOFF == '#' or cfg.CWDCUTOFF == 0:
                msg = ('ERROR: CWD cutoff distance is required for pinch point'
                        ' analyses.')
                lu.raise_error(msg)

            # Make a local grid copy of resistance raster-
            # will run faster than gdb.
            lu.delete_data(cfg.RESRAST)
            if not arcpy.Exists(cfg.RESRAST_IN):
                msg = ('ERROR: Resistance raster is required for pinch point'
                        ' analyses, but was not found.')
                lu.raise_error(msg)

            arcpy.env.extent = cfg.RESRAST_IN
            desc = arcpy.Describe(cfg.RESRAST_IN)
            if hasattr(desc, "catalogPath"):
                cfg.RESRAST_IN = arcpy.Describe(cfg.RESRAST_IN).catalogPath

            gprint('\nMaking local copy of resistance raster.')
            try:
                arcpy.CopyRaster_management(cfg.RESRAST_IN, cfg.RESRAST)
            except Exception:
                msg = ('ERROR: Could not make a copy of your resistance raster. ' +
                    'Try re-starting ArcMap to release the file lock.')
                lu.raise_error(msg)

            arcpy.env.snapRaster = cfg.RESRAST

        if cfg.DOCENTRALITY:
            gprint("Creating output folder: " + cfg.CENTRALITYBASEDIR)
            if path.exists(cfg.CENTRALITYBASEDIR):
                shutil.rmtree(cfg.CENTRALITYBASEDIR)
            lu.create_dir(cfg.CENTRALITYBASEDIR)
            arcpy.CreateFolder_management(cfg.CENTRALITYBASEDIR,
                                        cfg.CIRCUITOUTPUTDIR_NM)
            arcpy.CreateFolder_management(cfg.CENTRALITYBASEDIR,
                                        cfg.CIRCUITCONFIGDIR_NM)
            lu.clean_out_workspace(cfg.CORECENTRALITYGDB)

            s7.STEP7_calc_centrality()
            if not cfg.SAVECENTRALITYDIR:
                lu.delete_dir(cfg.CENTRALITYBASEDIR)

        if cfg.DOPINCH:
            if cfg.CWDCUTOFF > 0: # Negative values mean we're restarting
                gprint("Creating output folder: " + cfg.CIRCUITBASEDIR)
                lu.delete_dir(cfg.CIRCUITBASEDIR)
                lu.create_dir(cfg.CIRCUITBASEDIR)
                arcpy.CreateFolder_management(cfg.CIRCUITBASEDIR,
                                        cfg.CIRCUITOUTPUTDIR_NM)
                arcpy.CreateFolder_management(cfg.CIRCUITBASEDIR,
                                        cfg.CIRCUITCONFIGDIR_NM)

            s8.STEP8_calc_pinchpoints()

            if not cfg.SAVE_TEMP_CIRCUIT_FILES:
                lu.delete_dir(cfg.SCRATCHDIR)
            if not cfg.SAVECIRCUITDIR:
                lu.delete_dir(cfg.CIRCUITBASEDIR)

        gprint('\nDONE!\n')

    # Return GEOPROCESSING specific errors
    except arcpy.ExecuteError:
        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except Exception:
        lu.exit_with_python_error(_SCRIPT_NAME)


def get_cs_path():
    """Return path to Circuitscape installation."""
    env_list = ["ProgramW6432", "ProgramFiles", "ProgramFiles(x86)"]
    cs_app_path = None

    for env_path in [os.environ[key] for key in env_list if key in os.environ]:
        chk_path = path.join(env_path, "Circuitscape\\cs_run.exe")
        if path.exists(chk_path):
            cs_app_path = chk_path
            break

    return cs_app_path


if __name__ == "__main__":
    circuitscape_master()
