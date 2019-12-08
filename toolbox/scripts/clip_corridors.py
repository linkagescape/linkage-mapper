#!/usr/bin/env python2
# Author: Brad McRae

from os import path
import sys

import traceback

import arcpy


arcpy.CheckOutExtension("spatial")

gprint = arcpy.AddMessage

_SCRIPT_NAME = "clip_corridors"


def clip_corridor():
    """Truncates corridors at user-specified cutoff width in CWD units.

    """
    try:
        inRaster = sys.argv[1]
        cutoffVal = sys.argv[2]
        outputGDB = sys.argv[3]

        cutoffText = str(cutoffVal)
        if cutoffText[-6:] == '000000':
            cutoffText = cutoffText[0:-6]+'m'
        elif cutoffText[-3:] == '000':
            cutoffText = cutoffText[0:-3]+'k'

        inPath,FN = path.split(inRaster) # In case raster is in a group layer
        outRasterFN = FN + '_truncated_' + cutoffText
        outRaster = path.join(outputGDB,outRasterFN)
        delete_data(outRaster)

        desc = arcpy.Describe(inRaster)
        if hasattr(desc, "catalogPath"):
            inRaster = arcpy.Describe(inRaster).catalogPath
        arcpy.env.overwriteOutput = True
        arcpy.env.workspace = outputGDB
        arcpy.env.scratchWorkspace = outputGDB
        arcpy.env.extent = inRaster
        arcpy.env.cellSize = inRaster
        output = arcpy.sa.Con(
            arcpy.sa.Raster(inRaster) <= float(cutoffVal), inRaster)
        output.save(outRaster)

        gprint('Building output statistics for truncated raster')
        build_stats(outRaster)

        gprint('\nThe new truncated corridor raster "'+ outRasterFN + '" can '
                'be found in the corridor geodatabase:')
        gprint(outputGDB)

    # Return GEOPROCESSING specific errors
    except arcpy.ExecuteError:
        exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)

    return

def exit_with_geoproc_error(filename):
    """Handle geoprocessor errors and provide details to user"""
    tb = sys.exc_info()[2]  # get the traceback object
    # tbinfo contains the error's line number and the code
    tbinfo = traceback.format_tb(tb)[0]
    line = tbinfo.split(", ")[1]

    arcpy.AddError("Geoprocessing error on **" + line + "** of " + filename +
                " :")
    for msg in range(0, arcpy.GetMessageCount() - 1):
        if arcpy.GetSeverity(msg) == 2:
            arcpy.AddReturnMessage(msg)
        print(arcpy.AddReturnMessage(msg))
    exit(0)

def exit_with_python_error(filename):
    """Handle python errors and provide details to user"""
    tb = sys.exc_info()[2]  # get the traceback object
    # tbinfo contains the error's line number and the code
    tbinfo = traceback.format_tb(tb)[0]
    line = tbinfo.split(", ")[1]
    err = traceback.format_exc().splitlines()[-1]
    arcpy.AddError("Python error on **" + line + "** of " + filename)
    arcpy.AddError(err)
    exit(0)

def delete_data(dataset):
    try:
        if arcpy.Exists(dataset):
            arcpy.Delete_management(dataset)
    except Exception:
        pass


def build_stats(raster):
    """Builds statistics and pyramids for output rasters"""
    try:
        arcpy.CalculateStatistics_management(raster, "1", "1", "#")
    except Exception:
        gprint('Statistics failed. They can still be calculated manually.')
    try:
        arcpy.BuildPyramids_management(raster)
    except Exception:
        gprint('Pyramids failed. They can still be built manually.')
    return


if __name__ == "__main__":
    clip_corridor()
