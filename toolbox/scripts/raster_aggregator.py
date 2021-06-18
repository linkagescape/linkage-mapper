"""Aggregates resistance rasters to coarser cell size using average resistance.
    Written by Brad McRae
Reguired Software:
ArcGIS Desktop 10.3+ or ArcGIS Pro with Spatial Analyst extension

"""

__filename__ = "raster_aggregator.py"
__version__ = "2013_0610"

from os import path
import sys
import traceback

import arcpy

import lm_util_config as util


arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

gprint = arcpy.AddMessage


def raster_aggregator(argv=None):
    """Main function

    Called by ArcMap with parameters or run from command line with parameters
    entered in script below.

    """
    if argv is None:
        argv = sys.argv  # Get parameters from ArcGIS tool dialog

    try:

        OUTPUTDIR = argv[1]  # Output directory
        AG_FACTOR =  int(argv[2])
        METHOD = argv[3]
        SMOOTH = util.str2bool(argv[4])
        RESRAS = {}#list of resistance rasters
        RESRAS[1] = argv[5]
        RESRAS[2] = argv[6]
        RESRAS[3] = argv[7]
        RESRAS[4] = argv[8]
        RESRAS[5] = argv[9]

        if AG_FACTOR < 2 or AG_FACTOR > 99:
            msg = ('ERROR: Cell factor must be between 2 and 99.')
            arcpy.AddError(msg)
            exit(1)

        OUTPUTGDBNAME = 'agregated_rasters.gdb'
        OUTPUTGDB = path.join(OUTPUTDIR, path.basename(OUTPUTGDBNAME))
        if not arcpy.Exists(OUTPUTGDB):
            arcpy.CreateFileGDB_management(OUTPUTDIR, path.basename(OUTPUTGDBNAME))

        numRasters = 1
        if RESRAS[2] != '#':
            numRasters = 2
        if RESRAS[3] != '#':
            numRasters = 3
        if RESRAS[4] != '#':
            numRasters = 4
        if RESRAS[5] != '#':
            numRasters = 5

        gprint('\nThere are ' + str(numRasters) + ' rasters to aggregate.')
        gprint('\nCell factor is ' + str(AG_FACTOR))
        gprint('\nCell sizes will be multiplied by this amount')
        for rasterNum in range(1,numRasters+1):
            inputRaster = RESRAS[rasterNum]
            arcpy.env.snapRaster = inputRaster
            oldCellSize = arcpy.Describe(inputRaster).MeanCellHeight
            dir,fileName = path.split(inputRaster)

            if SMOOTH == True and METHOD == "MEAN":
                gprint('\nSmoothing cell values by taking mean of ' + str(AG_FACTOR) +'x' + str(AG_FACTOR) + ' neighborhood')

                InNeighborhood = ("Rectangle " + str(AG_FACTOR) + " " +
                                     str(AG_FACTOR) + " CELL")
                InNoDataOption = "NODATA"

                # FocalStatistics
                inputRaster = arcpy.sa.FocalStatistics(
                    inputRaster, InNeighborhood, "", InNoDataOption)

            gprint('\nAggregating raster "' + str(inputRaster) + '"')
            gprint('Old size was ' + str(oldCellSize))
            gprint('New cell size will be ' + str(oldCellSize * AG_FACTOR))
            gprint('Aggregation method: ' + METHOD)
            agg_rast = arcpy.sa.Aggregate(
                inputRaster, AG_FACTOR, METHOD, "TRUNCATE", "NODATA")

            if SMOOTH == True and METHOD == 'MEAN':
                finalRasterFN = fileName + '_CellFactor' + str(AG_FACTOR) + '_' + METHOD + '_sm'
            else:
                finalRasterFN = fileName + '_CellFactor' + str(AG_FACTOR) + '_' + METHOD

            gprint("\nAggregated raster will be saved to: " + finalRasterFN)
            gprint("in geodatabase " + OUTPUTGDB)
            agg_rast.save(path.join(OUTPUTGDB, finalRasterFN))

        gprint('Done.')

    # Return GEOPROCESSING specific errors
    except arcpy.ExecuteError:
        raise_geoproc_error(__filename__)

    # Return any PYTHON or system specific errors
    except Exception:
        raise_python_error(__filename__)


def raise_geoproc_error(filename):
    """Handle geoprocessor errors and provide details to user"""
    dashline(1)
    tb = sys.exc_info()[2]  # get the traceback object
    # tbinfo contains the error's line number and the code
    tbinfo = traceback.format_tb(tb)[0]
    line = tbinfo.split(", ")[1]

    arcpy.AddError("Geoprocessing error on **" + line + "** of " + filename +
                " :")
    dashline(1)
    for msg in range(0, arcpy.GetMessageCount() - 1):
        if arcpy.GetSeverity(msg) == 2:
            arcpy.AddReturnMessage(msg)
        print(arcpy.AddReturnMessage(msg))
    exit(0)


def raise_python_error(filename):

    """Handle python errors and provide details to user"""
    dashline(1)
    tb = sys.exc_info()[2]  # get the traceback object
    # tbinfo contains the error's line number and the code
    tbinfo = traceback.format_tb(tb)[0]
    line = tbinfo.split(", ")[1]

    err = traceback.format_exc().splitlines()[-1]

    arcpy.AddError("Python error on **" + line + "** of " + filename)
    arcpy.AddError(err)
    exit(0)


def dashline(lspace=0):
    """Output dashed line in tool output dialog.

       0 = No empty line
       1 = Empty line before
       2 = Empty line after

    """
    if lspace == 1:
        arcpy.AddMessage('\n')
    arcpy.AddMessage('---------------------------------')
    if lspace == 2:
        arcpy.AddMessage('\n')


if __name__ == "__main__":
    raster_aggregator()
