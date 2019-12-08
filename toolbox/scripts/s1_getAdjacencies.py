# Authors: Brad McRae and Darren Kavanagh

""" Step 1: Get adjacencies.

Determines adjacencies between core areas in either or both Euclidean and
cost-weighted distance space

"""

from os import path
import time

import numpy as npy
import arcpy

from lm_config import tool_env as cfg
import lm_util as lu


_SCRIPT_NAME = "s1_getAdjacencies.py"

gprint = lu.gprint


def STEP1_get_adjacencies():
    """Determines adjacencies between core areas in either or both
    Euclidean and cost-weighted distance space.

    """
    try:
        lu.dashline(1)
        gprint('Running script ' + _SCRIPT_NAME)

        # Default behavior is to use same cell size as resistance raster,
        # but this can be changed here.
        arcpy.env.scratchWorkspace = cfg.ARCSCRATCHDIR
        arcpy.env.workspace = cfg.PROJECTDIR

        #remove adj directory and files from previous runs
        lu.delete_dir(cfg.ADJACENCYDIR)
        lu.delete_data(cfg.CWDADJFILE)
        lu.delete_data(cfg.EUCADJFILE)

        if not cfg.S2ADJMETH_CW and not cfg.S2ADJMETH_EU:
            # Adjacency not needed
            return

        lu.create_dir(cfg.ADJACENCYDIR)

        gprint('Adjacency files will be written to ' +
                          cfg.ADJACENCYDIR)

        # ------------------------------------------------------------------
        # Create bounding circles to limit cwd and allocation calculations
        if cfg.BUFFERDIST is not None:
            gprint('Reducing processing area using bounding circle '
                              'plus buffer of ' +
                              str(float(cfg.BUFFERDIST)) + ' map units')

            extentBoxList = npy.zeros((0, 5), dtype='float32')
            boxCoords = lu.get_ext_box_coords(cfg.COREFC)
            extentBoxList = npy.append(extentBoxList, boxCoords, axis=0)
            extentBoxList[0, 0] = 0

            # cwd bounding circle- used to clip raster to limit cwd
            # calculations
            boundingCirclePointArray = npy.zeros((0, 5), dtype='float32')
            circlePointData = lu.get_bounding_circle_data(extentBoxList, 0, 0,
                                                          cfg.BUFFERDIST)
            lu.make_points(cfg.SCRATCHDIR, circlePointData,
                           path.basename(cfg.BNDCIRCEN))

            lu.delete_data(cfg.BNDCIR)
            arcpy.Buffer_analysis(cfg.BNDCIRCEN, cfg.BNDCIR, "radius")

            del boundingCirclePointArray

        arcpy.env.pyramid = "NONE"
        arcpy.env.rasterStatistics = "NONE"
        arcpy.env.workspace = cfg.SCRATCHDIR

        if cfg.S1ADJMETH_CW:
            cwadjacency()
        if cfg.S1ADJMETH_EU:
            euadjacency()

    # Return GEOPROCESSING specific errors
    except arcpy.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 1. Details follow.****')
        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except Exception:
        lu.dashline(1)
        gprint('****Failed in step 1. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)
    return


def cwadjacency():
    """Calculate cost-weighted adjacency."""
    try:
        ALLOC_RASFN = "CWD_alloc_ras"


        gprint('\nCalculating cost-weighted distance adjacency')
        outcsvfile = cfg.CWDADJFILE
        outcsvLogfile = path.join(cfg.LOGDIR, "cwdAdj_STEP1.csv")
        PREFIX = cfg.PREFIX

        # May need to set extent prior to core poly to raster conversion...
        # ----------------------------------------------
        # Cost-weighted allocation code
        arcpy.env.cellSize = arcpy.Describe(cfg.RESRAST).MeanCellHeight
        arcpy.env.extent = arcpy.Describe(cfg.RESRAST).extent
        if cfg.BUFFERDIST is not None:
            # Clip resistance raster using bounding circle
            start_time = time.clock()
            arcpy.env.cellSize = arcpy.Describe(cfg.RESRAST).MeanCellHeight
            arcpy.env.extent = arcpy.Describe(cfg.RESRAST).Extent
            bResistance = arcpy.sa.ExtractByMask(cfg.RESRAST, cfg.BNDCIR)
            gprint('\nReduced resistance raster extracted using '
                              'bounding circle.')
            start_time = lu.elapsed_time(start_time)
        else:
            bResistance = cfg.RESRAST

        start_time = time.clock()
        gprint('Starting cost-weighted distance allocation...')

        if cfg.TMAXCWDIST is not None:
            gprint('Maximum cost-weighted distance set to ' +
                              str(cfg.TMAXCWDIST))
        arcpy.env.cellSize = arcpy.Describe(bResistance).MeanCellHeight
        arcpy.env.extent = "MAXOF"
        gprint('Processing cell size: ' + arcpy.env.cellSize)

        arcpy.env.workspace = cfg.ADJACENCYDIR
        arcpy.env.scratchWorkspace = cfg.ARCSCRATCHDIR

        lu.delete_data(cfg.CWDGDB)
        if not arcpy.Exists(cfg.CWDGDB):
            arcpy.CreateFileGDB_management(cfg.OUTPUTDIR, path.basename(cfg.CWDGDB))
        outDistanceRaster = path.join(cfg.CWDGDB, PREFIX + "_cwd")
        alloc_ras = path.join(cfg.ADJACENCYDIR, ALLOC_RASFN)
        lu.delete_data(alloc_ras)
        lu.delete_data(outDistanceRaster)

        statement = ('costAllocOut = arcpy.sa.CostAllocation(cfg.CORERAS, '
                     'bResistance, cfg.TMAXCWDIST, cfg.CORERAS,"VALUE", '
                     'outDistanceRaster);'
                     'costAllocOut.save(alloc_ras)')
        count = 0
        while True:
            try:
                exec(statement)
            except Exception:
                count, tryAgain = lu.retry_arc_error(count, statement)
                if not tryAgain:
                    exec(statement)
            else:
                break

        gprint('\nBuilding output statistics and pyramids for CWD raster.')
        lu.build_stats(outDistanceRaster)
        arcpy.env.scratchWorkspace = cfg.ARCSCRATCHDIR
        gprint('Cost-weighted distance allocation done.')
        start_time = lu.elapsed_time(start_time)
        adjshiftwrite(alloc_ras, outcsvfile, outcsvLogfile)

    # Return GEOPROCESSING specific errors
    except arcpy.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 1. Details follow.****')
        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except Exception:
        lu.dashline(1)
        gprint('****Failed in step 1. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)


def euadjacency():
    """Calculate Euclidean adjacency."""
    try:
        ALLOC_RASFN = "Euc_alloc_ras"
        lu.dashline()
        gprint('Calculating Euclidean adjacency')
        outcsvfile = cfg.EUCADJFILE
        outcsvLogfile = path.join(cfg.LOGDIR, "eucAdj_STEP1.csv")

        # ----------------------------------------------
        # Euclidean allocation code
        arcpy.env.workspace = cfg.ADJACENCYDIR
        gprint('Starting Euclidean adjacency processing...')
        # Euclidean cell size
        cellSizeEuclidean = arcpy.Describe(cfg.RESRAST).MeanCellHeight

        oldextent = arcpy.env.extent
        if cfg.BUFFERDIST is not None:
            arcpy.env.extent = arcpy.Describe(cfg.BNDCIR).extent

        start_time = time.clock()

        arcpy.env.scratchWorkspace = cfg.ARCSCRATCHDIR
        outDistanceRaster = path.join(cfg.ADJACENCYDIR, "euc")
        alloc_ras = path.join(cfg.ADJACENCYDIR, ALLOC_RASFN)
        lu.delete_data(alloc_ras)
        lu.delete_data(outDistanceRaster)

        count = 0
        statement = ('alloc_raster = arcpy.sa.EucAllocation('
                     'cfg.CORERAS, "", "", '
                     'cellSizeEuclidean, "", outDistanceRaster, ""); '
                     'alloc_raster.save(alloc_ras)')
        while True:
            try:
                exec(statement)
            except Exception:
                count, tryAgain = lu.retry_arc_error(count, statement)
                if not tryAgain:
                    exec(statement)
            else:
                break

        arcpy.env.scratchWorkspace = cfg.ARCSCRATCHDIR
        gprint('\nEuclidean distance allocation done.')
        start_time = lu.elapsed_time(start_time)
        arcpy.env.extent = oldextent
        adjshiftwrite(alloc_ras, outcsvfile, outcsvLogfile)

        # Clean up
        lu.delete_data(outDistanceRaster)

     # Return GEOPROCESSING specific errors
    except arcpy.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 1. Details follow.****')

        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except Exception:
        lu.dashline(1)
        gprint('****Failed in step 1. Details follow.****')

        lu.exit_with_python_error(_SCRIPT_NAME)


def adjshiftwrite(araster, csvfile, logfile):
    """Get adjacencies using shift method and write to disk"""
    # To be replaced by getLeastCostDistsUsingShiftMethod if implemented
    adjTable = lu.get_adj_using_shift_method(araster)
    lu.write_adj_file(csvfile, adjTable)
    lu.write_adj_file(logfile, adjTable)
