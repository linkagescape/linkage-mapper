#!/usr/bin/env python2.5
# Authors: Brad McRae and Darren Kavanagh

""" Step 1: Get adjacencies.

Determines adjacencies between core areas in either or both Euclidean and
cost-weighted distance space

"""

import time
import os.path as path

import numpy as npy

import arcgisscripting

from lm_config import tool_env as cfg
import lm_util as lu

try:
    import arcpy
    from arcpy.sa import *
    arcpy.CheckOutExtension("spatial")
    gp = arcpy.gp
    arcgisscripting = arcpy
except:
    arcpy = False
    import arcgisscripting
    gp = cfg.gp


_SCRIPT_NAME = "s1_getAdjacencies.py"

gp = cfg.gp
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
        gp.scratchWorkspace = cfg.ARCSCRATCHDIR
        gp.workspace = cfg.PROJECTDIR

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

            gp.MakeFeatureLayer(cfg.COREFC, cfg.FCORES)

            extentBoxList = npy.zeros((0, 5), dtype='float32')
            boxCoords = lu.get_extent_box_coords()
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
            gp.buffer_analysis(cfg.BNDCIRCEN, cfg.BNDCIR, "radius")

            del boundingCirclePointArray

        gp.pyramid = "NONE"
        gp.rasterstatistics = "NONE"
        gp.workspace = cfg.SCRATCHDIR

        if cfg.S1ADJMETH_CW:
            cwadjacency()
        if cfg.S1ADJMETH_EU:
            euadjacency()

    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 1. Details follow.****')
        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 1. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)
    return


def cwadjacency():
    """Calculate cost-weighted adjacency

       Inputs: gp - geoprocessing object

    """
    try:
        ALLOC_RASFN = "CWD_alloc_ras"

        
        gprint('\nCalculating cost-weighted distance adjacency')
        outcsvfile = cfg.CWDADJFILE
        outcsvLogfile = path.join(cfg.LOGDIR, "cwdAdj_STEP1.csv")
        PREFIX = cfg.PREFIX

        # May need to set extent prior to core poly to raster conversion...
        # ----------------------------------------------
        # Cost-weighted allocation code
        gp.cellSize = gp.Describe(cfg.RESRAST).MeanCellHeight
        gp.extent = gp.Describe(cfg.RESRAST).extent
        if cfg.BUFFERDIST is not None:
            # Clip resistance raster using bounding circle
            start_time = time.clock()
            gp.cellSize = gp.Describe(cfg.RESRAST).MeanCellHeight#xxx
            gp.extent = gp.Describe(cfg.RESRAST).Extent#xxx
            bResistance = path.join(cfg.SCRATCHDIR, "bResistance")
            gp.ExtractByMask_sa(cfg.RESRAST, cfg.BNDCIR,
                                    bResistance)
            gprint('\nReduced resistance raster extracted using '
                              'bounding circle.')
            start_time = lu.elapsed_time(start_time)
        else:
            bResistance = cfg.RESRAST

        start_time = time.clock()
        gprint('Starting cost-weighted distance allocation...')

        # core_rastmp = 'core_rastmp'
        if cfg.TMAXCWDIST is not None:
            gprint('Maximum cost-weighted distance set to ' +
                              str(cfg.TMAXCWDIST))
        gp.CellSize = gp.Describe(bResistance).MeanCellHeight
        gp.extent = "MAXOF"
        gprint('Processing cell size: ' + gp.CellSize)

        gp.workspace = cfg.ADJACENCYDIR
        gp.scratchworkspace = cfg.ARCSCRATCHDIR

        lu.delete_data(cfg.CWDGDB)
        if not gp.exists(cfg.CWDGDB):
            gp.createfilegdb(cfg.OUTPUTDIR, path.basename(cfg.CWDGDB))
        outDistanceRaster = path.join(cfg.CWDGDB, PREFIX + "_cwd")
        alloc_ras = path.join(cfg.ADJACENCYDIR, ALLOC_RASFN)
        lu.delete_data(alloc_ras)
        lu.delete_data(outDistanceRaster)

        count = 0


        if arcpy:
            statement = ('costAllocOut = CostAllocation(cfg.CORERAS, '
                        'bResistance, cfg.TMAXCWDIST, cfg.CORERAS,"VALUE", '
                        'outDistanceRaster);'
                        'costAllocOut.save(alloc_ras)')
        else:
            statement = ('gp.Costallocation_sa(cfg.CORERAS, bResistance, '
                     'alloc_ras, cfg.TMAXCWDIST, cfg.CORERAS, "VALUE", '
                     'outDistanceRaster, "")')
        while True:
            try:
                exec statement
            except:
                count, tryAgain = lu.retry_arc_error(count, statement)
                if not tryAgain:
                    exec statement
            else:
                break
        gprint('\nBuilding output statistics and pyramids for CWD raster.')
        lu.build_stats(outDistanceRaster)
        gp.scratchworkspace = cfg.ARCSCRATCHDIR
        gprint('Cost-weighted distance allocation done.')
        start_time = lu.elapsed_time(start_time)
        adjshiftwrite(alloc_ras, outcsvfile, outcsvLogfile)

    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 1. Details follow.****')
        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 1. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)


def euadjacency():
    """Calculate Euclidean adjacency

       Inputs: gp - geoprocessing object

    """
    try:
        ALLOC_RASFN = "Euc_alloc_ras"
        lu.dashline()
        gprint('Calculating Euclidean adjacency')
        outcsvfile = cfg.EUCADJFILE
        outcsvLogfile = path.join(cfg.LOGDIR, "eucAdj_STEP1.csv")

        # ----------------------------------------------
        # Euclidean allocation code
        gp.workspace = cfg.ADJACENCYDIR
        gprint('Starting Euclidean adjacency processing...')
        # Euclidean cell size
        cellSizeEuclidean = gp.Describe(cfg.RESRAST).MeanCellHeight

        oldextent = gp.extent
        if cfg.BUFFERDIST is not None:
            gp.extent = gp.Describe(cfg.BNDCIR).extent

        start_time = time.clock()

        gp.scratchworkspace = cfg.ARCSCRATCHDIR
        outDistanceRaster = path.join(cfg.ADJACENCYDIR, "euc")
        alloc_ras = path.join(cfg.ADJACENCYDIR, ALLOC_RASFN)
        lu.delete_data(alloc_ras)
        lu.delete_data(outDistanceRaster)

        count = 0
        statement = ('gp.EucAllocation_sa(cfg.CORERAS, alloc_ras, "","", '
                     'cellSizeEuclidean, "", outDistanceRaster, "")')
        while True:
            try:
                exec statement
            except:
                count, tryAgain = lu.retry_arc_error(count, statement)
                if not tryAgain:
                    exec statement
            else:
                break

        gp.scratchworkspace = cfg.ARCSCRATCHDIR
        gprint('\nEuclidean distance allocation done.')
        start_time = lu.elapsed_time(start_time)
        gp.extent = oldextent
        adjshiftwrite(alloc_ras, outcsvfile, outcsvLogfile)

        # Clean up
        lu.delete_data(outDistanceRaster)

     # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 1. Details follow.****')

        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 1. Details follow.****')

        lu.exit_with_python_error(_SCRIPT_NAME)


def adjshiftwrite(araster, csvfile, logfile):
    """Get adjacencies using shift method and write to disk"""
    # To be replaced by getLeastCostDistsUsingShiftMethod if implemented
    adjTable = lu.get_adj_using_shift_method(araster)
    lu.write_adj_file(csvfile, adjTable)
    lu.write_adj_file(logfile, adjTable)
