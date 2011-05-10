#!/usr/bin/env python2.5

""" Step 1: Get adjacencies.

Determines adjacencies between core areas in either or both Euclidean and
cost-weighted distance space

"""

__filename__ = "s1_getAdjacencies.py"
__version__ = "0.6.2"

import shutil
import time
import os.path as path

import numpy as npy

import arcgisscripting

from lm_config import Config as Cfg
import lm_util as lu

BNDCIRCENWD = path.join(Cfg.SCRATCHDIR, Cfg.BNDCIRCEN)
BNDCIRWD = path.join(Cfg.SCRATCHDIR, Cfg.BNDCIR)
EUC_BNDCIRCEN = "euc" + Cfg.BNDCIRCEN
EUC_BNDCIRCENWD = path.join(Cfg.SCRATCHDIR, EUC_BNDCIRCEN)
EUC_BNDCIR = "euc" + Cfg.BNDCIR
EUC_BNDCIRWD = path.join(Cfg.SCRATCHDIR, EUC_BNDCIR)
S1CORE_RAS = "s1core_ras"


def STEP1_get_adjacencies():
    """Determines adjacencies between core areas in either or both
    Euclidean and cost-weighted distance space.

    """
    try:
        lu.dashline(1)
        Cfg.gp.addmessage('Running script ' + __filename__)

        # Default behavior is to use same cell size as resistance raster,
        # but this can be changed here.
        Cfg.gp.scratchWorkspace = Cfg.SCRATCHDIR
        Cfg.gp.workspace = Cfg.PROJECTDIR

        Cfg.gp.AddMessage('Adjacency files will be written to ' +
                          Cfg.ADJACENCYDIR)
        if path.exists(Cfg.ADJACENCYDIR):
            shutil.rmtree(Cfg.ADJACENCYDIR)
        Cfg.gp.CreateFolder_management(path.dirname(Cfg.ADJACENCYDIR),
                                       path.basename(Cfg.ADJACENCYDIR))

        # To set spatial reference for shapefiles we create later
        SR = Cfg.gp.describe(Cfg.COREFC).SpatialReference

        # ------------------------------------------------------------------
        # Create bounding circles to limit cwd and allocation calculations
        if Cfg.BUFFERDIST is not None:
            Cfg.gp.addmessage('Reducing processing area using bounding circle '
                              'plus buffer of ' +
                              str(float(Cfg.BUFFERDIST)) + ' map units')

            Cfg.gp.MakeFeatureLayer(Cfg.COREFC, Cfg.FCORES)

            extentBoxList = npy.zeros((0, 5), dtype='float32')
            boxCoords = lu.get_extent_box_coords()
            extentBoxList = npy.append(extentBoxList, boxCoords, axis=0)
            extentBoxList[0, 0] = 0

            # cwd bounding circle- used to clip raster to limit cwd
            # calculations
            boundingCirclePointArray = npy.zeros((0, 5), dtype='float32')
            circlePointData = lu.get_bounding_circle_data(extentBoxList, 0, 0,
                                                          Cfg.BUFFERDIST)
            lu.make_points(Cfg.SCRATCHDIR, circlePointData, Cfg.BNDCIRCEN)
            Cfg.gp.defineprojection(path.join(Cfg.SCRATCHDIR, Cfg.BNDCIRCEN),
                                    SR)
            if Cfg.gp.Exists(BNDCIRWD):
                Cfg.gp.delete_management(BNDCIRWD)
            Cfg.gp.buffer_analysis(path.join(Cfg.SCRATCHDIR, Cfg.BNDCIRCEN),
                                             BNDCIRWD, "radius")
            Cfg.gp.defineprojection(BNDCIRWD, SR)
            # boundingCirclePointArray  = npy.zeros((0, 5),dtype='float32')
            # circlePointData = lu.get_bounding_circle_data(extentBoxList, 0,
            #                                               0, Cfg.BUFFERDIST)

            # euc bounding circle- at the moment just limits extent of
            # euclidean allocation calculations
            boundingCirclePointArray = npy.zeros((0, 5), dtype='float32')
            # no buffer needed
            circlePointData = lu.get_bounding_circle_data(
                extentBoxList, 0, 0, 0)
            lu.make_points(Cfg.SCRATCHDIR, circlePointData, EUC_BNDCIRCEN)
            Cfg.gp.defineprojection(EUC_BNDCIRCENWD, SR)
            if Cfg.gp.Exists(EUC_BNDCIRWD):
                Cfg.gp.delete_management(EUC_BNDCIRWD)
            Cfg.gp.buffer_analysis(EUC_BNDCIRCENWD, EUC_BNDCIRWD, "radius")
            Cfg.gp.defineprojection(EUC_BNDCIRWD, SR)

            del boundingCirclePointArray

        Cfg.gp.pyramid = "NONE"
        Cfg.gp.rasterstatistics = "NONE"
        Cfg.gp.workspace = Cfg.SCRATCHDIR

        if Cfg.S1ADJMETH_CW:
            cwadjacency()
        if Cfg.S1ADJMETH_EU:
            euadjacency()

    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        Cfg.gp.addmessage('****Failed in step 1. Details follow.****')
        lu.raise_geoproc_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        Cfg.gp.addmessage('****Failed in step 1. Details follow.****')
        lu.raise_python_error(__filename__)
    return


def cwadjacency():
    """Calculate cost-weighted adjacency

       Inputs: Cfg.gp - geoprocessing object

    """
    try:
        alloc_rasFN = "Cwd_alloc_ras"
        Cfg.gp.addmessage('\nCalculating cost-weighted distance adjacency')
        outcsvfile = path.join(Cfg.DATAPASSDIR, "cwdAdj.csv")
        outcsvLogfile = path.join(Cfg.LOGDIR, "cwdAdj_STEP1.csv")
        # May need to set extent prior to core poly to raster conversion...
        # ----------------------------------------------
        # Cost-weighted allocation code
        if Cfg.BUFFERDIST is not None:
            # Clip resistance raster using bounding circle
            start_time = time.clock()
            bResistance = path.join(Cfg.SCRATCHDIR, "bResistance")
            Cfg.gp.ExtractByMask_sa(Cfg.RESRAST, BNDCIRWD,
                                    bResistance)
            Cfg.gp.addmessage('\nReduced resistance raster extracted using '
                              'bounding circle.')
            start_time = lu.elapsed_time(start_time)
        else:
            bResistance = Cfg.RESRAST

        start_time = time.clock()
        Cfg.gp.addmessage('Starting cost weighted distance allocation...')

        if Cfg.TMAXCWDIST is not None:
            Cfg.gp.addmessage('Maximum cost-weighted distance set to ' +
                              str(Cfg.TMAXCWDIST))
        Cfg.gp.CellSize = Cfg.gp.Describe(bResistance).MeanCellHeight
        Cfg.gp.extent = "MAXOF"
        Cfg.gp.AddMessage('Processing cell size: ' + Cfg.gp.CellSize)
        count = 0
        statement = ('Cfg.gp.FeatureToRaster_conversion(Cfg.COREFC, '
                     'Cfg.COREFN, S1CORE_RAS, Cfg.gp.Cellsize)')
        while True:
            try:
                exec statement
            except:
                count, tryAgain = lu.hiccup_test(count, statement)
                if not tryAgain:
                    exec statement
            else:
                break

        Cfg.gp.workspace = Cfg.ADJACENCYDIR
        Cfg.gp.scratchworkspace = Cfg.gp.workspace

        # fixme: put this in geodatabase instead
        outDistanceRaster = path.join(Cfg.OUTPUTDIR, "cwd")
        alloc_ras = path.join(Cfg.ADJACENCYDIR, alloc_rasFN)
        s1core_ras_path = path.join(Cfg.SCRATCHDIR, S1CORE_RAS)
        count = 0
        statement = ('Cfg.gp.Costallocation_sa(s1core_ras_path, bResistance, '
                     'alloc_ras, Cfg.TMAXCWDIST, s1core_ras_path, "VALUE", '
                     'outDistanceRaster, "")')
        while True:
            try:
                exec statement
            except:
                count, tryAgain = lu.hiccup_test(count, statement)
                if not tryAgain:
                    exec statement
            else:
                break

        Cfg.gp.scratchworkspace = Cfg.SCRATCHDIR
        Cfg.gp.addmessage('\nCost-weighted distance allocation done.')
        start_time = lu.elapsed_time(start_time)
        adjshiftwrite(alloc_ras, outcsvfile, outcsvLogfile)

    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        Cfg.gp.addmessage('****Failed in step 1. Details follow.****')
        lu.raise_python_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        Cfg.gp.addmessage('****Failed in step 1. Details follow.****')
        lu.raise_python_error(__filename__)


def euadjacency():
    """Calculate Euclidean adjacency

       Inputs: Cfg.gp - geoprocessing object

    """
    try:
        alloc_rasFN = "Euc_alloc_ras"
        lu.dashline()
        Cfg.gp.addmessage('Calculating Euclidean adjacency')
        outcsvfile = path.join(Cfg.DATAPASSDIR, "eucAdj.csv")
        outcsvLogfile = path.join(Cfg.LOGDIR, "eucAdj_STEP1.csv")

        # FIXME - would be good to have bounding circle affect euclidean calcs
        # too
        # ----------------------------------------------
        # Euclidean allocation code
        Cfg.gp.workspace = Cfg.ADJACENCYDIR
        Cfg.gp.addmessage('Starting Euclidean adjacency processing...')
        # Euclidean cell size
        cellSizeEuclidean = Cfg.gp.Describe(Cfg.RESRAST).MeanCellHeight
        # Cfg.gp.addmessage('Euclidean cell size set equal to resistance '
                      # 'raster cell size (' +
                      # str(cellSizeEuclidean) + ').')

        oldextent = Cfg.gp.extent
        if Cfg.BUFFERDIST is not None:
            Cfg.gp.extent = Cfg.gp.Describe(BNDCIRWD).extent
        count = 0
        statement = ('Cfg.gp.FeatureToRaster_conversion(Cfg.COREFC, '
                     'Cfg.COREFN, S1CORE_RAS, cellSizeEuclidean)')
        while True:
            try:
                exec statement
            except:
                count, tryAgain = lu.hiccup_test(count, statement)
                if not tryAgain:
                    exec statement
            else:
                break

        start_time = time.clock()

        Cfg.gp.scratchworkspace = Cfg.gp.workspace
        outDistanceRaster = path.join(Cfg.ADJACENCYDIR, "euc")
        alloc_ras = path.join(Cfg.ADJACENCYDIR, alloc_rasFN)
        count = 0
        statement = ('Cfg.gp.EucAllocation_sa(S1CORE_RAS, alloc_ras, "","", '
                     'cellSizeEuclidean, "", outDistanceRaster, "")')
        while True:
            try:
                exec statement
            except:
                count, tryAgain = lu.hiccup_test(count, statement)
                if not tryAgain:
                    exec statement
            else:
                break

        Cfg.gp.scratchworkspace = Cfg.SCRATCHDIR
        Cfg.gp.addmessage('\nEuclidean distance allocation done.')
        start_time = lu.elapsed_time(start_time)
        Cfg.gp.extent = oldextent
        adjshiftwrite(alloc_ras, outcsvfile, outcsvLogfile)

     # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        Cfg.gp.addmessage('****Failed in step 1. Details follow.****')

        lu.raise_python_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        Cfg.gp.addmessage('****Failed in step 1. Details follow.****')

        lu.raise_python_error(__filename__)


def adjshiftwrite(araster, csvfile, logfile):
    """Get adjacencies using shift method and write to disk"""
    # To be replaced by getLeastCostDistsUsingShiftMethod if implemented
    adjTable = lu.get_adj_using_shift_method(araster)
    lu.write_adj_file(csvfile, adjTable)
    lu.write_adj_file(logfile, adjTable)
