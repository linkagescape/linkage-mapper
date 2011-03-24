##*****************************************************************
## 2011_0128
## NAME: s1_getAdjacencies.py
##
## SUMMARY: Determines adjacencies between core areas in either or both
## Euclidean and cost-weighted distance space
##
## SOFTWARE: ArcGIS 9.3 (requires Spatial Analyst extension)
##           Python 2.5
##
##*****************************************************************

import shutil
import arcgisscripting, sys, time
from time import localtime, strftime
import os.path as path
import string,csv
from numpy import *
from string import split
from numpy import loadtxt, where, delete, arange

import lm_config
import lm_util as lu

PROJECTDIR = lm_config.PROJECTDIR
OUTPUTDIR = lm_config.OUTPUTDIR
LOGDIR = lm_config.LOGDIR
SCRATCHDIR = lm_config.SCRATCHDIR
DATAPASSDIR = lm_config.DATAPASSDIR
ADJACENCYDIR = lm_config.ADJACENCYDIR
RESRAST = lm_config.RESRAST
COREFC = lm_config.COREFC
COREFN = lm_config.COREFN
S1ADJMETH_CW = lm_config.S1ADJMETH_CW
S1ADJMETH_EU = lm_config.S1ADJMETH_EU
BUFFERDIST = lm_config.BUFFERDIST
TMAXCWDIST = lm_config.TMAXCWDIST
FCORES = lm_config.FCORES
BNDCIRCEN = lm_config.BNDCIRCEN
BNDCIR = lm_config.BNDCIR
GP = lm_config.GP
BNDCIRWD  = path.join(SCRATCHDIR, BNDCIR)
S1CORE_RAS = "s1core_ras"
ECBNDCIRCEN = "eucBoundingCircleCenter.shp"
ECBNDCIRCENWD = path.join(SCRATCHDIR, ECBNDCIRCEN)
ECBNDCIR = "eucBoundingCircle.shp"
ECBNDCIRWD = path.join(SCRATCHDIR, ECBNDCIR)

def step1_get_adjacencies():
    """Determines adjacencies between core areas in either or both
    Euclidean and cost-weighted distance space.

    """
    try:
        global GP
        # Default behavior is to use same cell size as resistance raster,
        # but this can be changed here.

        GP.scratchWorkspace = SCRATCHDIR
        GP.workspace = PROJECTDIR

        GP.AddMessage('Adjacency files will be written to ' + ADJACENCYDIR)
        if path.exists(ADJACENCYDIR):
            shutil.rmtree(ADJACENCYDIR)
        GP.CreateFolder_management(path.dirname(ADJACENCYDIR),
                                       path.basename(ADJACENCYDIR))

        # To set spatial reference for shapefiles we create later
        SR = GP.describe(COREFC).SpatialReference

        # ------------------------------------------------------------------
        # Create bounding circles to limit cwd and allocation calculations
        if BUFFERDIST is not None:
            GP.addmessage('Reducing processing area using bounding circle plus'
                          'buffer of ' + str(float(BUFFERDIST)/1000) +
                          ' km.')

            startTime = time.clock()
            GP.MakeFeatureLayer(COREFC, FCORES)

            extentBoxList = zeros((0,5),dtype='float32')
            boxCoords = lu.get_extent_box_coords()
            extentBoxList = append(extentBoxList,boxCoords,axis=0)
            extentBoxList[0,0] = 0

            # cwd bounding circle- used to clip raster to limit cwd
            # calculations
            boundingCirclePointArray  = zeros((0,5),dtype='float32')
            circlePointData = lu.get_bounding_circle_data(extentBoxList, 0, 0,
                                                          BUFFERDIST)
            lu.make_points(SCRATCHDIR, circlePointData, BNDCIRCEN)
            GP.defineprojection(path.join(SCRATCHDIR, BNDCIRCEN), SR)
            if GP.Exists(BNDCIRWD):
                GP.delete_management(BNDCIRWD)
            GP.buffer_analysis(path.join(SCRATCHDIR, BNDCIRCEN), BNDCIRWD,
                               "radius")
            GP.defineprojection(BNDCIRWD, SR)
            # boundingCirclePointArray  = zeros((0,5),dtype='float32')
            # circlePointData = lu.get_bounding_circle_data(extentBoxList, 0,
            #                                                    0, BUFFERDIST)

            # euc bounding circle- at the moment just limits extent of
            # euclidean allocation calculations
            boundingCirclePointArray  = zeros((0,5),dtype='float32')
            circlePointData = lu.get_bounding_circle_data(extentBoxList, 0, 0,
                                                          0) # no buffer needed
            lu.make_points(SCRATCHDIR, circlePointData, ECBNDCIRCEN)
            GP.defineprojection(ECBNDCIRCENWD, SR)
            if GP.Exists(ECBNDCIRWD):
                GP.delete_management(ECBNDCIRWD)
            GP.buffer_analysis(ECBNDCIRCENWD, ECBNDCIRWD, "radius")
            GP.defineprojection(ECBNDCIRWD, SR)

            del boundingCirclePointArray

        GP.pyramid = "NONE"
        GP.rasterstatistics = "NONE"
        GP.workspace = SCRATCHDIR

        if S1ADJMETH_CW:
            cwadjacency()
        if S1ADJMETH_EU:
            euadjacency()

    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
       lu.dashline(1)
       GP.addmessage('****Failed in step 1. Details follow.****')
       filename =  __file__
       lu.raise_geoproc_error(filename)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        GP.addmessage('****Failed in step 1. Details follow.****')
        filename =  __file__
        lu.raise_python_error(filename)

    del GP
    return

def cwadjacency():
    """Calculate cost-weighted adjacency

       Inputs: GP - geoprocessing object

    """
    try:
        alloc_rasFN = "Cwd_alloc_ras"
        GP.addmessage('\nCalculating cost-weighted distance adjacency')
        outcsvfile = path.join(DATAPASSDIR, "cwdAdj.csv")
        outcsvLogfile = path.join(LOGDIR, "cwdAdj_step1.csv")
        # May need to set extent prior to core poly to raster conversion...
        # ----------------------------------------------
        # Cost-weighted allocation code
        if BUFFERDIST is not None:
            # Clip resistance raster using bounding circle
            startTime = time.clock()
            bResistance = path.join(SCRATCHDIR, "bResistance")
            GP.ExtractByMask_sa(RESRAST, BNDCIRWD,
                                bResistance)
            GP.addmessage('\nReduced resistance raster extracted using '
                          'bounding circle.')
            startTime, hours, mins, secs = lu.elapsed_time(startTime)
        else:
            bResistance = RESRAST

        startTime = time.clock()
        GP.addmessage('Starting cost weighted distance allocation...')

        if TMAXCWDIST is not None:
            GP.addmessage('Maximum cost-weighted distance set to ' +
                          str(TMAXCWDIST))
        GP.CellSize = GP.Describe(bResistance).MeanCellHeight
        GP.extent = "MAXOF"
        GP.AddMessage('Processing cell size: ' + GP.CellSize)
        count = 0
        statement = ('GP.FeatureToRaster_conversion(COREFC, COREFN, '
                    'S1CORE_RAS, GP.Cellsize)')
        while True:
            try: exec statement
            except:
                count, tryAgain = lu.hiccup_test(count, statement)
                if not tryAgain: exec statement
            else: break

        GP.workspace = ADJACENCYDIR
        GP.scratchworkspace = GP.workspace

        # fixme: put this in geodatabase instead
        outDistanceRaster = path.join(OUTPUTDIR, "cwd")
        alloc_ras = path.join(ADJACENCYDIR, alloc_rasFN)
        s1core_ras_path = path.join(SCRATCHDIR, S1CORE_RAS)
        count = 0
        statement = ('GP.Costallocation_sa(s1core_ras_path, bResistance, '
                     'alloc_ras, TMAXCWDIST, s1core_ras_path, "VALUE", '
                     'outDistanceRaster, "")')
        while True:
            try: exec statement
            except:
                count, tryAgain = lu.hiccup_test(count, statement)
                if not tryAgain: exec statement
            else: break

        GP.scratchworkspace = SCRATCHDIR
        GP.addmessage('\nCost-weighted distance allocation done.')
        startTime, hours, mins, secs = lu.elapsed_time(startTime)
        adjshiftwrite(alloc_ras, outcsvfile, outcsvLogfile)

    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
       lu.dashline(1)
       GP.addmessage('****Failed in step 1. Details follow.****')
       filename =  __file__
       lu.raise_geoproc_error(filename)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        GP.addmessage('****Failed in step 1. Details follow.****')
        filename =  __file__
        lu.raise_python_error(filename)

def euadjacency():
    """Calculate Euclidean adjacency

       Inputs: GP - geoprocessing object

    """
    alloc_rasFN = "Euc_alloc_ras"
    GP.addmessage('\nCalculating Euclidean adjacency')
    outcsvfile = path.join(DATAPASSDIR, "eucAdj.csv")
    outcsvLogfile = path.join(LOGDIR, "eucAdj_step1.csv")

    # FIXME- would be good to have bounding circle affect euclidean calcs too
    # ----------------------------------------------
    # Euclidean allocation code
    GP.workspace = ADJACENCYDIR
    GP.addmessage ('Starting Euclidean adjacency processing...')
    # Euclidean cell size
    cellSizeEuclidean = GP.Describe(RESRAST).MeanCellHeight
    GP.addmessage('Euclidean cell size set equal to resistance '
                  'raster cell size (' +
                  str(cellSizeEuclidean) + ').')

    oldextent = GP.extent
    if BUFFERDIST is not None:
        GP.extent = GP.Describe(ECBNDCIRWD).extent
    count = 0
    statement = ('GP.FeatureToRaster_conversion(COREFC, COREFN, '
                 'S1CORE_RAS, cellSizeEuclidean)')
    while True:
        try: exec statement
        except:
            count,tryAgain = lu.hiccup_test(count,statement)
            if not tryAgain: exec statement
        else: break

    startTime=time.clock()

    GP.scratchworkspace = GP.workspace
    outDistanceRaster = path.join(ADJACENCYDIR, "euc")
    alloc_ras = path.join(ADJACENCYDIR, alloc_rasFN)
    count = 0
    statement = ('GP.EucAllocation_sa(S1CORE_RAS, alloc_ras, "","", '
                 'cellSizeEuclidean, "", outDistanceRaster, "")')
    while True:
        try: exec statement
        except:
            count,tryAgain = lu.hiccup_test(count,statement)
            if not tryAgain: exec statement
        else: break

    GP.scratchworkspace = SCRATCHDIR
    GP.addmessage('\nEuclidean distance allocation done.')
    startTime, hours, mins, secs = lu.elapsed_time(startTime)
    GP.extent = oldextent
    adjshiftwrite(alloc_ras, outcsvfile, outcsvLogfile)

def adjshiftwrite(araster, csvfile, logfile):
    """Get adjacencies using shift method and write to disk"""
    # To be replaced by getLeastCostDistsUsingShiftMethod if implemented
    adjTable = lu.get_adj_using_shift_method(araster)
    lu.write_adj_file(csvfile, adjTable)
    lu.write_adj_file(logfile, adjTable)