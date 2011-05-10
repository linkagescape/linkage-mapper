#!/usr/bin/env python2.5

"""Step 3: Calculate cost-weighted distances.

Calculates cost-weighted distances from each core area.
Uses bounding circles around source and target cores to limit
extent of cwd calculations and speed computation.

"""

__filename__ = "s3_calcCwds.py"
__version__ = "0.6.3"

import os.path as path
import shutil
import time

import arcgisscripting
import numpy as npy

from lm_config import Config as Cfg
import lm_util as lu

BNDCIRCENS = "boundingCircleCenters.shp"
BNDCIRS = "boundingCircles.shp"
BNDFC = "boundingFeature.shp"
ZNSTATS = path.join(Cfg.SCRATCHDIR, "zonestats.dbf")

def STEP3_calc_cwds():
    """Calculates cost-weighted distances from each core area.
    Uses bounding circles around source and target cores to limit
    extent of cwd calculations and speed computation.

    """
    try:
        lu.dashline(1)
        Cfg.gp.addmessage('Running script ' + __filename__)
        Cfg.gp.scratchWorkspace = Cfg.SCRATCHDIR
        linkTableFile = lu.get_prev_step_link_table(step=3)

        if Cfg.TMAXCWDIST is None:
           	Cfg.gp.AddMessage('NOT using a maximum cost-weighted distance.')
        else:
            Cfg.gp.addmessage('Max cost-weighted distance for CWD calcs set '
                              'to ' + str(Cfg.TMAXCWDIST) + '\n')

        # FIXME: because it's integer, it fills in 0 if not entered.
        if (Cfg.BUFFERDIST) is not None:
            Cfg.gp.addmessage('Bounding circles plus a buffer of ' +
                              str(float(Cfg.BUFFERDIST)) + ' map units will '
                              'be used \n to limit extent of cost distance '
                              'calculations.')
        else:
            Cfg.gp.addmessage('NOT using bounding circles in cost distance '
                              'calculations.')

        # set the analysis extent and cell size
        Cfg.gp.CellSize = Cfg.gp.Describe(Cfg.RESRAST).MeanCellHeight
        # So we don't extract rasters that go beyond extent of original raster
        Cfg.gp.Extent = "MINOF"
        Cfg.gp.mask = Cfg.RESRAST
        Cfg.gp.Workspace = Cfg.SCRATCHDIR
        # for later shapefiles
        SR = Cfg.gp.describe(Cfg.COREFC).SpatialReference

        # Load linkTable (created in previous script)
        linkTable = lu.load_link_table(linkTableFile)
        lu.report_links(linkTable)

        rows,cols = npy.where(
            linkTable[:, Cfg.LTB_LINKTYPE:Cfg.LTB_LINKTYPE + 1] == Cfg.LT_CORR)
        coresToProcess = npy.unique(linkTable[:, Cfg.LTB_CORE1:Cfg.LTB_CORE2 + 1])
        maxCoreNum = max(coresToProcess)
        del rows, cols, coresToProcess

        # Set up cwd directories.
        # To keep there from being > 100 grids in any one directory,
        # outputs are written to:
        # cwd\cw for cores 1-99
        # cwd\cw1 for cores 100-199
        # etc.
        if path.exists(Cfg.CWDBASEDIR):
            shutil.rmtree(Cfg.CWDBASEDIR)
        # lu.dashline(1)
        Cfg.gp.addmessage("\nCreating cost-weighted distance grid output folders"
                          ":")
        Cfg.gp.addmessage(path.join(Cfg.CWDBASEDIR, Cfg.CWDSUBDIR_NM))
        Cfg.gp.CreateFolder_management(path.dirname(Cfg.CWDBASEDIR),
                                       path.basename(Cfg.CWDBASEDIR))
        Cfg.gp.CreateFolder_management(Cfg.CWDBASEDIR, Cfg.CWDSUBDIR_NM)
        if maxCoreNum > 100:
            maxDirCount = int(maxCoreNum/100)
            for dirCount in range(1, maxDirCount + 1):
                ccwdir = Cfg.CWDSUBDIR_NM + str(dirCount)
                Cfg.gp.addmessage(ccwdir)
                Cfg.gp.CreateFolder_management(Cfg.CWDBASEDIR, ccwdir)
        # lu.dashline(2)

        # make a feature layer for input cores to select from
        Cfg.gp.MakeFeatureLayer(Cfg.COREFC, Cfg.FCORES)

        # Identify cores to map from LinkTable
        rows,cols = npy.where(linkTable[:,Cfg.LTB_LINKTYPE:Cfg.LTB_LINKTYPE+1] ==
                          Cfg.LT_CORR)
        coresToMap = npy.unique(linkTable[:,Cfg.LTB_CORE1:Cfg.LTB_CORE2+1])
        numCoresToMap = len(coresToMap)
        del rows,cols
        if numCoresToMap < 3:
            # No need to check for intermediate cores, because there aren't any
            Cfg.S3DROPLCCSic = False
        else:
            Cfg.S3DROPLCCSic = Cfg.S3DROPLCCS

        Cfg.gp.addmessage('Number of core areas to connect:' +
                          str(numCoresToMap))

        # Drop links that are too long
        Cfg.gp.addmessage('\nChecking for corridors that are too long to map.')
        disableLeastCostNoVal = False
        linkTable,numDroppedLinks = lu.drop_links(linkTable, Cfg.MAXEUCDIST, 0,
                                                  Cfg.MINEUCDIST, 0,
                                                  disableLeastCostNoVal)

        # ------------------------------------------------------------------
        # Bounding boxes
        if (Cfg.BUFFERDIST) is not None:
            # create bounding boxes around cores
            start_time = time.clock()
            # lu.dashline(1)
            Cfg.gp.addmessage('Calculating bounding boxes for core areas.')
            extentBoxList = npy.zeros((0,5), dtype='float32')
            for x in range(len(coresToMap)):
                core = coresToMap[x]
                if len(coresToMap) > 20:
                    lu.report_pct_done(x, len(coresToMap))
                boxCoords = lu.get_extent_box_coords(core)
                extentBoxList = npy.append(extentBoxList, boxCoords, axis=0)
            Cfg.gp.addmessage('\nDone calculating bounding boxes.')
            start_time = lu.elapsed_time(start_time)
            # lu.dashline()

        # Bounding circle code
        if Cfg.BUFFERDIST is not None:
            # Make a set of circles encompassing core areas we'll be connecting
            start_time = time.clock()
            Cfg.gp.addmessage('Calculating bounding circles around potential'
                          ' corridors.')

            # x y corex corey radius- stores data for bounding circle centroids
            boundingCirclePointArray  = npy.zeros((0,5), dtype='float32')

            circleList = npy.zeros((0,3), dtype='int32')

            numLinks = linkTable.shape[0]
            for x in range(0, numLinks):
                if numLinks > 20:
                    lu.report_pct_done(x, numLinks)
                if linkTable[x,Cfg.LTB_LINKTYPE] == Cfg.LT_CORR:
                    # if it's a valid corridor link
                    linkId = int(linkTable[x,Cfg.LTB_LINKID])
                    # fixme- this code is clumsy- can trim down
                    cores = npy.zeros((1,3), dtype='int32')
                    cores[0,:] = npy.sort([0, linkTable[x,Cfg.LTB_CORE1],
                                      linkTable[x,Cfg.LTB_CORE2]])
                    corex = cores[0,1]
                    corey = cores[0,2]
                    cores[0,0] = linkId

                    ###################
                    foundFlag = False
                    for y in range(0,len(circleList)):  # clumsy
                        if (circleList[y,1] == corex and
                            circleList[y,2] == corey):
                            foundFlag = True
                    if not foundFlag:
                        circlePointData = (
                            lu.get_bounding_circle_data(extentBoxList,
                            corex, corey, Cfg.BUFFERDIST))
                        boundingCirclePointArray = (
                            npy.append(boundingCirclePointArray,circlePointData,
                            axis=0))
                        # keep track of which cores we draw bounding circles
                        # around
                        circleList = npy.append(circleList, cores, axis=0)

            Cfg.gp.addmessage('\nCreating bounding circles using buffer '
                              'analysis.')
            lu.make_points(Cfg.SCRATCHDIR, boundingCirclePointArray,
                           BNDCIRCENS)
            Cfg.gp.defineprojection(BNDCIRCENS, SR)
            if Cfg.gp.Exists(BNDCIRS):
                Cfg.gp.delete_management(BNDCIRS)
            Cfg.gp.buffer_analysis(BNDCIRCENS, BNDCIRS, "radius")
            Cfg.gp.defineprojection(BNDCIRS, SR)
            Cfg.gp.deletefield (BNDCIRS, "BUFF_DIST")

            Cfg.gp.addmessage('Successfully created bounding circles around '
                              'potential corridors using \na buffer of ' +
                              str(float(Cfg.BUFFERDIST)) + ' map units.')
            start_time = lu.elapsed_time(start_time)

            Cfg.gp.addmessage('Reducing global processing area using bounding '
                              'circle plus buffer of ' +
                              str(float(Cfg.BUFFERDIST)) + ' map units.\n')
            start_time = time.clock()

            extentBoxList = npy.zeros((0,5),dtype='float32')
            boxCoords = lu.get_extent_box_coords()
            extentBoxList = npy.append(extentBoxList,boxCoords,axis=0)
            extentBoxList[0,0] = 0

            boundingCirclePointArray  = npy.zeros((0,5),dtype='float32')
            circlePointData=lu.get_bounding_circle_data(extentBoxList, 0,
                                                        0, Cfg.BUFFERDIST)

            lu.make_points(Cfg.SCRATCHDIR, circlePointData, Cfg.BNDCIRCEN)
            Cfg.gp.defineprojection(Cfg.BNDCIRCEN, SR)
            if Cfg.gp.Exists(Cfg.BNDCIR):
                Cfg.gp.delete_management(Cfg.BNDCIR)
            Cfg.gp.buffer_analysis(Cfg.BNDCIRCEN, Cfg.BNDCIR, "radius")
            Cfg.gp.defineprojection(Cfg.BNDCIR, SR)

            boundResis = "boundResis"
            Cfg.gp.addmessage('Extracting raster....')

             # FIXME: wishlist- would extract by circle be faster?
            count = 0
            statement = ('Cfg.gp.ExtractByMask_sa(Cfg.RESRAST, Cfg.BNDCIR, '
                         'boundResis)')
            while True:
                try: exec statement
                except:
                    count,tryAgain = lu.hiccup_test(count,statement)
                    if not tryAgain: exec statement
                else: break
            Cfg.gp.addmessage('\nReduced resistance raster extracted using '
                              'bounding circle.')
            start_time = lu.elapsed_time(start_time)

        else: #if not using bounding circles, just go with resistance raster.
            boundResis = Cfg.RESRAST

        # ---------------------------------------------------------------------
        # Rasterize core areas to speed cost distance calcs
        # lu.dashline(1)
        Cfg.gp.addmessage("Creating core area raster.")
        s3core_ras="s3core_ras"
        Cfg.gp.SelectLayerByAttribute(Cfg.FCORES, "CLEAR_SELECTION")
        Cfg.gp.CellSize = Cfg.gp.Describe(boundResis).MeanCellHeight
        if Cfg.gp.exists(s3core_ras):
            Cfg.gp.delete_management(s3core_ras)
        Cfg.gp.extent = Cfg.gp.Describe(boundResis).extent
        count = 0
        statement = ('Cfg.gp.FeatureToRaster_conversion(Cfg.FCORES, '
                     'Cfg.COREFN, s3core_ras, Cfg.gp.Cellsize)')
        while True:
            try: exec statement
            except:
                count,tryAgain = lu.hiccup_test(count,statement)
                if not tryAgain: exec statement
            else: break
        Cfg.gp.extent = "MINOF"

        #----------------------------------------------------------------------
        # Loop through cores, do cwd calcs for each
        Cfg.gp.addmessage("\nStarting cost distance calculations.\n")
        lcpLoop = 0
        for x in range(len(coresToMap)):
            startTime1 = time.clock()
            # This is the focal core we're running cwd out from
            sourceCore = int(coresToMap[x])

            # Get target cores based on linktable with reinstated links
            # (we temporarily disable them below by adding 100)
            linkTableTemp = linkTable.copy()
            # reinstate temporarily disabled links
            rows = npy.where(linkTableTemp[:,Cfg.LTB_LINKTYPE] > 100)
            linkTableTemp[rows,Cfg.LTB_LINKTYPE] = (
                linkTableTemp[rows,Cfg.LTB_LINKTYPE] - 100)
            # get core areas to be connected to focal core
            targetCores = lu.get_core_targets(sourceCore, linkTableTemp)
            del linkTableTemp

            if len(targetCores)>0:
                lu.dashline(0)
                Cfg.gp.addmessage('Target core areas for core area #' +
                                  str(sourceCore) + ' = ' + str(targetCores))

                # -------------------------------------------------------------
                # Create BOUNDING FEATURE to limit extent of cost distance
                # calculations-This is a set of circles encompassing core areas
                # we'll be connecting each core area to.
                if Cfg.BUFFERDIST is not None:
                    # fixme: move outside of loop   # new circle
                    Cfg.gp.MakeFeatureLayer(
                        path.join(Cfg.gp.workspace, BNDCIRS),
                        "fGlobalBoundingFeat")

                    start_time = time.clock()
                    # loop through targets and get bounding circles that
                    # contain focal core and target cores
                    Cfg.gp.AddMessage("\nAdding up bounding circles for source"
                                      " core " + str(sourceCore))
                    Cfg.gp.SelectLayerByAttribute("fGlobalBoundingFeat",
                                                  "CLEAR_SELECTION")
                    for i in range(len(targetCores)):
                        # run thru circleList, find link that core pair
                        # corresponds to.
                        if sourceCore < targetCores[i]:
                            corex = sourceCore
                            corey = targetCores[i]
                        else:
                            corey = sourceCore
                            corex = targetCores[i]

                        cores_x_y = str(int(corex))+'_'+str(int(corey))
                        field = "cores_x_y"
                        # fixme: need to check for case where link is not found.
                        Cfg.gp.SelectLayerByAttribute(
                            "fGlobalBoundingFeat", "ADD_TO_SELECTION", field +
                            " = '" + cores_x_y + "'")

                    if Cfg.gp.Exists(BNDFC):
                        Cfg.gp.delete_management(BNDFC) # fixme: necessary?
                    # fixme: may not be needed- can we just clip raster
                    # using selected?
                    Cfg.gp.CopyFeatures_management("fGlobalBoundingFeat",
                                                   BNDFC)

                    # Clip out bounded area of resistance raster for cwd
                    # calculations from focal core
                    bResistance = "bResistance"
                    count = 0
                    statement = ('Cfg.gp.ExtractByMask_sa(boundResis, '
                                 'BNDFC, bResistance)')
                    while True:
                        try: exec statement
                        except:
                            count,tryAgain = lu.hiccup_test(count, statement)
                            if not tryAgain: exec statement
                        else: break
                    Cfg.gp.addmessage('Successfully extracted a reduced '
                                      ' resistance raster using')
                    Cfg.gp.addmessage('bounding circles plus a buffer of ' +
                                      str(float(Cfg.BUFFERDIST)) + ' map '
                                      'units.')
                    start_time = lu.elapsed_time(start_time)
                else:
                    bResistance = boundResis

                # ---------------------------------------------------------
                # CWD Calculations
                outDistanceRaster = lu.get_cwd_path(sourceCore)
                start_time = time.clock()

                # Create raster that just has source core in it
                # Note: this seems faster than setnull with LI grid.
                expression = ("con(" + s3core_ras + "== " +
                              str(int(sourceCore)) + ",1)")
                SRCRASTER = 'source'
                count = 0
                statement = ('Cfg.gp.SingleOutputMapAlgebra_sa(expression, '
                             'SRCRASTER)')
                while True:
                    try: exec statement
                    except:
                        count, tryAgain = lu.hiccup_test(count, statement)
                        if not tryAgain: exec statement
                    else: break

                # Cost distance raster creation
                count = 0
                statement = ('Cfg.gp.CostDistance_sa(SRCRASTER, bResistance, '
                             'outDistanceRaster, Cfg.TMAXCWDIST, "BACK")')
                while True:
                    try: exec statement
                    except:
                        count, tryAgain = lu.hiccup_test(count, statement)
                        if not tryAgain: exec statement
                    else: break
                Cfg.gp.addmessage('Cost distances for source core ' +
                              str(int(sourceCore)) + ' calculated.')
                start_time = lu.elapsed_time(start_time)

                # Extract cost distances from source core to target cores
                # Fixme: there will be redundant calls to b-a when already
                # done a-b
                Cfg.gp.addmessage('Getting least cost distances from source '
                                  'core #' + str(int(sourceCore)) + ' to ' +
                                  str(len(targetCores)) + ' potential targets')
                count = 0
                statement = ('Cfg.gp.zonalstatisticsastable_sa(s3core_ras, '
                             '"VALUE", outDistanceRaster, ZNSTATS)')
                while True:
                    try: exec statement
                    except:
                        count,tryAgain = lu.hiccup_test(count,statement)
                        if not tryAgain: exec statement
                    else: break

                tableRows = Cfg.gp.searchcursor(ZNSTATS)
                tableRow = tableRows.Next()
                while tableRow:
                    if tableRow.Value > sourceCore:
                        link = lu.get_links_from_core_pairs(linkTable,
                                                            sourceCore,
                                                            tableRow.Value)
                        if linkTable[link,Cfg.LTB_LINKTYPE] > 0: # valid link
                            linkTable[link,Cfg.LTB_CWDIST] = tableRow.Min
                            if Cfg.MAXCOSTDIST is not None:
                                if tableRow.Min > Cfg.MAXCOSTDIST:
                                     # Disable link, it's too long
                                    linkTable[link,Cfg.LTB_LINKTYPE] = Cfg.LT_TLLC
                            if Cfg.MINCOSTDIST is not None:
                                if tableRow.Min < Cfg.MINCOSTDIST:
                                    # Disable link, it's too short
                                    linkTable[link,Cfg.LTB_LINKTYPE] = Cfg.LT_TSLC
                    tableRow = tableRows.next()
                del tableRow, tableRows
                start_time = lu.elapsed_time(start_time)

                # ---------------------------------------------------------
                # Check for intermediate cores AND map LCP lines
                for y in range(len(targetCores)):
                    targetCore = targetCores[y]
                    rows = lu.get_links_from_core_pairs(linkTable, sourceCore,
                                                        targetCore)
                    # Map all links for which above code successfully extracted
                    #  cwds in above code
                    if (linkTable[rows[0],Cfg.LTB_LINKTYPE] > 0 and
                        linkTable[rows[0],Cfg.LTB_LINKTYPE] < 100 and
                        linkTable[rows[0],Cfg.LTB_CWDIST] != -1):
                        # Flag so that we only evaluate this pair once
                        linkTable[rows,Cfg.LTB_LINKTYPE] = (linkTable
                                                       [rows,Cfg.LTB_LINKTYPE]
                                                       + 100)
                        # Create raster that just has target core in it
                        expression = ("con(" + s3core_ras + "== " +
                                      str(int(targetCore)) + ",1)")
                        TARGETRASTER = 'targ'
                        count = 0
                        statement = (
                            'Cfg.gp.SingleOutputMapAlgebra_sa(expression,'
                            ' TARGETRASTER)')
                        while True:
                            try: exec statement
                            except:
                                count, tryAgain = lu.hiccup_test(count,
                                                                 statement)
                                if not tryAgain: exec statement
                            else: break

                        try:
                            # Cost path allows us to map the least cost path
                            # between source and target
                            count = 0
                            statement = ('Cfg.gp.CostPath_sa(TARGETRASTER, '
                                        'outDistanceRaster, "BACK",  "lcp", '
                                        '"BEST_SINGLE", "")')
                            try:
                                exec statement
                            except:
                                exec statement
                        except:
                            link = lu.get_links_from_core_pairs(linkTable,
                                                                sourceCore,
                                                                targetCore)
                            # Not picked up by above cwd calc code
                            if (Cfg.MAXCOSTDIST is not None
                                and linkTable[link,Cfg.LTB_CWDIST] == -1):
                                Cfg.gp.addmessage('Cost path failed- should '
                                                  'not have gotten to this '
                                                  'point?')
                                continue
                            else:
                                lu.dashline(1)
                                msg = ("Error in COST PATH function for link "
                                       "#" + str(int(link)) +
                                       ".\nPlease report error.\n")
                                Cfg.gp.AddError(msg)
                                exit(0)

                        # fixme: may be fastest to not do selection, do
                        # EXTRACTBYMASK, getvaluelist, use code snippet at end
                        # of file to discard src and target values. Still this
                        # is fast- 13 sec for LI data...But I'm not very
                        # comfortable using failed coreMin as our test....
                        if Cfg.S3DROPLCCSic:
                            # -------------------------------------------------
                            # Drop links where lcp passes through intermediate
                            # core area. Method below is faster than valuelist
                            # method because of soma in valuelist method.
                            Cfg.gp.SelectLayerByAttribute(Cfg.FCORES,
                                                      "NEW_SELECTION",
                                                      Cfg.COREFN + ' <> ' +
                                                      str(int(targetCore)) +
                                                      ' AND ' + Cfg.COREFN +
                                                      ' <> ' +
                                                      str(int(sourceCore)))
                            count = 0
                            statement = ('Cfg.gp.zonalstatisticsastable('
                                         'Cfg.FCORES, Cfg.COREFN, "lcp", '
                                         'ZNSTATS, "DATA", "MINIMUM")')
                            while True:
                                try: exec statement
                                except:
                                    count, tryAgain = lu.hiccup_test(count,
                                                                     statement)
                                    if not tryAgain: exec statement
                                else: break
                            rows = lu.get_links_from_core_pairs(linkTable,
                                                                sourceCore,
                                                                targetCore)
                            coreMin = lu.get_zonal_minimum(ZNSTATS)
                            #  Found a valid value, indicating overlap
                            if coreMin:
                                lu.dashline()
                                Cfg.gp.addmessage(
                                    "Found an intermediate core in the "
                                    "least-cost path between cores " +
                                    str(int(sourceCore)) + " and " +
                                    str(int(targetCore)) + ".  The corridor "
                                    "will be removed.\n")
                                # disable link
                                linkTable[rows,Cfg.LTB_LINKTYPE] = Cfg.LT_INT
                            # -------------------------------------------------

                        # Create lcp shapefile.  lcploop just keeps track of
                        # whether this is first time function is called.
                        lcpLoop = lu.create_lcp_shapefile(linkTable,
                                                          sourceCore,
                                                          targetCore,
                                                          lcpLoop, SR)

                endTime = time.clock()
                processTime = round((endTime - start_time), 2)
                Cfg.gp.addmessage('Intermediate core checks and LCP shapefiles'
                                  'for core #' + str(sourceCore) + ' took ' +
                                  str(processTime) + ' seconds.')
                # -----------------------------------------------------------


                Cfg.gp.addmessage('Done with all calculations for core #' +
                                  str(sourceCore) + '.')
                start_time = lu.elapsed_time(startTime1)
            # -----------------------------------------------------------

        # reinstate temporarily disabled links
        rows = npy.where(linkTable[:,Cfg.LTB_LINKTYPE] > 100)
        linkTable[rows,Cfg.LTB_LINKTYPE] = (linkTable[rows,Cfg.LTB_LINKTYPE] -
                                            100)

        # Drop links that are too long
        disableLeastCostNoVal = True
        linkTable,numDroppedLinks = lu.drop_links(linkTable, Cfg.MAXEUCDIST,
                                               Cfg.MINEUCDIST, Cfg.MAXCOSTDIST,
                                               Cfg.MINCOSTDIST,
                                               disableLeastCostNoVal)

        # Write link table file
        outlinkTableFile = lu.get_this_step_link_table(step=3)
        Cfg.gp.addmessage('Updating ' + outlinkTableFile)
        lu.write_link_table(linkTable, outlinkTableFile)
        linkTableLogFile = path.join(Cfg.LOGDIR, "linkTable_s3.csv")
        lu.write_link_table(linkTable, linkTableLogFile)

        # lu.dashline()
        start_time = time.clock()
        Cfg.gp.addmessage('Creating shapefiles with linework for links...')
        lu.write_link_maps(outlinkTableFile, step=3)
        start_time = lu.elapsed_time(start_time)

        Cfg.gp.addmessage('\nIndividual cost-weighted distance layers written '
                          'to "cwd" directory. \n')
        Cfg.gp.addmessage(
            outlinkTableFile + '\n updated with cost-weighted distances between '
            'core areas.')

    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        Cfg.gp.addmessage('****Failed in step 3. Details follow.****')
        lu.raise_python_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        Cfg.gp.addmessage('****Failed in step 3. Details follow.****')
        lu.raise_python_error(__filename__)

    return