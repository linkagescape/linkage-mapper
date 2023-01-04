#!/usr/bin/env python2.5
# Authors: Brad McRae and Darren Kavanagh

"""Step 3: Calculate cost-weighted distances.

Calculates cost-weighted distances from each core area.
Uses bounding circles around source and target cores to limit
extent of cwd calculations and speed computation.

"""


import os.path as path
import time
import numpy as npy

from lm_config import tool_env as cfg
import lm_util as lu

_SCRIPT_NAME = "s3_calcCwds.py"

try:
    import arcpy
    from arcpy.sa import *
    arcpy.CheckOutExtension("spatial")
    gp = arcpy.gp
    arcgisscripting = arcpy
    tif = ''
except:
    arcpy = False
    import arcgisscripting
    gp = cfg.gp
    tif = ''

gprint = lu.gprint


def write_cores_to_map(x, coresToMap):
    """ Save core list at start of loop to allow a run to be re-started
        if it fails.

    """
    try:
        coreListFile = path.join(cfg.DATAPASSDIR, "temp_cores_to_map.csv")
        outFile = open(coreListFile, "w")
        outFile.write("#INDEX of last core being processed:\n")
        outFile.write(str(int(x)))
        outFile.write("\n")
        outFile.write("#Full list of IDs of cores to be processed:\n")
        for coreToMap in range(len(coresToMap)):
            outFile.write(str(int(coresToMap[coreToMap])))
            outFile.write("\n")
        outFile.close()

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        # gprint('****Failed in step 3. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)



def STEP3_calc_cwds():
    """Calculates cost-weighted distances from each core area.
    Uses bounding circles around source and target cores to limit
    extent of cwd calculations and speed computation.

    """
    try:
        lu.dashline(1)
        gprint('Running script ' + _SCRIPT_NAME)
        lu.dashline(0)

        # Super secret setting to re-start failed run.  Enter 'RESTART' as the
        # Name of the pairwise distance table in step 2, and uncheck step 2.
        # We can eventually place this in a .ini file.
        rerun = False
        if cfg.S2EUCDISTFILE != None:
            if cfg.S2EUCDISTFILE.lower() == "restart":
                rerun = True

        # if cfg.TMAXCWDIST is None:
           	# gprint('NOT using a maximum cost-weighted distance.')
        # else:
            # gprint('Max cost-weighted distance for CWD calcs set '
                              # 'to ' + str(cfg.TMAXCWDIST) + '\n')

                              
        if (cfg.BUFFERDIST) is not None:
            gprint('Bounding circles plus a buffer of ' +
                              str(float(cfg.BUFFERDIST)) + ' map units will '
                              'be used \n to limit extent of cost distance '
                              'calculations.')
        elif cfg.TOOL <> cfg.TOOL_CC:
            gprint('NOT using bounding circles in cost distance '
                              'calculations.')

        # set the analysis extent and cell size
        # So we don't extract rasters that go beyond extent of original raster
        if arcpy:
            arcpy.env.cellSize = cfg.RESRAST
            arcpy.env.extent="MINOF"
        else:
            gp.cellSize = gp.Describe(cfg.RESRAST).MeanCellHeight
            gp.Extent = "MINOF"
        gp.mask = cfg.RESRAST
        if arcpy:
            arcpy.env.overwriteOutput = True
            arcpy.env.workspace = cfg.SCRATCHDIR
            arcpy.env.scratchWorkspace = cfg.ARCSCRATCHDIR
        else:
            gp.OverwriteOutput = True
            gp.workspace = cfg.SCRATCHDIR
            gp.scratchWorkspace = cfg.ARCSCRATCHDIR

        # Load linkTable (created in previous script)
        linkTableFile = lu.get_prev_step_link_table(step=3)
        linkTable = lu.load_link_table(linkTableFile)
        lu.report_links(linkTable)

        # Identify cores to map from LinkTable
        coresToMap = npy.unique(linkTable[:, cfg.LTB_CORE1:cfg.LTB_CORE2 + 1])
        numCoresToMap = len(coresToMap)

        if numCoresToMap < 3:
            # No need to check for intermediate cores, because there aren't any
            cfg.S3DROPLCCSic = False
        else:
            cfg.S3DROPLCCSic = cfg.S3DROPLCCS
        gprint('\nNumber of core areas to connect: ' +
                          str(numCoresToMap))

        if rerun:
            # If picking up a failed run, make sure needed files are there
            lu.dashline(1)
            gprint ('\n****** RESTART MODE ENABLED ******\n')
            gprint ('**** NOTE: This mode picks up step 3 where a\n'
                    'previous run left off due to a crash or user\n'
                    'abort.  It assumes you are using the same input\n'
                    'data used in the terminated run.\n\n')
            lu.warn('IMPORTANT: Your LCP and stick feature classes\n'
                    'will LOSE LCPs that were already created, but\n'
                    'your final raster corridor map should be complete.\n')
                    
            lu.dashline(0)
            lu.snooze(10)
            savedLinkTableFile = path.join(cfg.DATAPASSDIR,
                                           "temp_linkTable_s3_partial.csv")
            coreListFile = path.join(cfg.DATAPASSDIR, "temp_cores_to_map.csv")

            if not path.exists(savedLinkTableFile) or not path.exists(
                                                          coreListFile):

                gprint('No partial results file found from previous '
                       'stopped run. Starting run from beginning.\n')
                lu.dashline(0)
                rerun = False

        # If picking up a failed run, use old folders
        if not rerun:
            startIndex = 0
            if cfg.TOOL <> cfg.TOOL_CC:
                lu.make_cwd_paths(max(coresToMap)) # Set up cwd directories

        # make a feature layer for input cores to select from
        gp.MakeFeatureLayer(cfg.COREFC, cfg.FCORES)

        # Drop links that are too long
        gprint('\nChecking for corridors that are too long to map.')
        DISABLE_LEAST_COST_NO_VAL = False
        linkTable,numDroppedLinks = lu.drop_links(linkTable, cfg.MAXEUCDIST, 0,
                                                  cfg.MAXCOSTDIST, 0,
                                                  DISABLE_LEAST_COST_NO_VAL)
        # ------------------------------------------------------------------
        # Bounding boxes
        if (cfg.BUFFERDIST) is not None:
            # create bounding boxes around cores
            start_time = time.clock()
            # lu.dashline(1)
            gprint('Calculating bounding boxes for core areas.')
            extentBoxList = npy.zeros((0,5), dtype='float32')
            for x in range(len(coresToMap)):
                core = coresToMap[x]
                boxCoords = lu.get_extent_box_coords(core)
                extentBoxList = npy.append(extentBoxList, boxCoords, axis=0)
            gprint('\nDone calculating bounding boxes.')
            start_time = lu.elapsed_time(start_time)
            # lu.dashline()

        # Bounding circle code
        if cfg.BUFFERDIST is not None:
            # Make a set of circles encompassing core areas we'll be connecting
            start_time = time.clock()
            gprint('Calculating bounding circles around potential'
                          ' corridors.')

            # x y corex corey radius- stores data for bounding circle centroids
            boundingCirclePointArray  = npy.zeros((0,5), dtype='float32')

            circleList = npy.zeros((0,3), dtype='int32')

            numLinks = linkTable.shape[0]
            for x in range(0, numLinks):
                if ((linkTable[x,cfg.LTB_LINKTYPE] == cfg.LT_CORR) or
                    (linkTable[x,cfg.LTB_LINKTYPE] == cfg.LT_KEEP)):
                    # if it's a valid corridor link
                    linkId = int(linkTable[x,cfg.LTB_LINKID])
                    # fixme- this code is clumsy- can trim down
                    cores = npy.zeros((1,3), dtype='int32')
                    cores[0,:] = npy.sort([0, linkTable[x,cfg.LTB_CORE1],
                                      linkTable[x,cfg.LTB_CORE2]])
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
                            corex, corey, cfg.BUFFERDIST))
                        boundingCirclePointArray = (
                            npy.append(boundingCirclePointArray,
                            circlePointData, axis=0))
                        # keep track of which cores we draw bounding circles
                        # around
                        circleList = npy.append(circleList, cores, axis=0)

            gprint('\nCreating bounding circles using buffer '
                              'analysis.')

            dir, BNDCIRCENS = path.split(cfg.BNDCIRCENS)
            lu.make_points(cfg.SCRATCHDIR, boundingCirclePointArray,
                           BNDCIRCENS)
            lu.delete_data(cfg.BNDCIRS)
            gp.buffer_analysis(cfg.BNDCIRCENS, cfg.BNDCIRS, "radius")
            gp.deletefield (cfg.BNDCIRS, "BUFF_DIST")

            gprint('Successfully created bounding circles around '
                              'potential corridors using \na buffer of ' +
                              str(float(cfg.BUFFERDIST)) + ' map units.')
            start_time = lu.elapsed_time(start_time)

            gprint('Reducing global processing area using bounding '
                              'circle plus buffer of ' +
                              str(float(cfg.BUFFERDIST)) + ' map units.\n')


            extentBoxList = npy.zeros((0,5),dtype='float32')
            boxCoords = lu.get_extent_box_coords()
            extentBoxList = npy.append(extentBoxList,boxCoords,axis=0)
            extentBoxList[0,0] = 0

            boundingCirclePointArray  = npy.zeros((0,5),dtype='float32')
            circlePointData=lu.get_bounding_circle_data(extentBoxList, 0,
                                                        0, cfg.BUFFERDIST)

            dir, BNDCIRCEN = path.split(cfg.BNDCIRCEN)
            lu.make_points(cfg.SCRATCHDIR, circlePointData, BNDCIRCEN)
            lu.delete_data(cfg.BNDCIR)
            gp.buffer_analysis(cfg.BNDCIRCEN, cfg.BNDCIR, "radius")

            gprint('Extracting raster....')
            cfg.BOUNDRESIS = cfg.BOUNDRESIS + tif
            lu.delete_data(cfg.BOUNDRESIS)
            count = 0
            statement = (
                'gp.ExtractByMask_sa(cfg.RESRAST, cfg.BNDCIR, cfg.BOUNDRESIS)')
            while True:
                try:
                    exec statement
                    randomerror()
                except:
                    count,tryAgain = lu.retry_arc_error(count,statement)
                    if not tryAgain: exec statement
                else: break
            gprint('\nReduced resistance raster extracted using '
                              'bounding circle.')

        else: #if not using bounding circles, just go with resistance raster.
            cfg.BOUNDRESIS = cfg.RESRAST

        # ---------------------------------------------------------------------
        # Rasterize core areas to speed cost distance calcs
        # lu.dashline(1)
        gprint("Creating core area raster.")

        gp.SelectLayerByAttribute(cfg.FCORES, "CLEAR_SELECTION")

        if arcpy:
            arcpy.env.cellSize = cfg.BOUNDRESIS
            arcpy.env.extent = cfg.BOUNDRESIS
        else:
            gp.cellSize = gp.Describe(cfg.BOUNDRESIS).MeanCellHeight
            gp.extent = gp.Describe(cfg.BOUNDRESIS).extent

        if rerun:
            # saved linktable replaces the one now in memory
            linkTable = lu.load_link_table(savedLinkTableFile)
            coresToMapSaved = npy.loadtxt(coreListFile, dtype='Float64',
                                          comments='#', delimiter=',')
            startIndex = coresToMapSaved[0] # Index of core where we left off
            del coresToMapSaved
            gprint ('\n****** Re-starting run at core area number '
                    + str(int(coresToMap[startIndex]))+ ' ******\n')
            lu.dashline(0)

        if arcpy:
            arcpy.env.extent = "MINOF"
        else:
            gp.extent = "MINOF"

        #----------------------------------------------------------------------
        # Loop through cores, do cwd calcs for each
        if cfg.TOOL == cfg.TOOL_CC:
            gprint("\nMapping least-cost paths.\n")
        else:
            gprint("\nStarting cost distance calculations.\n")
        lcpLoop = 0
        failures = 0
        x = startIndex
        endIndex = len(coresToMap)
        linkTableMod = linkTable.copy()
        while x < endIndex:
            startTime1 = time.clock()
            # Modification of linkTable in function was causing problems. so
            # make a copy:
            linkTablePassed = linkTableMod.copy()

            (linkTableReturned, failures, lcpLoop) = do_cwd_calcs(x,
                        linkTablePassed, coresToMap, lcpLoop, failures)
            if failures == 0:
                # If iteration was successful, continue with next core
                linkTableMod = linkTableReturned
                sourceCore = int(coresToMap[x])
                gprint('Done with all calculations for core ID #' +
                        str(sourceCore) + '. ' + str(int(x + 1)) + ' of ' +
                        str(endIndex) + ' cores have been processed.')
                start_time = lu.elapsed_time(startTime1)

                outlinkTableFile = path.join(cfg.DATAPASSDIR,
                                             "temp_linkTable_s3_partial.csv")
                lu.write_link_table(linkTableMod, outlinkTableFile)
                # Increment  loop counter
                x = x + 1
            else:
                # If iteration failed, try again after a wait period
                delay_restart(failures)
        #----------------------------------------------------------------------

        linkTable = linkTableMod

        # reinstate temporarily disabled links
        rows = npy.where(linkTable[:,cfg.LTB_LINKTYPE] > 1000)
        linkTable[rows,cfg.LTB_LINKTYPE] = (linkTable[rows,cfg.LTB_LINKTYPE] -
                                            1000)

        # Drop links that are too long
        DISABLE_LEAST_COST_NO_VAL = True
        linkTable,numDroppedLinks = lu.drop_links(linkTable, cfg.MAXEUCDIST,
                                               cfg.MINEUCDIST, cfg.MAXCOSTDIST,
                                               cfg.MINCOSTDIST,
                                               DISABLE_LEAST_COST_NO_VAL)

        # Write link table file
        outlinkTableFile = lu.get_this_step_link_table(step=3)
        gprint('Updating ' + outlinkTableFile)
        lu.write_link_table(linkTable, outlinkTableFile)
        linkTableLogFile = path.join(cfg.LOGDIR, "linkTable_s3.csv")
        lu.write_link_table(linkTable, linkTableLogFile)

        start_time = time.clock()
        gprint('Creating shapefiles with linework for links...')
        try:
            lu.write_link_maps(outlinkTableFile, step=3)
        except:
            lu.write_link_maps(outlinkTableFile, step=3)
        start_time = lu.elapsed_time(start_time)

        gprint('\nIndividual cost-weighted distance layers written '
                          'to "cwd" directory. \n')
        gprint(outlinkTableFile +
                '\n updated with cost-weighted distances between core areas.')

        #Clean up temporary files for restart code
        tempFile = path.join(cfg.DATAPASSDIR, "temp_cores_to_map.csv")
        lu.delete_file(tempFile)
        tempFile = path.join(cfg.DATAPASSDIR, "temp_linkTable_s3_partial.csv")
        lu.delete_file(tempFile)

        # Check if climate tool is calling linkage mapper
        if cfg.TOOL == cfg.TOOL_CC:
            coreList = npy.unique(linkTable[:, cfg.LTB_CORE1:cfg.LTB_CORE2 + 1])
            for core in coreList:
                cwdRaster = lu.get_cwd_path(int(core))
                back_rast = cwdRaster.replace("cwd_", "back_")        
                lu.delete_data(back_rast)
        

    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 3. Details follow.****')
        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 3. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)

    return



def do_cwd_calcs(x, linkTable, coresToMap, lcpLoop, failures):
    try:
        # This is the focal core area we're running cwd out from
        sourceCore = int(coresToMap[x])

        # Create temporary scratch directory just this focal core
        coreDir = path.join(cfg.SCRATCHDIR, 'core' + str(sourceCore))
        lu.delete_dir(coreDir)
        lu.create_dir(coreDir)

        if arcpy:
            gp = arcpy.gp
            arcpy.env.workspace = coreDir
            arcpy.env.scratchWorkspace = cfg.ARCSCRATCHDIR
            arcpy.env.overwriteOutput = True
            arcpy.env.extent = "MINOF"
        else:
            gp = cfg.gp
            gp.workspace = coreDir
            gp.scratchWorkspace = cfg.ARCSCRATCHDIR
            gp.OverwriteOutput = True
            gp.Extent = "MINOF"

        write_cores_to_map(x, coresToMap)

        # Get target cores based on linktable with reinstated links
        # (we temporarily disable them below by adding 1000)
        linkTableTemp = linkTable.copy()
        # reinstate temporarily disabled links
        rows = npy.where(linkTableTemp[:,cfg.LTB_LINKTYPE] > 1000)
        linkTableTemp[rows,cfg.LTB_LINKTYPE] = (
            linkTableTemp[rows,cfg.LTB_LINKTYPE] - 1000)
        del rows

        # get core areas to be connected to focal core
        targetCores = lu.get_core_targets(sourceCore, linkTableTemp)
        # gprint( str(sourceCore))
        # gprint(str(linkTableTemp.astype('int32')))
        # gprint('targets'+str(targetCores))
        del linkTableTemp

        if len(targetCores)==0:
            # Nothing to do, so reset failure count and return.
            failures = 0
            return linkTable, failures, lcpLoop

        lu.dashline(0)
        gprint('Target core areas for core area #' +
                          str(sourceCore) + ' = ' + str(targetCores))

        # -------------------------------------------------------------
        # Create BOUNDING FEATURE to limit extent of cost distance
        # calculations-This is a set of circles encompassing core areas
        # we'll be connecting each core area to.
        if cfg.BUFFERDIST is not None:
            # fixme: move outside of loop   # new circle
            gp.MakeFeatureLayer(cfg.BNDCIRS,"fGlobalBoundingFeat")

            start_time = time.clock()
            # loop through targets and get bounding circles that
            # contain focal core and target cores
            # gprint("\nAdding up bounding circles for source"
                              # " core " + str(sourceCore))
            gp.SelectLayerByAttribute("fGlobalBoundingFeat",
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
                # fixme: need to check for case where link is not found
                gp.SelectLayerByAttribute(
                    "fGlobalBoundingFeat", "ADD_TO_SELECTION", field +
                    " = '" + cores_x_y + "'")

            lu.delete_data(path.join(coreDir,cfg.BNDFC))
            # fixme: may not be needed- can we just clip raster
            # using selected?
            gp.CopyFeatures_management("fGlobalBoundingFeat",
                                           cfg.BNDFC)

            # Clip out bounded area of resistance raster for cwd
            # calculations from focal core
            bResistance = path.join(coreDir,"bResistance") # Can't be tif-
                                                           # need STA for CWD
            lu.delete_data(bResistance)
            statement = (
                'gp.ExtractByMask_sa(cfg.BOUNDRESIS, cfg.BNDFC, bResistance)')
            try:
                exec statement
                randomerror()
            except:
                failures = lu.print_arcgis_failures(statement, failures)
                if failures < 20:
                    return None,failures,lcpLoop
                else: exec statement

        else:
            bResistance = cfg.BOUNDRESIS
        # ---------------------------------------------------------
        # CWD Calculations
        outDistanceRaster = lu.get_cwd_path(sourceCore)
        # Check if climate tool is calling linkage mapper
        if cfg.TOOL == cfg.TOOL_CC:
            back_rast = outDistanceRaster.replace("cwd_", "back_")
        else:
            back_rast = "BACK"
            lu.delete_data(path.join(coreDir, back_rast))
            lu.delete_data(outDistanceRaster)
            start_time = time.clock()

            # Create raster that just has source core in it
            # Note: this seems faster than setnull with LI grid.
            SRCRASTER = 'source' + tif
            lu.delete_data(path.join(coreDir,SRCRASTER))
            if arcpy:
                statement = ('conRaster = '
                             'Con(Raster(cfg.CORERAS) == int(sourceCore), 1);'
                             'conRaster.save(SRCRASTER)')
            else:
                expression = ("con(" + cfg.CORERAS + " == " +
                               str(int(sourceCore)) + ", 1)")
                statement = ('gp.SingleOutputMapAlgebra_sa'
                            '(expression, SRCRASTER)')

            try:
                exec statement
                randomerror()
            except:
                failures = lu.print_arcgis_failures(statement, failures)
                if failures < 20:
                    return None, failures, lcpLoop
                else: exec statement

            # Cost distance raster creation
            if arcpy:
                arcpy.env.extent = "MINOF"
            else:
                gp.Extent = "MINOF"

            lu.delete_data(path.join(coreDir,"BACK"))
            
            if arcpy:
                statement = ('outCostDist = CostDistance(SRCRASTER, '
                             'bResistance, cfg.TMAXCWDIST, back_rast);'
                             'outCostDist.save(outDistanceRaster)')
            else:
                statement = ('gp.CostDistance_sa(SRCRASTER, bResistance, '
                             'outDistanceRaster, cfg.TMAXCWDIST, back_rast)')
            try:
                exec statement
                randomerror()
            except:
                failures = lu.print_arcgis_failures(statement, failures)
                if failures < 20:
                    return None, failures, lcpLoop
                else:
                    exec statement

        start_time = time.clock()
        # Extract cost distances from source core to target cores
        # Fixme: there will be redundant calls to b-a when already
        # done a-b
        ZNSTATS = path.join(coreDir, "zonestats.dbf")
        lu.delete_data(ZNSTATS)
        #Fixme: zonalstatistics is returning integer values for minimum. Why???
        #Extra zonalstatistics code implemented later in script to correct
        #values.
        if arcpy:
            statement = ('outZSaT = ZonalStatisticsAsTable(cfg.CORERAS, '
                    '"VALUE", outDistanceRaster,ZNSTATS, "DATA", "MINIMUM")')
        else:
            statement = ('gp.zonalstatisticsastable_sa('
                      'cfg.CORERAS, "VALUE", outDistanceRaster, ZNSTATS)')
                      
        try:
            exec statement
            randomerror()
        except:
            failures = lu.print_arcgis_failures(statement, failures)
            if failures < 20:
                return None,failures,lcpLoop
            else:
                if cfg.TOOL == cfg.TOOL_CC:
                    msg = ('ERROR in Zonal Stats. Please restart ArcMap '
                        'and try again.')                
                else:
                    msg = ('ERROR in Zonal Stats. Restarting ArcMap '
                        'then restarting Linkage Mapper at step 3 usually\n'
                        'solves this one so please restart and try again.')

                lu.raise_error(msg)
        tableRows = gp.searchcursor(ZNSTATS)
        tableRow = tableRows.Next()
        while tableRow:
            if tableRow.Value > sourceCore:
                link = lu.get_links_from_core_pairs(linkTable,
                                                    sourceCore,
                                                    tableRow.Value)
                if linkTable[link,cfg.LTB_LINKTYPE] > 0: # valid link
                    linkTable[link,cfg.LTB_CWDIST] = tableRow.Min
                    if cfg.MAXCOSTDIST is not None:
                        if ((tableRow.Min > cfg.MAXCOSTDIST) and
                           (linkTable[link,cfg.LTB_LINKTYPE] != cfg.LT_KEEP)):
                             # Disable link, it's too long
                            linkTable[link,cfg.LTB_LINKTYPE] = cfg.LT_TLLC
                    if cfg.MINCOSTDIST is not None:
                        if (tableRow.Min < cfg.MINCOSTDIST and
                           (linkTable[link,cfg.LTB_LINKTYPE] != cfg.LT_KEEP)):
                            # Disable link, it's too short
                            linkTable[link,cfg.LTB_LINKTYPE] = cfg.LT_TSLC
            tableRow = tableRows.next()
        del tableRow, tableRows
        #start_time = lu.elapsed_time(start_time)

        # ---------------------------------------------------------
        # Check for intermediate cores AND map LCP lines
        for y in range(0,len(targetCores)):
            targetCore = targetCores[y]
            rows = lu.get_links_from_core_pairs(linkTable, sourceCore,
                                                targetCore)
            # Map all links for which we successfully extracted
            #  cwds in above code
            link = rows[0]
            if (linkTable[link,cfg.LTB_LINKTYPE] > 0 and
                linkTable[link,cfg.LTB_LINKTYPE] < 1000 and
                linkTable[link,cfg.LTB_CWDIST] != -1):
                # Flag so that we only evaluate this pair once
                linkTable[rows,cfg.LTB_LINKTYPE] = (linkTable
                                               [rows,cfg.LTB_LINKTYPE]
                                               + 1000)

                # Create raster that just has target core in it
                TARGETRASTER = 'targ' + tif
                lu.delete_data(path.join(coreDir,TARGETRASTER))
                try:
                    if arcpy:
                        # For climate corridors, errors occur when core raster
                        # overlaps null values in cwd rasters
                        statement = ('conRaster = Con(IsNull(outDistanceRaster'
                            '), Int(outDistanceRaster), Con(Raster'
                            '(cfg.CORERAS) == int(targetCore), 1));'
                            'conRaster.save(TARGETRASTER)') 
                        # statement = ('conRaster = Con(Raster('
                                    # 'cfg.CORERAS) == int(targetCore), 1);'
                                    # 'conRaster.save(TARGETRASTER)')

                    else:
                        expression = ("con(" + cfg.CORERAS + " == " +
                        str(int(targetCore)) + ",1)")
                        statement = ('gp.SingleOutputMapAlgebra_sa(expression,'
                                     ' TARGETRASTER)')
                    exec statement
                    randomerror()
                except:
                    failures = lu.print_arcgis_failures(statement, failures)
                    if failures < 20:
                        return None,failures,lcpLoop
                    else: exec statement
                # Execute ZonalStatistics to get more precise cw distance if
                # arc rounded it earlier (not critical, hence the try/pass)
                if (linkTable[link,cfg.LTB_CWDIST] ==
                                int(linkTable[link,cfg.LTB_CWDIST])):
                    try:
                        zonalRas = path.join(coreDir,'zonal')
                        gp.ZonalStatistics_sa(TARGETRASTER, "VALUE",
                            outDistanceRaster, zonalRas, "MINIMUM", "DATA")
                        minObject = gp.GetRasterProperties_management(zonalRas,
                                                                "MINIMUM")
                        rasterMin = float(str(minObject.getOutput(0)))
                        linkTable[link,cfg.LTB_CWDIST] = rasterMin
                        lu.delete_data(zonalRas)
                    except:
                        pass
                # Cost path maps the least cost path
                # between source and target
                lcpRas = path.join(coreDir,"lcp" + tif)
                lu.delete_data(lcpRas)

                # Note: costpath (both gp and arcpy versions) uses GDAL.               
                if arcpy:
                    statement = ('outCostPath = CostPath(TARGETRASTER,'
                          'outDistanceRaster, back_rast, "BEST_SINGLE", ""); '
                          'outCostPath.save(lcpRas)')
                else:
                    statement = ('gp.CostPath_sa(TARGETRASTER, '
                                 'outDistanceRaster, back_rast, '
                                 'lcpRas, "BEST_SINGLE", "")')
                try:
                    exec statement                    
                    randomerror()
                except:
                    failures = lu.print_arcgis_failures(statement, failures)
                    if failures < 20:
                        return None,failures,lcpLoop
                    else:
                        lu.dashline(1)
                        gprint('\nCost path is failing for Link #'
                           + str(int(link)) + ' connecting core areas ' +
                            str(int(sourceCore)) + ' and ' +
                            str(int(targetCore)) + '\n.'
                            'Retrying one more time in 5 minutes.')
                        lu.snooze(300)
                        exec statement
                
                # fixme: may be fastest to not do selection, do
                # EXTRACTBYMASK, getvaluelist, use code snippet at end
                # of file to discard src and target values. Still this
                # is fast- 13 sec for LI data...But I'm not very
                # comfortable using failed coreMin as our test....
                if (cfg.S3DROPLCCSic and
                    (linkTable[link,cfg.LTB_LINKTYPE] != cfg.LT_KEEP) and
                    (linkTable[link,cfg.LTB_LINKTYPE] != cfg.LT_KEEP + 1000)):
                    # -------------------------------------------------
                    # Drop links where lcp passes through intermediate
                    # core area. Method below is faster than valuelist
                    # method because of soma in valuelist method.
                    # make a feature layer for input cores to select from
                    gp.MakeFeatureLayer(cfg.COREFC, cfg.FCORES)

                    gp.SelectLayerByAttribute(cfg.FCORES,
                                              "NEW_SELECTION",
                                              cfg.COREFN + ' <> ' +
                                              str(int(targetCore)) +
                                              ' AND ' + cfg.COREFN +
                                              ' <> ' +
                                              str(int(sourceCore)))

                    corePairRas = path.join(coreDir,"s3corepair"+ tif)
                    if arcpy:
                        arcpy.env.extent = cfg.BOUNDRESIS
                    else:
                        gp.extent = gp.Describe(cfg.BOUNDRESIS).extent


                    statement = ('gp.FeatureToRaster_conversion(cfg.FCORES, '
                                'cfg.COREFN, corePairRas, gp.cellSize)')
                    try:
                        exec statement
                        randomerror()
                    except:
                        failures = lu.print_arcgis_failures(statement,
                                                            failures)
                        if failures < 20:
                            return None,failures,lcpLoop
                        else: exec statement

                    #------------------------------------------
                    # Intermediate core test
                    try:
                        coreDetected = test_for_intermediate_core(coreDir,
                                                lcpRas, corePairRas)
                        randomerror()
                    except:
                        statement = 'test_for_intermediate_core'
                        failures = lu.print_arcgis_failures(statement,
                                                            failures)
                        if failures < 20:
                            return None,failures,lcpLoop
                        else:
                            coreDetected = test_for_intermediate_core(
                                        coreDir, lcpRas, corePairRas)

                    if coreDetected:
                        # lu.dashline()
                        gprint(
                            "Found an intermediate core in the "
                            "least-cost path between cores " +
                            str(int(sourceCore)) + " and " +
                            str(int(targetCore)) + ".  The corridor "
                            "will be removed.")
                        # disable link
                        rows = lu.get_links_from_core_pairs(linkTable,
                                                    sourceCore,
                                                    targetCore)
                        linkTable[rows,cfg.LTB_LINKTYPE] = cfg.LT_INT
                    #------------------------------------------

                # Create lcp shapefile.  lcploop just keeps track of
                # whether this is first time function is called.
                lcpLoop = lu.create_lcp_shapefile(coreDir, linkTable,
                                                  sourceCore, targetCore,
                                                  lcpLoop)

        # Made it through, so reset failure count and return.
        failures = 0
        lu.delete_dir(coreDir)
        # if cfg.TOOL == cfg.TOOL_CC:
            # lu.delete_data(back_rast)
        return linkTable, failures, lcpLoop

    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        # gprint('****Failed in step 3. Details follow.****')
        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        # gprint('****Failed in step 3. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)


def test_for_intermediate_core(workspace,lcpRas,corePairRas):
    """ Test if there is an intermediate core by seeing if least-cost
        path and remaining cores intersect

    """
    try:
        gp.workspace = workspace
        gp.OverwriteOutput = True
        if gp.exists("addRas"): #Can't use tif for getrasterproperties
            gp.delete_management("addRas")
        count = 0
        if arcpy:
            statement = ('outRas = Raster(lcpRas) + Raster(corePairRas); '
                         'outRas.save("addRas")')
        else:
            expression = (lcpRas + " + " + corePairRas)
            statement = ('gp.SingleOutputMapAlgebra_sa(expression, "addRas")')
        while True:
            try:
                exec statement
                randomerror()
            except:
                count,tryAgain = lu.retry_arc_error(count,statement)
                if not tryAgain: exec statement
            else: break

        # Test to see if raster has data
        if gp.GetRasterProperties("addRas", "ALLNODATA").getOutput(0) == "0":
            return True  # Data present and therefore overlap
        else:
            return False  # Empty and therefore no overlap

    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        # gprint('****Failed in step 3. Details follow.****')
        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        # gprint('****Failed in step 3. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)

def delay_restart(failures):
    gprint('That was try #' + str(failures) + ' of 20 for this core area.')
    if failures < 7:
        gprint('Restarting iteration in ' + str(10*failures) + ' seconds. ')
        lu.dashline(2)
        lu.snooze(10*failures)
    else:
        gprint('Restarting iteration in 5 minutes. ')
        lu.dashline(2)
        lu.snooze(300)


def test_for_intermediate_core_old_method(workspace,lcpRas,corePairRas):
    """ Zonal stats method to test for intermediate core test (discontinued)

    """
    ZNSTATS2 = path.join(cfg.SCRATCHDIR, "zonestats2.dbf")
    value = "VALUE"
    gp.ZonalStatisticsAsTable_sa(corePairRas, value, lcpRas, ZNSTATS2,
                                 "DATA", "MINIMUM")
    coreMin = lu.get_zonal_minimum(ZNSTATS2)
    if coreMin:
        return True
    else:
        return False


def randomerror():
    """ Used to test error recovery.

    """
    generateError = False # Set to True to create random errors
    if generateError:
        gprint('\n***Rolling dice for random error***')
        import random
        test = random.randrange(1, 30)
        if test == 2:
            gprint('Creating artificial ArcGIS error')
            gp.MosaicToNewRaster_management(
                            "rasterString","mosaicDir","mosFN", "",
                            "32_BIT_FLOAT", "gp.cellSize", "1", "MINIMUM",
                            "MATCH")
        elif test == 3:
            gprint('Creating artificial python error')
            artificialPythonError
    return

