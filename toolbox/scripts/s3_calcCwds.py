#!/usr/bin/env python2.5
# Authors: Brad McRae and Darren Kavanagh

"""Step 3: Calculate cost-weighted distances.

Calculates cost-weighted distances from each core area.
Uses bounding circles around source and target cores to limit
extent of cwd calculations and speed computation.

"""

__filename__ = "s3_calcCwds.py"
__version__ = "0.6.6"

import os.path as path
import shutil
import time

import arcgisscripting
import numpy as npy

from lm_config import Config as Cfg
import lm_util as lu
            
gp = Cfg.gp
gprint = gp.addmessage
   
BNDCIRCENS = "boundingCircleCenters.shp"
BNDCIRS = "boundingCircles.shp"
BNDFC = "boundingFeature.shp"
ZNSTATS = path.join(Cfg.SCRATCHDIR, "zonestats.dbf")

def write_cores_to_map(x, coresToMap):
    """ save core list at start of loop to allow a run to be re-started 
        if it fails.
        
    """
    try:
        coreListFile = path.join(Cfg.DATAPASSDIR, "temp_cores_to_map.csv")
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
        gprint('****Failed in step 3. Details follow.****')
        lu.raise_python_error(__filename__)

    

def STEP3_calc_cwds():
    """Calculates cost-weighted distances from each core area.
    Uses bounding circles around source and target cores to limit
    extent of cwd calculations and speed computation.

    """
    try:
            
        lu.dashline(1)
        gprint('Running script ' + __filename__)
        gp.scratchWorkspace = Cfg.SCRATCHDIR
        lu.dashline(0)

        # Super secret setting to re-start failed run.  Enter 'RESTART' as the
        # Name of the pairwise distance table in step 2, and uncheck step 2.
        # We can eventually place this in a .ini file.
        if Cfg.S2EUCDISTFILE == "RESTART": 
            rerun = True 
        else:
            rerun = False        
             
        if Cfg.TMAXCWDIST is None:
           	gprint('NOT using a maximum cost-weighted distance.')
        else:
            gprint('Max cost-weighted distance for CWD calcs set '
                              'to ' + str(Cfg.TMAXCWDIST) + '\n')

        # FIXME: because it's integer, it fills in 0 if not entered.
        if (Cfg.BUFFERDIST) is not None:
            gprint('Bounding circles plus a buffer of ' +
                              str(float(Cfg.BUFFERDIST)) + ' map units will '
                              'be used \n to limit extent of cost distance '
                              'calculations.')
        else:
            gprint('NOT using bounding circles in cost distance '
                              'calculations.')

        # set the analysis extent and cell size
        gp.CellSize = gp.Describe(Cfg.RESRAST).MeanCellHeight
        # So we don't extract rasters that go beyond extent of original raster
        gp.Extent = "MINOF"
        gp.mask = Cfg.RESRAST
        gp.Workspace = Cfg.SCRATCHDIR       

        # Load linkTable (created in previous script)
        linkTableFile = lu.get_prev_step_link_table(step=3)
        linkTable = lu.load_link_table(linkTableFile)
        lu.report_links(linkTable)

        rows,cols = npy.where(
            linkTable[:, Cfg.LTB_LINKTYPE:Cfg.LTB_LINKTYPE + 1] == Cfg.LT_CORR)
        coresToProcess = npy.unique(
            linkTable[:, Cfg.LTB_CORE1:Cfg.LTB_CORE2 + 1])
        maxCoreNum = max(coresToProcess)
        del rows, cols, coresToProcess

        
        # Identify cores to map from LinkTable
        rows,cols = npy.where(linkTable[:,Cfg.LTB_LINKTYPE:Cfg.LTB_LINKTYPE+1] 
                              == Cfg.LT_CORR)
        coresToMap = npy.unique(linkTable[:,Cfg.LTB_CORE1:Cfg.LTB_CORE2+1])
        numCoresToMap = len(coresToMap)
        del rows,cols
        if numCoresToMap < 3:
            # No need to check for intermediate cores, because there aren't any
            Cfg.S3DROPLCCSic = False
        else:
            Cfg.S3DROPLCCSic = Cfg.S3DROPLCCS

        gprint('Number of core areas to connect:' +
                          str(numCoresToMap))        
        
        if rerun == True: 
            # If picking up a failed run, make sure needed files are there
            lu.dashline(1)
            gprint ('\n**** RESTART MODE ENABLED ****\n')   
            gprint ('**** WARNING: This mode picks up step 3 where\n'
                    'a previous run left off due to a crash or\n' 
                    'user abort.  It assumes you are using exactly\n' 
                    'the same settings as the terminated run.****\n')
            lu.dashline(0)
            time.sleep(20)
            savedLinkTableFile = path.join(Cfg.DATAPASSDIR, 
                                           "temp_linkTable_s3_partial.csv")
            coreListFile = path.join(Cfg.DATAPASSDIR, "temp_cores_to_map.csv")

            if not path.exists(savedLinkTableFile) or not path.exists(
                                                          coreListFile):
                
                gprint('No partial results file found from previous '
                       'stopped run. Starting run from beginning.\n')
                lu.dashline(0)
                rerun = False
                

        if rerun == False: # If picking up a failed run, use old folders
            startIndex = 0
            
            #remove cwd directory
            lu.delete_dir(Cfg.CWDBASEDIR)
                
            # Set up cwd directories.
            # To keep there from being > 100 grids in any one directory,
            # outputs are written to:
            # cwd\cw for cores 1-99
            # cwd\cw1 for cores 100-199
            # etc.

            gprint("\nCreating cost-weighted distance output folders:")
            gprint(path.join(Cfg.CWDBASEDIR, Cfg.CWDSUBDIR_NM))
            gp.CreateFolder_management(path.dirname(Cfg.CWDBASEDIR),
                                           path.basename(Cfg.CWDBASEDIR))
            gp.CreateFolder_management(Cfg.CWDBASEDIR, Cfg.CWDSUBDIR_NM)
            if maxCoreNum > 100:
                maxDirCount = int(maxCoreNum/100)
                for dirCount in range(1, maxDirCount + 1):
                    ccwdir = Cfg.CWDSUBDIR_NM + str(dirCount)
                    gprint(ccwdir)
                    gp.CreateFolder_management(Cfg.CWDBASEDIR, ccwdir)
                        
        # make a feature layer for input cores to select from
        gp.MakeFeatureLayer(Cfg.COREFC, Cfg.FCORES)

        # Drop links that are too long
        gprint('\nChecking for corridors that are too long to map.')
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
            gprint('Calculating bounding boxes for core areas.')
            extentBoxList = npy.zeros((0,5), dtype='float32')
            pctDone = 0
            for x in range(len(coresToMap)):
                core = coresToMap[x]
                pctDone = lu.report_pct_done(x, len(coresToMap), pctDone)
                boxCoords = lu.get_extent_box_coords(core)
                extentBoxList = npy.append(extentBoxList, boxCoords, axis=0)
            gprint('\nDone calculating bounding boxes.')
            start_time = lu.elapsed_time(start_time)
            # lu.dashline()

        # Bounding circle code
        if Cfg.BUFFERDIST is not None:
            # Make a set of circles encompassing core areas we'll be connecting
            start_time = time.clock()
            gprint('Calculating bounding circles around potential'
                          ' corridors.')

            # x y corex corey radius- stores data for bounding circle centroids
            boundingCirclePointArray  = npy.zeros((0,5), dtype='float32')

            circleList = npy.zeros((0,3), dtype='int32')

            numLinks = linkTable.shape[0]
            pctDone = 0
            for x in range(0, numLinks):
                pctDone = lu.report_pct_done(x, numLinks, pctDone)
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
                            npy.append(boundingCirclePointArray,
                            circlePointData, axis=0))
                        # keep track of which cores we draw bounding circles
                        # around
                        circleList = npy.append(circleList, cores, axis=0)

            gprint('\nCreating bounding circles using buffer '
                              'analysis.')
            lu.make_points(Cfg.SCRATCHDIR, boundingCirclePointArray,
                           BNDCIRCENS)
            lu.delete_data(BNDCIRS)
            gp.buffer_analysis(BNDCIRCENS, BNDCIRS, "radius")
            gp.deletefield (BNDCIRS, "BUFF_DIST")

            gprint('Successfully created bounding circles around '
                              'potential corridors using \na buffer of ' +
                              str(float(Cfg.BUFFERDIST)) + ' map units.')
            start_time = lu.elapsed_time(start_time)

            gprint('Reducing global processing area using bounding '
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
            lu.delete_data(Cfg.BNDCIR)
            gp.buffer_analysis(Cfg.BNDCIRCEN, Cfg.BNDCIR, "radius")

            boundResis = "boundResis"
            gprint('Extracting raster....')

             # FIXME: wishlist- would extract by circle be faster?
            count = 0
            statement = ('gp.ExtractByMask_sa(Cfg.RESRAST, Cfg.BNDCIR, '
                         'boundResis)')
            while True:
                try: exec statement
                except:
                    count,tryAgain = lu.hiccup_test(count,statement)
                    if not tryAgain: exec statement
                else: break
            gprint('\nReduced resistance raster extracted using '
                              'bounding circle.')
            start_time = lu.elapsed_time(start_time)

        else: #if not using bounding circles, just go with resistance raster.
            boundResis = Cfg.RESRAST

        # ---------------------------------------------------------------------
        # Rasterize core areas to speed cost distance calcs
        # lu.dashline(1)
        gprint("Creating core area raster.")
        # core_rastmp="core_rastmp"
        s3core_ras="s3core_ras"
        gp.SelectLayerByAttribute(Cfg.FCORES, "CLEAR_SELECTION")
        gp.CellSize = gp.Describe(boundResis).MeanCellHeight
        lu.delete_data(s3core_ras)
        gp.extent = gp.Describe(boundResis).extent
        count = 0
        statement = ('gp.FeatureToRaster_conversion(Cfg.FCORES, '
                     'Cfg.COREFN, s3core_ras, gp.Cellsize)')
        while True:
            try: exec statement
            except:
                count,tryAgain = lu.hiccup_test(count,statement)
                if not tryAgain: exec statement
            else: break
        
        #Convert core raster to integer format  #Not implemented- we only allow integer inputs.
        # s3core_ras="s3core_ras"
        # if gp.exists(s3core_ras):
            # gp.delete_management(s3core_ras)           
        # gp.Int_sa(core_rastmp, s3core_ras)
        # gp.delete_management(core_rastmp)
        
                
        if rerun == True:
            # saved linktable replaces the one now in memory
            linkTable = lu.load_link_table(savedLinkTableFile) 
            coresToMapSaved = npy.loadtxt(coreListFile, dtype='Float64',
                                          comments='#', delimiter=',')        
            startIndex = coresToMapSaved[0] # Index of core where we left off
            del coresToMapSaved
            gprint ('\n****** Re-starting run at core area number ' 
                    + str(int(coresToMap[startIndex]))+ ' ******\n')
            lu.dashline(0)            
            
            
        gp.extent = "MINOF"

        #----------------------------------------------------------------------
        # Loop through cores, do cwd calcs for each
        gprint("\nStarting cost distance calculations.\n")
        lcpLoop = 0

        for x in range(startIndex, len(coresToMap)):              
            startTime1 = time.clock()    
            write_cores_to_map(x, coresToMap)             

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
                gprint('Target core areas for core area #' +
                                  str(sourceCore) + ' = ' + str(targetCores))

                # -------------------------------------------------------------
                # Create BOUNDING FEATURE to limit extent of cost distance
                # calculations-This is a set of circles encompassing core areas
                # we'll be connecting each core area to.
                if Cfg.BUFFERDIST is not None:
                    # fixme: move outside of loop   # new circle
                    gp.MakeFeatureLayer(
                        path.join(gp.workspace, BNDCIRS),
                        "fGlobalBoundingFeat")

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

                    lu.delete_data(BNDFC) # fixme: necessary?
                    # fixme: may not be needed- can we just clip raster
                    # using selected?
                    gp.CopyFeatures_management("fGlobalBoundingFeat",
                                                   BNDFC)

                    # Clip out bounded area of resistance raster for cwd
                    # calculations from focal core
                    bResistance = "bResistance"
                    count = 0
                    statement = ('gp.ExtractByMask_sa(boundResis, '
                                 'BNDFC, bResistance)')
                    while True:
                        try: exec statement
                        except:
                            count,tryAgain = lu.hiccup_test(count, statement)
                            if not tryAgain: exec statement
                        else: break
                    # gprint('Successfully extracted a reduced '
                                      # ' resistance raster using')
                    # gprint('bounding circles plus a buffer of ' +
                                      # str(float(Cfg.BUFFERDIST)) + ' map '
                                      # 'units.')
                    # start_time = lu.elapsed_time(start_time)
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
                statement = ('gp.SingleOutputMapAlgebra_sa(expression, '
                             'SRCRASTER)')
                while True:
                    try: exec statement
                    except:
                        count, tryAgain = lu.hiccup_test(count, statement)
                        if not tryAgain: exec statement
                    else: break
                # Cost distance raster creation

                count = 0
                statement = ('gp.CostDistance_sa(SRCRASTER, bResistance, '
                             'outDistanceRaster, Cfg.TMAXCWDIST, "BACK")')                            
                while True:
                    try: exec statement
                    except:
                        count, tryAgain = lu.hiccup_test(count, statement)
                        if not tryAgain: exec statement
                    else: break
                start_time = time.clock()
                # Extract cost distances from source core to target cores
                # Fixme: there will be redundant calls to b-a when already
                # done a-b
                # gprint('Getting least cost distances from source '
                            # 'core #' + str(int(sourceCore)) + ' to ' +
                            # str(len(targetCores)) + ' potential targets')
                count = 0
               
                statement = ('gp.zonalstatisticsastable_sa(s3core_ras, '
                             '"VALUE", outDistanceRaster, ZNSTATS)')
                try: exec statement
                except:

                    lu.dashline(1)
                    msg = ('ERROR in Zonal Stats.  This seems to be an ArcGIS'
                           ' bug or a problem with intermediate files being open'
                           ' in ArcMap.\n  Please make sure no intermediate files'
                           ' (e.g., CWD maps) are open, then restart ArcMap and try again.')
                    Cfg.gp.AddError(msg)
                    exit(1)                               
                
                tableRows = gp.searchcursor(ZNSTATS)
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
                #start_time = lu.elapsed_time(start_time)

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
                            'gp.SingleOutputMapAlgebra_sa(expression,'
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
                            statement = ('gp.CostPath_sa(TARGETRASTER, '
                                        'outDistanceRaster, "BACK",  "lcp", '
                                        '"BEST_SINGLE", "")')
                                        
                            #statement = ('gp.CostPath_sa(TARGETRASTER, '
                             #           'outDistanceRaster, backraster,  '
                              #          '"lcp", "BEST_SINGLE", "")')
                                        
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
                                msg = ("Cost path failed.  Should not have "
                                       "gotten to this point. Link "
                                       "#" + str(int(link)) +
                                       ".\nThis is mysterious.  "
                                       "Please report error.\n")

                            else:
                                lu.dashline(1)
                                msg = ("Error in COST PATH function for link "
                                       "#" + str(int(link)) +
                                       ".\nPlease report error.\n")
                                gp.AddError(msg)
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
                            gp.SelectLayerByAttribute(Cfg.FCORES,
                                                      "NEW_SELECTION",
                                                      Cfg.COREFN + ' <> ' +
                                                      str(int(targetCore)) +
                                                      ' AND ' + Cfg.COREFN +
                                                      ' <> ' +
                                                      str(int(sourceCore)))
                            count = 0
                            statement = ('gp.zonalstatisticsastable('
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
                                gprint(
                                    "Found an intermediate core in the "
                                    "least-cost path between cores " +
                                    str(int(sourceCore)) + " and " +
                                    str(int(targetCore)) + ".  The corridor "
                                    "will be removed.")
                                # disable link
                                linkTable[rows,Cfg.LTB_LINKTYPE] = Cfg.LT_INT
                            # -------------------------------------------------

                        # Create lcp shapefile.  lcploop just keeps track of
                        # whether this is first time function is called.
                        lcpLoop = lu.create_lcp_shapefile(linkTable,
                                                          sourceCore,
                                                          targetCore,
                                                          lcpLoop)

                # endTime = time.clock()
                # processTime = round((endTime - start_time), 2)
                # gprint('Intermediate core checks and LCP shapefiles'
                                  # ' for core #' + str(sourceCore) + ' took '
                                  # + str(processTime) + ' seconds.')
                # -----------------------------------------------------------


                gprint('Done with all calculations for core #' +
                                  str(sourceCore) + '.')
                start_time = lu.elapsed_time(startTime1)
  
                outlinkTableFile = path.join(Cfg.DATAPASSDIR, 
                                             "temp_linkTable_s3_partial.csv")
                lu.write_link_table(linkTable, outlinkTableFile)
                    
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
        gprint('Updating ' + outlinkTableFile)
        lu.write_link_table(linkTable, outlinkTableFile)
        linkTableLogFile = path.join(Cfg.LOGDIR, "linkTable_s3.csv")
        lu.write_link_table(linkTable, linkTableLogFile)

        # lu.dashline()
        start_time = time.clock()
        gprint('Creating shapefiles with linework for links...')
        lu.write_link_maps(outlinkTableFile, step=3)
        start_time = lu.elapsed_time(start_time)

        gprint('\nIndividual cost-weighted distance layers written '
                          'to "cwd" directory. \n')
        gprint(outlinkTableFile + 
                '\n updated with cost-weighted distances between core areas.')

        #Clean up temporary files for restart code
        tempFile = path.join(Cfg.DATAPASSDIR, "temp_cores_to_map.csv")
        lu.delete_file(tempFile)
        tempFile = path.join(Cfg.DATAPASSDIR, "temp_linkTable_s3_partial.csv")        
        lu.delete_file(tempFile)
    
    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 3. Details follow.****')
        lu.raise_geoproc_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 3. Details follow.****')
        lu.raise_python_error(__filename__)

    return
