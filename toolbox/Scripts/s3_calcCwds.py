##*****************************************************************
## 2011_0128
## NAME: s3_calcCwds.py
##
## SUMMARY: Calculates cost-weighted distances from each core area.
## Uses bounding circles around source and target cores to limit
## extent of cwd calculations and speed computation.
##
## SOFTWARE: ArcGIS 9.3 (requires Spatial Analyst extension)
##           Python 2.5
##
##*****************************************************************

# import required modules
import sys
import os.path as path
import shutil
import time

import arcgisscripting
from numpy import *

import lm_config
import lm_util as lu

LOGDIR = lm_config.LOGDIR
SCRATCHDIR = lm_config.SCRATCHDIR
CWDBASEDIR = lm_config.CWDBASEDIR
CWDSUBDIR_MAME = lm_config.CWDSUBDIR_MAME
COREFC = lm_config.COREFC
COREFN = lm_config.COREFN
RESRAST = lm_config.RESRAST
S3DROPLCCS = lm_config.S3DROPLCCS
TMAXCWDIST = lm_config.TMAXCWDIST
BUFFERDIST = lm_config.BUFFERDIST
MAXCOSTDIST = lm_config.MAXCOSTDIST
MINCOSTDIST = lm_config.MINCOSTDIST
MAXEUCDIST = lm_config.MAXEUCDIST
MINEUCDIST = lm_config.MINEUCDIST
FCORES = lm_config.FCORES
BNDCIRCEN = lm_config.BNDCIRCEN
BNDCIR = lm_config.BNDCIR
BNDCIRCENS = "boundingCircleCenters.shp"    
BNDCIRS = "boundingCircles.shp"    
GP = lm_config.GP
BNDFC = "boundingFeature.shp"
ZNSTATS = path.join(SCRATCHDIR, "zonestats.dbf")

def step3_calc_cwds():
    """Calculates cost-weighted distances from each core area.
    Uses bounding circles around source and target cores to limit
    extent of cwd calculations and speed computation.    

    """
    try:
        global GP
        GP.scratchWorkspace = SCRATCHDIR         
        linkTableFile = lu.get_prev_step_link_table(step=3)      
  
        # linkTable column numbers
        linkIdCol = 0 # Link ID
        core1Col = 1 # core ID of1st core area link connects
        core2Col = 2 # core ID of2nd core area link connects
        
        # linkTypeCol 
        # 0=no link. 2=corridor, 3=intermediate core area detected, 
        # 4=too long EucDist, 5=too long lcDist 6=too short EucDist, 
        # 7=too short cwDist, 10 = cluster link
        linkTypeCol = 5
        
        eucDistCol = 6
        cwDistCol = 7
        
        if TMAXCWDIST is None: 
           	GP.AddMessage('NOT using a maximum cost-weighted distance.')          
        else:
            GP.addmessage('Max cost-weighted distance for CWD calcs set to ' 
                          + str(TMAXCWDIST) + '\n')		            
		
        if (BUFFERDIST) is not None: # FIXME: because it's integer, it fills in 0 if not entered. 
            GP.addmessage('\nBounding circles plus a buffer of ' + str(float(BUFFERDIST)/1000) + ' km will be used \nto limit extent of cost distance calculations.')
        else:
            GP.addmessage('NOT using bounding circles in cost distance calculations.')

        # get the number of cores
        coreCount = int(GP.GetCount_management(COREFC).GetOutput(0))
        
        # set the analysis extent and cell size 
        GP.CellSize = GP.Describe(RESRAST).MeanCellHeight
        GP.Extent = "MINOF" # So we don't extract rasters that go beyond extent of original raster
        GP.mask = RESRAST
        GP.Workspace = SCRATCHDIR      
        SR = GP.describe(COREFC).SpatialReference # for later shapefiles
    
        # Load linkTable (created in previous script)
        linkTable = lu.load_link_table(linkTableFile)           
        lu.report_links(linkTable)  
          
        rows,cols = where(linkTable[:,linkTypeCol:linkTypeCol+1] == 2) 
        coresToProcess = unique(linkTable[:,core1Col:core2Col+1])
        maxCoreNum = max(coresToProcess)
        del rows,cols,coresToProcess
     
        # Set up cwd directories.
        # To keep there from being > 100 grids in any one directory,
        # outputs are written to:    
        # cwd\cw for cores 1-99
        # cwd\cw1 for cores 100-199
        # etc.        
        if path.exists(CWDBASEDIR):           
            shutil.rmtree(CWDBASEDIR)
        lu.dashline(1) 
        GP.addmessage("Creating cost-weighted distance grid output folders:")
        GP.addmessage(path.join(CWDBASEDIR, CWDSUBDIR_MAME))  
        GP.CreateFolder_management(path.dirname(CWDBASEDIR), 
                                       path.basename(CWDBASEDIR))                   
        GP.CreateFolder_management(CWDBASEDIR, CWDSUBDIR_MAME)        
        if maxCoreNum > 100:
            maxDirCount = int(maxCoreNum/100)
            for dirCount in range(1, maxDirCount + 1):
                ccwdir = CWDSUBDIR_MAME + str(dirCount)
                GP.addmessage(ccwdir)    
                GP.CreateFolder_management(CWDBASEDIR, ccwdir)
        lu.dashline(2)
                
        # make a feature layer for input cores to select from
        GP.MakeFeatureLayer(COREFC, FCORES)
    
        # Identify cores to map from LinkTable
        rows,cols = where(linkTable[:,linkTypeCol:linkTypeCol+1] == 2)
        coresToMap = unique(linkTable[:,core1Col:core2Col+1])
        numCoresToMap = len(coresToMap)
        del rows,cols
        if numCoresToMap < 3:
            S3DROPLCCSic = False # No need to check for intermediate cores, because there aren't any
        else:
            S3DROPLCCSic = S3DROPLCCS
    
        GP.addmessage('Number of core areas to connect:' + str(numCoresToMap))

        # Drop links that are too long 
        GP.addmessage('\nChecking for corridors that are too long to map.')
        disableLeastCostNoVal = False
        linkTable,numDroppedLinks = lu.drop_links(linkTable, MAXEUCDIST, 0, 
                                               MINEUCDIST, 0, 
                                               disableLeastCostNoVal)
        
        # ------------------------------------------------------------------
        # Bounding boxes
        if (BUFFERDIST) is not None:
            # create bounding boxes around cores
            startTime = time.clock()
            lu.dashline(1)
            GP.addmessage('Calculating bounding boxes for core areas....')
            extentBoxList = zeros((0,5), dtype='float32')
            for x in range(len(coresToMap)):
                core = coresToMap[x]
                if len(coresToMap) > 20:
                    report_pct_done(x, len(coresToMap))                
                boxCoords = lu.get_extent_box_coords(core)
                extentBoxList=append(extentBoxList, boxCoords, axis=0)
            GP.addmessage('\nDone calculating bounding boxes.')
            startTime, hours, mins, secs = lu.elapsed_time(startTime)
            lu.dashline()   
    
        # Bounding circle code
        if BUFFERDIST is not None:
            # Make a set of circles encompassing core areas we'll be connecting                 
            startTime = time.clock()
            GP.addmessage('\nCalculating bounding circles around potential corridors.')
    
            # x y corex corey radius- stores data for bounding circle centroids
            boundingCirclePointArray  = zeros((0,5),dtype='float32')
            
            circleList = zeros((0,3),dtype='int32')           
            
            numLinks = linkTable.shape[0]
            for x in range(0,numLinks):
                if numLinks > 20:
                    report_pct_done(x, numLinks)
                if linkTable[x,linkTypeCol]==2: # if it's a valid corridor link
                    linkId = int(linkTable[x,linkIdCol])        
                    cores = zeros((1,3),dtype='int32')#fixme- this code is clumsy- can trim down
                    cores[0,:]=sort([0,linkTable[x,core1Col],linkTable[x,core2Col]])
                    corex=cores[0,1]
                    corey=cores[0,2]
                    cores[0,0]=linkId
                    
                    ###################
                    foundFlag=False
                    for y in range(0,len(circleList)):#clumsy
                        if circleList[y,1]==corex and circleList[y,2]==corey:
                            foundFlag=True
                    if foundFlag==False:
                        circlePointData = lu.get_bounding_circle_data(extentBoxList,corex,corey,BUFFERDIST)
                        boundingCirclePointArray = append(boundingCirclePointArray,circlePointData,axis=0)
                        circleList=append(circleList,cores,axis=0) # keep track of which cores we draw bounding circles around
    
            GP.addmessage('\nCreating bounding circles using buffer analysis.')                    
            lu.make_points(SCRATCHDIR, boundingCirclePointArray, BNDCIRCENS)
            GP.defineprojection(BNDCIRCENS, SR)                                      
            if GP.Exists(BNDCIRS):
                GP.delete_management(BNDCIRS)
            GP.buffer_analysis(BNDCIRCENS, BNDCIRS, "radius")
            GP.defineprojection(BNDCIRS, SR)
            GP.deletefield (BNDCIRS, "BUFF_DIST")

            GP.addmessage('Successfully created bounding circles around potential corridors using \na buffer of ' + str(float(BUFFERDIST)/1000) + ' km.')
            startTime, hours, mins, secs = lu.elapsed_time(startTime)
                
            GP.addmessage('Reducing global processing area using bounding circle plus buffer of ' + str(float(BUFFERDIST)/1000) + ' km.\n')           
            startTime = time.clock()

            extentBoxList = zeros((0,5),dtype='float32')
            boxCoords = lu.get_extent_box_coords()
            extentBoxList=append(extentBoxList,boxCoords,axis=0)
            extentBoxList[0,0]=0
    
            boundingCirclePointArray  = zeros((0,5),dtype='float32')
            circlePointData=lu.get_bounding_circle_data(extentBoxList, 0, 
                                                        0, BUFFERDIST)
                
            lu.make_points(SCRATCHDIR, circlePointData, BNDCIRCEN)
            GP.defineprojection(BNDCIRCEN, SR)                                  
            if GP.Exists(BNDCIR):
                GP.delete_management(BNDCIR)
            GP.buffer_analysis(BNDCIRCEN, BNDCIR, "radius")
            GP.defineprojection(BNDCIR, SR)

            boundResis = "boundResis"
            GP.addmessage('Extracting raster....')

             # FIXME: wishlist- would extract by circle be faster?
            count = 0
            statement = ('GP.ExtractByMask_sa(RESRAST, BNDCIR, '
                         'boundResis)')
            while True:
                try: exec statement
                except:
                    count,tryAgain = lu.hiccup_test(count,statement)
                    if not tryAgain: exec statement
                else: break
            GP.addmessage('\nReduced resistance raster extracted using '
                          'bounding circle.')
            startTime, hours, mins, secs = lu.elapsed_time(startTime)

        else: #if not using bounding circles, just go with resistance raster.
            boundResis = RESRAST 
       
        # ---------------------------------------------------------------------       
        # Rasterize core areas to speed cost distance calcs
        lu.dashline(1)
        GP.addmessage("Creating core area raster.")    
        s3core_ras="s3core_ras"
        GP.SelectLayerByAttribute(FCORES, "CLEAR_SELECTION")
        GP.CellSize = GP.Describe(boundResis).MeanCellHeight
        if GP.exists(s3core_ras):
            GP.delete_management(s3core_ras)       
        GP.extent = GP.Describe(boundResis).extent
        count = 0 
        statement = ('GP.FeatureToRaster_conversion(FCORES, COREFN, '
                     's3core_ras, GP.Cellsize)')
        while True:
            try: exec statement
            except:
                count,tryAgain = lu.hiccup_test(count,statement)
                if not tryAgain: exec statement
            else: break
        GP.extent = "MINOF"
 
        #----------------------------------------------------------------------
        # Loop through cores, do cwd calcs for each
        GP.addmessage("Starting cost distance calculations.")   
        lcpLoop = 0
        for x in range(len(coresToMap)):
            startTime1 = time.clock()    
            # This is the focal core we're running cwd out from            
            sourceCore = int(coresToMap[x]) 
    
            # Get target cores based on linktable with reinstated links 
            # (we temporarily disable them below by adding 100)
            linkTableTemp = linkTable.copy()
            # reinstate temporarily disabled links
            rows = where(linkTableTemp[:,linkTypeCol] > 100)
            linkTableTemp[rows,linkTypeCol] = (linkTableTemp[rows,linkTypeCol] 
                                               - 100)
            # get core areas to be connected to focal core
            targetCores = lu.get_core_targets(sourceCore, linkTableTemp) 
            del linkTableTemp
            
            if len(targetCores)>0: 
                lu.dashline(1)
                GP.addmessage('Target core areas for core area #' 
                              + str(sourceCore) + ' = ' + str(targetCores))

                # ----------------------------------------------------------------------------------------
                # Create BOUNDING FEATURE to limit extent of cost distance calculations-
                # This is a set of circles encompassing core areas we'll be connecting each core area to.    
                if BUFFERDIST is not None:        
                    # fixme: move outside of loop   # new circle 
                    GP.MakeFeatureLayer(path.join(GP.workspace, 
                                        BNDCIRS), "fGlobalBoundingFeat")
                    startTime = time.clock()    
                    # loop through targets and get bounding circles that 
                    # contain focal core and target cores   
                    GP.AddMessage("\nAdding up bounding circles for source "
                                  "core " + str(sourceCore))
                    GP.SelectLayerByAttribute("fGlobalBoundingFeat", 
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
                        # fixme:need to check for case where link is not found.
                        GP.SelectLayerByAttribute("fGlobalBoundingFeat", 
                                                  "ADD_TO_SELECTION", field 
                                                  + " = '" + cores_x_y + "'")
                        
                    if GP.Exists(BNDFC):
                        GP.delete_management(BNDFC) # fixme: necessary?
                    # fixme: may not be needed- can we just clip raster 
                    # using selected?
                    GP.CopyFeatures_management ("fGlobalBoundingFeat", BNDFC)
    
                    # Clip out bounded area of resistance raster for cwd 
                    # calculations from focal core
                    bResistance = "bResistance"
                    count = 0
                    statement = ('GP.ExtractByMask_sa(boundResis, '
                                 'BNDFC, bResistance)')
                    while True:
                        try: exec statement
                        except:
                            count,tryAgain = lu.hiccup_test(count,statement)
                            if not tryAgain: exec statement
                        else: break
                    GP.addmessage('Successfully extracted a reduced resistance raster using')
                    GP.addmessage('bounding circles plus a buffer of ' 
                                  + str(float(BUFFERDIST)/1000) + ' km.')
                    startTime, hours, mins, secs = lu.elapsed_time(startTime)
                else:
                    bResistance=boundResis
 
                # ---------------------------------------------------------
                # CWD Calculations
                outDistanceRaster = lu.get_cwd_path(sourceCore)
                startTime = time.clock()

                # Create raster that just has source core in it                
                # Note: this seems faster than setnull with LI grid.
                expression = ("con(" + s3core_ras + "== " 
                              + str(int(sourceCore)) + ",1)")
                SRCRASTER = 'source'               
                count = 0
                statement = ('GP.SingleOutputMapAlgebra_sa(expression, '
                             'SRCRASTER)')
                while True:
                    try: exec statement
                    except:
                        count,tryAgain = lu.hiccup_test(count,statement)
                        if not tryAgain: exec statement
                    else: break

                # Cost distance raster creation
                count = 0 
                statement = ('GP.CostDistance_sa(SRCRASTER, bResistance, '
                             'outDistanceRaster, TMAXCWDIST, "BACK")')
                while True:
                    try: exec statement
                    except:
                        count, tryAgain = lu.hiccup_test(count, statement)
                        if not tryAgain: exec statement
                    else: break
                GP.addmessage('Cost distances for source core ' 
                              + str(int(sourceCore)) + ' calculated.')
                startTime, hours, mins, secs = lu.elapsed_time(startTime)

                # Extract cost distances from source core to target cores
                # Fixme: there will be redundant calls to b-a when already 
                # done a-b 
                GP.addmessage('Getting least cost distances from source core #'                
                              + str(int(sourceCore)) + ' to ' 
                              + str(len(targetCores)) + ' potential targets') 
                count = 0
                statement = ('GP.zonalstatisticsastable_sa(s3core_ras, '
                             '"VALUE", outDistanceRaster, ZNSTATS)')
                while True:
                    try: exec statement
                    except:
                        count,tryAgain = lu.hiccup_test(count,statement)
                        if not tryAgain: exec statement
                    else: break

                tableRows = GP.searchcursor(ZNSTATS)               
                tableRow = tableRows.Next()
                while tableRow:
                    if tableRow.Value > sourceCore:
                        link = lu.get_links_from_core_pairs(linkTable,sourceCore,tableRow.Value)
                        linkTable[link,cwDistCol]=tableRow.Min                       
                        if MAXCOSTDIST is not None:
                            if tableRow.Min > MAXCOSTDIST:
                                linkTable[link,linkTypeCol] = 5 # Disable link, it's too long 
                        if MINCOSTDIST is not None:
                            if tableRow.Min < MINCOSTDIST:
                                linkTable[link,linkTypeCol] = 7 # Disable link, it's too short
                    tableRow = tableRows.next()
                del tableRow, tableRows
                startTime, hours, mins, secs = lu.elapsed_time(startTime)

                # ---------------------------------------------------------                
                # Check for intermediate cores AND map LCP lines 
                for y in range(len(targetCores)):
                    targetCore = targetCores[y]    
                    rows = lu.get_links_from_core_pairs(linkTable, sourceCore, targetCore)
                    if linkTable[rows[0],linkTypeCol] < 100 and linkTable[rows[0],cwDistCol]!=-1: # Map all links for which above code successfully extracted cwds in above code
                        linkTable[rows,linkTypeCol]=linkTable[rows,linkTypeCol] + 100 # Flag so that we only evaluate this pair once
                        # Create raster that just has target core in it
                        expression = "con(" + s3core_ras + "== " + str(int(targetCore)) + ",1)"
                        TARGETRASTER = 'targ'
                        count,statement = 0, 'GP.SingleOutputMapAlgebra_sa(expression, TARGETRASTER)'
                        while True:
                            try: exec statement
                            except:
                                count,tryAgain = lu.hiccup_test(count,statement)
                                if not tryAgain: exec statement
                            else: break

                        try:
                            # Cost path allows us to map the least cost path between source and target
                            count,statement=0,'GP.CostPath_sa(TARGETRASTER,outDistanceRaster,"BACK", "lcp","BEST_SINGLE","")'
                            try:
                                exec statement
                            except:
                                exec statement                      
                        except:
                            link = lu.get_links_from_core_pairs(linkTable, sourceCore, targetCore)
                            # Not picked up by above cwd calc code
                            if (MAXCOSTDIST is not None 
                                and linkTable[link,cwDistCol] == -1): 
                                GP.addmessage('Cost path failed- should not '
                                              'have gotten to this point?')
                                continue
                            else:
                                lu.dashline(1)       
                                msg = "Error in COST PATH function for link # " + str(int(link)) + ".\nTell Brad you broke his new code.\n"
                                GP.AddError(msg)
                                exit(0)
    
                        # fixme: may be fastest to not do selection, do EXTRACTBYMASK, getvaluelist, use code snippet at end of file to discard src and target values.
                        # still this is fast- 13 sec for LI data...
                        # But I'm not very comfortable using failed coreMin as our test....
                        if S3DROPLCCSic:
                            # ---------------------------------------------------------
                            # Drop links where lcp passes through intermediate core area
                            # Method below is faster than valuelist method because of soma in valuelist method.
                            GP.SelectLayerByAttribute(FCORES, 
                                                      "NEW_SELECTION", 
                                                      COREFN + ' <> ' 
                                                      + str(int(targetCore)) 
                                                      + ' AND ' + COREFN 
                                                      + ' <> ' 
                                                      + str(int(sourceCore)))
                            count,statement = 0, 'GP.zonalstatisticsastable(FCORES, COREFN, "lcp", ZNSTATS, "DATA", "MINIMUM")'
                            while True:
                                try: exec statement
                                except:
                                    count,tryAgain = lu.hiccup_test(count,statement)
                                    if not tryAgain: exec statement
                                else: break        
                            rows = lu.get_links_from_core_pairs(linkTable,sourceCore,targetCore) 
                            coreMin = lu.get_zonal_minimum(ZNSTATS)
                            if coreMin != 'Failed': #  Found a valid value, indicating overlap
                                lu.dashline(1)
                                GP.addmessage("Found an intermediate core in the least-cost path between cores " + str(int(sourceCore)) + " and " + str(int(targetCore))+".  The corridor will be removed.\n")
                                linkTable[rows,linkTypeCol]=3 #  disable link
                            # ---------------------------------------------------------

                        # Create lcp shapefile.  lcploop just keeps track of whether this is first time function is called.                            
                        lcpLoop = lu.create_lcp_shapefile(linkTable,sourceCore,targetCore,lcpLoop,SR)
            
                endTime = time.clock()
                processTime = round((endTime - startTime), 2)
                GP.addmessage('Intermediate core checks and LCP shapefiles for core #' + str(sourceCore) + ' took '+str(processTime)+' seconds.')
                # -----------------------------------------------------------
                
                
                GP.addmessage('Done with all calculations for core #' + str(sourceCore) + '.')
                startTime,hours,mins,secs = lu.elapsed_time(startTime1)
            # -----------------------------------------------------------
         
        # reinstate temporarily disabled links
        rows = where(linkTable[:,linkTypeCol] > 100) 
        linkTable[rows,linkTypeCol] = linkTable[rows,linkTypeCol] - 100
        
        # Drop links that are too long
        disableLeastCostNoVal = True
        linkTable,numDroppedLinks = lu.drop_links(linkTable, MAXEUCDIST, 
                                               MINEUCDIST, MAXCOSTDIST, 
                                               MINCOSTDIST, 
                                               disableLeastCostNoVal)
        
        # Write link table file
        linkTableFile = lu.get_this_step_link_table(step=3)
        GP.addmessage('\nUpdating ' + linkTableFile)
        lu.write_link_table(linkTable, linkTableFile)
        linkTableLogFile = path.join(LOGDIR, "linkTable_step3.csv")
        lu.write_link_table(linkTable, linkTableLogFile)
        
        lu.dashline()
        GP.addmessage('\nCreating shapefiles with linework for links...')
        lu.write_link_maps(linkTableFile, step=3)

        startTime = time.clock()
        dummy = lu.update_lcp_shapefile(linkTable, lastStep=3, thisStep=3)
        startTime, hours, mins, secs = lu.elapsed_time(startTime)
            
        GP.addmessage('\nIndividual cost-weighted distance layers written to '
                      '"cwd" directory. \n')
        GP.addmessage(linkTableFile + 
                      '\n updated with cost-weighted distances between core '
                      'areas\n')
    
    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        GP.addmessage('****Failed in step 3. Details follow.****')        
        filename =  __file__
        lu.raise_geoproc_error(filename)
    
    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        GP.addmessage('****Failed in step 3. Details follow.****')        
        filename =  __file__
        lu.raise_python_error(filename)    
        
    return
