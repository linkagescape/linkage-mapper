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
import arcgisscripting, sys, time, os
from numpy import *
from watools_util import *

def step3_calc_cwds(gp,Version,options):
    """Calculates cost-weighted distances from each core area.
    Uses bounding circles around source and target cores to limit
    extent of cwd calculations and speed computation.    

    """
    try:
        # ------------------------------------------------------------------                    
        # Unpack options
        projectDir=options["projectDir"]
        coreShapefile=options["coreShapefile"]
        coreIds=options["coreIds"]
        resistanceRas=options["resistanceRas"]
        bufferDist=options["bufferDist"]
        useMinEucDist=options["useMinEucDist"]
        useMaxEucDist=options["useMaxEucDist"]
        maxEucDist=options["maxEucDist"]
        minEucDist=options["minEucDist"]
        useMinCostDist=options["useMinCostDist"]
        useMaxCostDist=options["useMaxCostDist"]
        maxCostDist=options["maxCostDist"]
        minCostDist=options["minCostDist"]
        maxCwDist=options["maxCwDist"]
        dropLccsWithIntermediateCores=options["dropLccsWithIntermediateCores"]
        defaultNullValue = options["defaultNullValue"]
        coreAreaIds = coreIds # We now set variable coreAreaIds field name to user input ID field name                
        # ------------------------------------------------------------------            

        # ------------------------------------------------------------------                
        # Folders, paths, filenames, and workspaces
        scratchDir = projectDir + "\\scratch\\"
        if os.path.exists(scratchDir)==False:
            gp.CreateFolder_management(projectDir,"scratch")
        gp.scratchWorkspace = scratchDir 
        dbfFile = scratchDir + "\\" + "zonestats.dbf"
        outputRootDir = projectDir + "\\output\\"
        dataPassDir = projectDir + "\\dataPass\\"
        linkTableFile = get_prev_step_link_table(projectDir,step=3)      
        # ------------------------------------------------------------------            
            
        # ------------------------------------------------------------------
        ## linkTable column numbers
        linkIdCol=0 # Link ID
        core1Col=1 # core ID of1st core area link connects
        core2Col=2 # core ID of2nd core area link connects
        linkTypeCol=5 # 0=no link. 2=corridor, 3=intermediate core area detected, 4=too long EucDist, 5=too long lcDist 6=too short EucDist, 7=too short cwDist, 10 = cluster link
        eucDistCol=6
        cwDistCol=7
        # ------------------------------------------------------------------
        
        if (maxCwDist != defaultNullValue): 
            gp.addmessage('Max cost-weighted distance for CWD calcs set to ' + str(maxCwDist) + '\n')			
            maxD=maxCwDist			
        else:
            gp.AddMessage('NOT using a maximum cost-weighted distance.')
            maxD="#" # "#" is for gp.CostDistance call
		
        if (bufferDist) != defaultNullValue: # FIXME: because it's integer, it fills in 0 if not entered. 
            gp.addmessage('\nBounding circles plus a buffer of ' + str(float(bufferDist)/1000) + ' km will be used \nto limit extent of cost distance calculations.')
        else:
            gp.addmessage('NOT using bounding circles in cost distance calculations.')

        # get the number of cores
        coreCount = int(gp.GetCount_management(coreShapefile).GetOutput(0))
        
        # set the analysis extent and cell size 
        gp.CellSize = gp.Describe(resistanceRas).MeanCellHeight
        gp.Extent = "MINOF" # So we don't extract rasters that go beyond extent of original raster
        gp.SnapRaster = resistanceRas
        gp.mask = resistanceRas
        gp.Workspace = scratchDir      
        SR = gp.describe(coreShapefile).SpatialReference # for later shapefiles
    
        # Load linkTable (created in previous script)
        linkTable = load_link_table(linkTableFile)           
        report_links(linkTable)  
          
        rows,cols = where(linkTable[:,linkTypeCol:linkTypeCol+1] == 2) 
        coresToProcess = unique(linkTable[:,core1Col:core2Col+1])
        maxCoreNum = max(coresToProcess)
        del rows,cols,coresToProcess
    
        # ------------------------------------------------------------------
        # Set up cwd directories.
        # To keep there from being > 100 grids in any one directory,
        # outputs are written to:    
        # cwd\cw for cores 1-99
        # cwd\cw1 for cores 100-199
        # etc.
        cwdBasePath = projectDir + "\\cwd\\"
        if os.path.exists(cwdBasePath)==True:
            import shutil
            shutil.rmtree(cwdBasePath)
        gp.addmessage('\n---------------------')       
        gp.addmessage("Creating cost-weighted distance grid output folders:")
        gp.addmessage("   cwd\\cw")            
        gp.CreateFolder_management(projectDir,"cwd\\")
        gp.CreateFolder_management(projectDir,"cwd\\cw\\")
        if maxCoreNum > 100:
            maxDirCount = int(maxCoreNum/100)
            for dirCount in range(1,maxDirCount+1):
                gp.addmessage("   cwd\\cw"+str(dirCount))    
                gp.CreateFolder_management(projectDir,"\\cwd\\cw"+str(dirCount))
        gp.addmessage('---------------------\n')
        # ------------------------------------------------------------------

        
        # make a feature layer for input cores to select from
        gp.MakeFeatureLayer(coreShapefile,"fcores")
    
        # Identify cores to map from LinkTable
        rows,cols = where(linkTable[:,linkTypeCol:linkTypeCol+1] == 2)
        coresToMap = unique(linkTable[:,core1Col:core2Col+1])
        numCoresToMap = len(coresToMap)
        del rows,cols
        if numCoresToMap < 3:
            dropLccsWithIntermediateCores = False # No need to check for intermediate cores, because there aren't any
    
        gp.addmessage('Number of core areas to connect:' + str(numCoresToMap))

        # Drop links that are too long 
        gp.addmessage('\nChecking for corridors that are too long to map.')
        disableLeastCostNoVal = False
        linkTable,numDroppedLinks = drop_links(linkTable,useMaxEucDist,maxEucDist,False,maxCostDist,useMinEucDist,minEucDist,False,minCostDist,disableLeastCostNoVal)
        
        # ------------------------------------------------------------------
        # Bounding boxes
        if (bufferDist) != defaultNullValue:
            # create bounding boxes around cores
            startTime = time.clock()
            gp.addmessage('\n---------------------')
            gp.addmessage('Calculating bounding boxes for core areas....')
            extentBoxList = zeros((0,5),dtype='float32')
            for x in range(len(coresToMap)):
                core = coresToMap[x]
                if len(coresToMap) > 20:
                    report_pct_done(x, len(coresToMap))                
                boxCoords = get_extent_box_coords(gp.workspace,"fcores",coreAreaIds,core)
                extentBoxList=append(extentBoxList,boxCoords,axis=0)
            gp.addmessage('\nDone calculating bounding boxes.')
            startTime,hours,mins,secs = elapsed_time(startTime)
            gp.addmessage('---------------------')    
        # ------------------------------------------------------------------


        # ------------------------------------------------------------------
        # Bounding circle code
        if bufferDist != defaultNullValue:
            # Make a set of circles encompassing core areas we'll be connecting                 
            startTime = time.clock()
            gp.addmessage('\nCalculating bounding circles around potential corridors.')
    
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
                        circlePointData=get_bounding_circle_data(extentBoxList,corex,corey,bufferDist)
                        boundingCirclePointArray = append(boundingCirclePointArray,circlePointData,axis=0)
                        circleList=append(circleList,cores,axis=0) # keep track of which cores we draw bounding circles around
    
            gp.addmessage('\nCreating bounding circles using buffer analysis.')        
            boundingCircleCenters = "boundingCircleCenters.shp"
            make_points(scratchDir,boundingCirclePointArray,boundingCircleCenters,coreIds)
            gp.defineprojection(boundingCircleCenters, SR)
            boundingCircles="boundingCircles.shp"                                  
            if gp.Exists(boundingCircles):
                gp.delete_management(boundingCircles)
            gp.buffer_analysis(boundingCircleCenters, boundingCircles, "radius")
            gp.defineprojection(boundingCircles, SR)
            gp.deletefield (boundingCircles, "BUFF_DIST")

            gp.addmessage('Successfully created bounding circles around potential corridors using \na buffer of ' + str(float(bufferDist)/1000) + ' km.')
            startTime,hours,mins,secs = elapsed_time(startTime)
                
            gp.addmessage('Reducing global processing area using bounding circle plus buffer of ' + str(float(bufferDist)/1000) + ' km.\n')           
            startTime = time.clock()

            extentBoxList = zeros((0,5),dtype='float32')
            boxCoords = get_extent_box_coords(gp.workspace,"fcores",coreIds,"#")
            extentBoxList=append(extentBoxList,boxCoords,axis=0)
            extentBoxList[0,0]=0
    
            boundingCirclePointArray  = zeros((0,5),dtype='float32')
            circlePointData=get_bounding_circle_data(extentBoxList,0,0,bufferDist)
    
            boundingCircleCenter = "boundingCircleCenter.shp"
            make_points(scratchDir,circlePointData,boundingCircleCenter,coreIds)
            gp.defineprojection(boundingCircleCenter, SR)
            boundingCircle="boundingCircle.shp"                                  
            if gp.Exists(boundingCircle):
                gp.delete_management(boundingCircle)
            gp.buffer_analysis(boundingCircleCenter, boundingCircle, "radius")
            gp.defineprojection(boundingCircle, SR)

            boundResis="boundResis"
            gp.addmessage('Extracting raster....')

             # FIXME: wishlist- would extract by circle be faster?
            count,statement = 0, 'gp.ExtractByMask_sa(resistanceRas, boundingCircle, boundResis)'
            while True:
                try: exec statement
                except:
                    count,tryAgain = hiccup_test(count,statement)
                    if not tryAgain: exec statement
                else: break
            gp.addmessage('\nReduced resistance raster extracted using bounding circle.')
            startTime,hours,mins,secs = elapsed_time(startTime)

        else: #if not using bounding circles, just go with resistance raster.
            boundResis=resistanceRas
        # ---------------------------------------------------------------------
       
        # ---------------------------------------------------------------------       
        # Rasterize core areas to speed cost distance calcs
        gp.addmessage('\n---------------------')
        gp.addmessage("Creating core area raster.")    
        core_ras="core_ras"
        gp.SelectLayerByAttribute("fcores", "CLEAR_SELECTION")
        gp.CellSize = gp.Describe(boundResis).MeanCellHeight
        if gp.exists(core_ras):
            gp.delete_management(core_ras)       
        gp.extent = gp.Describe(boundResis).extent
        count,statement = 0, 'gp.FeatureToRaster_conversion("fcores",coreAreaIds,core_ras,gp.Cellsize)'
        while True:
            try: exec statement
            except:
                count,tryAgain = hiccup_test(count,statement)
                if not tryAgain: exec statement
            else: break
        gp.extent = "MINOF"
        # ---------------------------------------------------------------------    
    
    
        #------------------------------------------------------------------------------
        # Loop through cores, do cwd calcs for each
        gp.addmessage("Starting cost distance calculations.")   
        lcpLoop=0
        for x in range(len(coresToMap)):
            startTime1 = time.clock()       
            sourceCore = int(coresToMap[x]) # This is the focal core we're running cwd out from
    
            # Get target cores based on linktable with reinstated links (we temporarily disable them below by adding 100)
            linkTableTemp = linkTable.copy()
            rows = where(linkTableTemp[:,linkTypeCol]>100) # reinstate temporarily disabled links
            linkTableTemp[rows,linkTypeCol]=linkTableTemp[rows,linkTypeCol]-100           
            targetCores = get_core_targets(sourceCore,linkTableTemp) # get core areas to be connected to focal core
            del linkTableTemp
            
            if len(targetCores)>0: 
                gp.addmessage('\n---------------------')
                gp.addmessage('Target core areas for core area #'+ str(int(sourceCore)) + ' = ' + str(targetCores))

                # ----------------------------------------------------------------------------------------
                # Create BOUNDING FEATURE to limit extent of cost distance calculations-
                # This is a set of circles encompassing core areas we'll be connecting each core area to.    
                if bufferDist != defaultNullValue:
                    gp.MakeFeatureLayer(gp.workspace+"\\boundingCircles.shp","fGlobalBoundingFeat") # fixme: move outside of loop   # new circle 
                    startTime = time.clock()    
                    # loop through targets and get bounding circles that contain focal core and target cores   
                    gp.AddMessage("\nAdding up bounding circles for source core " + str(int(sourceCore)))
                    gp.SelectLayerByAttribute("fGlobalBoundingFeat", "CLEAR_SELECTION") 
                    for i in range(len(targetCores)):
                        # run thru circleList, find link that core pair corresponds to. 
                        if sourceCore < targetCores[i]: 
                            corex=sourceCore
                            corey=targetCores[i] 
                        else:
                            corey=sourceCore
                            corex=targetCores[i]
                            
                        cores_x_y = str(int(corex))+'_'+str(int(corey))                      
                        field = "cores_x_y"
                        gp.SelectLayerByAttribute("fGlobalBoundingFeat", "ADD_TO_SELECTION", field + " = '" + cores_x_y + "'")     # fixme:need to check for case where link is not found.
    
                    boundingFeature="boundingFeature.shp"
                    if gp.Exists(boundingFeature):
                        gp.delete_management(boundingFeature) # fixme: necessary?
                    gp.CopyFeatures_management ("fGlobalBoundingFeat", boundingFeature) # fixme: may not be needed- can we just clip raster using selected?
    
                    # Clip out bounded area of resistance raster for cwd calculations from focal core.
                    bResistance="bResistance"
                    count,statement = 0, 'gp.ExtractByMask_sa(boundResis, boundingFeature, bResistance)'
                    while True:
                        try: exec statement
                        except:
                            count,tryAgain = hiccup_test(count,statement)
                            if not tryAgain: exec statement
                        else: break
                    gp.addmessage('Successfully extracted a reduced resistance raster using')
                    gp.addmessage('bounding circles plus a buffer of ' + str(float(bufferDist)/1000) + ' km.')
                    startTime,hours,mins,secs = elapsed_time(startTime)
                else:
                    bResistance=boundResis
                # ---------------------------------------------------------

                    
                # ---------------------------------------------------------
                # CWD Calculations
                outDistanceRaster = get_cwd_path(projectDir,sourceCore)
                startTime = time.clock()
                
                # Create raster that just has source core in it                
                expression = "con(" + core_ras + "== " + str(int(sourceCore)) + ",1)" # Note: this seems faster than setnull with LI grid.
                SRCRASTER = 'source'               
                count,statement = 0, 'gp.SingleOutputMapAlgebra_sa(expression, SRCRASTER)'
                while True:
                    try: exec statement
                    except:
                        count,tryAgain = hiccup_test(count,statement)
                        if not tryAgain: exec statement
                    else: break


###################################### INSERTING IN GRASS ASCII CONVERSIONS HERE
######################################
                # Take grass cwd and grass back asciis and write them as ARCINFO grids 
                gp.ASCIIToRaster(projectDir + "\\cwdascii\\cwdascii_" + str(int(sourceCore)) + ".asc", outDistanceRaster, "FLOAT")
                # I don't think that the Backlink rasters need to be kept, because they seem to be temporary files in the script - but I'm not sure. 
                gp.ASCIIToRaster(projectDir + "\\backascii\\backascii_" + str(int(sourceCore)) + ".asc", "grBACK", "FLOAT")
                #NEED TO RECLASSIFY from the directional degree output from GRASS to Arc's 1 to 8 directions format. Not clear exactly how to convert. Need to confirm.  
                gp.Reclassify_sa("grBACK", "Value", "0 5;45 4;90 3;135 2;180 1;225 8;270 7;315 6", "BACK", "DATA")
                print "completed grass ascii to arc raster conversion FOR CWD AND BACK RASTERS" 
######################################
######################################
 
### COMMENTING OUT THE ORIGINAL CWD CALCULATION 
                # Cost distance raster creation
##              ######ORIGINAL CWD CALCULATION ######  count,statement = 0, 'gp.CostDistance_sa(SRCRASTER, bResistance, outDistanceRaster, str(maxD), "BACK")'
##                while True:
##                    try: exec statement
##                    except:
##                        count,tryAgain = hiccup_test(count,statement)
##                        if not tryAgain: exec statement
##                    else: break
                gp.addmessage('Cost distances for source core ' + str(int(sourceCore)) + ' calculated.')
                startTime,hours,mins,secs = elapsed_time(startTime)

                # Extract cost distances from source core to target cores
                gp.addmessage('Getting least cost distances from source core #' + str(int(sourceCore)) + ' to ' + str(len(targetCores)) + ' potential targets') # Fixme: there will be redundant calls to b-a when already done a-b 
                count,statement = 0, 'gp.zonalstatisticsastable_sa(core_ras,"VALUE",outDistanceRaster,dbfFile)'
                while True:
                    try: exec statement
                    except:
                        count,tryAgain = hiccup_test(count,statement)
                        if not tryAgain: exec statement
                    else: break

                tableRows = gp.searchcursor(dbfFile)               
                tableRow = tableRows.Next()
                while tableRow:
                    if tableRow.Value > sourceCore:
                        link = get_links_from_core_pairs(linkTable,sourceCore,tableRow.Value)
                        linkTable[link,cwDistCol]=tableRow.Min                       
                        if useMaxCostDist == True:
                            if tableRow.Min > maxCostDist:
                                linkTable[link,linkTypeCol] = 5 # Disable link, it's too long 
                        if useMinCostDist == True:
                            if tableRow.Min < minCostDist:
                                linkTable[link,linkTypeCol] = 7 # Disable link, it's too short
                    tableRow = tableRows.next()
                del tableRow, tableRows
                startTime,hours,mins,secs = elapsed_time(startTime)

                # ---------------------------------------------------------                
                # Check for intermediate cores AND map LCP lines 
                for y in range(len(targetCores)):
                    targetCore = targetCores[y]    
                    rows = get_links_from_core_pairs(linkTable,sourceCore,targetCore)
                    if linkTable[rows[0],linkTypeCol] < 100 and linkTable[rows[0],cwDistCol]!=-1: # Map all links for which above code successfully extracted cwds in above code
                        linkTable[rows,linkTypeCol]=linkTable[rows,linkTypeCol] + 100 # Flag so that we only evaluate this pair once
                        # Create raster that just has target core in it
                        expression = "con(" + core_ras + "== " + str(int(targetCore)) + ",1)"
                        TARGETRASTER = 'targ'
                        count,statement = 0, 'gp.SingleOutputMapAlgebra_sa(expression, TARGETRASTER)'
                        while True:
                            try: exec statement
                            except:
                                count,tryAgain = hiccup_test(count,statement)
                                if not tryAgain: exec statement
                            else: break

                        try:
                            # Cost path allows us to map the least cost path between source and target

############################## HERE IS THE COST PATH CALCULATION; DOESNT CHANGE BECAUSE THE OUTPUTS
############################### FROM THE GRASS ASCII CONVERSION HAVE THE SAME NAMES AS THE INPUTS TO GP.COSTPATH_SA

                            count,statement=0,'gp.CostPath_sa(TARGETRASTER,outDistanceRaster,"BACK", "lcp","BEST_SINGLE","")'
                            try:
                                exec statement
                            except:
                                exec statement                      
                        except:
                            link = get_links_from_core_pairs(linkTable,sourceCore,targetCore)
                            if useMaxCostDist == True and linkTable[link,cwDistCol]==-1: # Not picked up by above cwd calc code
                                gp.addmessage('Cost path failed- should not have gotten to this point?') # Debug code
                                continue
                            else:
                                gp.AddMessage("\n--------------------------------------")        
                                msg = "Error in COST PATH function for link # " + str(int(link)) + ".\nTell Brad you broke his new code.\n"
                                gp.AddError(msg)
                                exit(0)
    
                        # fixme: may be fastest to not do selection, do EXTRACTBYMASK, getvaluelist, use code snippet at end of file to discard src and target values.
                        # still this is fast- 13 sec for LI data...
                        # But I'm not very comfortable using failed coreMin as our test....
                        if dropLccsWithIntermediateCores == True:
                            # ---------------------------------------------------------
                            # Drop links where lcp passes through intermediate core area
                            # Method below is faster than valuelist method because of soma in valuelist method.
                            gp.SelectLayerByAttribute("fcores", "NEW_SELECTION", coreAreaIds + ' <> ' + str(int(targetCore)) + ' AND ' + coreAreaIds + ' <> ' + str(int(sourceCore)) )###
                            count,statement = 0, 'gp.zonalstatisticsastable("fcores",coreAreaIds,"lcp",dbfFile,"DATA","MINIMUM")' #
                            while True:
                                try: exec statement
                                except:
                                    count,tryAgain = hiccup_test(count,statement)
                                    if not tryAgain: exec statement
                                else: break        
                            rows = get_links_from_core_pairs(linkTable,sourceCore,targetCore) 
                            coreMin = get_zonal_minimum(dbfFile)
                            if coreMin != 'Failed': #  Found a valid value, indicating overlap
                                gp.addmessage('\n---------------------------------')
                                gp.addmessage("Found an intermediate core in the least-cost path between cores " + str(int(sourceCore)) + " and " + str(int(targetCore))+".  The corridor will be removed.\n")
                                linkTable[rows,linkTypeCol]=3 #  disable link
                            # ---------------------------------------------------------

                        # Create lcp shapefile.  lcploop just keeps track of whether this is first time function is called.                            
                        lcpLoop = create_lcp_shapefile(projectDir,linkTable,sourceCore,targetCore,lcpLoop,SR)
            
                    endTime = time.clock()
                    processTime = round((endTime - startTime), 2)
                    gp.addmessage('Intermediate core checks and LCP shapefiles for core #' + str(sourceCore) + ' took '+str(processTime)+' seconds.')
                # -----------------------------------------------------------
                
                
                gp.addmessage('\nDone with all calculations for core #' + str(sourceCore) + '.')
                startTime,hours,mins,secs = elapsed_time(startTime1)
            # -----------------------------------------------------------
         
        # reinstate temporarily disabled links
        rows = where(linkTable[:,linkTypeCol]>100) 
        linkTable[rows,linkTypeCol]=linkTable[rows,linkTypeCol]-100
        
        # Drop links that are too long
        disableLeastCostNoVal = True
        linkTable,numDroppedLinks = drop_links(linkTable,useMaxEucDist,maxEucDist,useMaxCostDist,maxCostDist,useMinEucDist,minEucDist,useMinCostDist,minCostDist,disableLeastCostNoVal)
        
        # Write link table file
        linkTableFile = get_this_step_link_table(projectDir,step=3)
        gp.addmessage('\nUpdating ' + linkTableFile)
        write_link_table(linkTable,coreIds,linkTableFile)
        linkTableLogFile=projectDir+"\\log\\"+"linkTable_step3.csv"   
        write_link_table(linkTable,coreIds,linkTableLogFile)
        
        gp.addmessage('---------------------')
        gp.addmessage('\nCreating shapefiles with linework for links...')
        write_link_maps(outputRootDir,linkTableFile,coreShapefile,coreIds,coreAreaIds,step=3)

        startTime = time.clock()
        dummy = update_lcp_shapefile(projectDir,linkTable,lastStep=3,thisStep=3)
        startTime,hours,mins,secs = elapsed_time(startTime)
            
        gp.addmessage('\nIndividual cost-weighted distance layers written to "cwd" directory. \n')
        gp.addmessage(linkTableFile+'\n updated with cost-weighted distances between core areas\n')
    
    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        gp.addmessage("\n--------------------")
        gp.addmessage('****Failed in step 3. Details follow.****')        
        filename =  __file__
        raise_geoproc_error(filename)
    
    # Return any PYTHON or system specific errors
    except:
        gp.addmessage("\n--------------------")
        gp.addmessage('****Failed in step 3. Details follow.****')        
        filename =  __file__
        raise_python_error(filename)    

    del gp    
    return
