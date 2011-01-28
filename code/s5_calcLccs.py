##*****************************************************************
## 2011_0128
## NAME: s5_calcLccs.py
##
## SUMMARY: Creates and mosaics normalized least-cost corridors
## using connected core area pairs specified in linkTable and
## cwd layers
##
## SOFTWARE: ArcGIS 9.3 (requires Spatial Analyst extension)
##           Python 2.5
##
##*****************************************************************

# import required modules
import arcgisscripting, sys, time
from numpy import *
from watools_util import *
import shutil

def step5_calc_lccs(gp,Version,options):
    """Creates and mosaics normalized least-cost corridors
    using connected core area pairs specified in linkTable and
    cwd layers
    
    """
       
# Fixme: add option to saveRawLccs? Or mosaicLccs that already exist?
    try:
        gp.OverwriteOutput = 1
    
        # --------------------------------------------
        # Unpack options
        coreShapefile=options["coreShapefile"]
        coreIds=options["coreIds"]
        resistanceRas=options["resistanceRas"]
        projectDir=options["projectDir"]
        useMinEucDist=options["useMinEucDist"]
        useMaxEucDist=options["useMaxEucDist"]
        maxEucDist=options["maxEucDist"]
        minEucDist=options["minEucDist"]
        useMinCostDist=options["useMinCostDist"]
        useMaxCostDist=options["useMaxCostDist"]
        maxCostDist=options["maxCostDist"]
        minCostDist=options["minCostDist"]
        saveNormLccs=options["saveNormLccs"]
        step3=options["step3"]
        step4=options["step4"]
        defaultNullValue = options["defaultNullValue"]
        # --------------------------------------------

        # --------------------------------------------
        # Constants
        coreAreaIds = coreIds # We now set variable coreAreaIds field name to user input ID field name  
        threshold = 100000 # This is how "wide" corridors will be (measured in cost-weighted distances) in truncated raster         
        # --------------------------------------------
        
        # --------------------------------------------
        # Folders, paths, filenames, and workspaces        
        outputRootDir = projectDir + "\\output\\"
        dataPassDir = projectDir + "\\dataPass\\"
        linkTableFile = get_prev_step_link_table(projectDir,step=5)      
        scratchDir = projectDir + "\\scratch\\"
        if os.path.exists(scratchDir)==False:
            gp.CreateFolder_management(projectDir,"scratch")
        gp.workspace=scratchDir 
        # --------------------------------------------
        
        if not maxEucDist:     
            useMaxEucDist = False
        else:
            useMaxEucDist = True
            gp.addmessage('Max Euclidean distance between cores')
            gp.addmessage('for linkage mapping set to ' + str(maxEucDist))
        
        if not maxCostDist: # == '#':
            useMaxCostDist=False
        else:
            useMaxCostDist=True
            gp.addmessage('Max cost-weighted distance between cores')
            gp.addmessage('for linkage mapping set to ' + str(maxCostDist))
       
        # ------------------------------------------------------------------
        ## linkTable column numbers
        linkIdCol=0 # Link ID
        core1Col=1 # core ID of 1st core area link connects
        core2Col=2 # core ID of 2nd core area link connects
        linkTypeCol=5 # 0=no link. 1=core, 2=corridor, 3=intermediate core area detected, 4=too long EucDist, 5=too long lcDist
        eucDistCol=6
        cwDistCol=7
        # ------------------------------------------------------------------

        # ------------------------------------------------------------------
        # set the analysis extent and cell size to that of the resistanceRas surface
        gp.Extent = gp.Describe(resistanceRas).Extent
        gp.CellSize = gp.Describe(resistanceRas).MeanCellHeight
        gp.Extent = "MINOF"
        gp.mask = resistanceRas
        gp.snapraster = resistanceRas
        # ------------------------------------------------------------------
        
        linkTable = load_link_table(linkTableFile)
        numLinks = linkTable.shape[0]        
        report_links(linkTable)


        if step3==False and step4 == False:        
            # re-check for links that are too long or in case script run out of sequence with more stringent settings
            gp.addmessage('Double-checking for corridors that are too long to map.')
            disableLeastCostNoVal = True
            linkTable,numDroppedLinks = drop_links(linkTable,useMaxEucDist,maxEucDist,useMaxCostDist,maxCostDist,useMinEucDist,minEucDist,useMinCostDist,minCostDist,disableLeastCostNoVal)        
                            
        # Added to try to speed up:
        gp.pyramid = "NONE"
        gp.rasterstatistics = "NONE"

        # ------------------------------------------------------------------
        # set up directories for normalized lcc and mosaic grids 
        dirCount=0
        baseLccPath = projectDir + "\\nlcc"
        gp.addmessage("Creating output folder: " + str(baseLccPath))    
        if os.path.exists(baseLccPath)==True:
            shutil.rmtree(baseLccPath)
        gp.CreateFolder_management(projectDir,"nlcc")    
        gp.CreateFolder_management(projectDir+"\\nlcc","nlc")    
        lccPath = projectDir + "\\nlcc\\nlc"
        mosaicPath = baseLccPath + "\\mosaic"
        if os.path.exists(mosaicPath)==True:
            shutil.rmtree(mosaicPath)
        gp.CreateFolder_management(baseLccPath,"mosaic")    
        mosaicRaster = mosaicPath +"\\nlcc_mos"           
        gp.addmessage('\nNormalized Least-cost corridors will be written to ' + lccPath + '\n')   
        # ------------------------------------------------------------------

        
        # ------------------------------------------------------------------
        # Add CWD layers for core area pairs to produce NORMALIZED LCC layers
        numGridsWritten = 0
        coreList=linkTable[:,core1Col:core2Col+1]
        coreList=sort(coreList)
        for x in range(0,numLinks):
            linkId = str(int(linkTable[x,linkIdCol]))
            if linkTable[x,linkTypeCol]==2 or linkTable[x,linkTypeCol]==10: # if it's a valid corridor link  
                # source and target cores
                corex=str(int(coreList[x,0]))
                corey=str(int(coreList[x,1]))
    
                # Get cwd rasters for source and target cores
                cwdRaster1 = get_cwd_path(projectDir,corex)
                cwdRaster2 = get_cwd_path(projectDir,corey)
    
                lccNormRaster = lccPath + "\\nlc" + str(corex) + "_" + str(corey)
                gp.Extent = "MINOF"

                # FIXME: need to check for this?:
                #if not gp.Exists(lccRaster): # if exists already, don't re-create

                link = get_links_from_core_pairs(linkTable,int(corex),int(corey))
                lcDist = str(linkTable[link,cwDistCol])

                # Normalized lcc rasters are created by adding cwd rasters and subtracting the least cost distance between them.
                expression = cwdRaster1 + " + " + cwdRaster2 + " - " + lcDist
                count,statement = 0, 'gp.SingleOutputMapAlgebra_sa(expression, lccNormRaster)'
                startTime = time.clock()
                while True:
                    try: exec statement
                    except:
                        count,tryAgain = hiccup_test(count,statement)
                        if not tryAgain: exec statement
                    else: break

                gp.Extent = "MAXOF"            
                if numGridsWritten == 0 and dirCount == 0: 
                    #If this is the first grid then copy rather than mosaic
                    gp.CopyRaster_management(lccNormRaster,mosaicRaster)
                else:                    
                    # Note: cannot use SOMA to mosaic.  It is a different process entirely.                    
                    count,statement = 0, 'gp.Mosaic_management(lccNormRaster,mosaicRaster,"MINIMUM","MATCH")'
                    while True:
                        try: exec statement
                        except:
                            count,tryAgain = hiccup_test(count,statement)
                            if not tryAgain: exec statement
                        else: break

                endTime = time.clock()
                processTime = round((endTime - startTime), 2)        
                gp.addmessage("Normalized and mosaicked corridor for link #" + str(linkId) + " connecting core areas " + str(corex) + " and " + str(corey)+ " in " + str(processTime) + " seconds.")              
                
                if saveNormLccs == False:
                    gp.delete_management(lccNormRaster)

                # temporarily disable links in linktable- don't want to mosaic them twice 
                for y in range (x+1,numLinks):
                    corex1=int(coreList[y,0])
                    corey1=int(coreList[y,1])
                    if int(corex1)==int(corex) and int(corey1)==int(corey):
                        linkTable[y,linkTypeCol]=linkTable[y,linkTypeCol]+100 
                    elif corex1==corey and corey1==corex:
                        linkTable[y,linkTypeCol]=linkTable[y,linkTypeCol]+100 

                
                numGridsWritten = numGridsWritten + 1
                if saveNormLccs == True:
                    if numGridsWritten == 100:
                        # We only write up to 100 grids to any one folder because
                        # otherwise Arc slows to a crawl
                        dirCount=dirCount+1
                        numGridsWritten = 0
                        gp.addmessage("Creating output folder: " + str(lccPath))    
                        gp.CreateFolder_management(projectDir+"\\nlcc","nlc"+str(dirCount))
                        lccPath = projectDir + "\\nlcc\\nlc"+str(dirCount)                        

        rows = where(linkTable[:,linkTypeCol]>100)#rows that were temporarily disabled 
        linkTable[rows,linkTypeCol]=linkTable[rows,linkTypeCol]-100
        # ---------------------------------------------------------------------
        
        # Create output geodatabase
        outputGDBname = "linkages"      
        outputGdb = outputRootDir + outputGDBname + ".gdb"
        gp.workspace = outputRootDir
        gp.createfilegdb(gp.Workspace + "\\", outputGDBname)

        # Copy mosaic raster to output geodatabase
        mosRaster = outputGdb + "\\lcc_mos"
        count,statement = 0, 'gp.CopyRaster_management(mosaicRaster,mosRaster)'
        while True:
            try: exec statement
            except:
                count,tryAgain = hiccup_test(count,statement)
                if not tryAgain: exec statement
            else: break
        
        gp.workspace = outputGdb
        gp.pyramid = "PYRAMIDS 2"
        gp.rasterStatistics = "STATISTICS 10 10"

        # ---------------------------------------------------------------------
        # convert mosaic raster to integer, set anything beyond threshold to NODATA.       
        truncRaster = "lcc_mosaic_meters_100000_max"
        expression = "(" + mosaicRaster + " * (con(" + mosaicRaster + "<= " + str(threshold) + ",1)))"
        count,statement = 0, 'gp.SingleOutputMapAlgebra_sa(expression, truncRaster)'
        while True:
            try: exec statement
            except:
                count,tryAgain = hiccup_test(count,statement)
                if not tryAgain: exec statement
            else: break
        kmRaster = "lcc_mosaic_km_100_max"
        expression = "float(int((" + truncRaster + ") / 10)) / 100"
        count,statement = 0, 'gp.SingleOutputMapAlgebra_sa(expression, kmRaster)'
        while True:
            try: exec statement
            except:
                count,tryAgain = hiccup_test(count,statement)
                if not tryAgain: exec statement
            else: break
        try:
            gp.delete_management(truncRaster)
        except:
            pass
        # ---------------------------------------------------------------------
        
        
        startTime = time.clock()

        if step4 == True:
            finalLinkTable = update_lcp_shapefile(projectDir,linkTable,lastStep=4,thisStep=5)
            startTime,hours,mins,secs = elapsed_time(startTime)
        elif step3 == True:
            finalLinkTable = update_lcp_shapefile(projectDir,linkTable,lastStep=3,thisStep=5)
            startTime,hours,mins,secs = elapsed_time(startTime)
        else:
            # Don't know if step 4 was run, since this is started at step 5.  Will look for step 4 lcp file, then step 3.  FIXME: this is one reason to remove old LCP files- or make a copy of step 3 with step 4 filename.  Otherwise could retrieve a step 4 file when a new run superceded it.
            finalLinkTable = update_lcp_shapefile(projectDir,linkTable,lastStep=4,thisStep=5) 
            startTime,hours,mins,secs = elapsed_time(startTime)
            
        linkTableFile = get_this_step_link_table(projectDir,step=5)       
        gp.addmessage('\nUpdating ' + linkTableFile)
        write_link_table(linkTable,coreIds,linkTableFile)

        linkTableLogFile=projectDir+"\\log\\"+"linkTable_step5.csv"   
        write_link_table(linkTable,coreIds,linkTableLogFile)   

        linkTableFinalFile=projectDir+"\\output\\"+"linkTable_Final.csv"   
        write_link_table(finalLinkTable,coreIds,linkTableFinalFile)   
        gp.addmessage('Copy of final linkTable written to '+linkTableFinalFile)
       
        # Pull out active corridor and constellation links
        numLinks=finalLinkTable.shape[0]
        rows,cols = where(finalLinkTable[:,linkTypeCol:linkTypeCol+1]== 2)
        coreLinks=finalLinkTable[rows,:]
        rows,cols=where(finalLinkTable[:,linkTypeCol:linkTypeCol+1]== 10)
        componentLinks=finalLinkTable[rows,:]
        activeLinkTable=append(coreLinks,componentLinks,axis=0)        
        del componentLinks
        del coreLinks
        
        ind=argsort((activeLinkTable[:,linkIdCol])) # sort by linkIdCol
        activeLinkTable = activeLinkTable[ind]
        
        activeLinkTableFile=projectDir+"\\output\\"+"linkTable_Final_Active_Links_Only.csv"   
        write_link_table(activeLinkTable,coreIds,activeLinkTableFile)   
        gp.addmessage('Table of active links written to '+activeLinkTableFile)
       
        gp.addmessage('---------------------')
        gp.addmessage('\nCreating shapefiles with linework for links.')
        write_link_maps(outputRootDir,linkTableFile,coreShapefile,coreIds,coreAreaIds,step=5)

        # Create final linkmap files in output directory, and remove files from scratch.
        copy_final_link_maps(projectDir)
        shutil.rmtree(scratchDir)
		
    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        gp.addmessage("\n--------------------")
        gp.addmessage('****Failed in step 5. Details follow.****')        
        filename =  __file__
        raise_geoproc_error(filename)
    
    # Return any PYTHON or system specific errors
    except:
        gp.addmessage("\n--------------------")
        gp.addmessage('****Failed in step 5. Details follow.****')        
        filename =  __file__
        raise_python_error(filename)
        
    del gp
    return 