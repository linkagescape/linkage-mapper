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
import sys
import os.path as path
import time
import shutil

import arcgisscripting
from numpy import *

import lm_config
import lm_util as lu

OUTPUTDIR = lm_config.OUTPUTDIR
LOGDIR = lm_config.LOGDIR
SCRATCHDIR = lm_config.SCRATCHDIR
LCCMOSAICDIR = lm_config.LCCMOSAICDIR
LCCBASEDIR = lm_config.LCCBASEDIR
LCCNLCDIR_NAME = lm_config.LCCNLCDIR_NAME
RESRAST = lm_config.RESRAST
STEP3 = lm_config.STEP3
STEP4 = lm_config.STEP4
OUTPUTGDB = lm_config.OUTPUTGDB
MAXCOSTDIST = lm_config.MAXCOSTDIST
MINCOSTDIST = lm_config.MINCOSTDIST
MAXEUCDIST = lm_config.MAXEUCDIST   
MINEUCDIST = lm_config.MINEUCDIST
CWDTHRESH = lm_config.CWDTHRESH
SAVENORMLCCS = lm_config.SAVENORMLCCS
GP = lm_config.GP

def step5_calc_lccs():
    """Creates and mosaics normalized least-cost corridors
    using connected core area pairs specified in linkTable and
    cwd layers
    
    """

# Fixme: add option to saveRawLccs? Or mosaicLccs that already exist?
    try:
        global GP
                
        linkTableFile = lu.get_prev_step_link_table(step=5)      
        GP.workspace = SCRATCHDIR        
        
        if MAXEUCDIST is not None:              
            GP.addmessage('Max Euclidean distance between cores')
            GP.addmessage('for linkage mapping set to ' + str(MAXEUCDIST))
        
        if MAXCOSTDIST is not None:
            GP.addmessage('Max cost-weighted distance between cores')
            GP.addmessage('for linkage mapping set to ' + str(MAXCOSTDIST))
       
        ## linkTable column numbers
        linkIdCol = 0 # Link ID
        core1Col = 1 # core ID of 1st core area link connects
        core2Col = 2 # core ID of 2nd core area link connects
        linkTypeCol = 5 # 0=no link. 1=core, 2=corridor, 3=intermediate core area detected, 4=too long EucDist, 5=too long lcDist
        eucDistCol = 6
        cwDistCol = 7

        # set the analysis extent and cell size to that of the resistance 
        # surface
        GP.Extent = GP.Describe(RESRAST).Extent
        GP.CellSize = GP.Describe(RESRAST).MeanCellHeight
        GP.Extent = "MINOF"
        GP.mask = RESRAST
        GP.snapraster = RESRAST
        
        linkTable = lu.load_link_table(linkTableFile)
        numLinks = linkTable.shape[0]        
        lu.report_links(linkTable)

        if not STEP3 and not STEP4:
            # re-check for links that are too long or in case script run out of sequence with more stringent settings
            GP.addmessage('Double-checking for corridors that are too long to map.')
            disableLeastCostNoVal = True
            linkTable,numDroppedLinks = lu.drop_links(linkTable, MAXEUCDIST, 
                                               MINEUCDIST, MAXCOSTDIST, 
                                               MINCOSTDIST, 
                                               disableLeastCostNoVal)     
                            
        # Added to try to speed up:
        GP.pyramid = "NONE"
        GP.rasterstatistics = "NONE"

        # set up directories for normalized lcc and mosaic grids 
        dirCount = 0        
        GP.addmessage("Creating output folder: " + LCCBASEDIR)        
        if path.exists(LCCBASEDIR):
            shutil.rmtree(LCCBASEDIR)
        GP.CreateFolder_management(path.dirname(LCCBASEDIR), 
                                   path.basename(LCCBASEDIR))                 
        GP.CreateFolder_management(LCCBASEDIR, LCCNLCDIR_NAME)
        clccdir = path.join(LCCBASEDIR, LCCNLCDIR_NAME)              
        GP.CreateFolder_management(LCCBASEDIR, 
                                   path.basename(LCCMOSAICDIR))                                     
        mosaicRaster = path.join(LCCMOSAICDIR, "nlcc_mos")
        GP.addmessage('\nNormalized Least-cost corridors will be written to ' 
                      + clccdir + '\n')   
        
        # Add CWD layers for core area pairs to produce NORMALIZED LCC layers
        numGridsWritten = 0
        coreList = linkTable[:,core1Col:core2Col+1]
        coreList = sort(coreList)        
        
        for x in range(0,numLinks):
            linkId = str(int(linkTable[x,linkIdCol]))
            if linkTable[x,linkTypeCol]==2 or linkTable[x,linkTypeCol]==10: # if it's a valid corridor link  
                # source and target cores
                corex=int(coreList[x,0])
                corey=int(coreList[x,1])
    
                # Get cwd rasters for source and target cores
                cwdRaster1 = lu.get_cwd_path(corex)
                cwdRaster2 = lu.get_cwd_path(corey)
    
                lccNormRaster = path.join(clccdir, str(corex) + "_" 
                                          + str(corey))
                GP.Extent = "MINOF"

                # FIXME: need to check for this?:
                #if not GP.Exists(lccRaster): # if exists already, don't re-create

                link = lu.get_links_from_core_pairs(linkTable, corex, corey)
                lcDist = str(linkTable[link,cwDistCol])

                # Normalized lcc rasters are created by adding cwd rasters and subtracting the least cost distance between them.
                expression = cwdRaster1 + " + " + cwdRaster2 + " - " + lcDist
                count,statement = 0, 'GP.SingleOutputMapAlgebra_sa(expression, lccNormRaster)'
                startTime = time.clock()
                while True:
                    try: exec statement
                    except:
                        count,tryAgain = lu.hiccup_test(count,statement)
                        if not tryAgain: exec statement
                    else: break

                GP.Extent = "MAXOF"            
                if numGridsWritten == 0 and dirCount == 0: 
                    #If this is the first grid then copy rather than mosaic
                    GP.CopyRaster_management(lccNormRaster,mosaicRaster)
                else:                    
                    # Note: cannot use SOMA to mosaic. It is a different 
                    # process entirely.                    
                    count = 0
                    statement = ('GP.Mosaic_management(lccNormRaster, '
                                 'mosaicRaster, "MINIMUM", "MATCH")')
                    while True:
                        try: exec statement
                        except:
                            count,tryAgain = lu.hiccup_test(count,statement)
                            if not tryAgain: exec statement
                        else: break

                endTime = time.clock()
                processTime = round((endTime - startTime), 2)        
                GP.addmessage("Normalized and mosaicked corridor for link #" 
                              + str(linkId) + " connecting core areas " 
                              + str(corex) + " and " + str(corey)+ " in " 
                              + str(processTime) + " seconds.")              
                
                if not SAVENORMLCCS:
                    GP.delete_management(lccNormRaster)

                # temporarily disable links in linktable - don't want to mosaic 
                # them twice 
                for y in range (x+1,numLinks):
                    corex1=int(coreList[y,0])
                    corey1=int(coreList[y,1])
                    if corex1 == corex and corey1 == corey:
                        linkTable[y,linkTypeCol]=linkTable[y,linkTypeCol]+100 
                    elif corex1==corey and corey1==corex:
                        linkTable[y,linkTypeCol]=linkTable[y,linkTypeCol]+100 

                
                numGridsWritten = numGridsWritten + 1
                if SAVENORMLCCS:
                    if numGridsWritten == 100:
                        # We only write up to 100 grids to any one folder because
                        # otherwise Arc slows to a crawl
                        dirCount=dirCount+1
                        numGridsWritten = 0
                        clccdir = path.join(clccdir, str(dirCount))
                        GP.addmessage("Creating output folder: " + clccdir)                   
                        GP.CreateFolder_management(LCCBASEDIR, 
                                                   path.basename(clccdir))

        rows = where(linkTable[:,linkTypeCol]>100)#rows that were temporarily disabled 
        linkTable[rows,linkTypeCol]=linkTable[rows,linkTypeCol]-100
        # ---------------------------------------------------------------------
        
        # Create output geodatabase                        
        GP.createfilegdb(OUTPUTDIR, path.basename(OUTPUTGDB))
        GP.workspace = OUTPUTGDB
        
        # Copy mosaic raster to output geodatabase
        mosRaster = "lcc_mos"
        count,statement = 0, 'GP.CopyRaster_management(mosaicRaster,mosRaster)'
        while True:
            try: exec statement
            except:
                count,tryAgain = lu.hiccup_test(count,statement)
                if not tryAgain: exec statement
            else: break
                
        GP.pyramid = "PYRAMIDS 2"
        GP.rasterStatistics = "STATISTICS 10 10"

        # ---------------------------------------------------------------------
        # convert mosaic raster to integer, set anything beyond CWDTHRESH to NODATA.       
        truncRaster = "lcc_mosaic_meters_100000_max"
        expression = "(" + mosaicRaster + " * (con(" + mosaicRaster + "<= " + str(CWDTHRESH) + ",1)))"
        count,statement = 0, 'GP.SingleOutputMapAlgebra_sa(expression, truncRaster)'
        while True:
            try: exec statement
            except:
                count,tryAgain = lu.hiccup_test(count,statement)
                if not tryAgain: exec statement
            else: break
        kmRaster = "lcc_mosaic_km_100_max"
        expression = "float(int((" + truncRaster + ") / 10)) / 100"
        count,statement = 0, 'GP.SingleOutputMapAlgebra_sa(expression, kmRaster)'
        while True:
            try: exec statement
            except:
                count,tryAgain = lu.hiccup_test(count,statement)
                if not tryAgain: exec statement
            else: break
        try:
            GP.delete_management(truncRaster)
        except:
            pass
        # ---------------------------------------------------------------------
        
        
        startTime = time.clock()

        if STEP4:
            finalLinkTable = lu.update_lcp_shapefile(linkTable, lastStep=4, 
                                                  thisStep=5)
            startTime, hours, mins, secs = lu.elapsed_time(startTime)
        elif STEP3:
            finalLinkTable = lu.update_lcp_shapefile(linkTable, lastStep=3, 
                                                  thisStep=5)
            startTime, hours, mins, secs = lu.elapsed_time(startTime)
        else:
            # Don't know if step 4 was run, since this is started at step 5.  
            # Will look for step 4 lcp file, then step 3.  FIXME: this is one 
            # reason to remove old LCP files- or make a copy of step 3 with 
            # step 4 filename.  Otherwise could retrieve a step 4 file when a 
            # new run superceded it.
            finalLinkTable = lu.update_lcp_shapefile(linkTable, lastStep=4, 
                                                  thisStep=5) 
            startTime, hours, mins, secs = lu.elapsed_time(startTime)
            
        linkTableFile = lu.get_this_step_link_table(step=5)       
        GP.addmessage('\nUpdating ' + linkTableFile)
        lu.write_link_table(linkTable, linkTableFile)

        linkTableLogFile = path.join(LOGDIR, "linkTable_step5.csv")
        lu.write_link_table(linkTable, linkTableLogFile)   

        linkTableFinalFile = path.join(OUTPUTDIR, "linkTable_Final.csv")
        lu.write_link_table(finalLinkTable, linkTableFinalFile)   
        GP.addmessage('Copy of final linkTable written to ' 
                      + linkTableFinalFile)
       
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
        
        activeLinkTableFile = path.join(OUTPUTDIR, 
                              "linkTable_Final_Active_Links_Only.csv")
        lu.write_link_table(activeLinkTable, activeLinkTableFile)   
        GP.addmessage('Table of active links written to '+ activeLinkTableFile)
       
        lu.dashline()
        GP.addmessage('\nCreating shapefiles with linework for links.')
        lu.write_link_maps(linkTableFile, step=5)

        # Create final linkmap files in output directory, and remove files from scratch.
        lu.copy_final_link_maps()
        shutil.rmtree(SCRATCHDIR)
		
    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        GP.addmessage('****Failed in step 5. Details follow.****')        
        filename =  __file__
        lu.raise_geoproc_error(filename)
    
    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        GP.addmessage('****Failed in step 5. Details follow.****')        
        filename =  __file__
        lu.raise_python_error(filename)
        
    return 