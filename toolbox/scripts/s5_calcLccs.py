#!/usr/bin/env python2.5
# Authors: Brad McRae and Darren Kavanagh

"""Step 5: Calculate least cost corridors.

Creates and mosaics normalized least-cost corridors using connected core area
pairs specified in linkTable and cwd layers

"""

__filename__ = "s5_calcLccs.py"
__version__ = "0.7.0"

import os.path as path
import time
import shutil

import arcgisscripting
import numpy as npy

from lm_config import Config as Cfg
import lm_util as lu

gp = Cfg.gp
gprint = gp.addmessage


def STEP5_calc_lccs():
    """Creates and mosaics normalized least-cost corridors
    using connected core area pairs specified in linkTable and
    cwd layers

    """
    try:
        normalize = True
        calc_lccs(normalize)
        
        #Code to allow extra iteration to mosaic NON-normalized LCCs
        if Cfg.CALCNONNORMLCCS == True:
            normalize = False
            lu.dashline(1)
            gprint('\n**EXTRA STEP 5 RUN to mosaic NON-normalized corridors**')
            calc_lccs(normalize)
            
    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 5. Details follow.****')
        lu.raise_python_error(__filename__)

        
def calc_lccs(normalize):
    try:
        if normalize == True:
            mosaicBaseName = "_lcc_mosaic"
            writeTruncRaster = Cfg.WRITETRUNCRASTER
            outputGDB = Cfg.OUTPUTGDB
            if Cfg.CALCNONNORMLCCS:
                SAVENORMLCCS = False
            else:
                SAVENORMLCCS = Cfg.SAVENORMLCCS 
        else:
            mosaicBaseName = "_NON_NORMALIZED_lcc_mosaic"
            SAVENORMLCCS = False
            outputGDB = Cfg.EXTRAGDB
            writeTruncRaster = False

        lu.dashline(1)
        gprint('Running script' + __filename__)
        linkTableFile = lu.get_prev_step_link_table(step=5)
        gp.workspace = Cfg.SCRATCHDIR

        if Cfg.MAXEUCDIST is not None:
            gprint('Max Euclidean distance between cores')
            gprint('for linkage mapping set to ' +
                              str(Cfg.MAXEUCDIST))

        if Cfg.MAXCOSTDIST is not None:
            gprint('Max cost-weighted distance between cores')
            gprint('for linkage mapping set to ' +
                              str(Cfg.MAXCOSTDIST))


        # set the analysis extent and cell size to that of the resistance
        # surface
        gp.Extent = gp.Describe(Cfg.RESRAST).Extent
        gp.CellSize = gp.Describe(Cfg.RESRAST).MeanCellHeight
        gp.Extent = "MINOF"
        gp.mask = Cfg.RESRAST
        gp.snapraster = Cfg.RESRAST

        linkTable = lu.load_link_table(linkTableFile)
        numLinks = linkTable.shape[0]
        numCorridorLinks = lu.report_links(linkTable)
        if numCorridorLinks == 0:
            lu.dashline()
            gprint('\nThere are no corridors to map. Bailing.')
            time.sleep(5)
            return


        if not Cfg.STEP3 and not Cfg.STEP4:
            # re-check for links that are too long or in case script run out of
            # sequence with more stringent settings
            gprint('Double-checking for corridors that are too long'
                              ' to map.')
            disableLeastCostNoVal = True
            linkTable,numDroppedLinks = lu.drop_links(
                linkTable, Cfg.MAXEUCDIST, Cfg.MINEUCDIST, Cfg.MAXCOSTDIST,
                Cfg.MINCOSTDIST, disableLeastCostNoVal)

        # Added to try to speed up:
        gp.pyramid = "NONE"
        gp.rasterstatistics = "NONE"

        # set up directories for normalized lcc and mosaic grids
        dirCount = 0
        gprint("Creating output folder: " + Cfg.LCCBASEDIR)
        if path.exists(Cfg.LCCBASEDIR):
            shutil.rmtree(Cfg.LCCBASEDIR)
        gp.CreateFolder_management(path.dirname(Cfg.LCCBASEDIR),
                                       path.basename(Cfg.LCCBASEDIR))
        gp.CreateFolder_management(Cfg.LCCBASEDIR, Cfg.LCCNLCDIR_NM)
        clccdir = path.join(Cfg.LCCBASEDIR, Cfg.LCCNLCDIR_NM)
        gp.CreateFolder_management(Cfg.LCCBASEDIR,
                                       path.basename(Cfg.LCCMOSAICDIR))
        mosaicRaster = path.join(Cfg.LCCMOSAICDIR, "nlcc_mos")
        gprint("")
        if normalize == True:
            gprint('Normalized least-cost corridors will be written '
                          'to ' + clccdir + '\n')
        PREFIX = Cfg.PREFIX
        
        # Add CWD layers for core area pairs to produce NORMALIZED LCC layers
        numGridsWritten = 0
        coreList = linkTable[:,Cfg.LTB_CORE1:Cfg.LTB_CORE2+1]
        coreList = npy.sort(coreList)

        for x in range(0,numLinks):
            linkId = str(int(linkTable[x,Cfg.LTB_LINKID]))

            if (linkTable[x,Cfg.LTB_LINKTYPE] > 0):
                # source and target cores
                corex=int(coreList[x,0])
                corey=int(coreList[x,1])

                # Get cwd rasters for source and target cores
                cwdRaster1 = lu.get_cwd_path(corex)
                cwdRaster2 = lu.get_cwd_path(corey)

                lccNormRaster = path.join(clccdir, str(corex) + "_" +
                                          str(corey))
                gp.Extent = "MINOF"

                # FIXME: need to check for this?:
                # if exists already, don't re-create
                #if not gp.Exists(lccRaster):

                link = lu.get_links_from_core_pairs(linkTable, corex, corey)
                lcDist = str(linkTable[link,Cfg.LTB_CWDIST])

                # Normalized lcc rasters are created by adding cwd rasters and
                # subtracting the least cost distance between them.
                if normalize:
                    expression = cwdRaster1 + " + " + cwdRaster2 + " - " + lcDist
                else:
                    expression = cwdRaster1 + " + " + cwdRaster2 
                count = 0
                statement = ('gp.SingleOutputMapAlgebra_sa(expression, '
                            'lccNormRaster)')
                start_time = time.clock()
                while True:
                    try: exec statement
                    except:
                        count,tryAgain = lu.hiccup_test(count,statement)
                        if not tryAgain: exec statement
                    else: break

                gp.Extent = "MAXOF"
                if numGridsWritten == 0 and dirCount == 0:
                    #If this is the first grid then copy rather than mosaic
                    gp.CopyRaster_management(lccNormRaster, mosaicRaster)
                else:
                    # Note: cannot use SOMA to mosaic. It is a different
                    # process entirely.
                    count = 0
                    statement = ('gp.Mosaic_management(lccNormRaster, '
                                 'mosaicRaster, "MINIMUM", "MATCH")')
                    while True:
                        try: exec statement
                        except:
                            count,tryAgain = lu.hiccup_test(count,statement)
                            if not tryAgain: exec statement
                        else: break

                endTime = time.clock()
                processTime = round((endTime - start_time), 2)
                if normalize == True:
                    printText = "Normalized and mosaicked "
                else:
                    printText = "Mosaicked NON-normalized "
                gprint(printText + "corridor for link "
                                  "#" + str(linkId)
                                  + " connecting core areas " + str(corex) +
                                  " and " + str(corey)+ " in " +
                                  str(processTime) + " seconds.")

                if not SAVENORMLCCS:
                    lu.delete_data(lccNormRaster)

                # temporarily disable links in linktable - don't want to mosaic
                # them twice
                for y in range (x+1,numLinks):
                    corex1 = int(coreList[y,0])
                    corey1 = int(coreList[y,1])
                    if corex1 == corex and corey1 == corey:
                        linkTable[y,Cfg.LTB_LINKTYPE] = (
                            linkTable[y,Cfg.LTB_LINKTYPE] + 100)
                    elif corex1==corey and corey1==corex:
                        linkTable[y,Cfg.LTB_LINKTYPE] = (
                            linkTable[y,Cfg.LTB_LINKTYPE] + 100)


                numGridsWritten = numGridsWritten + 1
                if SAVENORMLCCS:
                    if numGridsWritten == 100:
                        # We only write up to 100 grids to any one folder
                        # because otherwise Arc slows to a crawl
                        dirCount = dirCount + 1
                        numGridsWritten = 0
                        clccdir = path.join(Cfg.LCCBASEDIR, Cfg.LCCNLCDIR_NM + str(dirCount))
                        #clccdir = path.join(clccdir, str(dirCount))
                        #clccdir = clccdir+str(dirCount)
                        gprint("Creating output folder: " + clccdir)
                        gp.CreateFolder_management(Cfg.LCCBASEDIR,
                                                       path.basename(clccdir))
        #rows that were temporarily disabled
        rows = npy.where(linkTable[:,Cfg.LTB_LINKTYPE]>100)
        linkTable[rows,Cfg.LTB_LINKTYPE] = (
            linkTable[rows,Cfg.LTB_LINKTYPE] - 100)
        # ---------------------------------------------------------------------

        # Create output geodatabase
        gp.createfilegdb(Cfg.OUTPUTDIR, path.basename(outputGDB))
        gp.workspace = outputGDB

        gp.pyramid = "NONE"
        gp.rasterstatistics = "NONE"
        
        # Copy mosaic raster to output geodatabase
        mosRaster = PREFIX + mosaicBaseName
        count = 0
        statement = 'gp.CopyRaster_management(mosaicRaster, mosRaster)'
        while True:
            try: exec statement
            except:
                count,tryAgain = lu.hiccup_test(count,statement)
                if not tryAgain:
                    exec statement
            else: break


        # ---------------------------------------------------------------------
        # convert mosaic raster to integer
        intRaster = PREFIX + mosaicBaseName + "_int"
        expression = "int(" + mosaicRaster + " + 0.5)"
        count = 0
        statement = 'gp.SingleOutputMapAlgebra_sa(expression, intRaster)'
        while True:
            try: exec statement
            except:
                count,tryAgain = lu.hiccup_test(count,statement)
                if not tryAgain: exec statement
            else: break
        # ---------------------------------------------------------------------

        # generate pyramids and statistics for final output 
        # DISABLE? Seems to build statistics too coarsely for clear
        # corridor display.
        try:
            gp.addmessage('\nBuilding output statistics and pyramids ' 
                              'for corridor raster\n')        
            gp.CalculateStatistics_management(intRaster, "1", "1", "#")
            gp.BuildPyramids_management(intRaster)    
        except:
            pass
        
        
        saveFloatRaster = False
        if saveFloatRaster == False:
            lu.delete_data(mosRaster)

        
        if writeTruncRaster == True:
            # ---------------------------------------------------------------------
            # Set anything beyond Cfg.CWDTHRESH to NODATA.
            truncRaster = PREFIX + mosaicBaseName + "_truncated_values"
            expression = ("(" + intRaster + " * (con(" + intRaster + "<= " +
                          str(Cfg.CWDTHRESH) + ",1)))")
            count = 0
            statement = 'gp.SingleOutputMapAlgebra_sa(expression, truncRaster)'
            while True:
                try: exec statement
                except:
                    count,tryAgain = lu.hiccup_test(count,statement)
                    if not tryAgain: exec statement
                else: break
        # ---------------------------------------------------------------------
            
            try:
                gp.addmessage('Building output statistics and pyramids ' 
                              'for truncated corridor raster\n')        
                gp.CalculateStatistics_management(intRaster, "1", "1", "#")
                gp.BuildPyramids_management(intRaster)    
            except:
                pass
            
        start_time = time.clock()
        gprint('Writing final LCP maps...')
        if Cfg.STEP4:
            finalLinkTable = lu.update_lcp_shapefile(linkTable, lastStep=4,
                                                     thisStep=5)
        elif Cfg.STEP3:
            finalLinkTable = lu.update_lcp_shapefile(linkTable, lastStep=3,
                                                     thisStep=5)
        else:
            # Don't know if step 4 was run, since this is started at step 5.
            # Use presence of previous linktable files to figure this out.
            # Linktable name includes step number.
            prevLinkTableFile = lu.get_prev_step_link_table(step=5)
            prevStepInd = len(prevLinkTableFile)-5 
            lastStep = prevLinkTableFile[prevStepInd]
        
            finalLinkTable = lu.update_lcp_shapefile(linkTable, lastStep,
                                                     thisStep=5)

        outlinkTableFile = lu.get_this_step_link_table(step=5)
        gprint('Updating ' + outlinkTableFile)
        lu.write_link_table(linkTable, outlinkTableFile)

        linkTableLogFile = path.join(Cfg.LOGDIR, "linkTable_s5.csv")
        lu.write_link_table(linkTable, linkTableLogFile)

        linkTableFinalFile = path.join(Cfg.OUTPUTDIR, PREFIX + "_linkTable_s5.csv")
        lu.write_link_table(finalLinkTable, linkTableFinalFile)
        gprint('Copy of final linkTable written to '+
                          linkTableFinalFile)

        gprint('Creating shapefiles with linework for links.')
        lu.write_link_maps(outlinkTableFile, step=5)

        # Create final linkmap files in output directory, and remove files from
        # scratch.
        lu.copy_final_link_maps(step=5)
       
        # Check for unreasonably low minimum NLCC values
        propertyType = "MINIMUM"
        minObject = gp.GetRasterProperties(mosaicRaster, propertyType)
        rasterMin = float(str(minObject.getoutput(0)))
        tolerance = (float(gp.CellSize) * -10)

        if not SAVENORMLCCS:
            lu.delete_dir(Cfg.LCCBASEDIR)
            
        if rasterMin < tolerance:
            lu.dashline(1)
            msg = ('WARNING: Minimum value of mosaicked corridor map is ' 
                   'much less than zero ('+str(rasterMin)+').'
                   '\nThis could mean that BOUNDING CIRCLE BUFFER DISTANCES '
                   'were too small and a corridor passed outside of a '
                   'bounding circle. All steps were completed. ')
            gp.AddError(msg)
            exit(1) 

        
    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 5. Details follow.****')
        lu.raise_geoproc_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 5. Details follow.****')
        lu.raise_python_error(__filename__)

    return
