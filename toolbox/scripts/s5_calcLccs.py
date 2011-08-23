#!/usr/bin/env python2.5

"""Step 5: Calculate least cost corridors.

Creates and mosaics normalized least-cost corridors using connected core area
pairs specified in linkTable and cwd layers

"""

__filename__ = "s5_calcLccs.py"
__version__ = "0.6.4"

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

# Fixme: add option to saveRawLccs? Or mosaicLccs that already exist?
    try:
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
        gprint('\nNormalized Least-cost corridors will be written '
                          'to ' + clccdir + '\n')
        prefix = path.basename(Cfg.PROJECTDIR)
        
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
                expression = cwdRaster1 + " + " + cwdRaster2 + " - " + lcDist
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
                gprint("Normalized and mosaicked corridor for link "
                                  "#" + str(linkId)
                                  + " connecting core areas " + str(corex) +
                                  " and " + str(corey)+ " in " +
                                  str(processTime) + " seconds.")

                if not Cfg.SAVENORMLCCS:
                    gp.delete_management(lccNormRaster)

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
                if Cfg.SAVENORMLCCS:
                    if numGridsWritten == 100:
                        # We only write up to 100 grids to any one folder
                        # because otherwise Arc slows to a crawl
                        dirCount = dirCount + 1
                        numGridsWritten = 0
                        clccdir = path.join(clccdir, str(dirCount))
                        gprint("Creating output folder: " + clccdir)
                        gp.CreateFolder_management(Cfg.LCCBASEDIR,
                                                       path.basename(clccdir))
        #rows that were temporarily disabled
        rows = npy.where(linkTable[:,Cfg.LTB_LINKTYPE]>100)
        linkTable[rows,Cfg.LTB_LINKTYPE] = (
            linkTable[rows,Cfg.LTB_LINKTYPE] - 100)
        # ---------------------------------------------------------------------

        # Create output geodatabase
        gp.createfilegdb(Cfg.OUTPUTDIR, path.basename(Cfg.OUTPUTGDB))
        gp.workspace = Cfg.OUTPUTGDB
        
        # Copy mosaic raster to output geodatabase
        mosRaster = prefix + "_lcc_mosaic"
        count = 0
        statement = 'gp.CopyRaster_management(mosaicRaster, mosRaster)'
        while True:
            try: exec statement
            except:
                count,tryAgain = lu.hiccup_test(count,statement)
                if not tryAgain:
                    exec statement
            else: break

        gp.pyramid = "PYRAMIDS 2"
        gp.rasterStatistics = "STATISTICS 10 10"

        # ---------------------------------------------------------------------
        # convert mosaic raster to integer, set anything beyond Cfg.CWDTHRESH
        # to NODATA.
        truncRaster = "lcc_mosaic_100k_max"
        expression = ("(" + mosaicRaster + " * (con(" + mosaicRaster + "<= " +
                      str(Cfg.CWDTHRESH) + ",1)))")
        count = 0
        statement = 'gp.SingleOutputMapAlgebra_sa(expression, truncRaster)'
        while True:
            try: exec statement
            except:
                count,tryAgain = lu.hiccup_test(count,statement)
                if not tryAgain: exec statement
            else: break
        intRaster = prefix + "_lcc_mosaic_100k_max_int"
        expression = "int(" + truncRaster + ")"
        count = 0
        statement = 'gp.SingleOutputMapAlgebra_sa(expression, intRaster)'
        while True:
            try: exec statement
            except:
                count,tryAgain = lu.hiccup_test(count,statement)
                if not tryAgain: exec statement
            else: break
        try:
            gp.delete_management(truncRaster)
        except:
            pass
        # ---------------------------------------------------------------------

        start_time = time.clock()

        if Cfg.STEP4:
            finalLinkTable = lu.update_lcp_shapefile(linkTable, lastStep=4,
                                                     thisStep=5)
            start_time = lu.elapsed_time(start_time)
        elif Cfg.STEP3:
            finalLinkTable = lu.update_lcp_shapefile(linkTable, lastStep=3,
                                                     thisStep=5)
            start_time = lu.elapsed_time(start_time)
        else:
            # Don't know if step 4 was run, since this is started at step 5.
            # Will look for step 4 lcp file, then step 3.  FIXME: this is one
            # reason to remove old LCP files- or make a copy of step 3 with
            # step 4 filename.  Otherwise could retrieve a step 4 file when a
            # new run superceded it.
            finalLinkTable = lu.update_lcp_shapefile(linkTable, lastStep=4,
                                                     thisStep=5)
            start_time = lu.elapsed_time(start_time)

        outlinkTableFile = lu.get_this_step_link_table(step=5)
        gprint('Updating ' + outlinkTableFile)
        lu.write_link_table(linkTable, outlinkTableFile)

        linkTableLogFile = path.join(Cfg.LOGDIR, "linkTable_s5.csv")
        lu.write_link_table(linkTable, linkTableLogFile)

        linkTableFinalFile = path.join(Cfg.OUTPUTDIR, prefix + "_linkTable_Final.csv")
        lu.write_link_table(finalLinkTable, linkTableFinalFile)
        gprint('Copy of final linkTable written to '+
                          linkTableFinalFile)

        # # Pull out active corridor and constellation links
        # numLinks = finalLinkTable.shape[0]
        # rows, cols = npy.where(
            # finalLinkTable[:,Cfg.LTB_LINKTYPE:Cfg.LTB_LINKTYPE+1] ==
            # Cfg.LT_CORR)
        # coreLinks = finalLinkTable[rows,:]
        # rows, cols = npy.where(
            # finalLinkTable[:,Cfg.LTB_LINKTYPE:Cfg.LTB_LINKTYPE+1] ==
            # Cfg.LT_CLU)
        # componentLinks = finalLinkTable[rows,:]
        # activeLinkTable = npy.append(coreLinks, componentLinks, axis=0)
        # del componentLinks
        # del coreLinks

        # # sort by Cfg.LTB_LINKID
        # ind = argsort((activeLinkTable[:,Cfg.LTB_LINKID]))
        # activeLinkTable = activeLinkTable[ind]

        # activeLinkTableFile = path.join(
            # Cfg.OUTPUTDIR, "linkTable_.csv")
        # lu.write_link_table(activeLinkTable, activeLinkTableFile)
        # gprint('Table of active links written to ' +
                          # activeLinkTableFile)

        # lu.dashline(1)
        gprint('Creating shapefiles with linework for links.')
        lu.write_link_maps(outlinkTableFile, step=5)

        # Create final linkmap files in output directory, and remove files from
        # scratch.
        lu.copy_final_link_maps()

        # Clean up
        lu.delete_dir(Cfg.SCRATCHDIR)
        if not Cfg.SAVENORMLCCS:
            lu.delete_dir(Cfg.LCCBASEDIR)

            
    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 5. Details follow.****')
        lu.raise_python_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 5. Details follow.****')
        lu.raise_python_error(__filename__)

    return