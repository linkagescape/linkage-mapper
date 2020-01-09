# Authors: Brad McRae and Darren Kavanagh

"""Step 5: Calculate least cost corridors.

Creates and mosaics normalized least-cost corridors using connected core area
pairs specified in linkTable and cwd layers

"""

from os import path
import time

import numpy as npy
import arcpy

from lm_config import tool_env as cfg
import lm_util as lu

_SCRIPT_NAME = "s5_calcLccs.py"

gprint = lu.gprint


def STEP5_calc_lccs():
    """Creates and mosaics normalized least-cost corridors
    using connected core area pairs specified in linkTable and
    cwd layers

    """
    try:

        normalize = True
        calc_lccs(normalize)

        #Code to allow extra iteration to mosaic NON-normalized LCCs
        if cfg.CALCNONNORMLCCS:
            normalize = False
            lu.dashline(1)
            gprint('\n**EXTRA STEP 5 RUN to mosaic NON-normalized corridors**')
            calc_lccs(normalize)

    # Return any PYTHON or system specific errors
    except Exception:
        lu.dashline(1)
        gprint('****Failed in step 5. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)


def calc_lccs(normalize):
    try:
        if normalize:
            mosaicBaseName = "_corridors"
            writeTruncRaster = cfg.WRITETRUNCRASTER
            outputGDB = cfg.OUTPUTGDB
            SAVENORMLCCS = cfg.SAVENORMLCCS
        else:
            mosaicBaseName = "_NON_NORMALIZED_corridors"
            SAVENORMLCCS = False
            outputGDB = cfg.EXTRAGDB
            writeTruncRaster = False

        lu.dashline(1)
        gprint('Running script ' + _SCRIPT_NAME)
        linkTableFile = lu.get_prev_step_link_table(step=5)
        arcpy.env.workspace = cfg.SCRATCHDIR
        arcpy.env.scratchWorkspace = cfg.ARCSCRATCHDIR
        arcpy.env.compression = "NONE"

        if cfg.MAXEUCDIST is not None:
            gprint('Max Euclidean distance between cores')
            gprint('for linkage mapping set to ' +
                              str(cfg.MAXEUCDIST))

        if cfg.MAXCOSTDIST is not None:
            gprint('Max cost-weighted distance between cores')
            gprint('for linkage mapping set to ' +
                              str(cfg.MAXCOSTDIST))


        # set the analysis extent and cell size to that of the resistance
        # surface
        arcpy.env.extent = cfg.RESRAST
        arcpy.env.cellSize = cfg.RESRAST
        arcpy.env.snapRaster = cfg.RESRAST
        arcpy.env.mask = cfg.RESRAST

        linkTable = lu.load_link_table(linkTableFile)
        numLinks = linkTable.shape[0]
        numCorridorLinks = lu.report_links(linkTable)
        if numCorridorLinks == 0:
            lu.dashline(1)
            msg =('\nThere are no corridors to map. Bailing.')
            lu.raise_error(msg)


        if not cfg.STEP3 and not cfg.STEP4:
            # re-check for links that are too long or in case script run out of
            # sequence with more stringent settings
            gprint('Double-checking for corridors that are too long to map.')
            DISABLE_LEAST_COST_NO_VAL = True
            linkTable,numDroppedLinks = lu.drop_links(
                linkTable, cfg.MAXEUCDIST, cfg.MINEUCDIST, cfg.MAXCOSTDIST,
                cfg.MINCOSTDIST, DISABLE_LEAST_COST_NO_VAL)

        # Added to try to speed up:
        arcpy.env.pyramid = "NONE"
        arcpy.env.rasterStatistics = "NONE"

        # set up directories for normalized lcc and mosaic grids
        dirCount = 0
        gprint("Creating output folder: " + cfg.LCCBASEDIR)
        lu.delete_dir(cfg.LCCBASEDIR)
        arcpy.CreateFolder_management(path.dirname(cfg.LCCBASEDIR),
                                       path.basename(cfg.LCCBASEDIR))
        arcpy.CreateFolder_management(cfg.LCCBASEDIR, cfg.LCCNLCDIR_NM)
        clccdir = path.join(cfg.LCCBASEDIR, cfg.LCCNLCDIR_NM)
        gprint("")
        if normalize:
            gprint('Normalized least-cost corridors will be written '
                          'to ' + clccdir + '\n')
        PREFIX = cfg.PREFIX

        # Add CWD layers for core area pairs to produce NORMALIZED LCC layers
        numGridsWritten = 0
        coreList = linkTable[:,cfg.LTB_CORE1:cfg.LTB_CORE2+1]
        coreList = npy.sort(coreList)

        x = 0
        linkCount = 0
        endIndex = numLinks
        while x < endIndex:
            if (linkTable[x, cfg.LTB_LINKTYPE] < 1): # If not a valid link
                x = x + 1
                continue

            linkCount = linkCount + 1
            start_time = time.clock()

            linkId = str(int(linkTable[x, cfg.LTB_LINKID]))

            # source and target cores
            corex=int(coreList[x,0])
            corey=int(coreList[x,1])

            # Get cwd rasters for source and target cores
            cwdRaster1 = lu.get_cwd_path(corex)
            cwdRaster2 = lu.get_cwd_path(corey)

            if not arcpy.Exists(cwdRaster1):
                msg =('\nError: cannot find cwd raster:\n' + cwdRaster1)
            if not arcpy.Exists(cwdRaster2):
                msg =('\nError: cannot find cwd raster:\n' + cwdRaster2)
                lu.raise_error(msg)


            lccNormRaster = path.join(clccdir, str(corex) + "_" +
                                      str(corey))# + ".tif")
            arcpy.env.extent = "MINOF"

            link = lu.get_links_from_core_pairs(linkTable, corex, corey)

            offset = 10000

            # Normalized lcc rasters are created by adding cwd rasters and
            # subtracting the least cost distance between them.
            lcDist = (float(linkTable[link,cfg.LTB_CWDIST]) - offset)

            if normalize:
                statement = ('outras = arcpy.sa.Raster(cwdRaster1) '
                             '+ arcpy.sa.Raster(cwdRaster2) - lcDist; '
                             'outras.save(lccNormRaster)')
            else:
                statement = ('outras = arcpy.sa.Raster(cwdRaster1) '
                             '+ arcpy.sa.Raster(cwdRaster2); '
                             'outras.save(lccNormRaster)')

            count = 0
            while True:
                try:
                    exec(statement)
                except Exception:
                    count,tryAgain = lu.retry_arc_error(count,statement)
                    if not tryAgain:
                        exec(statement)
                else: break

            if normalize:
                try:
                    minObject = arcpy.GetRasterProperties_management(lccNormRaster, "MINIMUM")
                    rasterMin = float(str(minObject.getOutput(0)))
                except Exception:
                    lu.warn('\n------------------------------------------------')
                    lu.warn('WARNING: Raster minimum check failed in step 5. \n'
                        'This may mean the output rasters are corrupted. Please \n'
                        'be sure to check for valid rasters in '+ outputGDB)
                    rasterMin = 0
                tolerance = (float(arcpy.env.cellSize) * -10)
                if rasterMin < tolerance:
                    lu.dashline(1)
                    msg = ('WARNING: Minimum value of a corridor #' + str(x+1)
                           + ' is much less than zero ('+str(rasterMin)+').'
                           '\nThis could mean that BOUNDING CIRCLE BUFFER DISTANCES '
                           'were too small and a corridor passed outside of a '
                           'bounding circle, or that a corridor passed outside of the '
                           'resistance map. \n')
                    lu.warn(msg)

            arcpy.env.extent = cfg.RESRAST

            mosaicDir = path.join(cfg.LCCBASEDIR,'mos'+str(x+1))
            lu.create_dir(mosaicDir)
            mosFN = 'mos'#.tif' change and move
            mosaicRaster = path.join(mosaicDir,mosFN)

            if numGridsWritten == 0 and dirCount == 0:
                #If this is the first grid then copy rather than mosaic
                arcpy.CopyRaster_management(lccNormRaster, mosaicRaster)
            else:
                statement = (
                    'arcpy.MosaicToNewRaster_management('
                    'input_rasters=";".join([lccNormRaster, '
                    'lastMosaicRaster]), output_location=mosaicDir, '
                    'raster_dataset_name_with_extension=mosFN, '
                    'pixel_type="32_BIT_FLOAT", cellsize=arcpy.env.cellSize, '
                    'number_of_bands="1", mosaic_method="MINIMUM", '
                    'mosaic_colormap_mode="MATCH")')

                count = 0
                while True:
                    try:
                        lu.write_log('Executing mosaic for link #'+str(linkId))
                        exec(statement)
                        lu.write_log('Done with mosaic.')
                    except Exception:
                        count,tryAgain = lu.retry_arc_error(count,statement)
                        lu.delete_data(mosaicRaster)
                        lu.delete_dir(mosaicDir)
                        # Try a new directory
                        mosaicDir = path.join(cfg.LCCBASEDIR,'mos'+str(x+1)+ '_' + str(count))
                        lu.create_dir(mosaicDir)
                        mosaicRaster = path.join(mosaicDir,mosFN)
                        if not tryAgain:
                            exec(statement)
                    else: break
            endTime = time.clock()
            processTime = round((endTime - start_time), 2)

            if normalize == True:
                printText = "Normalized and mosaicked "
            else:
                printText = "Mosaicked NON-normalized "
            gprint(printText + "corridor for link ID #" + str(linkId) +
                    " connecting core areas " + str(corex) +
                    " and " + str(corey)+ " in " +
                    str(processTime) + " seconds. " + str(int(linkCount)) +
                    " out of " + str(int(numCorridorLinks)) + " links have been "
                    "processed.")

            # temporarily disable links in linktable - don't want to mosaic
            # them twice
            for y in range (x+1,numLinks):
                corex1 = int(coreList[y,0])
                corey1 = int(coreList[y,1])
                if corex1 == corex and corey1 == corey:
                    linkTable[y,cfg.LTB_LINKTYPE] = (
                        linkTable[y,cfg.LTB_LINKTYPE] + 1000)
                elif corex1==corey and corey1==corex:
                    linkTable[y,cfg.LTB_LINKTYPE] = (
                            linkTable[y,cfg.LTB_LINKTYPE] + 1000)

            numGridsWritten = numGridsWritten + 1
            if not SAVENORMLCCS:
                lu.delete_data(lccNormRaster)
                lu.delete_dir(clccdir)
                lu.create_dir(clccdir)
            else:
                if numGridsWritten == 100:
                    # We only write up to 100 grids to any one folder
                    # because otherwise Arc slows to a crawl
                    dirCount = dirCount + 1
                    numGridsWritten = 0
                    clccdir = path.join(cfg.LCCBASEDIR,
                                        cfg.LCCNLCDIR_NM + str(dirCount))
                    gprint("Creating output folder: " + clccdir)
                    arcpy.CreateFolder_management(cfg.LCCBASEDIR,
                                               path.basename(clccdir))

            if numGridsWritten > 1 or dirCount > 0:
                lu.delete_data(lastMosaicRaster)
                lu.delete_dir(path.dirname(lastMosaicRaster))

            lastMosaicRaster = mosaicRaster
            x = x + 1

        #rows that were temporarily disabled
        rows = npy.where(linkTable[:,cfg.LTB_LINKTYPE]>1000)
        linkTable[rows,cfg.LTB_LINKTYPE] = (
            linkTable[rows,cfg.LTB_LINKTYPE] - 1000)
        # ---------------------------------------------------------------------

        # Create output geodatabase
        if not arcpy.Exists(outputGDB):
            arcpy.CreateFileGDB_management(cfg.OUTPUTDIR, path.basename(outputGDB))

        arcpy.env.workspace = outputGDB

        arcpy.env.pyramid = "NONE"
        arcpy.env.rasterStatistics = "NONE"

        # ---------------------------------------------------------------------
        # convert mosaic raster to integer
        intRaster = path.join(outputGDB,PREFIX + mosaicBaseName)
        statement = ('outras = arcpy.sa.Int(arcpy.sa.Raster(mosaicRaster) '
                     '- offset + 0.5); '
                     'outras.save(intRaster)')
        count = 0
        while True:
            try:
                exec(statement)
            except Exception:
                count,tryAgain = lu.retry_arc_error(count,statement)
                if not tryAgain: exec(statement)
            else: break
        # ---------------------------------------------------------------------


        if writeTruncRaster:
            # -----------------------------------------------------------------
            # Set anything beyond cfg.CWDTHRESH to NODATA.
            truncRaster = (outputGDB + '\\' + PREFIX + mosaicBaseName +
                           '_truncated_at_' + lu.cwd_cutoff_str(cfg.CWDTHRESH))

            statement = ('outRas = arcpy.sa.Raster(intRaster)'
                         '* (arcpy.sa.Con(arcpy.sa.Raster(intRaster) '
                         '<= cfg.CWDTHRESH, 1)); '
                         'outRas.save(truncRaster)')

            count = 0
            while True:
                try:
                    exec(statement)
                except Exception:
                    count,tryAgain = lu.retry_arc_error(count,statement)
                    if not tryAgain: exec(statement)
                else: break
        # ---------------------------------------------------------------------
        # Check for unreasonably low minimum NLCC values
        try:
            mosaicGrid = path.join(cfg.LCCBASEDIR,'mos')
            # Copy to grid to test
            arcpy.CopyRaster_management(mosaicRaster, mosaicGrid)
            minObject = arcpy.GetRasterProperties_management(mosaicGrid, "MINIMUM")
            rasterMin = float(str(minObject.getOutput(0)))
        except Exception:
            lu.warn('\n------------------------------------------------')
            lu.warn('WARNING: Raster minimum check failed in step 5. \n'
                'This may mean the output rasters are corrupted. Please \n'
                'be sure to check for valid rasters in '+ outputGDB)
            rasterMin = 0
        tolerance = (float(arcpy.env.cellSize) * -10)
        if rasterMin < tolerance:
            lu.dashline(1)
            msg = ('WARNING: Minimum value of mosaicked corridor map is '
                   'much less than zero ('+str(rasterMin)+').'
                   '\nThis could mean that BOUNDING CIRCLE BUFFER DISTANCES '
                   'were too small and a corridor passed outside of a '
                   'bounding circle, or that a corridor passed outside of the '
                   'resistance map. \n')
            lu.warn(msg)


        gprint('\nWriting final LCP maps...')
        if cfg.STEP4:
            finalLinkTable = lu.update_lcp_shapefile(linkTable, lastStep=4,
                                                     thisStep=5)
        elif cfg.STEP3:
            finalLinkTable = lu.update_lcp_shapefile(linkTable, lastStep=3,
                                                     thisStep=5)
        else:
            # Don't know if step 4 was run, since this is started at step 5.
            # Use presence of previous linktable files to figure this out.
            # Linktable name includes step number.
            prevLinkTableFile = lu.get_prev_step_link_table(step=5)
            prevStepInd = len(prevLinkTableFile) - 5
            lastStep = prevLinkTableFile[prevStepInd]

            finalLinkTable = lu.update_lcp_shapefile(linkTable, lastStep,
                                                     thisStep=5)

        outlinkTableFile = lu.get_this_step_link_table(step=5)
        gprint('Updating ' + outlinkTableFile)
        lu.write_link_table(linkTable, outlinkTableFile)

        linkTableLogFile = path.join(cfg.LOGDIR, "linkTable_s5.csv")
        lu.write_link_table(linkTable, linkTableLogFile)

        linkTableFinalFile = path.join(cfg.OUTPUTDIR, PREFIX +
                                       "_linkTable_s5.csv")
        lu.write_link_table(finalLinkTable, linkTableFinalFile)
        gprint('Copy of final linkTable written to '+
                          linkTableFinalFile)

        gprint('Creating shapefiles with linework for links.')
        try:
            lu.write_link_maps(outlinkTableFile, step=5)
        except Exception:
            lu.write_link_maps(outlinkTableFile, step=5)

        # Create final linkmap files in output directory, and remove files from
        # scratch.
        lu.copy_final_link_maps(step=5)

        if not SAVENORMLCCS:
            lu.delete_dir(cfg.LCCBASEDIR)

        # Build statistics for corridor rasters
        arcpy.AddMessage('\nBuilding output statistics and pyramids '
                          'for corridor raster')
        lu.build_stats(intRaster)

        if writeTruncRaster:
            arcpy.AddMessage('Building output statistics '
                              'for truncated corridor raster')
            lu.build_stats(truncRaster)

        if cfg.OUTPUTFORMODELBUILDER:
            arcpy.CopyFeatures_management(cfg.COREFC, cfg.OUTPUTFORMODELBUILDER)

    # Return GEOPROCESSING specific errors
    except arcpy.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 5. Details follow.****')
        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except Exception:
        lu.dashline(1)
        gprint('****Failed in step 5. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)

    return
