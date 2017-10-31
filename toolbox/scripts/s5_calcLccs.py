#!/usr/bin/env python2.5
# Authors: Brad McRae and Darren Kavanagh

"""Step 5: Calculate least cost corridors.

Creates and mosaics normalized least-cost corridors using connected core area
pairs specified in linkTable and cwd layers

"""

import os.path as path
import time

import numpy as npy

from lm_config import tool_env as cfg
import lm_util as lu

_SCRIPT_NAME = "s5_calcLccs.py"

# try:
    # import arcpy
    # gp = arcpy.gp
    # from arcpy.sa import *
    # arcgisscripting = arcpy
    # arcpy.CheckOutExtension("spatial")
    # arcObj = arcpy
# except:
try:
    import arcpy # want to know if arcpy is available
    from arcpy.sa import *
    arcpyAvailable = True
except:
    arcpyAvailable = False

# # Arcpy not working in some cases unless s3 executed in same run.
# # For nso3 dataset, map algebra adding rasters results in empty rasters. 
# # for now there seems to be no speed improvement with arcpy.

cfg.useArcpy = False
gp = cfg.gp
import arcgisscripting
arcObj = cfg.gp
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
    except:
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
        if cfg.useArcpy:
            arcpy.env.workspace = cfg.SCRATCHDIR
            arcpy.env.scratchWorkspace = cfg.ARCSCRATCHDIR
            arcpy.env.overwriteOutput = True  
            arcpy.env.compression = "NONE"
        else:        
            gp.workspace = cfg.SCRATCHDIR
            gp.scratchWorkspace = cfg.ARCSCRATCHDIR
            gp.OverwriteOutput = True            

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
        if cfg.useArcpy:
            arcpy.env.Extent = cfg.RESRAST
            arcpy.env.cellSize = cfg.RESRAST
            arcpy.env.snapRaster = cfg.RESRAST
            arcpy.env.mask = cfg.RESRAST
        else:
            gp.Extent = (gp.Describe(cfg.RESRAST)).Extent 
            gp.cellSize = gp.Describe(cfg.RESRAST).MeanCellHeight
            gp.mask = cfg.RESRAST
            gp.snapraster = cfg.RESRAST

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
        gp.pyramid = "NONE"
        gp.rasterstatistics = "NONE"

        # set up directories for normalized lcc and mosaic grids
        dirCount = 0
        gprint("Creating output folder: " + cfg.LCCBASEDIR)
        lu.delete_dir(cfg.LCCBASEDIR)
        gp.CreateFolder_management(path.dirname(cfg.LCCBASEDIR),
                                       path.basename(cfg.LCCBASEDIR))
        gp.CreateFolder_management(cfg.LCCBASEDIR, cfg.LCCNLCDIR_NM)
        clccdir = path.join(cfg.LCCBASEDIR, cfg.LCCNLCDIR_NM)
        # mosaicGDB = path.join(cfg.LCCBASEDIR, "mosaic.gdb")
        # gp.createfilegdb(cfg.LCCBASEDIR, "mosaic.gdb")
        #mosaicRaster = mosaicGDB + '\\' + "nlcc_mos" # Full path
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

            if not gp.Exists(cwdRaster1):
                msg =('\nError: cannot find cwd raster:\n' + cwdRaster1) 
            if not gp.Exists(cwdRaster2):
                msg =('\nError: cannot find cwd raster:\n' + cwdRaster2) 
                lu.raise_error(msg)

            
            lccNormRaster = path.join(clccdir, str(corex) + "_" +
                                      str(corey))# + ".tif")
            if cfg.useArcpy: 
                arcpy.env.Extent = "MINOF"
            else:
                gp.Extent = "MINOF"

            # FIXME: need to check for this?:
            # if exists already, don't re-create
            #if not gp.Exists(lccRaster):

            link = lu.get_links_from_core_pairs(linkTable, corex, corey)

            offset = 10000 

            # Normalized lcc rasters are created by adding cwd rasters and
            # subtracting the least cost distance between them.
            count = 0
            if arcpyAvailable:
                cfg.useArcpy = True # Fixes Canran Liu's bug with lcDist
            if cfg.useArcpy:
                
                lcDist = (float(linkTable[link,cfg.LTB_CWDIST]) - offset) 
                
                if normalize:
                    statement = ('outras = Raster(cwdRaster1) + Raster('
                        'cwdRaster2) - lcDist; outras.save(lccNormRaster)') 
                                                
                else:
                    statement = ('outras =Raster(cwdRaster1) + Raster('
                                'cwdRaster2); outras.save(lccNormRaster)')
            else:
                if normalize:
                    lcDist = str(linkTable[link,cfg.LTB_CWDIST] - offset)  
                    expression = (cwdRaster1 + " + " + cwdRaster2 + " - " 
                                  + lcDist)
                else:
                    expression = (cwdRaster1 + " + " + cwdRaster2) 
                statement = ('gp.SingleOutputMapAlgebra_sa(expression, '
                     'lccNormRaster)')
            count = 0
            while True:
                try: 
                    exec statement                    
                    randomerror()
                except:
                    count,tryAgain = lu.retry_arc_error(count,statement)
                    if not tryAgain:    
                        exec statement
                else: break
            cfg.useArcpy = False # End fix for Conran Liu's bug with lcDist
            
            if normalize and cfg.useArcpy: 
                try: 
                    minObject = gp.GetRasterProperties(lccNormRaster, "MINIMUM") 
                    rasterMin = float(str(minObject.getoutput(0)))
                except:
                    lu.warn('\n------------------------------------------------')
                    lu.warn('WARNING: Raster minimum check failed in step 5. \n'
                        'This may mean the output rasters are corrupted. Please \n'
                        'be sure to check for valid rasters in '+ outputGDB)
                    rasterMin = 0
                if rasterMin < tolerance:
                    lu.dashline(1)
                    msg = ('WARNING: Minimum value of a corridor #' + str(x+1) 
                           + ' is much less than zero ('+str(rasterMin)+').'
                           '\nThis could mean that BOUNDING CIRCLE BUFFER DISTANCES '
                           'were too small and a corridor passed outside of a '
                           'bounding circle, or that a corridor passed outside of the '
                           'resistance map. \n')
                    lu.warn(msg)

            
            if cfg.useArcpy: 
                arcpy.env.Extent = cfg.RESRAST
            else:
                gp.Extent = (gp.Describe(cfg.RESRAST)).Extent

            mosaicDir = path.join(cfg.LCCBASEDIR,'mos'+str(x+1))  
            lu.create_dir(mosaicDir) 
            mosFN = 'mos'#.tif' change and move
            mosaicRaster = path.join(mosaicDir,mosFN) 
                       
            if numGridsWritten == 0 and dirCount == 0:
                #If this is the first grid then copy rather than mosaic
                arcObj.CopyRaster_management(lccNormRaster, mosaicRaster) 
            else:
                
                rasterString = '"'+lccNormRaster+";"+lastMosaicRaster+'"'
                statement = ('arcObj.MosaicToNewRaster_management('
                            'rasterString,mosaicDir,mosFN, "", '
                            '"32_BIT_FLOAT", gp.cellSize, "1", "MINIMUM", '
                            '"MATCH")') 
                # statement = ('arcpy.Mosaic_management(lccNormRaster, '
                                 # 'mosaicRaster, "MINIMUM", "MATCH")') 
                
                count = 0
                while True:
                    try:
                        lu.write_log('Executing mosaic for link #'+str(linkId))
                        exec statement
                        lu.write_log('Done with mosaic.')
                        randomerror()
                    except:
                        count,tryAgain = lu.retry_arc_error(count,statement)
                        lu.delete_data(mosaicRaster)
                        lu.delete_dir(mosaicDir)
                        # Try a new directory
                        mosaicDir = path.join(cfg.LCCBASEDIR,'mos'+str(x+1)+ '_' + str(count))
                        lu.create_dir(mosaicDir)
                        mosaicRaster = path.join(mosaicDir,mosFN)                        
                        if not tryAgain:    
                            exec statement
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
                    gp.CreateFolder_management(cfg.LCCBASEDIR,
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
        if not gp.exists(outputGDB):
            gp.createfilegdb(cfg.OUTPUTDIR, path.basename(outputGDB))

        if cfg.useArcpy:
            arcpy.env.workspace = outputGDB
        else:        
            gp.workspace = outputGDB

        gp.pyramid = "NONE"
        gp.rasterstatistics = "NONE"

        
        # Copy mosaic raster to output geodatabase
        saveFloatRaster = False
        if saveFloatRaster == True:
            floatRaster = outputGDB + '\\' + PREFIX + mosaicBaseName + '_flt' # Full path 
            statement = 'arcObj.CopyRaster_management(mosaicRaster, floatRaster)'
            try:
                exec statement
            except:
                pass
                

        # ---------------------------------------------------------------------
        # convert mosaic raster to integer
        intRaster = path.join(outputGDB,PREFIX + mosaicBaseName)
        if cfg.useArcpy:
            statement = ('outras = Int(Raster(mosaicRaster) - offset + 0.5); ' 
                        'outras.save(intRaster)')
        else:
            expression = "int(" + mosaicRaster + " - " + str(offset) + " + 0.5)"
            statement = 'gp.SingleOutputMapAlgebra_sa(expression, intRaster)'
        count = 0
        while True:
            try: 
                exec statement
                randomerror()
            except:
                count,tryAgain = lu.retry_arc_error(count,statement)
                if not tryAgain: exec statement
            else: break
        # ---------------------------------------------------------------------       
        

        if writeTruncRaster:
            # -----------------------------------------------------------------
            # Set anything beyond cfg.CWDTHRESH to NODATA.
            if arcpyAvailable:
                cfg.useArcpy = True # For Alissa Pump's error with 10.1
            cutoffText = str(cfg.CWDTHRESH)
            if cutoffText[-6:] == '000000':
                cutoffText = cutoffText[0:-6]+'m' 
            elif cutoffText[-3:] == '000':
                cutoffText = cutoffText[0:-3]+'k' 
            
            truncRaster = (outputGDB + '\\' + PREFIX + mosaicBaseName + 
                           '_truncated_at_' + cutoffText)

            count = 0
            if cfg.useArcpy:
                statement = ('outRas = Raster(intRaster) * '
                            '(Con(Raster(intRaster) <= cfg.CWDTHRESH,1)); '
                            'outRas.save(truncRaster)')
            else:
                expression = ("(" + intRaster + " * (con(" + intRaster + "<= " 
                              + str(cfg.CWDTHRESH) + ",1)))")
                statement = ('gp.SingleOutputMapAlgebra_sa(expression, '
                                                          'truncRaster)')
            count = 0
            while True:
                try: 
                    exec statement
                    randomerror()
                except:
                    count,tryAgain = lu.retry_arc_error(count,statement)
                    if not tryAgain: exec statement
                else: break
            cfg.useArcpy = False # End fix for Alissa Pump's error with 10.1                
        # ---------------------------------------------------------------------
        # Check for unreasonably low minimum NLCC values    
        try:
            mosaicGrid = path.join(cfg.LCCBASEDIR,'mos') 
            # Copy to grid to test
            arcObj.CopyRaster_management(mosaicRaster, mosaicGrid)
            minObject = gp.GetRasterProperties(mosaicGrid, "MINIMUM") 
            rasterMin = float(str(minObject.getoutput(0)))
        except:
            lu.warn('\n------------------------------------------------')
            lu.warn('WARNING: Raster minimum check failed in step 5. \n'
                'This may mean the output rasters are corrupted. Please \n'
                'be sure to check for valid rasters in '+ outputGDB)
            rasterMin = 0
        tolerance = (float(gp.cellSize) * -10)
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
        except:
            lu.write_link_maps(outlinkTableFile, step=5)

        # Create final linkmap files in output directory, and remove files from
        # scratch.
        lu.copy_final_link_maps(step=5)

        if not SAVENORMLCCS:
            lu.delete_dir(cfg.LCCBASEDIR)

        # Build statistics for corridor rasters
        gp.addmessage('\nBuilding output statistics and pyramids '
                          'for corridor raster')
        lu.build_stats(intRaster)

        if writeTruncRaster:
            gp.addmessage('Building output statistics '
                              'for truncated corridor raster')
            lu.build_stats(truncRaster)
        
        if cfg.OUTPUTFORMODELBUILDER:
            arcpy.CopyFeatures_management(cfg.COREFC, cfg.OUTPUTFORMODELBUILDER)

    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 5. Details follow.****')
        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 5. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)

    return
       
    
def randomerror():
    """ Used to test error recovery.

    """
    generateError = False # Set to True to create random errors
    if generateError:
        gprint('\n***Rolling dice for random error***')
        import random
        test = random.randrange(1, 8)
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
