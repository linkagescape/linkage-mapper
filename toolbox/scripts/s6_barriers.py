#add options to do trm, max, etc
#remove datapass/barrier

#!/usr/bin/env python2.5
# Author: Brad McRae

"""Detects influential barriers given CWD calculations from
    linkage mapper step 3.
"""

import os.path as path

import arcpy
from arcpy.sa import *
import numpy as npy

from lm_config import tool_env as cfg
import lm_util as lu


_SCRIPT_NAME = "s6_barriers.py"


def STEP6_calc_barriers():
    """Detects influential barriers given CWD calculations from
       s3_calcCwds.py.

    """

# Fixme: add option to save individual barrier files?
    gprint = lu.gprint
    try:
        arcpy.CheckOutExtension("spatial")

        lu.dashline(0)
        gprint('Running script ' + _SCRIPT_NAME)

        startRadius = float(cfg.STARTRADIUS)
        endRadius = float(cfg.ENDRADIUS)
        radiusStep = float(cfg.RADIUSSTEP)
        if radiusStep == 0:
            endRadius = startRadius # Calculate at just one radius value
            radiusStep = 1
        linkTableFile = lu.get_prev_step_link_table(step=6)
        arcpy.env.workspace = cfg.SCRATCHDIR
        arcpy.env.scratchWorkspace = cfg.ARCSCRATCHDIR
        arcpy.RefreshCatalog(cfg.PROJECTDIR)
        PREFIX = path.basename(cfg.PROJECTDIR)
        # For speed:
        arcpy.env.pyramid = "NONE"
        arcpy.env.rasterStatistics = "NONE"

        # set the analysis extent and cell size to that of the resistance
        # surface
        arcpy.OverWriteOutput = True
        arcpy.env.extent = cfg.RESRAST
        arcpy.env.cellSize = cfg.RESRAST
        arcpy.env.snapRaster = cfg.RESRAST
        if float(arcpy.env.cellSize) > startRadius or startRadius > endRadius:
            msg = ('Error: minimum detection radius must be greater than '
                    'cell size (' + str(arcpy.env.cellSize) +
                    ') \nand less than or equal to maximum detection radius.')
            lu.raise_error(msg)

        linkTable = lu.load_link_table(linkTableFile)
        numLinks = linkTable.shape[0]
        numCorridorLinks = lu.report_links(linkTable)
        if numCorridorLinks == 0:
            lu.dashline(1)
            msg =('\nThere are no linkages. Bailing.')
            lu.raise_error(msg)

        # set up directories for barrier and barrier mosaic grids
        dirCount = 0
        gprint("Creating intermediate output folder: " + cfg.BARRIERBASEDIR)
        lu.delete_dir(cfg.BARRIERBASEDIR)
        arcpy.CreateFolder_management(path.dirname(cfg.BARRIERBASEDIR),
                                       path.basename(cfg.BARRIERBASEDIR))
        arcpy.CreateFolder_management(cfg.BARRIERBASEDIR, cfg.BARRIERDIR_NM)
        cbarrierdir = path.join(cfg.BARRIERBASEDIR, cfg.BARRIERDIR_NM)

        coresToProcess = npy.unique(linkTable[:, cfg.LTB_CORE1:cfg.LTB_CORE2 + 1])
        maxCoreNum = max(coresToProcess)


        # Set up focal directories.
        # To keep there from being > 100 grids in any one directory,
        # outputs are written to:
        # barrier\focalX_ for cores 1-99 at radius X
        # barrier\focalX_1 for cores 100-199
        # etc.
        lu.dashline(0)

        for radius in range(startRadius, endRadius + 1, radiusStep):
            core1path = lu.get_focal_path(1,radius)
            path1, dir1 = path.split(core1path)
            path2, dir2 = path.split(path1)
            arcpy.CreateFolder_management(path.dirname(path2),
                                       path.basename(path2))
            arcpy.CreateFolder_management(path.dirname(path1),
                                       path.basename(path1))

            # if maxCoreNum > 100:
                # maxDirCount = int(maxCoreNum/100)
                # focalDirBaseName = dir2
                # for dir in range(1, maxDirCount + 1):
                    # focalDir = focalDirBaseName + str(dir)
                    # arcpy.CreateFolder_management(path2, focalDir)


            if maxCoreNum > 99:
                gprint('Creating subdirectories')
                maxDirCount = int(maxCoreNum/100)
                focalDirBaseName = dir2
                for dirCount in range(1, maxDirCount + 1):
                    floor = dirCount * 100
                    ceiling = 99 + dirCount * 100
                    ind = npy.where(coresToProcess < ceiling)
                    lowCores = coresToProcess[ind]
                    if lowCores.shape[0] > 0:
                        ind = npy.where(lowCores > floor)
                        medCores = lowCores[ind]
                        if medCores.shape[0] > 0:
                            focalDir = focalDirBaseName + str(dirCount)
                            gprint('...' + focalDir)
                            arcpy.CreateFolder_management(path2, focalDir)



        del coresToProcess
        # Create resistance raster with filled-in Nodata values for later use
        arcpy.env.extent = cfg.RESRAST
        resistFillRaster = path.join(cfg.SCRATCHDIR, "resist_fill")
        output = arcpy.sa.Con(IsNull(cfg.RESRAST), 1000000000, cfg.RESRAST)
        output.save(resistFillRaster)

        numGridsWritten = 0
        coreList = linkTable[:,cfg.LTB_CORE1:cfg.LTB_CORE2+1]
        coreList = npy.sort(coreList)

        # Loop through each search radius to calculate barriers in each link
        for radius in range (startRadius, endRadius + 1, radiusStep):
            linkLoop = 0
            pctDone = 0
            gprint('\nMapping barriers at a radius of ' + str(radius) +
                   ' map units...')
            if numCorridorLinks > 1:
                gprint('0 percent done')
            for x in range(0,numLinks):
                pctDone = lu.report_pct_done(linkLoop, numCorridorLinks,
                                            pctDone)
                linkId = str(int(linkTable[x,cfg.LTB_LINKID]))
                if ((linkTable[x,LTB_LINKTYPE] > 0) and
                   (linkTable[x,LTB_LINKTYPE] < 1000)):
                    linkLoop = linkLoop + 1
                    # source and target cores
                    corex=int(coreList[x,0])
                    corey=int(coreList[x,1])

                    # Get cwd rasters for source and target cores
                    cwdRaster1 = lu.get_cwd_path(corex)
                    cwdRaster2 = lu.get_cwd_path(corey)
                    focalRaster1 = lu.get_focal_path(corex,radius)
                    focalRaster2 = lu.get_focal_path(corey,radius)
                    barrierRaster = path.join(cbarrierdir, "b" + str(radius)
                                              + "_" + str(corex) + "_" +
                                              str(corey))
                    arcpy.env.extent = cfg.RESRAST
                    arcpy.env.extent = "MINOF"

                    link = lu.get_links_from_core_pairs(linkTable,
                                                        corex, corey)
                    lcDist = float(linkTable[link,cfg.LTB_CWDIST])

                    # Detect barriers at radius using neighborhood stats
                    # Create the Neighborhood Object
                    innerRadius = radius - 1
                    outerRadius = radius

                    dia = 2 * radius
                    InNeighborhood = ("ANNULUS " + str(innerRadius) + " " +
                                     str(outerRadius) + " MAP")

                    # Execute FocalStatistics
                    if not path.exists(focalRaster1):
                        outFocalStats = arcpy.sa.FocalStatistics(cwdRaster1,
                                            InNeighborhood, "MINIMUM","DATA")
                        outFocalStats.save(focalRaster1)

                    if not path.exists(focalRaster2):
                        outFocalStats = arcpy.sa.FocalStatistics(cwdRaster2,
                                        InNeighborhood, "MINIMUM","DATA")
                        outFocalStats.save(focalRaster2)
                    #Calculate potential benefit per pixel restored
                    statement = ('outras = ((lcDist - Raster(focalRaster1) - '
                                 'Raster(focalRaster2) - dia) / dia); '
                                 'outras.save(barrierRaster)')

                    count = 0
                    while True:
                        try:
                            exec statement
                        except:
                            count,tryAgain = lu.hiccup_test(count,statement)
                            if not tryAgain:
                                exec statement
                        else: break

                    arcpy.env.workspace = cfg.SCRATCHDIR
                    arcpy.env.workspace = cfg.SCRATCHDIR
                    tempMosaicRaster = "mos_temp"
                    if linkLoop == 1:
                    #If this is the first grid then copy rather than mosaic
                        arcpy.CopyRaster_management(barrierRaster,
                                                 tempMosaicRaster)

                    else:
                        # How to combine across linkages?  For now take max....
                        tempMosaicRaster = path.join(cfg.SCRATCHDIR,"mos_temp")
                        arcpy.Mosaic_management(barrierRaster,
                                         tempMosaicRaster, "MAXIMUM", "MATCH")

                    if not cfg.SAVEBARRIERRASTERS:
                        lu.delete_data(barrierRaster)

                    # Temporarily disable links in linktable -
                    # don't want to mosaic them twice
                    for y in range (x+1,numLinks):
                        corex1 = int(coreList[y,0])
                        corey1 = int(coreList[y,1])
                        if corex1 == corex and corey1 == corey:
                            linkTable[y,LTB_LINKTYPE] = (
                                linkTable[y,LTB_LINKTYPE] + 1000)
                        elif corex1==corey and corey1==corex:
                            linkTable[y,LTB_LINKTYPE] = (
                                linkTable[y,LTB_LINKTYPE] + 1000)

                    if cfg.SAVEBARRIERRASTERS:
                        numGridsWritten = numGridsWritten + 1
                    if numGridsWritten == 100:
                        # We only write up to 100 grids to any one folder
                        # because otherwise Arc slows to a crawl
                        dirCount = dirCount + 1
                        numGridsWritten = 0
                        cbarrierdir = (path.join(cfg.BARRIERBASEDIR,
                                      cfg.BARRIERDIR_NM) + str(dirCount))
                        arcpy.CreateFolder_management(cfg.BARRIERBASEDIR,
                                                   path.basename(cbarrierdir))
            #rows that were temporarily disabled
            rows = npy.where(linkTable[:,LTB_LINKTYPE]>1000)
            linkTable[rows,LTB_LINKTYPE] = (
                linkTable[rows,LTB_LINKTYPE] - 1000)

            # -----------------------------------------------------------------

            mosaicFN = "barriers" + str(radius)
            #fixme- write final to geodatabase instead
            mosaicFN = PREFIX + "_barriers" + str(radius)
            arcpy.env.extent = cfg.RESRAST
            outSetNull = arcpy.sa.SetNull(tempMosaicRaster, tempMosaicRaster,
                                          "VALUE < 0")
            mosaicRaster = path.join(cfg.BARRIERGDB, mosaicFN)
            outSetNull.save(mosaicRaster)

            lu.delete_data(tempMosaicRaster)

            # 'Grow out' maximum restoration gain to
            # neighborhood size for display
            InNeighborhood = "CIRCLE " + str(outerRadius) + " MAP"
            # Execute FocalStatistics
            maxRasterFN = "bar_max" + str(outerRadius)
            maxRaster = path.join(cfg.BARRIERBASEDIR, maxRasterFN)
            outFocalStats = arcpy.sa.FocalStatistics(mosaicRaster,
                                            InNeighborhood, "MAXIMUM","DATA")
            outFocalStats.save(maxRaster)


            #Place a copy in output geodatabase
            arcpy.env.workspace = cfg.BARRIERGDB
            maxRasterFN = PREFIX + "_bar_max" + str(outerRadius)
            arcpy.CopyRaster_management(maxRaster, maxRasterFN)

            # Create pared-down version of maximum- remove pixels that
            # don't need restoring by allowing a pixel to only contribute its
            # resistance value to restoration gain
            outRasterFN = "bar_trm" + str(outerRadius)
            outRaster = path.join(cfg.BARRIERBASEDIR,outRasterFN)
            rasterList = [maxRaster, resistFillRaster]
            outCellStatistics = arcpy.sa.CellStatistics(rasterList, "MINIMUM")
            outCellStatistics.save(outRaster)

            #SECOND ROUND TO CLIP BY DATA VALUES IN BARRIER RASTER
            outRaster2 = outRaster + "2"

            output = arcpy.sa.Con(IsNull(maxRaster),maxRaster,outRaster)
            output.save(outRaster2)

            outRasterFN = PREFIX + "_bar_trm" + str(outerRadius)
            outRasterPath= path.join(cfg.BARRIERGDB, outRasterFN)
            arcpy.CopyRaster_management(outRaster2, outRasterFN)


        # Combine rasters across radii
        gprint('\nCreating summary rasters...')
        if startRadius != endRadius:
            mosaicFN = "bar_radii"
            arcpy.env.workspace = cfg.BARRIERBASEDIR
            for radius in range (startRadius, endRadius + 1, radiusStep):
                #radiusFN = "barriers" + str(radius)
                #radiusRaster = path.join(cfg.BARRIERBASEDIR, radiusFN)
                #Fixme: run speed test with gdb mosaicking above and here
                radiusFN = PREFIX + "_barriers" + str(radius)
                radiusRaster = path.join(cfg.BARRIERGDB, radiusFN)

                if radius == startRadius:
                #If this is the first grid then copy rather than mosaic
                    arcpy.CopyRaster_management(radiusRaster, mosaicFN)
                else:
                    mosaicRaster = path.join(cfg.BARRIERBASEDIR,mosaicFN)

                    arcpy.Mosaic_management(radiusRaster, mosaicRaster,
                                         "MAXIMUM", "MATCH")
            # Copy result to output geodatabase
            arcpy.env.workspace = cfg.BARRIERGDB
            mosaicFN = PREFIX + "_bar_radii"
            arcpy.CopyRaster_management(mosaicRaster, mosaicFN)

            #GROWN OUT rasters
            maxMosaicFN = "bar_radii_max"
            maxMosaicRaster = path.join(cfg.BARRIERBASEDIR,maxMosaicFN)
            arcpy.env.workspace = cfg.BARRIERBASEDIR
            for radius in range (startRadius, endRadius + 1, radiusStep):
                gprint(str(radius))
                radiusFN = "bar_max" + str(radius)
                #fixme- do this when only a single radius too
                radiusRaster = path.join(cfg.BARRIERBASEDIR, radiusFN)
                if radius == startRadius:
                #If this is the first grid then copy rather than mosaic
                    arcpy.CopyRaster_management(radiusRaster, maxMosaicFN)
                else:
                    arcpy.Mosaic_management(radiusRaster, maxMosaicRaster,
                                         "MAXIMUM", "MATCH")
            # Copy result to output geodatabase
            arcpy.env.workspace = cfg.BARRIERGDB
            maxMosaicFN = PREFIX + "_bar_radii_max"
            arcpy.CopyRaster_management(maxMosaicRaster, maxMosaicFN)

            #GROWN OUT AND TRIMMED rasters
            trimMosaicFN = "bar_radii_trm"
            arcpy.env.workspace = cfg.BARRIERBASEDIR
            trimMosaicRaster = path.join(cfg.BARRIERBASEDIR,trimMosaicFN)
            for radius in range (startRadius, endRadius + 1, radiusStep):
                radiusFN = PREFIX + "_bar_trm" + str(radius)
                #fixme- do this when only a single radius too
                radiusRaster = path.join(cfg.BARRIERGDB, radiusFN)
                gprint(radiusRaster)

                if radius == startRadius:
                #If this is the first grid then copy rather than mosaic
                    arcpy.CopyRaster_management(radiusRaster, trimMosaicFN)
                else:
                    arcpy.Mosaic_management(radiusRaster, trimMosaicRaster,
                                         "MAXIMUM", "MATCH")
            # Copy result to output geodatabase
            arcpy.env.workspace = cfg.BARRIERGDB
            trimMosaicFN = PREFIX + "_bar_radii_trm"
            arcpy.CopyRaster_management(trimMosaicRaster, trimMosaicFN)

        arcpy.env.workspace = cfg.BARRIERGDB
        rasters = arcpy.ListRasters()
        for raster in rasters:
            gprint('\nBuilding output statistics and pyramids\n '
                        'for raster ' + raster)
            lu.build_stats(raster)

        #Clean up temporary files and directories
        if not cfg.SAVEBARRIERRASTERS:
            cbarrierdir = path.join(cfg.BARRIERBASEDIR, cfg.BARRIERDIR_NM)
            lu.delete_dir(cbarrierdir)
            for dir in range(1,dirCount+1):
                cbarrierdir = path.join((cfg.BARRIERBASEDIR, cfg.BARRIERDIR_NM)
                                        + str(dir))
                lu.delete_dir(cbarrierdir)

        if not cfg.SAVEFOCALRASTERS:
            for radius in range(startRadius, endRadius + 1, radiusStep):
                core1path = lu.get_focal_path(1,radius)
                path1, dir1 = path.split(core1path)
                path2, dir2 = path.split(path1)
                lu.delete_dir(path2)


    # Return GEOPROCESSING specific errors
    except arcpy.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 6. Details follow.****')
        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 6. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)

    return
