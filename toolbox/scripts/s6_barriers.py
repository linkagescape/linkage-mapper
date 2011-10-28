#!/usr/bin/env python2.5

"""Detects influential barriers given CWD calculations from
    linkage mapper step 3.
"""

__filename__ = "s6_barriers.py"
__version__ = "BARRIER TEST"

import os.path as path
import time
import shutil

import arcgisscripting
import numpy as npy

from lm_config import Config as Cfg
import lm_util as lu

# Set local references to objects and constants from bar_config
gp = Cfg.gp
gprint = gp.addmessage
LTB_CORE1 = Cfg.LTB_CORE1
LTB_CORE2 = Cfg.LTB_CORE2
LTB_LINKID = Cfg.LTB_LINKID
LTB_LINKTYPE = Cfg.LTB_LINKTYPE
LTB_CWDIST = Cfg.LTB_CWDIST

def STEP6_calc_barriers():
    """Detects influential barriers given CWD calculations from
       s3_calcCwds.py.

    """

# Fixme: add option to save individual barrier files?
    try:        
        lu.dashline(0)
        gprint('Running script ' + __filename__)

        startRadius = float(Cfg.STARTRADIUS) 
        endRadius = float(Cfg.ENDRADIUS) 
        radiusStep = float(Cfg.RADIUSSTEP)
        if radiusStep == 0: 
            endRadius = startRadius # Calculate at just one radius value
            radiusStep = 1
        linkTableFile = lu.get_prev_step_link_table(step=6)
        gp.workspace = Cfg.SCRATCHDIR
        prefix = path.basename(Cfg.PROJECTDIR)        

        # set the analysis extent and cell size to that of the resistance
        # surface
        gp.Extent = gp.Describe(Cfg.RESRAST).Extent
        gp.CellSize = gp.Describe(Cfg.RESRAST).MeanCellHeight
        
        gp.Extent = "MINOF"
        #gp.mask = Cfg.RESRAST #BHM- need nodata values as barriers
        gp.snapraster = Cfg.RESRAST

        linkTable = lu.load_link_table(linkTableFile)
        numLinks = linkTable.shape[0]
        numCorridorLinks = lu.report_links(linkTable)
        if numCorridorLinks == 0:
            lu.dashline()
            gprint('\nThere are no linkages. Bailing.')
            time.sleep(5)
            return
            
        # For speed:
        gp.pyramid = "NONE"
        gp.rasterstatistics = "NONE"

                    
        # set up directories for barrier and barrier mosaic grids
        dirCount = 0
        gprint("Creating output folder: " + Cfg.BARRIERBASEDIR)
        if path.exists(Cfg.BARRIERBASEDIR):
            shutil.rmtree(Cfg.BARRIERBASEDIR)
        gp.CreateFolder_management(path.dirname(Cfg.BARRIERBASEDIR),
                                       path.basename(Cfg.BARRIERBASEDIR))
        gp.CreateFolder_management(Cfg.BARRIERBASEDIR, Cfg.BARRIERDIR_NM)
        cbarrierdir = path.join(Cfg.BARRIERBASEDIR, Cfg.BARRIERDIR_NM)
        # Create output geodatabase
        Cfg.gp.createfilegdb(Cfg.OUTPUTDIR, path.basename(Cfg.BARRIERGDB))
        
        coresToProcess = npy.unique(linkTable[:, LTB_CORE1:LTB_CORE2 + 1])
        maxCoreNum = max(coresToProcess)
        del coresToProcess

        # Set up focal directories.
        # To keep there from being > 100 grids in any one directory,
        # outputs are written to:
        # barrier\focalX_ for cores 1-99 at radius X
        # barrier\focalX_1 for cores 100-199
        # etc.
#code here:        
        lu.dashline(0)

        for radius in range(startRadius, endRadius + 1, radiusStep):
            core1path = lu.get_focal_path(1,radius)
            path1, dir1 = path.split(core1path)
            path2, dir2 = path.split(path1)
            gp.CreateFolder_management(path.dirname(path2),
                                       path.basename(path2))
            gp.CreateFolder_management(path.dirname(path1),
                                       path.basename(path1))            
            
            if maxCoreNum > 100:
                maxDirCount = int(maxCoreNum/100)
                focalDirBaseName = dir2
                for dir in range(1, maxDirCount + 1):
                    focalDir = focalDirBaseName + str(dir)
                    gp.CreateFolder_management(path2, focalDir)

        # Create resistance raster with filled-in Nodata values for later use
        isNullExpression = ("con(isnull(" + Cfg.RESRAST + "), " 
            + str(1000000000) + ", " + Cfg.RESRAST + ")")
        resistFillRaster = path.join(Cfg.SCRATCHDIR, "resist_fill")
        gp.SingleOutputMapAlgebra_sa(isNullExpression, resistFillRaster)
       
        numGridsWritten = 0
        coreList = linkTable[:,LTB_CORE1:LTB_CORE2+1]
        coreList = npy.sort(coreList)

        # Loop through each search radius to calculate barriers in each link
        for radius in range (startRadius, endRadius + 1, radiusStep):
            linkLoop = 0
            pctDone = 0
            gprint('\nMapping barriers at a radius of ' + str(radius) + 
                   ' map units...')
            gprint('0 percent done')
            for x in range(0,numLinks):
                pctDone = lu.report_pct_done(linkLoop, numCorridorLinks, 
                                            pctDone)
                linkId = str(int(linkTable[x,LTB_LINKID]))
                if ((linkTable[x,LTB_LINKTYPE] > 0) and 
                   (linkTable[x,LTB_LINKTYPE] < 100)):
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
                    gp.Extent = "MINOF"

                    link = lu.get_links_from_core_pairs(linkTable, 
                                                        corex, corey)
                    lcDist = str(float(linkTable[link,LTB_CWDIST]))

                    # Detect barriers at radius using neighborhood stats
                    # Create the Neighborhood Object
                    innerRadius = radius - 1
                    outerRadius = radius
                    
                    dia = 2 * radius 
                    InNeighborhood = ("ANNULUS " + str(innerRadius) + " " + 
                                     str(outerRadius) + " MAP")                              

                    # Execute FocalStatistics
                    if not path.exists(focalRaster1):
                        gp.FocalStatistics_sa(cwdRaster1, focalRaster1, 
                                              InNeighborhood, "MINIMUM","DATA")
                    if not path.exists(focalRaster2):
                        gp.FocalStatistics_sa(cwdRaster2, focalRaster2, 
                                              InNeighborhood, "MINIMUM","DATA")
                    
                    #Calculate potential benefit per pixel restored
                    deltaExpression = ("(" + lcDist + " - " + focalRaster1 
                        + " - " + focalRaster2 + " - " + str(dia) 
                        + ") / " + str(dia))
                    gp.SingleOutputMapAlgebra_sa(deltaExpression, 
                                                 barrierRaster)
                    
                    gp.workspace = Cfg.SCRATCHDIR
                    tempMosaicRaster = "mos_temp"                                    
                    if linkLoop == 1:
                    #If this is the first grid then copy rather than mosaic
                        gp.CopyRaster_management(barrierRaster, 
                                                 tempMosaicRaster)
                        
                    else:
                        # How to combine across linkages?  For now take max....                                                   
                        tempMosaicRaster = path.join(Cfg.SCRATCHDIR,"mos_temp")
                        gp.Mosaic_management(barrierRaster, tempMosaicRaster,
                                             "MAXIMUM", "MATCH")

                    if not Cfg.SAVEBARRIERRASTERS:
                        lu.delete_data(barrierRaster)

                    # Temporarily disable links in linktable - 
                    # don't want to mosaic them twice
                    for y in range (x+1,numLinks):
                        corex1 = int(coreList[y,0])
                        corey1 = int(coreList[y,1])
                        if corex1 == corex and corey1 == corey:
                            linkTable[y,LTB_LINKTYPE] = (
                                linkTable[y,LTB_LINKTYPE] + 100)
                        elif corex1==corey and corey1==corex:
                            linkTable[y,LTB_LINKTYPE] = (
                                linkTable[y,LTB_LINKTYPE] + 100)
                    
                    if Cfg.SAVEBARRIERRASTERS: 
                        numGridsWritten = numGridsWritten + 1
                    if numGridsWritten == 100:
                        # We only write up to 100 grids to any one folder
                        # because otherwise Arc slows to a crawl
                        dirCount = dirCount + 1
                        numGridsWritten = 0
                        cbarrierdir = (path.join(Cfg.BARRIERBASEDIR, 
                                      Cfg.BARRIERDIR_NM) + str(dirCount))
                        gp.CreateFolder_management(Cfg.BARRIERBASEDIR,
                                                   path.basename(cbarrierdir))
            #rows that were temporarily disabled
            rows = npy.where(linkTable[:,LTB_LINKTYPE]>100)
            linkTable[rows,LTB_LINKTYPE] = (
                linkTable[rows,LTB_LINKTYPE] - 100)
                        
            # -----------------------------------------------------------------
            mosaicFN = "barriers" + str(radius)
            #fixme- write final to geodatabase instead
            mosaicRaster = path.join(Cfg.BARRIERBASEDIR, mosaicFN)
            gp.SetNull_sa(tempMosaicRaster, tempMosaicRaster, mosaicRaster, 
                          "VALUE < 0")
            
            lu.delete_data(tempMosaicRaster)

            # Place copy of result in output geodatabase
            Cfg.gp.workspace = Cfg.BARRIERGDB
            mosaicFN = prefix + "_barriers" + str(radius)
            Cfg.gp.CopyRaster_management(mosaicRaster, mosaicFN)

            # 'Grow out' maximum restoration gain to 
            # neighborhood size for display
            InNeighborhood = "CIRCLE " + str(outerRadius) + " MAP"                               
            # Execute FocalStatistics
            maxRasterFN = "bar_max" + str(outerRadius)
            maxRaster = path.join(Cfg.BARRIERBASEDIR, maxRasterFN)
            gp.FocalStatistics_sa(mosaicRaster, maxRaster, InNeighborhood, 
                                  "MAXIMUM","DATA")
            #Place a copy in output geodatabase
            Cfg.gp.workspace = Cfg.BARRIERGDB
            maxRasterFN = prefix + "_bar_max" + str(outerRadius)
            Cfg.gp.CopyRaster_management(maxRaster, maxRasterFN)
            
            # Create pared-down version of maximum- remove pixels that
            # don't need restoring by allowing a pixel to only contribute its
            # resistance value to restoration gain
            outRasterFN = "bar_trm" + str(outerRadius)
            outRaster = path.join(Cfg.BARRIERBASEDIR,outRasterFN)
#            rasterList = "'" + resistFillRaster + "'; '" + maxRaster + "'"
            rasterList = maxRaster  + ";" + resistFillRaster 
            #gp.mask = maxRaster # 
            gp.CellStatistics_sa(rasterList, outRaster, "MINIMUM")
            #Place a copy in output geodatabase
            Cfg.gp.workspace = Cfg.BARRIERGDB
            outRasterFN = prefix + "_bar_trm" + str(outerRadius)
            Cfg.gp.CopyRaster_management(outRaster, outRasterFN)


            
        # Combine rasters across radii
        gprint('\nCreating summary rasters...')
        if startRadius != endRadius: 
            mosaicFN = "bar_radii"               
            gp.workspace = Cfg.BARRIERBASEDIR
            for radius in range (startRadius, endRadius + 1, radiusStep):
                radiusFN = "barriers" + str(radius)
                radiusRaster = path.join(Cfg.BARRIERBASEDIR, radiusFN)            
                if radius == startRadius:
                #If this is the first grid then copy rather than mosaic
                    gp.CopyRaster_management(radiusRaster, mosaicFN)
                else:
                    mosaicRaster = path.join(Cfg.BARRIERBASEDIR,mosaicFN)
                    gp.Mosaic_management(radiusRaster, mosaicRaster, 
                                         "MAXIMUM", "MATCH")
            # Copy result to output geodatabase
            Cfg.gp.workspace = Cfg.BARRIERGDB
            mosaicFN = prefix + "_bar_radii"
            Cfg.gp.CopyRaster_management(mosaicRaster, mosaicFN)
        
            #GROWN OUT rasters
            maxMosaicFN = "bar_radii_max"                                    
            gp.workspace = Cfg.BARRIERBASEDIR
            for radius in range (startRadius, endRadius + 1, radiusStep):
                radiusFN = "bar_max" + str(radius)
                #fixme- do this when only a single radius too
                radiusRaster = path.join(Cfg.BARRIERBASEDIR, radiusFN)            
                if radius == startRadius:
                #If this is the first grid then copy rather than mosaic
                    gp.CopyRaster_management(radiusRaster, maxMosaicFN)
                else:
                    maxMosaicRaster = path.join(Cfg.BARRIERBASEDIR,maxMosaicFN)
                    gp.Mosaic_management(radiusRaster, maxMosaicRaster, 
                                         "MAXIMUM", "MATCH")
            # Copy result to output geodatabase
            Cfg.gp.workspace = Cfg.BARRIERGDB
            maxMosaicFN = prefix + "_bar_radii_max"                                    
            Cfg.gp.CopyRaster_management(maxMosaicRaster, maxMosaicFN)
            
            #GROWN OUT AND TRIMMED rasters
            trimMosaicFN = "bar_radii_trm" 
            gp.workspace = Cfg.BARRIERBASEDIR            
            for radius in range (startRadius, endRadius + 1, radiusStep):
                radiusFN = "bar_trm" + str(radius)
                #fixme- do this when only a single radius too
                radiusRaster = path.join(Cfg.BARRIERBASEDIR, radiusFN)            
                if radius == startRadius:
                #If this is the first grid then copy rather than mosaic
                    gp.CopyRaster_management(radiusRaster, trimMosaicFN)
                else:
                    trimMosaicRaster = path.join(Cfg.BARRIERBASEDIR,
                                                 trimMosaicFN)
                    gp.Mosaic_management(radiusRaster, trimMosaicRaster, 
                                         "MAXIMUM", "MATCH")
            # Copy result to output geodatabase
            Cfg.gp.workspace = Cfg.BARRIERGDB
            trimMosaicFN = prefix + "_bar_radii_trm" 
            Cfg.gp.CopyRaster_management(trimMosaicRaster, trimMosaicFN)

                    
        #Clean up temporary files and directories
        try:
            shutil.rmtree(Cfg.SCRATCHDIR)
        except:
            pass
        
        if not Cfg.SAVEBARRIERRASTERS:
            cbarrierdir = path.join(Cfg.BARRIERBASEDIR, Cfg.BARRIERDIR_NM) 
            try:
                shutil.rmtree(cbarrierdir)
            except:
                pass
            for dir in range(1,dirCount+1):
                cbarrierdir = path.join((Cfg.BARRIERBASEDIR, Cfg.BARRIERDIR_NM) 
                                        + str(dir))
                try:
                    shutil.rmtree(cbarrierdir)
                except:
                    pass
                    
        if not Cfg.SAVEFOCALRASTERS:
            for radius in range(startRadius, endRadius + 1, radiusStep):
                core1path = lu.get_focal_path(1,radius)
                path1, dir1 = path.split(core1path)
                path2, dir2 = path.split(path1)
                try:
                    shutil.rmtree(path2)                
                except:
                    pass
                    
    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 6. Details follow.****')
        lu.raise_geoproc_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 6. Details follow.****')
        lu.raise_python_error(__filename__)

    return
    
    
    