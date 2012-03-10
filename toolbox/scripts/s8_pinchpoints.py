#!/usr/bin/env python2.5
# Author: Brad McRae

"""Maps pinch points using Circuitscape given CWD calculations from
       s3_calcCwds.py.
Reguired Software:
ArcGIS 10 with Spatial Analyst extension
Python 2.6
Numpy
"""

import os.path as path
import time
import subprocess
import gc

import numpy as npy
import arcpy
from arcpy.sa import *

from lm_config import tool_env as cfg
import lm_util as lu

_SCRIPT_NAME = "s8_pinchpoints.py"

arcpy.CheckOutExtension("spatial")

SETCORESTONULL = True
gprint = lu.gprint


def STEP8_calc_pinchpoints():
    """ Experimental code map pinch points using Circuitscape
        given CWD calculations from s3_calcCwds.py.

    """
    try:
        gc.collect()
        lu.dashline(0)
        gprint('Running script ' + _SCRIPT_NAME)

        CSPATH = lu.get_cs_path()
        if CSPATH == None:
            msg = ('Cannot find an installation of Circuitscape 3.5.5'
                    '\nor greater in your Program Files directory.')
            arcpy.AddError(msg)
            lu.write_log(msg)
            exit(1)

        arcpy.OverWriteOutput = True
        arcpy.env.workspace = cfg.SCRATCHDIR
        arcpy.env.scratchWorkspace = cfg.ARCSCRATCHDIR
        # For speed:
        arcpy.env.pyramid = "NONE"
        arcpy.env.rasterstatistics = "NONE"

        # set the analysis extent and cell size to that of the resistance
        # surface
        gprint (cfg.RESRAST)
        arcpy.env.extent = cfg.RESRAST
        arcpy.env.cellSize = cfg.RESRAST

        arcpy.env.extent = "MINOF"

        arcpy.snapraster = cfg.RESRAST
        resRaster = cfg.RESRAST

        if cfg.DO_ADJACENTPAIRS:
            prevLcpShapefile = lu.get_lcp_shapefile(None, thisStep = 8)
            if not arcpy.Exists(prevLcpShapefile):
                msg = ('Cannot find an LCP shapefile from step 5.  Please '
                        'rerun that step and any previous ones if necessary.')
                lu.raise_error(msg)

            # Remove lcp shapefile
            lcpShapefile = path.join(cfg.DATAPASSDIR, "lcpLines_s8.shp")
            lu.delete_data(lcpShapefile)

        inLinkTableFile = lu.get_prev_step_link_table(step=8)
        linkTable = lu.load_link_table(inLinkTableFile)
        numLinks = linkTable.shape[0]
        numCorridorLinks = lu.report_links(linkTable)
        if numCorridorLinks == 0:
            lu.dashline(1)
            msg =('\nThere are no linkages. Bailing.')
            lu.raise_error(msg)

        if linkTable.shape[1] < 16: # If linktable has no entries from prior
                                    # centrality or pinchpint analyses
            extraCols = npy.zeros((numLinks, 6), dtype="float64")
            linkTable = linkTable[:,0:10]
            linkTable = npy.append(linkTable, extraCols, axis=1)
            linkTable[:, cfg.LTB_LCPLEN] = -1
            linkTable[:, cfg.LTB_CWDEUCR] = -1
            linkTable[:, cfg.LTB_CWDPATHR] = -1
            linkTable[:, cfg.LTB_EFFRESIST] = -1
            linkTable[:, cfg.LTB_CWDTORR] = -1
            linkTable[:, cfg.LTB_CURRENT] = -1
            del extraCols

        # set up directories for circuit and circuit mosaic grids
        PREFIX = cfg.PREFIX

        # Create output geodatabase
        if not arcpy.Exists(cfg.PINCHGDB):
            arcpy.CreateFileGDB_management(cfg.OUTPUTDIR,
                                            path.basename(cfg.PINCHGDB))

        mosaicRaster = path.join(cfg.CIRCUITBASEDIR, "current_mos")
        coresToProcess = npy.unique(linkTable[:, cfg.LTB_CORE1:cfg.LTB_CORE2 + 1])
        maxCoreNum = max(coresToProcess)
        del coresToProcess

        lu.dashline(0)
        coreList = linkTable[:,cfg.LTB_CORE1:cfg.LTB_CORE2+1]
        coreList = npy.sort(coreList)
        #gprint('There are ' + str(len(npy.unique(coreList))) ' core areas.')

        INCIRCUITDIR = cfg.CIRCUITBASEDIR
        OUTCIRCUITDIR = path.join(cfg.CIRCUITBASEDIR,
                                  cfg.CIRCUITOUTPUTDIR_NM)
        CONFIGDIR = path.join(INCIRCUITDIR, cfg.CIRCUITCONFIGDIR_NM)

        if cfg.SQUARERESISTANCES:
            # Square resistance values
            squaredRaster = path.join(cfg.SCRATCHDIR,'res_sqr')
            arcpy.env.workspace = cfg.SCRATCHDIR
            arcpy.env.scratchWorkspace = cfg.ARCSCRATCHDIR
            outRas = Raster(resRaster) * Raster(resRaster)
            outRas.save(squaredRaster)
            resRaster = squaredRaster

        if cfg.DO_ADJACENTPAIRS:
            pctDone = 0
            linkLoop = 0
            lu.dashline(1)
            gprint('Mapping pinch points in individual corridors \n'
                    'using Circuitscape')
            lu.dashline(0)
            gprint('0 percent done')
            for x in range(0,numLinks):
                linkId = str(int(linkTable[x,cfg.LTB_LINKID]))
                if not (linkTable[x,cfg.LTB_LINKTYPE] > 0):
                    continue
                linkDir = path.join(cfg.SCRATCHDIR, 'link' + linkId)
                lu.create_dir(linkDir)
                start_time1 = time.clock()
                pctDone = lu.report_pct_done(linkLoop, numCorridorLinks,
                                             pctDone)
                linkLoop = linkLoop + 1

                # source and target cores
                corex=int(coreList[x,0])
                corey=int(coreList[x,1])

                # Get cwd rasters for source and target cores
                cwdRaster1 = lu.get_cwd_path(corex)
                cwdRaster2 = lu.get_cwd_path(corey)

                lccNormRaster = path.join(linkDir, 'lcc_norm')
                arcpy.env.extent = "MINOF"

                link = lu.get_links_from_core_pairs(linkTable, corex,
                                                    corey)
                lcDist = float(linkTable[link,cfg.LTB_CWDIST])

                # Normalized lcc rasters are created by adding cwd rasters
                # and subtracting the least cost distance between them.
                # expression = (cwdRaster1 + " + " + cwdRaster2 + " - "
                             # + lcDist)
                # gp.singleOutputMapAlgebra_sa(expression, lccNormRaster)

                outRas = Raster(cwdRaster1) + Raster(cwdRaster2) - lcDist
                outRas.save(lccNormRaster)

                #create raster mask
                resMaskRaster = path.join(linkDir, 'res_mask')

                #create raster mask
                # expression = ("(con(" + lccNormRaster + " <= " +
                              # str(cfg.CWDCUTOFF) + ", 1))")
                outCon = arcpy.sa.Con(Raster(lccNormRaster) <= cfg.CWDCUTOFF, 1)
                outCon.save(resMaskRaster)
                # gp.singleOutputMapAlgebra_sa(expression, resMaskRaster)

                # Convert to poly.  Use as mask to clip resistance raster.
                resMaskPoly = path.join(linkDir,
                                        'res_mask_poly.shp')
                arcpy.RasterToPolygon_conversion(resMaskRaster, resMaskPoly,
                                              "NO_SIMPLIFY")
                arcpy.env.extent = resMaskPoly

                resClipRasterMasked = path.join(linkDir,
                                                'res_clip_m')
                outRas = arcpy.sa.ExtractByMask(resRaster, resMaskPoly)

                outRas.save(resClipRasterMasked)

                resNpyFN = 'resistances_link_' + linkId + '.npy'
                resNpyFile = path.join(INCIRCUITDIR, resNpyFN)
                numElements = export_ras_to_npy(resClipRasterMasked,resNpyFile)
                corePairRaster = path.join(linkDir, 'core_pairs')

                arcpy.env.extent = resClipRasterMasked
                # expression = ("con(" + cwdRaster1 + " == 0, " + str(corex)
                    # + ", con(" + cwdRaster2 + " == 0, "
                    # + str(corey) + " + 0.0))")

                # Next result needs to be floating pt for numpy export
                outCon = arcpy.sa.Con(Raster(cwdRaster1) == 0, corex,
                            arcpy.sa.Con(Raster(cwdRaster2) == 0, corey + 0.0))
                outCon.save(corePairRaster)

                coreNpyFN = 'cores_link_' + linkId + '.npy'
                coreNpyFile = path.join(INCIRCUITDIR, coreNpyFN)
                numElements = export_ras_to_npy(corePairRaster,coreNpyFile)

                arcpy.env.extent = "MINOF"

                # Set circuitscape options and call
                options = lu.setCircuitscapeOptions()
                options['habitat_file'] = resNpyFile
                options['point_file'] = coreNpyFile
                options['set_focal_node_currents_to_zero']=True
                outputFN = 'Circuitscape_link' + linkId + '.out'
                options['output_file'] = path.join(OUTCIRCUITDIR, outputFN)
                if numElements > 250000:
                    options['print_timings']=True
                configFN = 'pinchpoint_config' + linkId + '.ini'

                outConfigFile = path.join(CONFIGDIR, configFN)
                lu.writeCircuitscapeConfigFile(outConfigFile, options)

                gprint('Calling Circuitscape...')
                if numElements > 250000:
                    test = subprocess.call([CSPATH, outConfigFile],
                                creationflags = subprocess.CREATE_NEW_CONSOLE)
                else:
                    subprocess.call([CSPATH, outConfigFile], shell=True)
                # Read in results                
                currentFN = ('Circuitscape_link' + linkId 
                            + '_cum_curmap.npy')
                currentMap = path.join(OUTCIRCUITDIR, currentFN)
                
                if not arcpy.Exists(currentMap):
                    gprint('Circuitscape failed.  Trying again in 10 seconds.')
                    lu.snooze(10)
                    test = subprocess.call([CSPATH, outConfigFile], 
                                creationflags = subprocess.CREATE_NEW_CONSOLE)
                

                    currentFN = ('Circuitscape_link' + linkId 
                                + '_cum_curmap.npy')
                    currentMap = path.join(OUTCIRCUITDIR, currentFN)
                
                if not arcpy.Exists(currentMap):
                    lu.dashline(1)
                    msg = ('ERROR: It looks like Circuitscape failed.  '
                         '\nResistance raster values that vary by > 6 orders of'
                         '\nmagnitude  can cause this, as can a mismatch in '
                         '\ncore area and resistance raster extents.')

                    arcpy.AddError(msg)
                    lu.write_log(msg)
                    exit(1)

                # Either set core areas to nodata in current map or
                # divide each by its radius
                currentRaster = path.join(linkDir, "current")
                import_npy_to_ras(currentMap,corePairRaster,currentRaster)

                arcpy.env.extent = currentRaster

                if SETCORESTONULL:
                    # Set core areas to NoData in current map for color ramping
                    currentRaster2 = currentRaster + '2'
                    outCon = arcpy.sa.Con(arcpy.sa.IsNull(Raster(corePairRaster)), Raster(currentRaster))
                    outCon.save(currentRaster2)
                    currentRaster = currentRaster2
                arcpy.env.extent = "MAXOF"
                if linkLoop == 1:
                    lu.delete_data(mosaicRaster)
                    arcpy.CopyRaster_management(currentRaster,
                                                 mosaicRaster)
                else:
                    arcpy.Mosaic_management(currentRaster,
                                         mosaicRaster, "MAXIMUM", "MATCH")

                resistancesFN = ('Circuitscape_link' + linkId
                            + '_resistances_3columns.out')

                resistancesFile = path.join(OUTCIRCUITDIR,resistancesFN)
                resistances = npy.loadtxt(resistancesFile,
                                          dtype = 'Float64', comments='#')

                resistance = float(str(arcpy.env.cellSize)) * resistances[2]
                linkTable[link,cfg.LTB_EFFRESIST] = resistance

                # Ratio
                if not cfg.SQUARERESISTANCES:
                    linkTable[link,cfg.LTB_CWDTORR] = (linkTable[link,
                           cfg.LTB_CWDIST] / linkTable[link,cfg.LTB_EFFRESIST])
                # Clean up
                lu.delete_file(coreNpyFile)
                lu.delete_file(resNpyFile)
                lu.delete_data(currentRaster) 
                lu.delete_dir(linkDir) 
                gprint('Finished with link #' + str(linkId))
                start_time1 = lu.elapsed_time(start_time1)

            outputGDB = path.join(cfg.OUTPUTDIR, path.basename(cfg.PINCHGDB))
            outputRaster = path.join(outputGDB,
                                     PREFIX + "_current_adjacent_pairs")
            lu.delete_data(outputRaster)
            statement = 'arcpy.CopyRaster_management(mosaicRaster, outputRaster)'
            count = 0
            while True:
                try: exec statement
                except:
                    count,tryAgain = lu.retry_arc_error(count,statement)
                    if not tryAgain: exec statement
                else: break


            gprint('Building output statistics and pyramids '
                                  'for corridor pinch point raster\n')
            lu.build_stats(outputRaster)

            finalLinkTable = lu.update_lcp_shapefile(linkTable, lastStep=5,
                                                      thisStep=8)

            linkTableFile = path.join(cfg.DATAPASSDIR, "linkTable_s8.csv")
            lu.write_link_table(finalLinkTable, linkTableFile, inLinkTableFile)
            linkTableFinalFile = path.join(cfg.OUTPUTDIR,
                                           PREFIX + "_linkTable_s8.csv")
            lu.write_link_table(finalLinkTable,
                                linkTableFinalFile, inLinkTableFile)
            gprint('Copy of linkTable written to '+
                              linkTableFinalFile)
            #fixme: update sticks?

            gprint('Creating shapefiles with linework for links.')
            lu.write_link_maps(linkTableFinalFile, step=8)

            # Copy final link maps to gdb.
            lu.copy_final_link_maps(step=8)

            lu.delete_data(mosaicRaster)

        if not cfg.DO_ALLPAIRS:
            # Clean up temporary files
            if not cfg.SAVECURRENTMAPS:
                lu.delete_dir(OUTCIRCUITDIR)
            return

        lu.dashline(1)
        gprint('Mapping global pinch points among all\n'
                'core area pairs using Circuitscape')
        lu.dashline(0)
        arcpy.env.workspace = cfg.SCRATCHDIR
        arcpy.env.scratchWorkspace = cfg.ARCSCRATCHDIR
        arcpy.env.extent = cfg.RESRAST
        arcpy.env.cellSize = cfg.RESRAST

        S8CORE_RAS = "s8core_ras"
        s8CoreRasPath = path.join(cfg.SCRATCHDIR,S8CORE_RAS)

        arcpy.FeatureToRaster_conversion(cfg.COREFC, cfg.COREFN,
                                         s8CoreRasPath, arcpy.env.cellSize)
        binaryCoreRaster = path.join(cfg.SCRATCHDIR,"core_ras_bin")

        # The following commands cause file lock problems on save.  using gp
        # instead.
        # outCon = arcpy.sa.Con(S8CORE_RAS, 1, "#", "VALUE > 0")
        # outCon.save(binaryCoreRaster)
        #gp.Con_sa(s8CoreRasPath, 1, binaryCoreRaster, "#", "VALUE > 0")
        outCon = arcpy.sa.Con(Raster(s8CoreRasPath) > 0, 1)
        outCon.save(binaryCoreRaster)
        s5corridorRas = path.join(cfg.OUTPUTGDB,PREFIX + "_lcc_mosaic_int")

        outCon = arcpy.sa.Con(Raster(s5corridorRas) <= cfg.CWDCUTOFF, Raster(
                              resRaster), arcpy.sa.Con(Raster(
                              binaryCoreRaster) > 0, Raster(resRaster)))

        resRasClipPath = path.join(cfg.SCRATCHDIR,'res_ras_clip')
        outCon.save(resRasClipPath)

        arcpy.env.cellSize = resRasClipPath
        arcpy.env.extent = resRasClipPath
        s8CoreRasClipped = s8CoreRasPath + '_c'

        # Produce core raster with same extent as clipped resistance raster
#bbb        outRas = Raster(s8CoreRasPath) * 1
#bbb        outRas.save(s8CoreRasClipped)
        outCon = arcpy.sa.Con(arcpy.sa.IsNull(Raster(s8CoreRasPath)), -9999, Raster(s8CoreRasPath))  #bbb added to ensure correct data type- nodata values were positive for cores otherwise
        outCon.save(s8CoreRasClipped)

        resNpyFN = 'resistances.npy'
        resNpyFile = path.join(INCIRCUITDIR, resNpyFN)
        numElements = export_ras_to_npy(resRasClipPath,resNpyFile)

        coreNpyFN = 'cores.npy'
        coreNpyFile = path.join(INCIRCUITDIR, coreNpyFN)
        numElements = export_ras_to_npy(s8CoreRasClipped,coreNpyFile)

        arcpy.env.extent = "MINOF"

        options = lu.setCircuitscapeOptions()
        options['scenario']='all-to-one'
        options['scenario']='pairwise'
        options['habitat_file'] = resNpyFile
        options['point_file'] = coreNpyFile
        options['set_focal_node_currents_to_zero']=True
        outputFN = 'Circuitscape.out'
        options['output_file'] = path.join(OUTCIRCUITDIR, outputFN)
        options['print_timings']=True
        configFN = 'pinchpoint_allpair_config.ini'
        outConfigFile = path.join(CONFIGDIR, configFN)
        lu.writeCircuitscapeConfigFile(outConfigFile, options)

        gprint('Calling Circuitscape...')
        test = subprocess.call([CSPATH, outConfigFile],
                               creationflags = subprocess.CREATE_NEW_CONSOLE)

        currentFNs = ['Circuitscape_cum_curmap.npy',
                      'Circuitscape_max_curmap.npy']
        if options['scenario']=='pairwise':
            rasterSuffixes =  ["_cum_current_all_pairs","_max_current_all_pairs"]
        else:
            rasterSuffixes =  ["_cum_current_all_to_one","_max_current_all_to_one"]
        for i in range(0,2):
            currentFN = currentFNs[i]
            currentMap = path.join(OUTCIRCUITDIR, currentFN)
            outputGDB = path.join(cfg.OUTPUTDIR, path.basename(cfg.PINCHGDB))
            outputRaster = path.join(outputGDB, PREFIX + rasterSuffixes[i])
            currentRaster = path.join(cfg.SCRATCHDIR, "current")

            try:
                import_npy_to_ras(currentMap,resRasClipPath,outputRaster)
            except:
                lu.dashline(1)
                msg = ('ERROR: Could not open Circuitscape output.  If '
                     '\nresistance raster values vary by > 6 orders of'
                     '\nmagnitude, that may have caused Circuitscape to fail.')
                arcpy.AddError(msg)
                lu.write_log(msg)
                exit(1)

            #set core areas to nodata 
            if SETCORESTONULL:                  
                # Set core areas to NoData in current map for color ramping
                outputRasterND = outputRaster + '_nodata_cores' 
                outCon = arcpy.sa.SetNull(Raster(s8CoreRasClipped) > 0, Raster(outputRaster))   
                outCon.save(outputRasterND)                

            gprint('\nBuilding output statistics and pyramids ' 
                    'for pinch point raster ' + str(i + 1))        
            lu.build_stats(outputRaster)
            lu.build_stats(outputRasterND)

        # Clean up temporary files
        if not cfg.SAVECURRENTMAPS:
            lu.delete_dir(OUTCIRCUITDIR)


    # Return GEOPROCESSING specific errors
    except arcpy.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 8. Details follow.****')
        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 8. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)



def export_ras_to_npy(raster,npyFile):
    try:
        descData=arcpy.Describe(raster)
        cellSize=descData.meanCellHeight
        extent=descData.Extent
        spatialReference=descData.spatialReference
        pnt=arcpy.Point(extent.XMin,extent.YMin)
        outData = arcpy.RasterToNumPyArray(raster,"#","#","#",-9999)

        if npy.array_equiv(outData, outData.astype('int32')):
            outData = outData.astype('int32')
        npy.save(npyFile, outData)
        write_header(raster,outData,npyFile)
        numElements = (outData.shape[0] * outData.shape[1])
        del outData
        return numElements
    # Return GEOPROCESSING specific errors
    except arcpy.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 8. Details follow.****')
        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 8. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)

def import_npy_to_ras(npyFile,baseRaster,outRasterPath):
    try:
        npyArray = npy.load(npyFile, mmap_mode=None)
        npyArray=npyArray.astype('float32')
        descData=arcpy.Describe(baseRaster)
        cellSize=descData.meanCellHeight
        extent=descData.Extent
        spatialReference=descData.spatialReference
        pnt=arcpy.Point(extent.XMin,extent.YMin)
        newRaster = arcpy.NumPyArrayToRaster(npyArray,pnt,
                                             cellSize,cellSize,-9999)
        newRaster.save(outRasterPath)
        return
    # Return GEOPROCESSING specific errors
    except arcpy.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 8. Details follow.****')
        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 8. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)

def write_header(raster,numpyArray,numpyFile):
        ncols=numpyArray.shape[1]
        nrows=numpyArray.shape[0]
        descData=arcpy.Describe(raster)
        cellSize=descData.meanCellHeight
        extent=descData.Extent
        xllcorner = extent.XMin
        yllcorner = extent.YMin
        nodata = -9999
        fileBase, fileExtension = path.splitext(numpyFile)
        headerFile = fileBase + '.hdr'
        f = False
        f = open(headerFile, 'w')

        f.write('ncols         ' + str(ncols) + '\n')
        f.write('nrows         ' + str(nrows) + '\n')
        f.write('xllcorner     ' + str(xllcorner) + '\n')
        f.write('yllcorner     ' + str(yllcorner) + '\n')
        f.write('cellsize      ' + str(cellSize) + '\n')
        f.write('NODATA_value  ' + str(nodata) + '\n')

        f.close()
