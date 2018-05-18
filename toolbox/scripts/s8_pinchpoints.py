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

import numpy as npy
from lm_retry_decorator import retry
import arcpy
from arcpy.sa import *

from lm_config import tool_env as cfg
import lm_util as lu

_SCRIPT_NAME = "s8_pinchpoints.py"

arcpy.CheckOutExtension("spatial")

SETCORESTONULL = True
gprint = lu.gprint
# gwarn = arcpy.AddWarning
tif = ".tif"
#tif = ""


@retry(2)
def STEP8_calc_pinchpoints():
    """ Maps pinch points in Linkage Mapper corridors using Circuitscape
        given CWD calculations from s3_calcCwds.py.

    """
    try:
        lu.dashline(0)
        gprint('Running script ' + _SCRIPT_NAME)
        
        restartFlag = False
        if cfg.CWDCUTOFF < 0:
            cfg.CWDCUTOFF = cfg.CWDCUTOFF * -1
            restartFlag = True # Restart code in progress
                
        outputGDB = path.join(cfg.OUTPUTDIR, path.basename(cfg.PINCHGDB))
        
        arcpy.OverWriteOutput = True
        arcpy.env.workspace = cfg.SCRATCHDIR
        arcpy.env.scratchWorkspace = cfg.ARCSCRATCHDIR
        arcpy.env.pyramid = "NONE"
        arcpy.env.rasterstatistics = "NONE"

        # set the analysis extent and cell size to that of the resistance
        # surface
        arcpy.env.extent = cfg.RESRAST
        arcpy.env.cellSize = cfg.RESRAST
        arcpy.snapraster = cfg.RESRAST

        resRaster = cfg.RESRAST
        arcpy.env.extent = "MINOF"

        
        minObject = arcpy.GetRasterProperties_management(resRaster, "MINIMUM") 
        rasterMin = float(str(minObject.getOutput(0)))
        if rasterMin <= 0:
            msg = ('Error: resistance raster cannot have 0 or negative values.')
            lu.raise_error(msg)
                
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
        # Create output geodatabase
        if not arcpy.Exists(cfg.PINCHGDB):
            arcpy.CreateFileGDB_management(cfg.OUTPUTDIR,
                                            path.basename(cfg.PINCHGDB))

        mosaicRaster = path.join(cfg.CIRCUITBASEDIR, "current_mos" + tif)
        coresToProcess = npy.unique(
                                linkTable[:, cfg.LTB_CORE1:cfg.LTB_CORE2 + 1])
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

        # Cutoff value text to append to filenames
        cutoffText = str(cfg.CWDCUTOFF)
        if cutoffText[-6:] == '000000':
            cutoffText = cutoffText[0:-6]+'m' 
        elif cutoffText[-3:] == '000':
            cutoffText = cutoffText[0:-3]+'k' 

        if cfg.SQUARERESISTANCES:
            # Square resistance values
            squaredRaster = path.join(cfg.SCRATCHDIR,'res_sqr')
            arcpy.env.workspace = cfg.SCRATCHDIR
            arcpy.env.scratchWorkspace = cfg.ARCSCRATCHDIR
            outRas = Raster(resRaster) * Raster(resRaster)
            outRas.save(squaredRaster)
            resRaster = squaredRaster

        if cfg.DO_ADJACENTPAIRS:
            linkLoop = 0
            lu.dashline(1)
            gprint('Mapping pinch points in individual corridors \n'
                    'using Circuitscape.')
            lu.dashline(1)
            gprint('If you try to cancel your run and the Arc dialog hangs, ')
            gprint('you can kill Circuitscape by opening Windows Task Manager')
            gprint('and ending the cs_run.exe process.')                    
            lu.dashline(2)

            for x in range(0,numLinks):            
                linkId = str(int(linkTable[x,cfg.LTB_LINKID]))
                if not (linkTable[x,cfg.LTB_LINKTYPE] > 0):
                    continue
                linkLoop = linkLoop + 1
                linkDir = path.join(cfg.SCRATCHDIR, 'link' + linkId)
                if restartFlag == True and path.exists(linkDir):
                    gprint('continuing')
                    continue
                restartFlag = False
                lu.create_dir(linkDir)
                start_time1 = time.clock()

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
                outRas = Raster(cwdRaster1) + Raster(cwdRaster2) - lcDist
                outRas.save(lccNormRaster)

                #create raster mask
                resMaskRaster = path.join(linkDir, 'res_mask'+tif)

                #create raster mask
                outCon = arcpy.sa.Con(Raster(lccNormRaster) <= cfg.CWDCUTOFF, 1)
                outCon.save(resMaskRaster)

                # Convert to poly.  Use as mask to clip resistance raster.
                resMaskPoly = path.join(linkDir,
                                        'res_mask_poly.shp')
                arcpy.RasterToPolygon_conversion(resMaskRaster, resMaskPoly,
                                              "NO_SIMPLIFY")
                arcpy.env.extent = resMaskPoly

                # Includes 0 values in some cases with CP LI model if tif
                # so using ESRI Grid format
                resClipRasterMasked = path.join(linkDir,
                                                'res_clip_m') 
                # Extract masked resistance raster.  
                # Needs to be float to get export to npy to work.
                outRas = arcpy.sa.ExtractByMask(resRaster, resMaskPoly) + 0.0 
                outRas.save(resClipRasterMasked)
               
                resNpyFN = 'resistances_link_' + linkId + '.npy'
                resNpyFile = path.join(INCIRCUITDIR, resNpyFN)
                numElements, numResistanceNodes = export_ras_to_npy(resClipRasterMasked,
                                                          resNpyFile)
                
                totMem, availMem = lu.get_mem()
                # gprint('Total memory: str(totMem))
                if numResistanceNodes / availMem > 2000000:
                    lu.dashline(1)
                    lu.warn('Warning:')
                    lu.warn('Circuitscape can only solve 2-3 million nodes')
                    lu.warn('per gigabyte of available RAM. \nTotal physical RAM'
                            ' on your machine is ~' + str(totMem) 
                            + ' GB. \nAvailable memory is ~'+ str(availMem) 
                            + ' GB. \nYour resistance raster has '
                            + str(numResistanceNodes) + ' nodes.')                                                          
                    lu.dashline(2)
                corePairRaster = path.join(linkDir, 'core_pairs'+tif)
                arcpy.env.extent = resClipRasterMasked

                # Next result needs to be floating pt for numpy export
                outCon = arcpy.sa.Con(Raster(cwdRaster1) == 0, corex,
                            arcpy.sa.Con(Raster(cwdRaster2) == 0, corey + 0.0))
                outCon.save(corePairRaster)

                coreNpyFN = 'cores_link_' + linkId + '.npy'
                coreNpyFile = path.join(INCIRCUITDIR, coreNpyFN)
                numElements, numNodes = export_ras_to_npy(corePairRaster,
                                                          coreNpyFile)

                arcpy.env.extent = "MINOF"

                # Set circuitscape options and call
                options = lu.setCircuitscapeOptions()
                if cfg.WRITE_VOLT_MAPS == True:
                    options['write_volt_maps']=True
                options['habitat_file'] = resNpyFile
                
                # if int(linkId) > 2:
                    # options['habitat_file'] = 'c:\\test.dummy'
                                
                options['point_file'] = coreNpyFile
                options['set_focal_node_currents_to_zero']=True
                outputFN = 'Circuitscape_link' + linkId + '.out'
                options['output_file'] = path.join(OUTCIRCUITDIR, outputFN)
                if numElements > 250000:
                    options['print_timings']=True
                configFN = 'pinchpoint_config' + linkId + '.ini'

                outConfigFile = path.join(CONFIGDIR, configFN)
                lu.writeCircuitscapeConfigFile(outConfigFile, options)                    
                gprint('Processing link ID #' + str(linkId) + '. Resistance map'
                        ' has ' + str(int(numResistanceNodes)) + ' nodes.') 

                memFlag = lu.call_circuitscape(cfg.CSPATH, outConfigFile)

                currentFN = ('Circuitscape_link' + linkId 
                            + '_cum_curmap.npy')
                currentMap = path.join(OUTCIRCUITDIR, currentFN)
                
                if not arcpy.Exists(currentMap):
                    print_failure(numResistanceNodes, memFlag, 10)
                    numElements, numNodes = export_ras_to_npy(
                                                resClipRasterMasked,resNpyFile)
                    memFlag = lu.call_circuitscape(cfg.CSPATH, outConfigFile)

                    currentFN = ('Circuitscape_link' + linkId 
                                + '_cum_curmap.npy')
                    currentMap = path.join(OUTCIRCUITDIR, currentFN)
                
                if not arcpy.Exists(currentMap):                
                    msg = ('\nCircuitscape failed. See error information above.')
                    arcpy.AddError(msg)
                    lu.write_log(msg)
                    exit(1)

                # Either set core areas to nodata in current map or
                # divide each by its radius
                currentRaster = path.join(linkDir, "current" + tif)
                import_npy_to_ras(currentMap,corePairRaster,currentRaster)
                
                if cfg.WRITE_VOLT_MAPS == True:
                    voltFN = ('Circuitscape_link' + linkId + '_voltmap_'
                           + str(corex) + '_'+str(corey) + '.npy')
                    voltMap = path.join(OUTCIRCUITDIR, voltFN)
                    voltRaster = path.join(outputGDB,
                             cfg.PREFIX + "_voltMap_"+ str(corex) + '_'+str(corey))                
                    import_npy_to_ras(voltMap,corePairRaster,voltRaster)
                    gprint('Building output statistics and pyramids '
                                   'for voltage raster\n')
                    lu.build_stats(voltRaster) 
                    
                arcpy.env.extent = currentRaster

                if SETCORESTONULL:
                    # Set core areas to NoData in current map for color ramping
                    currentRaster2 = currentRaster + '2' + tif
                    outCon = arcpy.sa.Con(arcpy.sa.IsNull(Raster
                                      (corePairRaster)), Raster(currentRaster))
                    outCon.save(currentRaster2)
                    currentRaster = currentRaster2
                arcpy.env.extent = "MAXOF"
                if linkLoop == 1:
                    lu.delete_data(mosaicRaster)
                    @retry(10)
                    def copyRas2():
                        arcpy.CopyRaster_management(currentRaster,
                                                    mosaicRaster)
                    copyRas2()
                else:
                    @retry(10)
                    def mosaicRas():                
                        arcpy.Mosaic_management(currentRaster,
                                         mosaicRaster, "MAXIMUM", "MATCH")
                    mosaicRas()
                    
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
                if cfg.SAVE_TEMP_CIRCUIT_FILES == False:
                    lu.delete_file(coreNpyFile)
                    coreNpyBase, extension = path.splitext(coreNpyFile)
                    lu.delete_data(coreNpyBase + '.hdr')                    
                    lu.delete_file(resNpyFile)
                    resNpyBase, extension = path.splitext(resNpyFile)
                    lu.delete_data(resNpyBase + '.hdr')                    
                    lu.delete_file(currentMap)
                    curMapBase, extension = path.splitext(currentMap)
                    lu.delete_data(curMapBase + '.hdr')
                    lu.delete_data(currentRaster) 
                    lu.clean_out_workspace(linkDir)
                    lu.delete_dir(linkDir) 
                gprint('Finished with link ID #' + str(linkId) + '. ' + 
                        str(linkLoop) + ' out of ' + str(numCorridorLinks) + 
                        ' links have been processed.')
                start_time1 = lu.elapsed_time(start_time1)
                
            outputRaster = path.join(outputGDB, cfg.PREFIX + 
                                     "_current_adjacentPairs_" + cutoffText)
            lu.delete_data(outputRaster)
            
            @retry(10)
            def copyRas():
                arcpy.CopyRaster_management(mosaicRaster, outputRaster)
            copyRas()

            gprint('Building output statistics and pyramids '
                                  'for corridor pinch point raster\n')
            lu.build_stats(outputRaster)
            
            finalLinkTable = lu.update_lcp_shapefile(linkTable, lastStep=5,
                                                      thisStep=8)

            linkTableFile = path.join(cfg.DATAPASSDIR, "linkTable_s5_plus.csv")
            lu.write_link_table(finalLinkTable, linkTableFile, inLinkTableFile)
            linkTableFinalFile = path.join(cfg.OUTPUTDIR, cfg.PREFIX + 
                                           "_linkTable_s5_plus.csv")
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
                'core area pairs using Circuitscape.')                   
                
        if cfg.ALL_PAIR_SCENARIO=='pairwise':
            gprint('Circuitscape will be run in PAIRWISE mode.')
                        
        else:
            gprint('Circuitscape will be run in ALL-TO-ONE mode.')     
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
        # gp.Con_sa(s8CoreRasPath, 1, binaryCoreRaster, "#", "VALUE > 0")
        outCon = arcpy.sa.Con(Raster(s8CoreRasPath) > 0, 1)
        outCon.save(binaryCoreRaster)
        s5corridorRas = path.join(cfg.OUTPUTGDB,cfg.PREFIX + "_corridors")
        
        if not arcpy.Exists(s5corridorRas):
            s5corridorRas = path.join(cfg.OUTPUTGDB,cfg.PREFIX + 
                                      "_lcc_mosaic_int")

        outCon = arcpy.sa.Con(Raster(s5corridorRas) <= cfg.CWDCUTOFF, Raster(
                              resRaster), arcpy.sa.Con(Raster(
                              binaryCoreRaster) > 0, Raster(resRaster)))

        resRasClipPath = path.join(cfg.SCRATCHDIR,'res_ras_clip')
        outCon.save(resRasClipPath)

        arcpy.env.cellSize = resRasClipPath
        arcpy.env.extent = resRasClipPath
        s8CoreRasClipped = s8CoreRasPath + '_c'

        # Produce core raster with same extent as clipped resistance raster
        # added to ensure correct data type- nodata values were positive for 
        # cores otherwise
        outCon = arcpy.sa.Con(arcpy.sa.IsNull(Raster(s8CoreRasPath)), 
                              -9999, Raster(s8CoreRasPath))  
        outCon.save(s8CoreRasClipped)

        resNpyFN = 'resistances.npy'
        resNpyFile = path.join(INCIRCUITDIR, resNpyFN)
        numElements, numResistanceNodes = export_ras_to_npy(resRasClipPath,resNpyFile)

        totMem, availMem = lu.get_mem()
        # gprint('Total memory: str(totMem))
        if numResistanceNodes / availMem > 2000000:
            lu.dashline(1)
            lu.warn('Warning:')
            lu.warn('Circuitscape can only solve 2-3 million nodes')
            lu.warn('per gigabyte of available RAM. \nTotal physical RAM '
                    'on your machine is ~' + str(totMem)
                    + ' GB. \nAvailable memory is ~'+ str(availMem)
                    + ' GB. \nYour resistance raster has '
                    + str(numResistanceNodes) + ' nodes.')   
            lu.dashline(0)

        coreNpyFN = 'cores.npy'
        coreNpyFile = path.join(INCIRCUITDIR, coreNpyFN)
        numElements, numNodes = export_ras_to_npy(s8CoreRasClipped,coreNpyFile)

        arcpy.env.extent = "MINOF"

        options = lu.setCircuitscapeOptions()
        options['scenario']=cfg.ALL_PAIR_SCENARIO
        options['habitat_file'] = resNpyFile
        options['point_file'] = coreNpyFile
        options['set_focal_node_currents_to_zero']=True
        outputFN = 'Circuitscape.out'
        options['output_file'] = path.join(OUTCIRCUITDIR, outputFN)
        options['print_timings']=True
        configFN = 'pinchpoint_allpair_config.ini'
        outConfigFile = path.join(CONFIGDIR, configFN)
        lu.writeCircuitscapeConfigFile(outConfigFile, options)
        gprint('\nResistance map has ' + str(int(numResistanceNodes)) + ' nodes.') 
        lu.dashline(1)
        gprint('If you try to cancel your run and the Arc dialog hangs, ')
        gprint('you can kill Circuitscape by opening Windows Task Manager')
        gprint('and ending the cs_run.exe process.')             
        lu.dashline(0)

        lu.call_circuitscape(cfg.CSPATH, outConfigFile)

        if options['scenario']=='pairwise':
            rasterSuffix =  "_current_allPairs_" + cutoffText
                        
        else:
            rasterSuffix =  "_current_allToOne_" + cutoffText

        currentFN = 'Circuitscape_cum_curmap.npy'
        currentMap = path.join(OUTCIRCUITDIR, currentFN)
        outputRaster = path.join(outputGDB, cfg.PREFIX + rasterSuffix)
        currentRaster = path.join(cfg.SCRATCHDIR, "current")

        try:
            import_npy_to_ras(currentMap,resRasClipPath,outputRaster)
        except:
            lu.dashline(1)
            msg = ('ERROR: Circuitscape failed. \n'
                  'Note: Circuitscape can only solve 2-3 million nodes'
                  '\nper gigabyte of available RAM. The resistance '
                  '\nraster for the last corridor had '
                  + str(numResistanceNodes) + ' nodes.\n\nResistance '
                  'raster values that vary by >6 orders of \nmagnitude'
                  ' can also cause failures, as can a mismatch in '
                  '\ncore area and resistance raster extents.')
            arcpy.AddError(msg)
            lu.write_log(msg)
            exit(1)

        #set core areas to nodata 
        if SETCORESTONULL:                  
            # Set core areas to NoData in current map for color ramping
            outputRasterND = outputRaster + '_noDataCores' 
            outCon = arcpy.sa.SetNull(Raster(s8CoreRasClipped) > 0, 
                                      Raster(outputRaster))   
            outCon.save(outputRasterND)                

        gprint('\nBuilding output statistics and pyramids ' 
                'for centrality raster.')        
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


@retry(10)
def export_ras_to_npy(raster,npyFile):
    descData=arcpy.Describe(raster)
    cellSize=descData.meanCellHeight
    extent=descData.Extent
    spatialReference=descData.spatialReference
    
    pnt=arcpy.Point(extent.XMin,extent.YMin)
    outData = arcpy.RasterToNumPyArray(raster,"#","#","#",-9999)
    #outData = npy.where(outData==noDataVal,-9999,outData)
    if npy.array_equiv(outData, outData.astype('int32')):
        outData = outData.astype('int32')
    npy.save(npyFile, outData)
    write_header(raster,outData,npyFile)
            
    numElements = (outData.shape[0] * outData.shape[1])
    #rows,cols = npy.where(outData != -9999)
    numNodes = (npy.where(outData != -9999, 1, 0)).sum() 
    #numZeros = (npy.where(outData != -9999, 1, 0)).sum() 
    #del rows
    
    del outData
    return numElements, numNodes

@retry(10)
def import_npy_to_ras(npyFile,baseRaster,outRasterPath):
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
        
    # # Return GEOPROCESSING specific errors
    # except arcpy.ExecuteError:
        # lu.dashline(1)
        # gprint('****Failed in step 8. Details follow.****')
        # lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # # Return any PYTHON or system specific errors
    # except:
        # lu.dashline(1)
        # gprint('****Failed in step 8. Details follow.****')
        # lu.exit_with_python_error(_SCRIPT_NAME)

@retry(10)
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

def print_failure(numResistanceNodes, memFlag, sleepTime):
    gprint('\nCircuitscape failed. See error information above.')
    if memFlag == True:
        totMem, availMem = lu.get_mem()                    
        gprint('Note: Circuitscape can only solve 2-3 million nodes')
        gprint('per gigabyte of available RAM. Your resistance raster had ')
        gprint(str(int(numResistanceNodes)) + ' nodes.\n')  
        gprint('Total physical RAM on your machine is ~' 
               + str(totMem) 
               + ' GB. \nAvailable memory is ~'
               + str(availMem) + ' GB. \n')
    gprint('Trying again in ' + str(sleepTime) + ' seconds.')
    lu.snooze(sleepTime)
