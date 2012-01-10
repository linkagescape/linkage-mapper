#!/usr/bin/env python2.5
# Author: Brad McRae 

"""Maps pinch points using Circuitscape given CWD calculations from
       s3_calcCwds.py.
"""

__filename__ = "s8_pinchpoints.py"
__version__ = "0.7.6"

import os.path as path
import os
import time
import shutil
import ConfigParser
import numpy as npy
import math
import subprocess
import gc

import arcpy
gprint = arcpy.AddMessage
gp = arcpy.gp

from lm_config import Config as Cfg
import lm_util as lu


# Set local references to objects and constants from config file
LTB_CORE1 = Cfg.LTB_CORE1
LTB_CORE2 = Cfg.LTB_CORE2
LTB_LINKID = Cfg.LTB_LINKID
LTB_LINKTYPE = Cfg.LTB_LINKTYPE
LTB_CWDIST = Cfg.LTB_CWDIST
LTB_EFFRESIST = Cfg.LTB_EFFRESIST
LTB_CWDTORR = Cfg.LTB_CWDTORR

DO_ALLPAIRS = Cfg.DO_ALLPAIRS
DO_ADJACENTPAIRS = True
NORMALIZECORECURRENTS = False
SETCORESTONULL = True

    
def STEP8_calc_pinchpoints():
    """ Experimental code map pinch points using Circuitscape 
        given CWD calculations from s3_calcCwds.py.

    """
    try:     
        gc.collect()
        lu.dashline(0)
        gprint('Running script ' + __filename__)        

        CSPATH = lu.get_cs_path()
        if CSPATH == None:
            msg = ('Cannot find an installation of Circuitscape 3.5.5' 
                    '\nor greater in your Program Files directory.') 
            arcpy.AddError(msg)
            gprint(arcpy.GetMessages(2))
            exit(1)

        arcpy.OverWriteOutput = True            
        gp.OverWriteOutput = True #sometimes still use gp object
        arcpy.env.workspace = Cfg.SCRATCHDIR
        arcpy.env.scratchWorkspace = Cfg.SCRATCHDIR
        # For speed:
        arcpy.env.pyramid = "NONE"
        arcpy.env.rasterstatistics = "NONE"

        # set the analysis extent and cell size to that of the resistance
        # surface
        gprint (Cfg.RESRAST)
        arcpy.env.extent = Cfg.RESRAST
        arcpy.env.cellSize = Cfg.RESRAST

        arcpy.env.extent = "MINOF"

        arcpy.snapraster = Cfg.RESRAST
        resRaster = Cfg.RESRAST 

        if DO_ADJACENTPAIRS == True:
            prevLcpShapefile = lu.get_lcp_shapefile(None, thisStep = 8)
            if not arcpy.Exists(prevLcpShapefile):
                msg = ('Cannot find an LCP shapefile from step 5.  Please '
                        'rerun that step and any previous ones if necessary.') 
                gp.AddError(msg)
                gp.AddMessage(gp.GetMessages(2))
                exit(1)

            # Remove lcp shapefile
            lcpShapefile = os.path.join(Cfg.DATAPASSDIR, "lcpLines_s8.shp")        
            lu.delete_data(lcpShapefile)
        
        
        
        
        inLinkTableFile = lu.get_prev_step_link_table(step=8)
        linkTable = lu.load_link_table(inLinkTableFile)
        numLinks = linkTable.shape[0]
        numCorridorLinks = lu.report_links(linkTable)
        if numCorridorLinks == 0:
            dashline(1)
            msg =('\nThere are no linkages. Bailing.')
            gp.AddError(msg)
            exit(1)
        if linkTable.shape[1] < 16: # If linktable has no entries from prior
                                    # centrality or pinchpint analyses
            extraCols = npy.zeros((numLinks, 6), dtype="float64")
            linkTable = linkTable[:,0:10]
            linkTable = npy.append(linkTable, extraCols, axis=1)
            linkTable[:, Cfg.LTB_LCPLEN] = -1
            linkTable[:, Cfg.LTB_CWDEUCR] = -1
            linkTable[:, Cfg.LTB_CWDPATHR] = -1
            linkTable[:, LTB_EFFRESIST] = -1
            linkTable[:, Cfg.LTB_CWDTORR] = -1
            linkTable[:, Cfg.LTB_CURRENT] = -1
            del extraCols
        
        # set up directories for circuit and circuit mosaic grids
        PREFIX = Cfg.PREFIX       

        # Create output geodatabase
        if not arcpy.Exists(Cfg.PINCHGDB):
            arcpy.CreateFileGDB_management(Cfg.OUTPUTDIR, 
                                            path.basename(Cfg.PINCHGDB))

        mosaicRaster = path.join(Cfg.CIRCUITBASEDIR, "current_mos")        
        coresToProcess = npy.unique(linkTable[:, LTB_CORE1:LTB_CORE2 + 1])
        maxCoreNum = max(coresToProcess)
        del coresToProcess

        lu.dashline(0)
        coreList = linkTable[:,LTB_CORE1:LTB_CORE2+1]
        coreList = npy.sort(coreList)
        #gprint('There are ' + str(len(npy.unique(coreList))) ' core areas.')

        INCIRCUITDIR = Cfg.CIRCUITBASEDIR
        OUTCIRCUITDIR = path.join(Cfg.CIRCUITBASEDIR, 
                                  Cfg.CIRCUITOUTPUTDIR_NM)
        CONFIGDIR = path.join(INCIRCUITDIR, Cfg.CIRCUITCONFIGDIR_NM)             
        
        if Cfg.SQUARERESISTANCES == True:
            # Square resistance values
            expression = (resRaster + " * " + resRaster)
            squaredRaster = path.join(Cfg.SCRATCHDIR,'res_sqr')
            arcpy.env.workspace = Cfg.SCRATCHDIR
            arcpy.env.scratchWorkspace = Cfg.SCRATCHDIR
            gp.singleOutputMapAlgebra_sa(expression, squaredRaster) 
            resRaster = squaredRaster    
            
        if DO_ADJACENTPAIRS == True:
            pctDone = 0
            linkLoop = 0
            lu.dashline(1)            
            gprint('Mapping pinch points in individual corridors \n'
                    'using Circuitscape')
            lu.dashline(0)            
            gprint('0 percent done')        
            for x in range(0,numLinks):
                linkId = str(int(linkTable[x,LTB_LINKID]))
                if not (linkTable[x,LTB_LINKTYPE] > 0):
                    continue
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

                lccNormRaster = path.join(Cfg.SCRATCHDIR, 'lcc_norm')
                arcpy.env.extent = "MINOF"

                link = lu.get_links_from_core_pairs(linkTable, corex, 
                                                    corey)
                lcDist = str(linkTable[link,LTB_CWDIST])



                # Normalized lcc rasters are created by adding cwd rasters 
                # and subtracting the least cost distance between them.
                expression = (cwdRaster1 + " + " + cwdRaster2 + " - " 
                             + lcDist)                
                gp.singleOutputMapAlgebra_sa(expression, lccNormRaster)

                #create raster mask 
                resMaskRaster = path.join(Cfg.SCRATCHDIR, 'res_mask')
                
                #create raster mask 
                expression = ("(con(" + lccNormRaster + " <= " +
                              str(Cfg.CWDCUTOFF) + ", 1))")                

                gp.singleOutputMapAlgebra_sa(expression, resMaskRaster)
                    
                # Convert to poly.  Use as mask to clip resistance raster.
                resMaskPoly = path.join(Cfg.SCRATCHDIR, 
                                        'res_mask_poly.shp')
                arcpy.RasterToPolygon_conversion(resMaskRaster, resMaskPoly, 
                                              "NO_SIMPLIFY")
                arcpy.env.extent = resMaskPoly
                
                resClipRasterMasked = path.join(Cfg.SCRATCHDIR, 
                                                'res_clip_m')
                arcpy.gp.ExtractByMask_sa(resRaster, resMaskPoly, 
                                    resClipRasterMasked)

                                    
                resNpyFN = 'resistances_link_' + linkId + '.npy'
                resNpyFile = path.join(INCIRCUITDIR, resNpyFN)
                numElements = export_ras_to_npy(resClipRasterMasked,resNpyFile)
                corePairRaster = path.join(Cfg.SCRATCHDIR, 'core_pairs') 
                
                arcpy.env.extent = resClipRasterMasked 
                expression = ("con(" + cwdRaster1 + " == 0, " + str(corex)  
                    + ", con(" + cwdRaster2 + " == 0, " 
                    + str(corey) + " + 0.0))")  # Needs to be floating pt 
                                                # for numpy export
                
                #gp.OverWriteOutput = True #shouldn't have to re-set...
                gp.SingleOutputMapAlgebra_sa(expression, corePairRaster)
                
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
                    lu.dashline(1)
                    msg = ('ERROR: It looks like Circuitscape failed.  If '
                         '\nresistance raster values vary by > 6 orders of'
                         '\nmagnitude, that may have caused this.')
                    arcpy.AddError(msg)
                    exit(1)
               
                # Either set core areas to nodata in current map or
                # divide each by its radius
                currentRaster = path.join(Cfg.SCRATCHDIR, "current")
                import_npy_to_ras(currentMap,corePairRaster,currentRaster)

                arcpy.env.extent = currentRaster

                if SETCORESTONULL == True:                  
                    # Set core areas to NoData in current map for color ramping
                    isNullExpression = ("con(isnull(" + corePairRaster + "), " 
                        + currentRaster + ")")   
                    currentRaster2 = currentRaster + '2'
                    gp.SingleOutputMapAlgebra_sa(isNullExpression, 
                                                 currentRaster2)
                    currentRaster = currentRaster2                  
                    
                arcpy.env.extent = "MAXOF"
                if linkLoop == 1: 
                    arcpy.CopyRaster_management(currentRaster, 
                                                 mosaicRaster)
                else:
                    gp.Mosaic_management(currentRaster, mosaicRaster, 
                                         "MAXIMUM", "MATCH")            
                    
                resistancesFN = ('Circuitscape_link' + linkId 
                            + '_resistances_3columns.out')
                        
                resistancesFile = path.join(OUTCIRCUITDIR,resistancesFN)
                resistances = npy.loadtxt(resistancesFile, 
                                          dtype = 'Float64', comments='#')

                resistance = float(str(arcpy.env.cellSize)) * resistances[2]
                linkTable[link,LTB_EFFRESIST] = resistance

                # Ratio
                if Cfg.SQUARERESISTANCES == True:
                    linkTable[link,LTB_CWDTORR] = -1
                else:
                    linkTable[link,LTB_CWDTORR] = (linkTable[link,LTB_CWDIST] /
                                                 linkTable[link,LTB_EFFRESIST])
                
                # Clean up 
                lu.delete_file(coreNpyFile)
                lu.delete_file(resNpyFile)
                gprint('Finished with link #' + str(linkId))
                start_time1 = lu.elapsed_time(start_time1)
                
            outputGDB = path.join(Cfg.OUTPUTDIR, path.basename(Cfg.PINCHGDB))
            outputRaster = path.join(outputGDB, 
                                     PREFIX + "_current_adjacent_pairs")            
            arcpy.CopyRaster_management(mosaicRaster, outputRaster)

            gprint('Building output statistics and pyramids ' 
                                  'for corridor pinch point raster\n')        
            lu.build_stats(outputRaster)
            
            finalLinkTable = lu.update_lcp_shapefile(linkTable, lastStep=5,
                                                      thisStep=8) 
                                                      
            linkTableFile = path.join(Cfg.DATAPASSDIR, "linkTable_s8.csv")
            lu.write_link_table(finalLinkTable, linkTableFile, inLinkTableFile)
            linkTableFinalFile = path.join(Cfg.OUTPUTDIR, 
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
            
        if DO_ALLPAIRS == False:
            # Clean up temporary files
            if not Cfg.SAVECURRENTMAPS:
                lu.delete_dir(OUTCIRCUITDIR)
            return

        lu.dashline(1)
        gprint('Mapping global pinch points among all\n'
                'core area pairs using Circuitscape')
        lu.dashline(0)
        arcpy.env.workspace = Cfg.SCRATCHDIR
        arcpy.env.scratchWorkspace = Cfg.SCRATCHDIR        
        S8CORE_RAS = "s8core_ras"
        s8CoreRasPath = os.path.join(Cfg.SCRATCHDIR,S8CORE_RAS)
        
        arcpy.FeatureToRaster_conversion(Cfg.COREFC, Cfg.COREFN, 
                                         s8CoreRasPath, arcpy.env.cellSize)
        binaryCoreRaster = "core_ras_bin"

        # The following commands cause file lock problems on save.  using gp
        # instead.
        # outCon = arcpy.sa.Con(S8CORE_RAS, 1, "#", "VALUE > 0")            
        # outCon.save(binaryCoreRaster)
        gp.Con_sa(s8CoreRasPath, 1, binaryCoreRaster, "#", "VALUE > 0")  

        s5corridorRas = os.path.join(Cfg.OUTPUTGDB,PREFIX + "_lcc_mosaic_int")
        expression = ("(con(" + s5corridorRas + "<= " +
                      str(Cfg.CWDCUTOFF) + ", "+ resRaster + ", con(" + 
                      binaryCoreRaster + " > 0, " + resRaster +  ")))")            
        
        resRasClip = 'res_ras_clip'
        resRasClipPath = os.path.join(Cfg.SCRATCHDIR,resRasClip)
        gp.singleOutputMapAlgebra_sa(expression, resRasClipPath)        
        arcpy.env.cellSize = resRasClipPath
        arcpy.env.extent = resRasClipPath  
        s8CoreRasClipped = s8CoreRasPath + '_c'
        expression = (s8CoreRasPath + " * 1")
        gp.singleOutputMapAlgebra_sa(expression, s8CoreRasClipped)

        resNpyFN = 'resistances.npy'
        resNpyFile = path.join(INCIRCUITDIR, resNpyFN)
        numElements = export_ras_to_npy(resRasClipPath,resNpyFile)
        
        coreNpyFN = 'cores.npy'
        coreNpyFile = path.join(INCIRCUITDIR, coreNpyFN)
        numElements = export_ras_to_npy(s8CoreRasClipped,coreNpyFile)

        arcpy.env.extent = "MINOF"
        
        options = lu.setCircuitscapeOptions() 
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
        rasterSuffixes =  ["_cum_current_all_pairs","_max_current_all_pairs"]
        for i in range(0,2):
            currentFN = currentFNs[i]
            currentMap = path.join(OUTCIRCUITDIR, currentFN)
            outputGDB = path.join(Cfg.OUTPUTDIR, path.basename(Cfg.PINCHGDB))
            outputRaster = path.join(outputGDB, PREFIX + rasterSuffixes[i])
            currentRaster = path.join(Cfg.SCRATCHDIR, "current")

            try:
                import_npy_to_ras(currentMap,resRasClipPath,outputRaster)
            except:
                lu.dashline(1)
                msg = ('ERROR: Could not open Circuitscape output.  If '
                     '\nresistance raster values vary by > 6 orders of'
                     '\nmagnitude, that may have caused Circuitscape to fail.')
                arcpy.AddError(msg)
                exit(1)
            
            gprint('\nBuilding output statistics and pyramids ' 
                    'for pinch point raster ' + str(i + 1))        
            lu.build_stats(outputRaster)
        
        # Clean up temporary files
        if not Cfg.SAVECURRENTMAPS:
            lu.delete_dir(OUTCIRCUITDIR)
    
        
    # Return GEOPROCESSING specific errors
    except arcpy.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 8. Details follow.****')
        lu.raise_geoproc_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 8. Details follow.****')
        lu.raise_python_error(__filename__)

        
        
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
        lu.raise_geoproc_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 8. Details follow.****')
        lu.raise_python_error(__filename__)
    
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
        lu.raise_geoproc_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 8. Details follow.****')
        lu.raise_python_error(__filename__)
    
def write_header(raster,numpyArray,numpyFile):
        ncols=numpyArray.shape[1]
        nrows=numpyArray.shape[0]
        descData=arcpy.Describe(raster)
        cellSize=descData.meanCellHeight
        extent=descData.Extent
        xllcorner = extent.XMin
        yllcorner = extent.YMin
        nodata = -9999
        fileBase, fileExtension = os.path.splitext(numpyFile) 
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
    
    
        ### REMOVED FOR NOW- no speed improvement found with LI data
        # cwdRaster = path.join(Cfg.OUTPUTDIR,"cwd")
        # if arcpy.Exists(lccMosaicRaster):
            # resRasterLccMask = path.join(Cfg.SCRATCHDIR, "res_lcc_mask")
            # expression = ("(" + resRaster + " * (con(" + lccMosaicRaster +
                           # "<= " + str(CWDCUTOFF) + ",1)))")
            # arcpy.SingleOutputMapAlgebra_sa(expression, resRasterLccMask)
            # resRaster = resRasterLccMask
            # gprint('masked')
        # elif arcpy.Exists(cwdRaster): # If linkages doesn't exist use cwd raster
            # resRasterCwdMask = path.join(Cfg.SCRATCHDIR, "res_cwd_mask")
            # expression = ("(" + resRaster + " * (con(" + cwdRaster + "<= " +
                                # str(CWDCUTOFF) + ",1)))")
            # arcpy.SingleOutputMapAlgebra_sa(expression, resRasterCwdMask)
            # resRaster = resRasterCwdMask
            # gprint('masked cwd')    
            
                #to detect architecture:
                # import platform
                # gprint(str(platform.architecture()))
                # gprint('done with prep')
                # start_time = lu.elapsed_time(start_time)            
                
#time trial.  0 secs vs 7 in demo data.                       
                # gprint('npy')
                # start_time = time.clock()
                # for i in range(1,20):
                    # data = arcpy.RasterToNumPyArray(currentRaster)
                    # outfile = "c://temp//test.npy"
                    # npy.save(outfile, data)
                # start_time = lu.elapsed_time(start_time)
                # gprint('aasc')
                # start_time = time.clock()
                # for i in range(1,20):                
                    # outfile2 = "c://temp//test.asc"
                    # arcpy.RasterToASCII_conversion(currentRaster,)
                # start_time = lu.elapsed_time(start_time)                
                
                
                # ##################experimental- insert in place of reg code
                    # resClipRasterMasked = path.join(Cfg.SCRATCHDIR, 
                                                    # 'res_clip_m')
                    # expression = (cwdRaster1 + " + " + cwdRaster2)                
                    # gp.singleOutputMapAlgebra_sa(expression, lccNormRaster)
                                        
                    # gprint('lcdist='+str(lcDist))                   
                    # expression = (resRaster + " + " + "( " + lccNormRaster + " / " + str(lcDist) + " )")
                    # #expression = (resRaster + " * " + lccNormRaster)
                    # gprint(expression)
                    # squaredRaster = path.join(Cfg.SCRATCHDIR,'res_sqr')
                    # arcpy.env.workspace = Cfg.SCRATCHDIR
                    # arcpy.env.scratchWorkspace = Cfg.SCRATCHDIR
                    # gp.singleOutputMapAlgebra_sa(expression, squaredRaster) 
                    # resClipRasterMasked = squaredRaster    
                # ###################

                # if NORMALIZECORECURRENTS == True:
                    # # Experimental- spread current over core areas 
                    # # dividing by their 'diameter'.
                    # corePairArray = arcpy.RasterToNumPyArray(corePairRaster)
                    # # Get areas of cores:
                    # numCellsX = (corePairArray==corex).sum()
                    # numCellsY = (corePairArray==corey).sum()
                    
                    # diaX = 2 * npy.sqrt(numCellsX / math.pi)
                    # diaY = 2 * npy.sqrt(numCellsY / math.pi)
                    # expression = ("con(isnull("+corePairRaster+"), 1, con(" + corePairRaster + " == " + str(corex) + ", " + str(diaX) + 
                        # ", con(" + corePairRaster + " == " + str(corey) + ", " + str(diaY) + ")))")
                    # divRaster = path.join(Cfg.SCRATCHDIR, 'divRaster')
                    # gp.SingleOutputMapAlgebra_sa(expression, divRaster)
                    
                    # expression = (currentRaster + " / " + divRaster)
                    # gp.SingleOutputMapAlgebra_sa(expression, currentRaster) 