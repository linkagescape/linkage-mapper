#NEED cores and r's in same # cells!
#cs call failing but works otherwise

#print number of nodata cells being exported
#option to coarsen resistance grids

#!/usr/bin/env python2.5

"""Maps pinch points using Circuitscape given CWD calculations from
       s3_calcCwds.py.
"""

__filename__ = "s7_pinchpoints.py"
__version__ = "0.6.4"

import os.path as path
import os
import time
import shutil
import ConfigParser
import arcgisscripting
import numpy as npy
import math
import subprocess

try:
    import arcpy
except:
    arcpy = False

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
LTB_EFFRESIST = Cfg.LTB_EFFRESIST
LTB_CWDTORR = Cfg.LTB_CWDTORR

ALLPAIRS = True
NORMALIZECORECURRENTS = False
SETCORESTONULL = True
    
def STEP7_calc_pinchpoints():
    """ Experimental code map pinch points using Circuitscape 
        given CWD calculations from s3_calcCwds.py.

    """
    try:     

    
        lu.dashline(0)
        gprint('Running script ' + __filename__)
            
        gp.workspace = Cfg.SCRATCHDIR

        # set the analysis extent and cell size to that of the resistance
        # surface
        gp.Extent = gp.Describe(Cfg.RESRAST).Extent
        gp.CellSize = gp.Describe(Cfg.RESRAST).MeanCellHeight

        gp.Extent = "MINOF"

        gp.snapraster = Cfg.RESRAST
        resRaster = Cfg.RESRAST 

        inLinkTableFile = lu.get_prev_step_link_table(step=7)
        linkTable = lu.load_link_table(inLinkTableFile)
        numLinks = linkTable.shape[0]
        numCorridorLinks = lu.report_links(linkTable)
        if numCorridorLinks == 0:
            lu.dashline()
            gprint('\nThere are no linkages. Bailing.')
            time.sleep(5)
            return
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
        
        # For speed:
        gp.pyramid = "NONE"
        gp.rasterstatistics = "NONE"

        # set up directories for circuit and circuit mosaic grids
        PREFIX = Cfg.PREFIX       

        # Create output geodatabase
        if not gp.exists(Cfg.PINCHGDB):
            gp.createfilegdb(Cfg.OUTPUTDIR, path.basename(Cfg.PINCHGDB))
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
            # expression = resRaster + " * " + resRaster
            cwdRaster = path.join(Cfg.CWDGDB, PREFIX + "_cwd")
            expression = resRaster + " * " + cwdRaster + " + 1 "
            
            squaredRaster = path.join(Cfg.SCRATCHDIR,'res_squared')
            gp.SingleOutputMapAlgebra_sa(expression, squaredRaster) 
            resRaster = squaredRaster    
            
# do global centrality?  either multiply indiv current rasters by link centrality
# or do all pairs (would give current in core areas...).        
        
        if ALLPAIRS == True:
            gp.extent = "MAXOF"
            S7CORE_RAS = "s7core_ras"
            gp.FeatureToRaster_conversion(Cfg.COREFC, Cfg.COREFN, S7CORE_RAS, gp.Cellsize)
            binaryCoreRaster = "core_ras_bin"
            gp.Con_sa(S7CORE_RAS, 1, binaryCoreRaster, "#", "VALUE > 0")

            
            
            s5corridorRas = os.path.join(Cfg.OUTPUTGDB,PREFIX + "_lcc_mosaic_int")
            expression = ("(con(" + s5corridorRas + "<= " +
                          str(Cfg.CWDCUTOFF) + ", "+ resRaster + ", con(" + binaryCoreRaster + " > 0, " + resRaster +  ")))")            
            resRasClip = 'resRasClip'
            gp.SingleOutputMapAlgebra_sa(expression, resRasClip)
            
        
            # export resistance map cut with linkage map and hcas
            # export hcas
            # do cs all pairs.
            s7CoreRasPath = os.path.join(Cfg.SCRATCHDIR,S7CORE_RAS)
            resRasClipPath = os.path.join(Cfg.SCRATCHDIR,resRasClip)
            resAsciiFN = 'resistances.asc'
            coreAsciiFN = 'cores.asc'
            gp.Extent = resClipRasterMasked 
            resAsciiFile = path.join(INCIRCUITDIR, resAsciiFN)
            gp.RasterToASCII_conversion(resRasClipPath, resAsciiFile)
            coreAsciiFile = path.join(INCIRCUITDIR, coreAsciiFN)
            gp.RasterToASCII_conversion(s7CoreRasPath, coreAsciiFile)
            gp.extent = "MINOF"


            options = lu.setCircuitscapeOptions() 
            options['habitat_file'] = resAsciiFile
            options['point_file'] = coreAsciiFile
            outputFN = 'Circuitscape.out'
            options['output_file'] = path.join(OUTCIRCUITDIR, outputFN)
            configFN = 'pinchpoint_allpair_config.ini'
            outConfigFile = path.join(CONFIGDIR, configFN)
            lu.writeCircuitscapeConfigFile(outConfigFile, options)

            csPath = lu.getCsPath()
            if csPath == None:
                msg = ('Cannot find an installation of Circuitscape 3.5.5' 
                        '\nor greater in your Program Files directory.') 
                gp.AddError(msg)
                gp.AddMessage(gp.GetMessages(2))
                exit(1)
            gprint('Calling Circuitscape...')
            #subprocess.check_call(systemCall, shell=True)
            #subprocess.call([csPath, outConfigFile], shell=True)                     
            test = subprocess.check_call([csPath, outConfigFile], shell=True)                     
            gprint(test)
            
            currentFN = ('Circuitscape_cum_curmap.asc')
            currentMap = path.join(OUTCIRCUITDIR, currentFN)

            outputGDB = path.join(Cfg.OUTPUTDIR, path.basename(Cfg.PINCHGDB))
            outputRaster = path.join(outputGDB, PREFIX + "_current_all_pairs")
            Cfg.gp.CopyRaster_management(currentMap, outputRaster)

            try:
                gp.addmessage('\nBuilding output statistics and pyramids ' 
                                  'for corridor raster\n')        
                gp.CalculateStatistics_management(outputRaster, "1", "1", "#")
                gp.BuildPyramids_management(outputRaster)    
            except:
                pass
            return
            



            
        pctDone = 0
        linkLoop = 0
        gprint('Mapping pinch points in corridors using Circuitscape....')
        gprint('0 percent done')        
        for x in range(0,numLinks):
            linkId = str(int(linkTable[x,LTB_LINKID]))
            if (linkTable[x,LTB_LINKTYPE] > 0):
                start_time = time.clock()
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
                gp.Extent = "MINOF"

                link = lu.get_links_from_core_pairs(linkTable, corex, 
                                                    corey)
                lcDist = str(linkTable[link,LTB_CWDIST])

                # Normalized lcc rasters are created by adding cwd rasters 
                # and subtracting the least cost distance between them.
                expression = (cwdRaster1 + " + " + cwdRaster2 + " - " 
                             + lcDist)
                gp.SingleOutputMapAlgebra_sa(expression, lccNormRaster)

                #create raster mask 
                resMaskRaster = path.join(Cfg.SCRATCHDIR, 'res_mask')
#                expression = ("(con(" + lccNormRaster + "<= " +
#                              str(Cfg.CWDCUTOFF) + ",1))")
                
                #create raster mask 
                #include core areas
                expression = ("(con(" + lccNormRaster + " <= " +
                              str(Cfg.CWDCUTOFF) + ", 1, con(" + cwdRaster1 + " == 0, 1, con(" + cwdRaster2 + " == 0, 1))))")                
                
                gp.SingleOutputMapAlgebra_sa(expression, resMaskRaster)
                    
                # Convert to poly.  Use as mask to clip resistance raster.
                resMaskPoly = path.join(Cfg.SCRATCHDIR, 
                                        'res_mask_poly.shp')
                gp.RasterToPolygon_conversion(resMaskRaster, resMaskPoly, 
                                              "NO_SIMPLIFY")
                                                             
                resClipRasterMasked = path.join(Cfg.SCRATCHDIR, 
                                                'res_clip_m')
                gp.ExtractByMask_sa(resRaster, resMaskPoly, 
                                    resClipRasterMasked)

                # Raster to ascii
                resAsciiFN = 'resistances_link_' + linkId + '.asc'
                resAsciiFile = path.join(INCIRCUITDIR, resAsciiFN)
                gp.RasterToASCII_conversion(resClipRasterMasked, 
                                            resAsciiFile)

                corePairRaster = path.join(Cfg.SCRATCHDIR, 'core_pairs') 
                
                gp.Extent = resClipRasterMasked 
                expression = ("con(" + cwdRaster1 + " == 0, " + str(corex)  
                    + ", con(" + cwdRaster2 + " == 0, " 
                    + str(corey) + "))")
                gp.SingleOutputMapAlgebra_sa(expression, corePairRaster)
                
                coreAsciiFN = 'cores_link_' + linkId + '.asc'
                coreAsciiFile = path.join(INCIRCUITDIR, coreAsciiFN)
                gp.RasterToASCII_conversion(corePairRaster, coreAsciiFile)
                        
                gp.Extent = "MINOF"
                
                options = lu.setCircuitscapeOptions() 
                options['habitat_file'] = resAsciiFile
                options['point_file'] = coreAsciiFile
                outputFN = 'Circuitscape_link' + linkId + '.out'
                options['output_file'] = path.join(OUTCIRCUITDIR, outputFN)
                configFN = 'pinchpoint_config' + linkId + '.ini'
                outConfigFile = path.join(CONFIGDIR, configFN)
                lu.writeCircuitscapeConfigFile(outConfigFile, options)

                #to detect architecture:
                # import platform
                # gprint(str(platform.architecture()))
                # gprint('done with prep')
                # start_time = lu.elapsed_time(start_time)

                csPath = lu.getCsPath()
                if csPath == None:
                    msg = ('Cannot find an installation of Circuitscape 3.5.5' 
                            '\nor greater in your Program Files directory.') 
                    gp.AddError(msg)
                    gp.AddMessage(gp.GetMessages(2))
                    exit(1)
                #gprint('Calling Circuitscape...')
                #subprocess.check_call(systemCall, shell=True)
                subprocess.call([csPath, outConfigFile], shell=True)                     
                
                # start_time = lu.elapsed_time(start_time)
                currentFN = ('Circuitscape_link' + linkId 
                            + '_cum_curmap.asc')
                currentMap = path.join(OUTCIRCUITDIR, currentFN)
        
                
                # Either set core areas to nodata in current map or
                # divide each by its radius
                currentRaster = path.join(Cfg.SCRATCHDIR, "current")
                gp.Extent = gp.Describe(currentMap).Extent
                SETCORESTONULL = True
                if not arcpy:
                    NORMALIZECORECURRENTS = False
                if SETCORESTONULL == True:                  
                    # Set core areas to NoData in current map for color ramping
                    isNullExpression = ("con(isnull(" + corePairRaster + "), " 
                        + currentMap + ")")   
                    gp.SingleOutputMapAlgebra_sa(isNullExpression, 
                                                 currentRaster)
                elif NORMALIZECORECURRENTS == True:
                    # Experimental- spread current over core areas 
                    # dividing by their 'diameter'.
                    corePairArray = arcpy.RasterToNumPyArray(corePairRaster)
                    # Get areas of cores:
                    numCellsX = (corePairArray==corex).sum()
                    numCellsY = (corePairArray==corey).sum()
                    
                    diaX = 2 * npy.sqrt(numCellsX / math.pi)
                    diaY = 2 * npy.sqrt(numCellsY / math.pi)
                    expression = ("con(isnull("+corePairRaster+"), 1, con(" + corePairRaster + " == " + str(corex) + ", " + str(diaX) + 
                        ", con(" + corePairRaster + " == " + str(corey) + ", " + str(diaY) + ")))")
                    divRaster = path.join(Cfg.SCRATCHDIR, 'divRaster')
                    gp.SingleOutputMapAlgebra_sa(expression, divRaster)
                    
                    expression = (currentMap + " / " + divRaster)
                    gp.SingleOutputMapAlgebra_sa(expression, currentRaster)                    
                
                gp.extent = "MAXOF"
                if linkLoop == 1:
                    Cfg.gp.CopyRaster_management(currentRaster, 
                                                 mosaicRaster)
                else:
                    # Note: cannot use SOMA to mosaic. It is a different
                    # process entirely.
                    count = 0
                    statement = ('Cfg.gp.Mosaic_management(currentRaster, '
                                 'mosaicRaster, "MAXIMUM", "MATCH")')
                    while True:
                        try: exec statement
                        except:
                            count,tryAgain = lu.hiccup_test(count,
                                                            statement)
                            if not tryAgain: exec statement
                        else: break                    

                resistancesFN = ('Circuitscape_link' + linkId 
                            + '_resistances_3columns.out')
                        
                resistancesFile = path.join(OUTCIRCUITDIR,resistancesFN)
                resistances = npy.loadtxt(resistancesFile, 
                                          dtype = 'Float64', comments='#')

                resistance = float(str(gp.CellSize)) * resistances[2]
                linkTable[link,LTB_EFFRESIST] = resistance
                # Ratio
                if Cfg.SQUARERESISTANCES == True:
                    linkTable[link,LTB_CWDTORR] = -1
                else:
                    linkTable[link,LTB_CWDTORR] = (linkTable[link,LTB_CWDIST] /
                                                 linkTable[link,LTB_EFFRESIST])
                
                # Clean up (fixme: put these in function)
                lu.delete_file(coreAsciiFile)
                lu.delete_file(resAsciiFile)
                if not Cfg.SAVECURRENTMAPS:
                    lu.delete_file(currentMap)
                gprint('Finished with link #' + str(linkId))
                start_time = lu.elapsed_time(start_time)
        outputGDB = path.join(Cfg.OUTPUTDIR, path.basename(Cfg.PINCHGDB))
        outputRaster = path.join(outputGDB, PREFIX + "_current_mos")
        Cfg.gp.CopyRaster_management(mosaicRaster, outputRaster)

        try:
            gp.addmessage('\nBuilding output statistics and pyramids ' 
                              'for corridor raster\n')        
            gp.CalculateStatistics_management(outputRaster, "1", "1", "#")
            gp.BuildPyramids_management(outputRaster)    
        except:
            pass
        
        finalLinkTable = lu.update_lcp_shapefile(linkTable, lastStep=5,
                                                  thisStep=7) 
                                                  
        linkTableFile = path.join(Cfg.DATAPASSDIR, "linkTable_s7.csv")
        lu.write_link_table(finalLinkTable, linkTableFile, inLinkTableFile)
        linkTableFinalFile = path.join(Cfg.OUTPUTDIR, 
                                       PREFIX + "_linkTable_s7.csv")
        lu.write_link_table(finalLinkTable, 
                            linkTableFinalFile, inLinkTableFile)
        gprint('Copy of final linkTable written to '+
                          linkTableFinalFile)
        #fixme: update sticks?
        
        gprint('Creating shapefiles with linework for links.')
        lu.write_link_maps(linkTableFinalFile, step=7) 
        
        # Copy final link maps to gdb.  
        lu.copy_final_link_maps(step=7)
        
        # Clean up temporary files
        # if not Cfg.SAVECURRENTMAPS:
            # lu.delete_dir(OUTCIRCUITDIR)
        lu.delete_dir(Cfg.SCRATCHDIR)
        lu.delete_data(mosaicRaster)


    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 7. Details follow.****')
        lu.raise_geoproc_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 7. Details follow.****')
        lu.raise_python_error(__filename__)
    return

    
    
        # Use linkage map mask out all portions of resrast that 
        # are above CWDCUTOFF
        ### REMOVED FOR NOW- no speed improvement found with LI data
        # cwdRaster = path.join(Cfg.OUTPUTDIR,"cwd")
        # if gp.exists(lccMosaicRaster):
            # resRasterLccMask = path.join(Cfg.SCRATCHDIR, "res_lcc_mask")
            # expression = ("(" + resRaster + " * (con(" + lccMosaicRaster +
                           # "<= " + str(CWDCUTOFF) + ",1)))")
            # gp.SingleOutputMapAlgebra_sa(expression, resRasterLccMask)
            # resRaster = resRasterLccMask
            # gprint('masked')
        # elif gp.exists(cwdRaster): # If linkages doesn't exist use cwd raster
            # resRasterCwdMask = path.join(Cfg.SCRATCHDIR, "res_cwd_mask")
            # expression = ("(" + resRaster + " * (con(" + cwdRaster + "<= " +
                                # str(CWDCUTOFF) + ",1)))")
            # gp.SingleOutputMapAlgebra_sa(expression, resRasterCwdMask)
            # resRaster = resRasterCwdMask
            # gprint('masked cwd')    