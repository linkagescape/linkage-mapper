#!/usr/bin/env python2.6

"""Script to iteratively run linkage mapper and barrier mapper tools restoring for max ROI"""

import sys
import os
import shutil
import time

import arcpy
from arcpy.sa import *
import numpy as npy

import lm_master
import barrier_master
from lm_config import tool_env as cfg
import lm_util as lu

gprint = lu.gprint

_SCRIPT_NAME = os.path.basename(__file__)

arcpy.CheckOutExtension("spatial")

def main():
    """Iterates over LM, BM, and restoration tasks"""

    ## USER SETTINGS ######################################################
    ## Restoration Settings
    ## ALL input data must be in the same projection
    start_time = time.clock()
    restoreMaxROI = False # Set to True to restore highest ROI
                         # Set to False to restore strongest barrier
    restoredResistanceVal = 1 # Resistance value of restored habitat.  Must be 1 or greater.
    restorationDataGDB = "C:\\barrierClassAnalysis\\RestorationINPUTS_July2013.gdb" # No spaces or special chars in paths or gdb names
    outputDir = "C:\\barrierClassAnalysis\\output" # No spaces in path, avoid using dropbox or network drive
                                                   # Project directories will be created in this (iter1, iter2...)
                                                   # as will an output geodatabase
    resistanceRaster = "URWA_resis"# Resistance raster.  Should be in input GDB
    coreFC = 'URWA_HCAs_Doug_Grant'# Core area feature class. Should be in input GDB 'URWA_HCAs_Doug_Grant'
    coreFN = 'HCA_ID' # Core area field name
    
    radius = 450 # restoration radius in meters
    iterations = 13 # number of restorations to perform
    minAgThreshold = 0.75 # if less than this proportion of ag in circle, don't consider restoring circle
    minImprovementVal = 0 # Don't consider barriers below this improvement score (average improvement per meter diameter restored)
    parcelCostRaster = 'DougGrantParcelCost_m2_projected_90m' # Average per-m2 parcel cost per pixel. Snapped to resistance raster.
    restorationCostRaster = 'restCostPer_m2' # Right now this is just a raster with all pixels set to 0.113174
    agRaster = "ARESmaskp_projected" # 1=Ag, 0 = not Ag
    barrierCombineMethod = 'Maximum' # Some restorations benefit multiple corridors. 
                                     # 'Maximum' takes the greatest improvement across core area pairs
                                     # 'Sum' adds improvement scores acreoss all pairs. 
    cwdThresh = None # Use cwdThresh = None for no threshold. Use cwdThresh = X to not consider 
                      # restorations more than X map units away from each core area.
    ## END USER SETTINGS ######################################################
    try:
        # Setup path and create directories
        gprint('Hey! Make sure everything is in the same projection!\n')
        gprint('Setting up paths and creating directories')
        sys.path.append('..\\toolbox\\scripts')
        resRast = os.path.join(restorationDataGDB, resistanceRaster)   
        coreFCPath = os.path.join(restorationDataGDB, coreFC)

        # Set up a NEW output gdb (leave previous ones on drive)
        for i in range (1,200):
            outputGDB = 'restorationOutput'+str(i)+'.gdb'
            if not arcpy.Exists(os.path.join(outputDir,outputGDB)):
                break
            gprint('Previous output GDB '+ outputGDB +' exists.  Delete to save disk space.')    
        arcpy.CreateFileGDB_management(outputDir,outputGDB)
        outputGDB = os.path.join(outputDir,outputGDB)
        logFile = os.path.join(outputGDB,'Iterate Barriers'+str(i)+'.py')
        shutil.copyfile(__file__, logFile) #write a copy of this file to output dir as a record of settings
        
        arcpy.env.cellSize = resRast
        arcpy.env.extent = resRast
        arcpy.env.snapRaster = resRast
        arcpy.env.overwriteOutput = True
        arcpy.env.scratchWorkspace = outputGDB
        arcpy.env.workspace = outputGDB
        
        spatialref = arcpy.Describe(resRast).spatialReference
        mapunits = spatialref.linearUnitName
        gprint('Cell size = ' + str(arcpy.env.cellSize) + ' ' + mapunits +'s')    
        
        
        # Calculate fraction of ag within radius of each pixel
        gprint('Calculating purchase cost, fraction of ag, etc within radius of each pixel.')
        agRaster = os.path.join(restorationDataGDB, agRaster)
        inNeighborhood = NbrCircle(radius, "MAP")
        arcpy.env.extent = agRaster
        outFocalStats = arcpy.sa.FocalStatistics(agRaster,
                                        inNeighborhood, "MEAN","NODATA")
        proportionAgRaster = os.path.join(outputGDB,'proportionAgRas')
        outFocalStats.save(proportionAgRaster)    
        arcpy.env.extent = resRast

        # Calculate purchase cost of circles
        parcelCostRaster = os.path.join(restorationDataGDB, parcelCostRaster)
        arcpy.env.extent = parcelCostRaster
        outFocalStats = arcpy.sa.FocalStatistics(parcelCostRaster,inNeighborhood, "MEAN","DATA")
        costFocalStatsRaster = os.path.join(outputGDB,'costFocalStatsRaster')
        outFocalStats.save(costFocalStatsRaster)
        arcpy.env.extent = resRast
        
        circleArea = float(npy.pi * radius * radius)
        outras = (Raster(costFocalStatsRaster) * circleArea)
        purchCostRaster = os.path.join(outputGDB,'purchaseCostRaster')
        outras.save(purchCostRaster)
        lu.delete_data(costFocalStatsRaster)
        
        # restCost = npy.pi * radius * radius * restCostPer_m2
        restorationCostRaster = os.path.join(restorationDataGDB, restorationCostRaster)
        outras = Raster(purchCostRaster) + (Raster(restorationCostRaster) * radius * radius * npy.pi)
        totalCostRaster = os.path.join(outputGDB,'totalCostRaster')
        outras.save(totalCostRaster)
        # lu.build_stats(totalCostRaster)
        
        # Create mask to remove areas without cost data
        arcpy.env.extent = totalCostRaster
        costMaskRaster = os.path.join(outputGDB,'costMaskRaster')
        costThresh = 0
        outCon = arcpy.sa.Con((Raster(totalCostRaster) > float(costThresh)), 1)
        outCon.save(costMaskRaster)
        arcpy.env.extent = resRast
        
        # Create mask to remove areas below ag threshold
        outCon = arcpy.sa.Con((Raster(proportionAgRaster) > float(minAgThreshold)), 1)
        agMaskRaster = os.path.join(outputGDB, 'agMaskRaster')
        outCon.save(agMaskRaster)       
        
        doStep1 = 'true'
        doStep2 = 'true'
        doStep5 = 'false'
        for iter in range(1,iterations+1): #xxx
            start_time1 = time.clock()
            arcpy.env.cellSize = resRast # Some env settings get changed by linkage mapper and must be reset here
            arcpy.env.extent = resRast
            arcpy.env.snapRaster = resRast
            arcpy.env.overwriteOutput = True
            arcpy.env.scratchWorkspace = outputGDB
            arcpy.env.workspace = outputGDB

            lu.dashline(1)
            gprint('Running iteration number '+str(iter))
            projDir = os.path.join(outputDir,'iter' + str(iter)+'Proj')    
            lu.create_dir(outputDir)
            lu.delete_dir(projDir) #xxx
            lu.create_dir(projDir)
            if iter > 1: # Copy previous s2 linktable to new project directory
                datapassDir = os.path.join(projDir,'datapass')
                lu.create_dir(datapassDir)
                projDir1 = os.path.join(outputDir,'iter1Proj')
                datapassDirIter1 = os.path.join(projDir1,'datapass')
                s2LinktableIter1 = os.path.join(datapassDirIter1 ,'linkTable_s2.csv')
                s2LinkTable = os.path.join(datapassDir ,'linkTable_s2.csv')
                shutil.copyfile(s2LinktableIter1, s2LinkTable)

            
            # Run Linkage Mapper
            distFile = os.path.join(outputDir, coreFC + '_dists.txt') # Copy distances text file from earlier LM run to the output directory- speeds things up!
            if not os.path.exists(distFile):
                if iter == 1:
                    gprint('Will calculate distance file.')
                    distFile = '#'
                else:
                    projDir1 = os.path.join(outputDir,'iter1Proj')
                    distFile1 = os.path.join(projDir1, coreFC + '_dists.txt')
                    shutil.copyfile(distFile1,distFile) # Put a copy here for future runs
                    
            arcpy.env.overwriteOutput = True
            arcpy.env.scratchWorkspace = outputGDB
            arcpy.env.workspace = outputGDB

            argv = ('lm_master.py', projDir, coreFCPath, coreFN, resRast,
                    doStep1, doStep2, 'Cost-Weighted & Euclidean', distFile,
                    'true', 'true', 'false', '4', 'Cost-Weighted', 'true',
                    doStep5, 'true', '200000', '10000', '#', '#', '#', '#') 
            gprint('Running ' + str(argv))
            cfg.lm_configured = False  # Insures lm_master uses current argv
            lm_master.lm_master(argv)    #xxx
            doStep1 = 'false' # Can skip for future iterations
            doStep2 = 'false' # Can skip for future iterations        
            doStep5 = 'false' # Skipping for future iterations
            
            startRadius = str(radius)
            endRadius = str(radius)
            radiusStep = '0'
            saveRadiusRasters= 'false'
            writePctRasters = 'false'

            argv = ('barrier_master.py', projDir, resRast, startRadius, endRadius, radiusStep, barrierCombineMethod,
                    saveRadiusRasters, writePctRasters, cwdThresh)
            gprint('Running ' + str(argv))
            barrier_master.bar_master(argv) #xxx

            arcpy.env.cellSize = resRast # Some env settings get changed by linkage mapper and must be reset here
            arcpy.env.extent = resRast
            arcpy.env.snapRaster = resRast
            arcpy.env.overwriteOutput = True
            arcpy.env.scratchWorkspace = outputGDB
            arcpy.env.workspace = outputGDB
            
            gprint('Finding restoration circles with max barrier score / ROI')
            # Find points with max ROI
            PREFIX = os.path.basename(projDir)
            if barrierCombineMethod == 'Sum':
                sumSuffix = 'Sum'
            else:
                sumSuffix = ''
            barrierFN = (PREFIX + "_BarrierCenters" + sumSuffix + "_Rad" + str(radius))
            barrierRaster = os.path.join(projDir,'output','barriers.gdb',barrierFN)
            if not arcpy.Exists(barrierRaster):
                msg = ('Error: cannot find barrier output: '+barrierRaster)
                lu.raise_error(msg)

            # arcpy.env.cellSize = agMaskRaster
            # arcpy.env.extent = agMaskRaster

            if iter > 1:
                gprint('Creating mask for previously restored areas')
                inNeighborhood = NbrCircle(radius, "MAP")
                arcpy.env.extent = allRestoredAreasRaster
                outFocalStats = arcpy.sa.FocalStatistics(allRestoredAreasRaster,inNeighborhood, "MEAN","DATA")
                allRestoredFocalRaster = os.path.join(outputGDB,'allRestFocRas_iter'+str(iter))
                outFocalStats.save(allRestoredFocalRaster) # Anything > 0 would include a restored area and 
                arcpy.env.extent = resRast
                restMaskRaster = os.path.join(outputGDB,'restMaskRaster_iter'+str(iter))
                minval = 0
                outCon = arcpy.sa.Con((Raster(allRestoredFocalRaster) == float(minval)), 1)
                outCon.save(restMaskRaster)
                
            # Candidate areas have not been restored, have cost data, meet
            # minimum improvement score criteria, and have enough ag in them
            candidateBarrierRaster = os.path.join(outputGDB, 'candidateBarrierRaster' + '_iter'+str(iter))
            if iter > 1:
                gprint('Creating candidate restoration raster using barrier results, previous restorations, and selection criteria')
                outCalc = (Raster(costMaskRaster) * Raster(agMaskRaster) * Raster(barrierRaster) * Raster(restMaskRaster) * (radius * 2)) # ROI scores will be in terms of total improvement (= score * diameter)
            else:
                outCalc = (Raster(costMaskRaster) * Raster(agMaskRaster) * Raster(barrierRaster) * radius * 2)
            
            minBarrierScore = minImprovementVal * radius * 2
            if restoredResistanceVal != 1:
                outCalc2 = (outCalc - (2 * radius * (restoredResistanceVal - 1)))
                outCon = arcpy.sa.Con((outCalc2 >= float(minBarrierScore)), outCalc2)
            else:
                outCon = arcpy.sa.Con((outCalc >= float(minBarrierScore)), outCalc)
            outCon.save(candidateBarrierRaster)
            lu.build_stats(candidateBarrierRaster)
            
            purchaseRoiRaster = os.path.join(outputGDB, 'purchaseRoiRaster' + '_iter'+str(iter))
            outCalc = Raster(candidateBarrierRaster) / Raster(purchCostRaster) 
            outCalc.save(purchaseRoiRaster)
            lu.build_stats(purchaseRoiRaster)
            
            totalRoiRaster = os.path.join(outputGDB, 'purchaseRestRoiRaster' + '_iter'+str(iter))
            outCalc = Raster(candidateBarrierRaster) / Raster(totalCostRaster)
            outCalc.save(totalRoiRaster)
            lu.build_stats(totalRoiRaster)

            maxBarrier = arcpy.GetRasterProperties_management(candidateBarrierRaster,"MAXIMUM")
            gprint('Maximum barrier improvement score: '+str(maxBarrier.getOutput(0)))
            if maxBarrier < 0:
                arcpy.AddWarning("\nNo barriers found that meet CWD or Ag threshold criteria.")
            
            maxPurchROI = arcpy.GetRasterProperties_management(purchaseRoiRaster,"MAXIMUM")
            gprint('Maximum purchase ROI score: '+str(maxPurchROI.getOutput(0)))

            maxROI = arcpy.GetRasterProperties_management(totalRoiRaster,"MAXIMUM")
            gprint('Maximum total ROI score: '+str(maxROI.getOutput(0)))

            if restoreMaxROI:
                outPoint = os.path.join(outputGDB, 'maxRoiPoint'+'_iter'+str(iter))
                gprint('Choosing circle with maximum ROI to restore')
                outCon = arcpy.sa.Con((Raster(totalRoiRaster) >= float(maxROI.getOutput(0))), totalRoiRaster)
                maxRoiRaster = os.path.join(outputGDB, 'maxRoiRaster')
                outCon.save(maxRoiRaster)    
                # Save max ROI to point
                try:
                    arcpy.RasterToPoint_conversion(maxRoiRaster, outPoint)
                except:
                    msg = ('Error: it looks like there are no viable restoration candidates.')
                    lu.raise_error(msg)
        
            else: #Restoring strongest barrier instead
                outPoint = os.path.join(outputGDB, 'maxBarrierPoint'+'_iter'+str(iter))
                gprint('Choosing circle with maximum BARRIER IMPROVEMENT SCORE to restore')
                outCon = arcpy.sa.Con((Raster(candidateBarrierRaster) >= float(maxBarrier.getOutput(0))), candidateBarrierRaster)
                maxBarrierRaster = os.path.join(outputGDB, 'maxBarrierRaster')
                outCon.save(maxBarrierRaster)            
                # Save max barrier to point
                try:
                    arcpy.RasterToPoint_conversion(maxBarrierRaster, outPoint)
                except:
                    msg = ('Error: it looks like there are no viable restoration candidates.')
                    lu.raise_error(msg)            
            
            gprint('Done evaluating candidate restorations')
            result = int(arcpy.GetCount_management(outPoint).getOutput(0)) 
            if result > 1:
                arcpy.AddWarning('Deleting points with identical ROI/improvement score values') # Would be better to retain point with max barrier score when we have multiple points with same ROI
                arcpy.DeleteIdentical_management(outPoint, "grid_code", 0.1, 0.1)            
            arcpy.sa.ExtractMultiValuesToPoints(outPoint, 
                [[candidateBarrierRaster, "barrierScore"],[purchCostRaster, "purchCost"],
                [totalCostRaster, "totalCost"],[purchaseRoiRaster, "purchaseROI"],
                [totalRoiRaster, "totalROI"]], "NONE")
            arcpy.AddField_management(outPoint, "restorationNumber", "SHORT")
            arcpy.CalculateField_management(outPoint, "restorationNumber", iter)        
            arcpy.AddField_management(outPoint, "radius", "DOUBLE")
            arcpy.CalculateField_management(outPoint, "radius", radius)        
            arcpy.AddField_management(outPoint, "barrierScore_per_m", "DOUBLE")
            arcpy.CalculateField_management(outPoint, "barrierScore_per_m", "(float(!barrierScore!) / (!radius! * 2))", "PYTHON")        

            gprint('\nCreating restoration circles')
            if restoreMaxROI:
                circleFC = os.path.join(outputGDB, 'maxRoiCircle'+'_iter'+str(iter))
            else:
                circleFC = os.path.join(outputGDB, 'maxBarrierCircle'+'_iter'+str(iter))
            arcpy.Buffer_analysis(outPoint, circleFC, radius)
            gprint('Rasterizing restoration circles')
            if restoreMaxROI:
                circleRas = os.path.join(outputGDB, 'maxRoiCircleRas'+'_iter'+str(iter))
            else:
                circleRas = os.path.join(outputGDB, 'maxBarrierCircleRas'+'_iter'+str(iter))
            arcpy.FeatureToRaster_conversion(circleFC, 'totalROI', circleRas, arcpy.env.cellSize)    

            # restore raster
            gprint('Digitally restoring resistance raster')
            resRastRestored = os.path.join(outputGDB, 'resRastRestored'+'_iter'+str(iter))
            outCon = arcpy.sa.Con(IsNull(circleRas), resRast, restoredResistanceVal)
            outCon.save(resRastRestored)

            allRestoredAreasRaster = os.path.join(outputGDB, 'allRestoredAreas_iter'+str(iter))
            PrevRestoredAreasRaster= os.path.join(outputGDB, 'allRestoredAreas_iter'+str(iter-1))
            if iter == 1:
                outCon = arcpy.sa.Con(IsNull(circleRas), 0, 1)
            else:
                outCon = arcpy.sa.Con(IsNull(circleRas), PrevRestoredAreasRaster, 1) # Add this restoration to areas restored
            outCon.save(allRestoredAreasRaster)
            
            lu.delete_data(circleRas)
            resRast = resRastRestored # Use for next iteration resistance raster
            
            #Add circle into feature class with all circles
            if restoreMaxROI:
                allCirclesFC = os.path.join(outputGDB,"allCirclesMaxROI")
            else:
                allCirclesFC = os.path.join(outputGDB,"allCirclesMaxBarriers")
            if iter == 1:
                arcpy.CopyFeatures_management(circleFC, allCirclesFC)
            else: 
                arcpy.Append_management(circleFC, allCirclesFC, "TEST")           
            gprint('Finished iteration #'+str(iter))
            start_time1 = lu.elapsed_time(start_time1)    

        gprint('\nDone with iterations.')
        start_time = lu.elapsed_time(start_time)    
        gprint('Outputs saved in: '+outputGDB)
        gprint('Back up your project directories if you want to save corridor/barrier results.')

    # Return GEOPROCESSING specific errors
    except arcpy.ExecuteError:
        lu.dashline(1)
        gprint('****Iteration script failed. Details follow.****')
        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Iteration script failed. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)        
if __name__ == "__main__":
    main()
    

