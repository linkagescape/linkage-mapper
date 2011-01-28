##*****************************************************************
## 2011_0128
## NAME: s1_getAdjacencies.py
##
## SUMMARY: Determines adjacencies between core areas in either or both
## Euclidean and cost-weighted distance space
##
## SOFTWARE: ArcGIS 9.3 (requires Spatial Analyst extension)
##           Python 2.5
##
##*****************************************************************

import arcgisscripting, sys, time
from time import localtime, strftime
import os, string,csv
from numpy import *
from string import split
from numpy import loadtxt, where, delete, arange

from watools_util import *

def step1_get_adjacencies(gp,Version,options):
    """Determines adjacencies between core areas in either or both
    Euclidean and cost-weighted distance space.
    
    """
    try:
        # ------------------------------------------------------------------                
        # Unpack options
        projectDir=options["projectDir"]
        inputAdjacencyMethod=options["inputAdjacencyMethod"]
        coreShapefile=options["coreShapefile"]
        coreIds=options["coreIds"]
        resistanceRas=options["resistanceRas"]
        bufferDist=options["bufferDist"]
        maxCwDist=options["maxCwDist"]
        defaultNullValue = options["defaultNullValue"]
        # ------------------------------------------------------------------            
        
        cellSizeEuclidean=defaultNullValue # Default behavior is to use same cell size as resistance raster, but this can be changed here.
       
        # ------------------------------------------------------------------                    
        # Folders and paths and workspace
        ouputFolder = "adj"            
        dataPassDir = projectDir + "\\dataPass\\"
        if os.path.exists(dataPassDir)==False:
            gp.CreateFolder_management(projectDir,"dataPass")    
        outputRootDir = projectDir + "\\output\\"
        if os.path.exists(outputRootDir)==False:
            gp.CreateFolder_management(projectDir,"output")    
        scratchDir = projectDir + "\\scratch"
        if os.path.exists(scratchDir)==False:
            gp.CreateFolder_management(projectDir,"scratch")
        gp.scratchWorkspace = scratchDir 
        gp.addmessage('Creating output folder: ' + ouputFolder)
        gp.workspace = projectDir        
        gp.AddMessage('Adjacency files will be written to ' + gp.workspace + "\\" + ouputFolder)
        if os.path.exists(gp.workspace + "\\" + ouputFolder)==True:
            import shutil
            shutil.rmtree(gp.workspace + "\\" + ouputFolder)
        gp.CreateFolder_management(gp.workspace, ouputFolder)    
        
        # Set output dataset names
        core_ras = "core_ras"
        focalvariety_ras = "focalvar_ras"
        expand_ras = "expand_ras"
        expand2_ras = "expand2_ras"
        edges_ras = "edges_ras"
        combine_ras = "combine_ras"
        # ------------------------------------------------------------------            
        
        if (resistanceRas != defaultNullValue):
            gp.SnapRaster = resistanceRas
        SR = gp.describe(coreShapefile).SpatialReference # To set spatial reference for shapefiles we create later
        
        # Euclidean cell size 
        if inputAdjacencyMethod != "CostWeightedDistance" and cellSizeEuclidean == defaultNullValue:
            if (resistanceRas != defaultNullValue):
                cellSizeEuclidean = gp.Describe(resistanceRas).MeanCellHeight 
                gp.addmessage('Euclidean cell size set equal to resistance raster cell size (' + str(cellSizeEuclidean) + ').')
            else:
                msg ="Error: Cell size required for Euclidean adjacency processing\n"
                gp.AddMessage(" ")
                gp.AddError(msg)
                exit(0)    
        
        # Check that resistance raster exists (if it is needed)
        if inputAdjacencyMethod !=  "EuclideanDistance" and not gp.exists(resistanceRas):
            msg = "Error: Resistance raster required for cost-weighted adjacency processing\n"
            gp.AddMessage(" ")
            gp.AddError(msg)
            exit(0)
        
        # ------------------------------------------------------------------                    
        # Create bounding circles to limit cwd and allocation calculations
        if bufferDist != defaultNullValue:
            gp.addmessage('Reducing processing area using bounding circle plus buffer of ' + str(float(bufferDist)/1000) + ' km.')
            
            startTime = time.clock()
            gp.MakeFeatureLayer(coreShapefile,"fcores")
    
            extentBoxList = zeros((0,5),dtype='float32')
            boxCoords = get_extent_box_coords(gp.workspace,"fcores",coreIds,"#")
            extentBoxList=append(extentBoxList,boxCoords,axis=0)
            extentBoxList[0,0]=0
    
            # cwd bounding circle- used to clip raster to limit cwd calculations
            boundingCirclePointArray  = zeros((0,5),dtype='float32')
            circlePointData=get_bounding_circle_data(extentBoxList,0,0,bufferDist)   
            boundingCircleCenter = scratchDir + "\\boundingCircleCenter.shp"
            make_points(scratchDir,circlePointData,"boundingCircleCenter.shp",coreIds)
            gp.defineprojection(boundingCircleCenter, SR)
            boundingCircle=scratchDir + "\\boundingCircle.shp"                                  
            if gp.Exists(boundingCircle):
                gp.delete_management(boundingCircle)
            gp.buffer_analysis(boundingCircleCenter, boundingCircle, "radius")
            gp.defineprojection(boundingCircle, SR)
            # boundingCirclePointArray  = zeros((0,5),dtype='float32')
            # circlePointData=get_bounding_circle_data(extentBoxList,0,0,bufferDist)
    
            # euc bounding circle- at the moment just limits extent of euclidean allocation calculations
            boundingCirclePointArray  = zeros((0,5),dtype='float32')
            circlePointData=get_bounding_circle_data(extentBoxList,0,0,0) # no buffer needed
            eucBoundingCircleCenter = scratchDir + "\\eucBoundingCircleCenter.shp"
            make_points(scratchDir,circlePointData,"eucBoundingCircleCenter.shp",coreIds)
            gp.defineprojection(eucBoundingCircleCenter, SR)
            eucBoundingCircle=scratchDir + "\\eucBoundingCircle.shp"                                  
            if gp.Exists(eucBoundingCircle):
                gp.delete_management(eucBoundingCircle)
            gp.buffer_analysis(eucBoundingCircleCenter, eucBoundingCircle, "radius")
            gp.defineprojection(eucBoundingCircle, SR)        

            del boundingCirclePointArray

        # ------------------------------------------------------------------                                
        # Loop through code twice if user wants to calculate BOTH Euclidean and cost-weighted adjacency.
        loop=0    
        if inputAdjacencyMethod=="Both":
            numLoops=2
        else:
            numLoops=1
            optAdjacencyMethod=inputAdjacencyMethod # Just doing one or the other
            
        while loop < numLoops:
            if inputAdjacencyMethod=="Both" and loop==0:
                optAdjacencyMethod="EuclideanDistance"
            elif inputAdjacencyMethod=="Both" and loop==1:
                optAdjacencyMethod="CostWeightedDistance"
            
            gp.OverwriteOutput = 1           
            gp.pyramid = "NONE"
            gp.rasterstatistics = "NONE"
            
            if optAdjacencyMethod == "EuclideanDistance":
                alloc_rasFN = "Euc_alloc_ras"        
            else:    
                alloc_rasFN = "Cwd_alloc_ras"        
            
            if optAdjacencyMethod ==  "CostWeightedDistance":
                gp.addmessage('\nCalculating cost-weighted distance adjacency')
                outcsvfile = dataPassDir + "cwdAdj.csv"
                outcsvLogfile = projectDir + "\\log\\cwdAdj_step1.csv"
            else:
                gp.addmessage('\nCalculating Euclidean adjacency')            
                outcsvfile = dataPassDir  + "eucAdj.csv"
                outcsvLogfile = projectDir + "\\log\\eucAdj_step1.csv"
                
            # --------------------------------------------
            
            gp.workspace = scratchDir 
            
            if optAdjacencyMethod == "CostWeightedDistance": # May need to set extent prior to core poly to raster conversion...
                # ----------------------------------------------
                # Cost-weighted allocation code
                if bufferDist != defaultNullValue:               
                    # Clip resistance raster using bounding circle
                    startTime = time.clock()
                    bResistance=scratchDir + "\\bResistance"
                    gp.ExtractByMask_sa(resistanceRas, boundingCircle, bResistance)
                    gp.addmessage('\nReduced resistance raster extracted using bounding circle.')
                    startTime,hours,mins,secs = elapsed_time(startTime)
                    
                else:
                    bResistance=resistanceRas
                
                startTime = time.clock()
                gp.addmessage('Starting cost weighted distance allocation...')
                
                if maxCwDist!=defaultNullValue:
                    gp.addmessage('Maximum cost-weighted distance set to ' + str(maxCwDist))
                gp.CellSize = gp.Describe(bResistance).MeanCellHeight 
                gp.extent = "MAXOF"
                gp.AddMessage('Processing cell size: ' + gp.CellSize)
                count,statement = 0, 'gp.FeatureToRaster_conversion(coreShapefile,coreIds,core_ras,gp.Cellsize)'
                while True:
                    try: exec statement
                    except:
                        count,tryAgain = hiccup_test(count,statement)
                        if not tryAgain: exec statement
                    else: break
                
                gp.workspace = projectDir + "\\" + ouputFolder        
                gp.scratchworkspace = gp.workspace
                outDistanceRaster = projectDir + "\\output\\cwd" #fixme: put this in geodatabase instead
                alloc_ras = projectDir + "\\" + ouputFolder + "\\" + alloc_rasFN
                core_ras_path = scratchDir + "\\" + core_ras
                count, statement = 0, 'gp.Costallocation_sa(core_ras_path,bResistance,alloc_ras, maxCwDist, core_ras_path, "VALUE",outDistanceRaster,"")'
                while True:
                    try: exec statement
                    except:
                        count,tryAgain = hiccup_test(count,statement)
                        if not tryAgain: exec statement
                    else: break                    
                    
                gp.scratchworkspace = scratchDir
                gp.addmessage('\nCost-weighted distance allocation done.')
                startTime,hours,mins,secs = elapsed_time(startTime)
                # ----------------------------------------------

                
            elif optAdjacencyMethod == "EuclideanDistance": # FIXME- would be good to have bounding circle affect euclidean calcs too
                # ----------------------------------------------
                # Euclidean allocation code
                gp.workspace = projectDir + "\\" + ouputFolder            
                gp.addmessage ('Starting Euclidean adjacency processing...')
                gp.AddMessage('Processing cell size: ',gp.CellSize)
                oldextent=gp.extent
                if bufferDist != defaultNullValue:    
                    gp.extent = gp.Describe(eucBoundingCircle).extent
                count,statement = 0, 'gp.FeatureToRaster_conversion(coreShapefile,coreIds,core_ras,cellSizeEuclidean)'
                while True:
                    try: exec statement
                    except:
                        count,tryAgain = hiccup_test(count,statement)
                        if not tryAgain: exec statement
                    else: break
                
                startTime=time.clock()
    
                gp.scratchworkspace = gp.workspace
                outDistanceRaster = projectDir + "\\" + ouputFolder + "\\euc"
                alloc_ras = projectDir + "\\" + ouputFolder + "\\" + alloc_rasFN
                count,statement = 0, 'gp.EucAllocation_sa(core_ras,alloc_ras,"","",cellSizeEuclidean,"",outDistanceRaster,"")'
                while True:
                    try: exec statement
                    except:
                        count,tryAgain = hiccup_test(count,statement)
                        if not tryAgain: exec statement
                    else: break

                gp.scratchworkspace = scratchDir
                gp.addmessage('\nEuclidean distance allocation done.')
                startTime,hours,mins,secs = elapsed_time(startTime)                   
                gp.extent=oldextent
                # ----------------------------------------------

            else:
                msg = "Invalid adjacency distance option. Must be Cost Distance Weighted, Euclidean, or Both"
                print msg
                gp.AddMessage(msg)
                exit(1)
                
    
            # ------------------------------------------------------------------            
            # Get adjacencies using shift method and write to disk.
            gp.addmessage('\nProcessing rasters...') 
            adjTable = get_adj_using_shift_method(scratchDir, alloc_ras) #to be replaced by getLeastCostDistsUsingShiftMethod if implemented
            write_adj_file(outcsvfile,coreIds,adjTable)
            write_adj_file(outcsvLogfile,coreIds,adjTable)                        
            # ------------------------------------------------------------------            
    
            loop=loop+1 # go to 2nd loop or quit
    
    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
       gp.addmessage("\n--------------------")
       gp.addmessage('****Failed in step 1. Details follow.****')        

       filename =  __file__
       raise_geoproc_error(filename)
       
    # Return any PYTHON or system specific errors
    except:
        gp.addmessage("\n--------------------")
        gp.addmessage('****Failed in step 1. Details follow.****')        

        filename =  __file__
        raise_python_error(filename)
    
    del gp
    return