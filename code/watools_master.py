##*****************************************************************
## 2011_0128
## NAME: watools_master.py
##
## SUMMARY: Master script for linkage mapper.  Called by ArcMap with
## parameters or run from command line with parameters entered in
## script below.  Calls functions in dedicated scripts for each of 5
## processing steps.
##
## SOFTWARE: ArcGIS 9.3 (requires Spatial Analyst extension)
##           Python 2.5
##
##*****************************************************************

import arcgisscripting, sys, time
from time import localtime, strftime
import os
import string
import csv    
from string import split

Version = '2011.01.14'

try:
    from numpy import *
except:    
    gp.AddError("Numpy is not installed.  Please get the Python 2.5-compatible version from http://sourceforge.net/projects/numpy/files/NumPy/1.4.1/numpy-1.4.1-win32-superpack-python2.5.exe/download ") 
    sys.exit(1) # end script process

from watools_util import *

gp = arcgisscripting.create(9.3)
gp.OverwriteOutput = 1
try:
    gp.CheckOutExtension("Spatial")
except RuntimeError:
    gp.AddError("You need to have a Spatial Analyst extension for this tool")
    sys.exit(1) # end script process

try:    
    #########################################################################
    #Users can manually set some options here
    ######################
    saveNormLccs = False # Set to True to save individual normalized LCC grids
    # Set minimum corridor lengths below
    useMinEucDist=False
    useMinCostDist=False        
    minEucDist="" # use "" for no minimum
    minCostDist="" # use "" for no minimum
    #########################################################################    
    if minCostDist:
        minCostDist = float(minCostDist)
    if minEucDist:
        minEucDist = float(minEucDist)

    defaultNullValue = "#"
    
    #####################################
    # Set input parameters or take from Toolbox
    narg = len(sys.argv)
    
    # If run from command line. Mainly for testing and debugging:
    if narg == 1: 
        step1 = False
        step2 = False
        step3 = False
        step4 = False
        step5 = True
        projectDir = "C:\\WATOOLS\\Demo\\project"
        logDir = projectDir + "\\log"
        coreShapefile = "C:\\WATOOLS\\Demo\\Data\\Cores.shp" 
        resistanceRas = "C:\\WATOOLS\\Demo\\Data\\Resistances"
        coreIds = "HCA_ID"             
        dropLccsWithIntermediateCores = True
        inputAdjacencyMethod = "Both"
        useBoundingFeatures = True
        bufferDist=10000

        #useMaxDist
        useMaxEucDist=True
        maxEucDist=50000
        useMaxCostDist=True
        maxCostDist=1000000
        
        connectNearestOnly = True
        connectComponents=True
        maxNumNearestNeighbors=1
        nearestNbrDistType = "CostWeightedDistance"
    
    else:  # If called from Arc Toolbox
        argList = {}
        for x in range(narg-1):
            argList[x]=gp.GetParameterAsText(x)
            if argList[x]=='true':
                argList[x]=True
            elif argList[x]=='false':
                argList[x]=False
   
        step1 = argList[0]
        step2 = argList[1]
        step3 = argList[2]
        step4 = argList[3]  
        step5 = argList[4]

        projectDir = argList[5]
        check_project_dir(projectDir)
        coreShapefile = argList[6]
        coreIds = argList[7]
        resistanceRas = argList[8]
        
        inputAdjacencyMethod = argList[9]
        
        useBoundingFeatures = argList[10]
        if useBoundingFeatures == True:
            bufferDist = argList[11]
            
        useMaxDist = argList[12]
        if useMaxDist:
            maxEucDist = argList[13]
            if maxEucDist: 
                maxEucDist = float(maxEucDist)
                useMaxEucDist = True
            else:
                useMaxEucDist = False
                maxEucDist=""
            
            maxCostDist = argList[14]
            if maxCostDist:
                useMaxCostDist = True
                maxCostDist = float(maxCostDist)
            else:
                useMaxCostDist = False
                maxCostDist=""
        else:
            useMaxEucDist = False
            useMaxCostDist = False
            maxCostDist=""
            maxEucDist=""
        
        ######################################################
        ######################
        ## Removing minimum corridor lengths from GUI for now.  
        ## Leaving in code in case there is demand.
        ## Users can manually set minimums AT TOP of this script. 
        ######################
        #useMinDist = argList[18]   
        #if useMinDist:
        #    minEucDist = argList[19]
        #    if minEucDist:
        #        minEucDist = float(minEucDist)
        #        useMinEucDist=True
        #    else:
        #        useMinEucDist=False
        #        minEucDist=""
        #    
        #    minCostDist = argList[20]
        #    if minCostDist:
        #        minCostDist = float(minCostDist)
        #        useMinCostDist=True
        #    else:
        #        useMinCostDist=False
        #        minCostDist=""
        #else:
        #    useMinEucDist=False
        #    useMinCostDist=False
        #    minEucDist=""
        #    minCostDist=""
        ######################
        
        dropLccsWithIntermediateCores = argList[15] 

        maxNumNearestNeighbors = argList[16]
        if maxNumNearestNeighbors == 'Connect all adjacent cores':
            connectNearestOnly = False
            step4 = False # Fixme- This is temporary- will eventually update linktable with nearest neighbor info regardless.  For now, skipping step 4.
            gp.addmessage('Skipping step 4. All valid links between adjacent cores will be retained.')
        else:
            maxNumNearestNeighbors = int(maxNumNearestNeighbors)            
            connectNearestOnly = True
        
        nearestNbrDistType = argList[17]        
        connectComponents = argList[18]
        if connectComponents=="":
            connectComponents=False

    logFile="watools_start.log"
    write_start_log(Version,projectDir,logFile, argList, step=0)
    
    if set("# ").intersection(projectDir):
        gp.AddError("\nSpaces are not allowed in project directory path.  Please avoid using any special characters in directory names.")
        sys.exit(1) # end script process
    
    outputRootDir = projectDir + "\\output\\"
    if os.path.exists(outputRootDir)==False:
        gp.CreateFolder_management(projectDir,"output")
    
    if useBoundingFeatures == True:
        if bufferDist:
            bufferDist = float(bufferDist)
        else:
            bufferDist = defaultNullValue
    else:
        bufferDist = defaultNullValue

    # Check to make sure no missing steps, and identify first step checked
    check_steps(step1,step2,step3,step4,step5) 
    if step5 == True: firstStep = 5
    if step4 == True:
        firstStep = 4
        #move_stick_maps(projectDir,firstStep)
    if step3 == True:
        firstStep = 3
        #move_stick_maps(projectDir,firstStep)
    if step2 == True:
        check_dist_file(projectDir,coreShapefile)
        firstStep = 2
        #move_stick_maps(projectDir,firstStep)
    if step1 == True:
        firstStep = 1
        #move_stick_maps(projectDir,firstStep)
    try:
        test=firstStep
    except:
        gp.addmessage('\n---------------------------')
        msg = 'ERROR: Please check at least one step.'
        gp.AddError(msg)
        exit(1)       

    # Make a backup copy of datapass directory.  This will have lcp maps and link tables from previous run.
    archive_datapass(projectDir)   

    clean_up_link_tables(projectDir,firstStep)    

    
    if step5 == True:    
        outputGDBname = "linkages"      
        outputRootDir = projectDir + '\\output\\'
        outputGdb = outputRootDir + outputGDBname + ".gdb" 
        gp.workspace = outputRootDir
        if gp.Exists(outputGdb):
            gp.addmessage('Deleting geodatabase ' + outputGdb)
            try:
                gp.delete_management(outputGdb)
            except:
                gp.addmessage('\n---------------------------')
                msg = 'ERROR: Could not remove geodatabase ' + outputGdb + '. Is it open in ArcMap?'
                gp.AddError(msg)
                exit(1)
        #move_stick_maps(projectDir,firstStep)
    
    if useMaxCostDist == True:
        maxCwDist=maxCostDist + 100000 # this will limit cw calcs.  Value of 100,000 assumes corridors can be no more than 100km 'wide' in cost-distance units.
        if step1==True or step3==True:
            gp.addmessage('\n---------------------------------')
            gp.addmessage('Note: maximum cost-weighted distance processing set to ' + str(maxCwDist) + ', \nwhich is 100km farther than the minimum specified corridor length.')
            gp.addmessage('This option can be edited in watools_master.py.')
            gp.addmessage('---------------------------------\n')
    else:
        maxCwDist=defaultNullValue
        
    if inputAdjacencyMethod=="Both":
        keepLinkCriteria = "Keep all adjacent links"
    elif inputAdjacencyMethod == "EuclideanDistance":
        keepLinkCriteria = "Keep Euclidean adjacent links" 
    else:
        keepLinkCriteria = "Keep cost-weighted adjacent links"

    check_project_dir(projectDir)
    coreShapefile = str(coreShapefile)
    coreIds = str(coreIds)

    if resistanceRas == "":
        resistanceRas = defaultNullValue
    resistanceRas = str(resistanceRas)   
    if connectNearestOnly == False:
        maxNumNearestNeighbors = defaultNullValue 

    # Pack options into dictionary to pass to main functions
    options={}
    options["projectDir"]=projectDir
    options["inputAdjacencyMethod"]=inputAdjacencyMethod
    options["coreShapefile"]=coreShapefile
    options["coreIds"]=coreIds
    options["resistanceRas"]=resistanceRas
    options["bufferDist"]=bufferDist
    options["useMinEucDist"]=useMinEucDist
    options["useMaxEucDist"]=useMaxEucDist
    options["useMinCostDist"]=useMinCostDist
    options["useMaxCostDist"]=useMaxCostDist
    options["maxEucDist"]=maxEucDist
    options["minEucDist"]=minEucDist
    options["maxCwDist"]=maxCwDist
    options["maxCostDist"]=maxCostDist
    options["minCostDist"]=minCostDist
    options["keepLinkCriteria"]=keepLinkCriteria
    options["dropLccsWithIntermediateCores"]=dropLccsWithIntermediateCores
    options["connectNearestOnly"]=connectNearestOnly
    options["connectComponents"]=connectComponents
    options["maxNumNearestNeighbors"]=maxNumNearestNeighbors
    options["nearestNbrDistType"]=nearestNbrDistType
    options["saveNormLccs"]=saveNormLccs
    options["step1"]=step1
    options["step2"]=step2
    options["step3"]=step3
    options["step4"]=step4
    options["step5"]=step5
    options["defaultNullValue"]=defaultNullValue
    
    # Step 1    
    if step1==True:
        gp.addmessage('\n---------------------------------')
        gp.addmessage('Running script s1_getAdjacencies.py')
        from s1_getAdjacencies import *       
        logFile="step1_start.log"
        write_start_log(Version,projectDir,logFile, options, step=1)
        step1_get_adjacencies(gp,Version,options)  #  Run step 1
        logFile="step1.log"
        write_log(Version,projectDir,logFile, options, step=1)

    # Step 2
    if step2 == True:
        gp.addmessage('\n---------------------------------')
        gp.addmessage('Running script s2_buildNetwork.py')
        from s2_buildNetwork import *
        logFile="step2_start.log"
        write_start_log(Version, projectDir,logFile, options, step=2)
        step2_build_network(gp,Version,options)  # Run step 2
        logFile="step2.log"
        write_log(Version,projectDir,logFile, options, step=2)
    
    # Step 3
    if step3 == True:
        gp.addmessage('\n---------------------------------')        
        gp.addmessage('Running script s3_calcCwds.py')
        from s3_calcCwds import *       
        logFile="step3_start.log"
        write_start_log(Version,projectDir,logFile, options, step=3)
        step3_calc_cwds(gp,Version,options)  # Run step 3
        logFile = "step3.log"
        write_log(Version,projectDir,logFile, options, step=3)                

    # Step 4    
    if step4 == True:
        gp.addmessage('\n---------------------------------')        
        gp.addmessage('Running script s4_refineNetwork.py')
        from s4_refineNetwork import *
        logFile="step4_start.log"
        write_start_log(Version,projectDir,logFile, options, step=4)
        step4_refine_network(gp,Version,options)  # Run step 4
        logFile="step4.log"
        write_log(Version,projectDir,logFile, options, step=4)
    
    # Step 5    
    if step5 == True:
        gp.addmessage('\n---------------------------------')        
        gp.addmessage('Running script s5_calcLccs.py')
        from s5_calcLccs import *
        logFile="step5_start.log"
        write_start_log(Version,projectDir,logFile, options, step=5)
        step5_calc_lccs(gp,Version,options) # Run step 5
        logFile="step5.log"
        write_log(Version,projectDir,logFile, options, step=5)
    
    
    gp.addmessage('\nDONE!\n')
    
# Return GEOPROCESSING specific errors
except arcgisscripting.ExecuteError:
    filename =  'watools_master.py'
    raise_geoproc_error(filename)

# Return any PYTHON or system specific errors
except:
    filename =  'watools_master.py'
    raise_python_error(filename)
            
del gp

                
