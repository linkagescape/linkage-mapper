


######################################################################   
      ###Currently Unused Functions- some are in development
      ##  KEEP  ##
######################################################################

#settings history code- disabled for now.

From write_link_table:
#Disabled settings history for now.        
        #if settingsHistory != {}:
        #    outFile.write ("# Settings History:\n")
        #    for item in settingsHistory:
        #        outFile.write ("# " + str(item) + "," + str(settingsHistory[item])+ "\n")
        #    outFile.write ("# End History\n")
        #    outFile.write ("#\n")
        
def read_settings_history(filename):
    """Reads header in link table (if exists) specifying history of settings in previous steps
    that generated that link table.  Currently disabled but will likely be re-enabled.
    
    """
    try:
        if os.path.isfile(filename)==False:
            raise RuntimeError('File "'  + filename + '" does not exist')
        else:
            headerLen = 0
            reader = csv.reader(open(filename))
            for row in reader:
                if (headerLen == 0) and (row[0] != '# Settings History:'):
                        defaultSettingsHistory = get_default_settings_history()
                        return defaultSettingsHistory               
                headerLen = headerLen+1
                if (row[0] == "# End History") or row[0] == "#link":
                    break
                if (row[0] != '# Settings History:'):                
                    setting = row[0]
                    setting=string.replace(setting,'# ','')
                    settingsHistory[setting] = row[1]
            if headerLen > 4:
                return settingsHistory
            else:
                defaultSettingsHistory = get_default_settings_history()
                return defaultSettingsHistory
    except:
        raise_python_error('watools_util')

def get_default_settings_history():
    """Returns a default placeholder settings history in case none available from link table"""
    try:
        gp.Addmessage('\nRetrieving default settings history.\n')
        settingsHistory = {}
        settingsHistory["projectDir"] = -1
        settingsHistory["keepLinkCriteria"] = -1
        settingsHistory["maxEucCorridorLength"] = -1
        settingsHistory["maxCwdCorridorLength"] = -1          
        #maxCwdForAdjacency
        #maxCwdForCostDistance
        #
        #coreShapefile = "C:\\litest\\project\\lica1050_150_eco.shp" 
        #resistanceRas = "C:\\litest\\Project\\resist_med"
        #coreIds = "HCA_ID"
        #      
        #dropLccsWithIntermediateCores = True
        #
        #inputAdjacencyMethod = "Both"
        #
        #adjBufferDist
        #cwdBufferDist
        #
        #connectNearestOnly = True
        #connectComponents=True
        #maxNumNearestNeighbors=1
        #nearestNbrDistType = "CostWeightedDistance"
        
        return settingsHistory
    except:
        raise_python_error('watools_util')        

        
        
        
# For removal?  Backing up stick maps isn't a priority.  Lcp maps are backed up in datapass archive.
def move_stick_maps(projectDir,step):
    """Move stick maps from previous runs to a backup directory."""
    try:
        outputDir = projectDir + '\\output'
        oldLinkMapsDir = outputDir + '\\old_results'
        if os.path.exists(oldLinkMapsDir) == False:
            gp.CreateFolder_management(outputDir,"old_results")
        coreLinksShapefile= outputDir + '\\Active_Sticks.shp' 
        if gp.exists(coreLinksShapefile):
            backupGroupLinksShapefile = outputDir + '\\old_results\\Active_Sticks.shp'
            gp.CopyFeatures_management (coreLinksShapefile, backupGroupLinksShapefile) 
            while 1:
                try:
                    gp.delete_management(coreLinksShapefile)
                    break
                except:
                    gp.addmessage('\n---------------------------')
                    msg = 'ERROR: Could not remove final sticks shapefile ' + coreLinksShapefile + '. Is it open in ArcMap?'
                    gp.AddError(msg)
                    exit(1)

        coreLinksShapefile=outputDir + '\\sticks_step' + str(step) + '.shp'
        if gp.exists(coreLinksShapefile):
            backupGroupLinksShapefile = outputDir + '\\old_results\\sticks_step' + str(step) + '.shp'
            gp.CopyFeatures_management (coreLinksShapefile, backupGroupLinksShapefile) 
            try:
                gp.delete_management(coreLinksShapefile)
            except:
                gp.addmessage('\n---------------------------')
                msg = 'ERROR: Could not remove links shapefile ' + coreLinksShapefile + '. Is it open in ArcMap?'
                gp.AddError(msg)
                exit(1)
        
    except arcgisscripting.ExecuteError:
        raise_geoproc_error('watools_util')
    except:
        raise_python_error('watools_util') 


def getlcDistTable(workspace,dbfFile):
    try:
        lcDistTable = zeros((0,2),dtype="float32")
        appendRow = zeros((1,2),dtype="float32")
        
        rows = gp.searchcursor(dbfFile)
        row = rows.Next()
        while row:
            appendRow[0,0] = row.Value
            appendRow[0,1] = row.Min
            lcDistTable = append(lcDistTable, appendRow, axis=0)
            row = rows.next()
        del row
        del rows
    
        return lcDistTable
    except arcgisscripting.ExecuteError:
        raise_geoproc_error('watools_util')
    except:
        raise_python_error('watools_util')   


def joinTables(lcDistTable,allocLookupTable):
    print 'lcd'
    print lcDistTable
    print 'allocLookupTable'
    print allocLookupTable
    try:
        lcDistJoined = zeros((len(lcDistTable),4),dtype="float32")
        lcDistJoined[:,0] = lcDistTable[:,0]
        lcDistJoined[:,1] = allocLookupTable[:,1]
        lcDistJoined[:,2] = allocLookupTable[:,2]    
        lcDistJoined[:,3] = lcDistTable[:,1]
        return lcDistJoined    
    except arcgisscripting.ExecuteError:
        raise_geoproc_error('watools_util')
    except:
        raise_python_error('watools_util')   


def getLeastCostDistancesFromShift(workspace,alloc,alloc_sh,lcdist_sh):
    try:
        combine_ras = gp.workspace + "\\combine"
        count,statement = 0, 'gp.SingleOutputMapAlgebra_sa("combine("+alloc+","+ alloc_sh+")", combine_ras, alloc, alloc_sh)'
        while True:
            try: exec statement
            except:
                count,tryAgain = hiccup_test(count,statement)
                if not tryAgain: exec statement
            else: break
        
        allocLookupTable = get_alloc_lookup_table(gp.workspace,combine_ras)
        
        dbfFile = workspace + "\\" + "zonestats.dbf"
        count,statement = 0, 'gp.zonalstatisticsastable_sa(combine_ras,"Value",lcdist_sh,dbfFile,"DATA")' # Try this in script 4 too... and anything else with minras (7)...
        while True:
            try: exec statement
            except:
                count,tryAgain = hiccup_test(count,statement)
                if not tryAgain: exec statement
            else: break
        
        lcDistTable1 = getlcDistTable(gp.workspace,dbfFile)
        lcDistTable = joinTables(lcDistTable1,allocLookupTable)
        return lcDistTable,allocLookupTable[:,1:3]
    except arcgisscripting.ExecuteError:
        raise_geoproc_error('watools_util')
    except:
        raise_python_error('watools_util')   


def writeConfig(dir, settings):
    try:
        if os.path.exists(dir + "\\log")==False:
            gp.CreateFolder_management(dir,"log")
        configFile=dir + "\\log\\config.ini"
        outFile = open(configFile,"w")
        for item in settings:
            outFile.write(str(item) + "," + str(settings[item]) + "\n")
        outFile.close()    
    except:
        raise_python_error('watools_util')
    return


def loadConfig(dir):
    try:
        configFile=dir + "\\log\\config.ini"
        if os.path.isfile(configFile)==False:
            raise RuntimeError('File "'  + configFile + '" does not exist')
        loadedSettings={}
        reader = csv.reader(open(configFile))
        for setting, value in reader:
            loadedSettings[setting]=value
    except:
        raise_python_error('watools_util')
    return loadedSettings


def zipdir(basedir, archivename): #Fixme: check for licensing of this recipe.
    from zipfile import ZipFile, ZIP_DEFLATED

#usage:
    #basedir = 'C:\WATEST\PROJECT1\cwd'
    #archivename = basedir + ".zip" # archive in the basedir   
    #zipdir(basedir, archivename)    

    assert os.path.isdir(basedir)
    with closing(ZipFile(archivename, "w", ZIP_DEFLATED)) as z:
        for root, dirs, files in os.walk(basedir):
            #NOTE: ignore empty directories
            for fn in files:
                absfn = os.path.join(root, fn)
                zfn = absfn[len(basedir)+len(os.sep):] #XXX: relative path
                z.write(absfn, zfn)


def getLeastCostDistsUsingShiftMethod(scratchDir, res, cwd, alloc): #in progress, NOT implemented
    #Least-cost distances calculated using this method can be biased upward
    #if least-cost path goes through an intermediate allocation zone.

    cellSize = gp.Describe(res).MeanCellHeight
    gp.CellSize = cellSize
    
    print 'Cell height is ',str(gp.Cellsize)
    posShift = gp.CellSize
    negShift = -1*float(gp.CellSize)
    
    gp.workspace = scratchDir
    
    #Least-cost distances calculated using this method can be biased upward
    #if least-cost path goes through an intermediate allocation zone.
    gp.addmessage('Calculating least-cost distances crossing horizontal allocation boundaries...')
    startTime=time.clock()
    gp.Shift_management(alloc, "alloc_r", posShift, "0")
    gp.Shift_management(cwd, "cwd_r", posShift, "0")
    gp.Shift_management(res, "res_r", posShift, "0")
    
    expression1 = "(con(" + alloc + " - alloc_r <> 0, 1))"
    expression = expression1 + " * ( " + cwd + " + cwd_r + (( " + res + " + res_r ) * " + str(cellSize) + " / 2 ))"
    count,statement = 0, 'gp.SingleOutputMapAlgebra_sa(expression, "lcdist_r")'
    while True:
        try: exec statement
        except:
            count,tryAgain = hiccup_test(count,statement)
            if not tryAgain: exec statement
        else: break
    
    
    alloc_r = "alloc_r"
    lcDistTable_r,adjTable_r = getLeastCostDistancesFromShift(gp.workspace,alloc,alloc_r,"lcdist_r")
    startTime,hours,mins,secs = elapsedTime(startTime)
    
    gp.addmessage('Calculating least-cost distances crossing upper-left diagonal allocation boundaries...')
    gp.Shift_management(alloc, "alloc_ul", negShift, posShift)
    gp.Shift_management(cwd, "cwd_ul", negShift, posShift)
    gp.Shift_management(res, "res_ul", negShift, posShift)
    
    expression1 = "(con(" + alloc + " - alloc_ul <> 0, 1))"
    expression = expression1 + " * ( " + cwd + " + cwd_ul + (( " + res + " + res_ul ) * " + str(cellSize) + " * Sqrt(2) / 2 ))"
    count,statement = 0, 'gp.SingleOutputMapAlgebra_sa(expression, "lcdist_ul")'
    while True:
        try: exec statement
        except:
            count,tryAgain = hiccup_test(count,statement)
            if not tryAgain: exec statement
        else: break
    
    
    alloc_ul = "alloc_ul"
    lcDistTable_ul,adjTable_ul = getLeastCostDistancesFromShift(gp.workspace,alloc,alloc_ul,"lcdist_ul")
    startTime,hours,mins,secs = elapsedTime(startTime)
    
    
    gp.addmessage('Calculating least-cost distances crossing upper-right diagonal allocation boundaries...')
    gp.Shift_management(alloc, "alloc_ur", posShift, posShift)
    gp.Shift_management(cwd, "cwd_ur", posShift, posShift)
    gp.Shift_management(res, "res_ur", posShift, posShift)
    
    expression1 = "(con(" + alloc + " - alloc_ur <> 0, 1))"
    expression = expression1 + " * ( " + cwd + " + cwd_ur + (( " + res + " + res_ur ) * " + str(cellSize) + " * Sqrt(2) / 2 ))"
    count,statement = 0, 'gp.SingleOutputMapAlgebra_sa(expression, "lcdist_ur")'
    while True:
        try: exec statement
        except:
            count,tryAgain = hiccup_test(count,statement)
            if not tryAgain: exec statement
        else: break
    
    
    alloc_ur = "alloc_ur"
    lcDistTable_ur,adjTable_ur = getLeastCostDistancesFromShift(gp.workspace,alloc,alloc_ur,"lcdist_ur")
    startTime,hours,mins,secs = elapsedTime(startTime)
    
    gp.addmessage('Calculating least-cost distances crossing vertical allocation boundaries...')
    gp.Shift_management(alloc, "alloc_u", "0", posShift)
    gp.Shift_management(cwd, "cwd_u", "0", posShift)
    gp.Shift_management(res, "res_u", "0", posShift)
    
    expression1 = "(con(" + alloc + " - alloc_u <> 0, 1))"
    expression = expression1 + " * ( " + cwd + " + cwd_u + (( " + res + " + res_u ) * " + str(cellSize) + " / 2 ))"
    count,statement = 0, 'gp.SingleOutputMapAlgebra_sa(expression, "lcdist_u")'
    while True:
        try: exec statement
        except:
            count,tryAgain = hiccup_test(count,statement)
            if not tryAgain: exec statement
        else: break
    
    
    alloc_u = "alloc_u"
    lcDistTable_u,adjTable_u = getLeastCostDistancesFromShift(gp.workspace,alloc,alloc_u,"lcdist_u")
    startTime,hours,mins,secs = elapsedTime(startTime)
    
    lcDistTable = combineLeastCostDistances(lcDistTable_r,lcDistTable_u,lcDistTable_ur,lcDistTable_ul)
    adjTable = combine_adjacency_tables(adjTable_r,adjTable_u,adjTable_ur,adjTable_ul)
    
    outcsvfile = scratchDir+'TEMP.CSV' #this will replace adj files...

    outfile = open(outcsvfile,"w")
    coreIds = 'CORE_ID' #fixme temp
    outfile.write ("#Edge" + "," + str(coreIds) + "," + str(coreIds) + "_1" + "\n")
    for x in range(0,len(adjTable)):
        outfile.write ( str(x) + "," + str(adjTable[x,0]) + "," + str(adjTable[x,1]) + "\n" )
    outfile.close()
    
    outcsvLogfile = scratchDir+'TEMP.CSV'
    outfile = open(outcsvLogfile,"w") #log file
    for x in range(0,len(adjTable)):
        outfile.write ( str(x) + "," + str(adjTable[x,0]) + "," + str(adjTable[x,1]) + "\n" )
    outfile.close()
        
    print 'final tables:'
    print lcDistTable
    print 'adj'
    print adjTable
    return adjTable,lcDistTable    


def combineLeastCostDistances(lcDistTable_r,lcDistTable_u,lcDistTable_ur,lcDistTable_ul):
    try:
        lcDistTable = append(lcDistTable_r,lcDistTable_u,axis=0)
        lcDistTable = append(lcDistTable,lcDistTable_ur,axis=0)
        lcDistTable = append(lcDistTable,lcDistTable_ul,axis=0)
        
        pairs = sort(lcDistTable[:,1:3])
        lcDistTable[:,1:3] = pairs
        
        ind=lexsort((lcDistTable[:,3],lcDistTable[:,2],lcDistTable[:,1])) #sort by 1st coreId then by 2nd coreId
        lcDistTable = lcDistTable[ind]    
        
        numDists=len(lcDistTable)
        x=1
        while x<numDists:
            if lcDistTable[x,1]==lcDistTable[x-1,1] and lcDistTable[x,2]==lcDistTable[x-1,2]:
                lcDistTable[x,0]=0 #mark for deletion
            x=x+1
        
        if numDists>0:
            delRows=asarray(where(lcDistTable[:,0]==0))
            delRowsVector = zeros((delRows.shape[1]),dtype="int32")
            delRowsVector[:] = delRows[0,:]
            lcDistTable=deleterow(lcDistTable,delRowsVector)        
            del delRows
            del delRowsVector       
        return lcDistTable[:,1:4]

    except arcgisscripting.ExecuteError:
        raise_geoproc_error('watools_util')
    except:
        raise_python_error('watools_util')   


def getLccGridPath(rootDir,baseDirName,gridName,maxNum):  #In progress
    try:
        returnedGridPath = "Null"
        dir=rootDir+"\\"+baseDirName #new version of watools directory structure
        if os.path.exists(dir):
            for num in range (maxNum+1):
                if num == 0:
                    dirCount=''
                else:
                    dirCount=str(num)
                dir=rootDir+"\\"+baseDirName+str(dirCount)#inputrootdir\\lcX
                gridPath = dir+"\\"+gridName#inputrootdir\\lcX\\rawLccgridName
                #if os.path.exists(dir)==True: #inputrootdir\\lcX
                if os.path.exists(gridPath)==True: #bhm 7/22/10
                    returnedGridPath = gridPath
                    break
            return returnedGridPath
        else: #return path created by old version (pre-7/20/2010) of WATOOLS
            fullName= os.path.split(rootDir)
            projectDir = fullName[0]
            gp.addmessage(str(projectDir))
            outputRootDir = projectDir + "\\output\\"
            rawlccGDBname = "rawlccs"
            lccPath=outputRootDir + rawlccGDBname + ".gdb\\"+gridName
            gp.addmessage(lccPath)
            return lccPath
    except:
        raise_python_error('watools_util')    


##########################################################################
###Currently defunct functions############################################
##########################################################################
def count_sub_dirs(rootDir):
    """Counts the number of subdirectories in a path."""
    if os.path.exists(rootDir):
        return len(os.walk(rootDir).next()[1])
    return 0

def copy_map(oldMap,newMap):
    """Copies a map to a new location """
    if gp.exists(oldMap):
        if gp.exists(newMap):
            try:
                gp.delete_management(newMap)
            except:
                pass
        try:
            gp.CopyFeatures_management (oldMap, newMap)
        except:
            pass
    return

def groupCores(coreList,groupCutoff,dists):
    try:
        groupCutoff = int(groupCutoff)
        justDists=(dists[:,2])
        distBelowCutoff = justDists !=-1
        distBelowCutoff[:] = where(justDists[:] > groupCutoff,False,distBelowCutoff)
        pairsToGroup=dists[distBelowCutoff,0:2]
        
        pairsToGroupCopy=pairsToGroup
        
        #Update group IDs for core area pairs that are combined
        numPairs = len(pairsToGroup)
        for x in range(numPairs):
            coreList[:,1]=where(coreList[:,1] == pairsToGroup[x,1],pairsToGroup[x,0],coreList[:,1])
            pairsToGroup=where(pairsToGroup == pairsToGroup[x,1],pairsToGroup[x,0],pairsToGroup)
    except arcgisscripting.ExecuteError:
        raise_geoproc_error('watools_util')
    except:
        raise_python_error('watools_util')
    return coreList,pairsToGroupCopy


def CreateValueList(raster):
    try:
        Values = []
        try:
            rows = gp.SearchCursor(raster)
        except:
            return 'Failed'
        row = rows.next()
        while row:
            Values.append(row.Value)
            row = rows.next()
        del row, rows
        Values.sort() 
        return Values 
    except arcgisscripting.ExecuteError:
        raise_geoproc_error('watools_util')
    except:
        raise_python_error('watools_util') 



def movePreV5LinkTable(projectDir,step):
    try:
        datapassDir = projectDir + '\\datapass'
        oldLinkTableDir = projectDir + '\\datapass\\old_linktables'
        if os.path.exists(oldLinkTableDir) == False:
            gp.CreateFolder_management(datapassDir,"old_linktables")
            
        if step == 3: oldPrevScripts = [3]            
        if step == 4: oldPrevScripts = [4,5]
        if step == 5: oldPrevScripts = [4,5,7]
        if step == 6: oldPrevScripts = [7,8] #8's not previous, but is defunct once step 6 is run.
        for script in oldPrevScripts:
            preV5LinkTableFile = datapassDir + '\\linkTable_script'+str(script)+'.csv'
            if os.path.exists(preV5LinkTableFile) == True:    
                newfilename = datapassDir + '\\old_linktables\\linkTable_script'+str(script)+'.csv'
                try:
                    shutil.move(preV5LinkTableFile , newfilename)
                except:
                    gp.addmessage('Unable to move old linktable file.')
                    pass
        
        preV5LinkTableFile = datapassDir + '\\linkTable.csv'
        if os.path.exists(preV5LinkTableFile) == True:
            newFilename = datapassDir + '\\old_linktables\\linkTable_preV5.csv' 
            os.move(preV5LinkTableFile,newFilename)            

    except arcgisscripting.ExecuteError:
        raise_geoproc_error('watools_util')
    except:
        raise_python_error('watools_util') 

   
   
def createExtentBox(workspace,inFC,outFC,field,fieldValue): #FIXME: Extent selected code above is more than 10x faster.
    try:
        gp.workspace=workspace

        if fieldValue == "#": #get all features, not just where field=fieldValue
            fieldValue = 1
            desc = gp.Describe
            extent=desc(inFC).extent
            lr = extent.lowerright
            ul = extent.upperleft
                
            ulx=str(ul.x)
            uly=str(ul.y)
            lrx=str(lr.x)
            lry=str(lr.y)
            
        else:        
            ulx,lry,lrx,uly =  newExtent(inFC,field,fieldValue)

        makeBox(workspace,inFC,outFC,field,fieldValue,ulx,lry,lrx,uly)    

    except arcgisscripting.ExecuteError:
        raise_geoproc_error('watools_util')
    except:
        raise_python_error('watools_util')
   
    return outFC
   

def makeBox(workspace,inFC,outFC,field,fieldValue,ulx,lry,lrx,uly):
    try:
        gp.workspace=workspace
        gp.OverwriteOutput = 1
        if gp.Exists(outFC):
            gp.delete_management(outFC)
        # Create an Array object.
        vertex_array = gp.createobject("Array")
        gp.CreateFeatureclass_Management(gp.workspace, outFC, "POLYGON", inFC) 
        #gp.MakeFeatureLayer(gp.workspace+outFC+".shp"....  #fixme: trying to do layer
         # clear the array ready for next time
        vertex_array.RemoveAll() 
        
        # List of coordinates.
        #
        coordList = [str(ulx)+";"+str(lry),str(ulx)+";"+str(uly),str(lrx)+";"+str(uly),str(lrx)+";"+str(lry)]

        # For each coordinate set, create a point object and add the x- and
        # y-coordinates to the point object, then add the point object
        # to the array object.
        #
        for coordPair in coordList:
            pnt = gp.createobject("Point")
            x, y = coordPair.split(";")
            pnt.x = x
            pnt.y = y
            vertex_array.add(pnt)
            # Create a polygon geometry object using the array object
            # created from the coordinate list above.   
            polyGeom = gp.createobject("geometry", "polygon", vertex_array)
    
        # add new cursor to the shapefile
        insert_cur = gp.InsertCursor(outFC) #outFC must be shp?  would shp work in memory?
        new_feature = insert_cur.newRow()
        new_feature.shape = polyGeom

        new_feature.setValue (field, fieldValue)
        insert_cur.InsertRow(new_feature)  
    
    except arcgisscripting.ExecuteError:
        raise_geoproc_error('watools_util')
    except:
        raise_python_error('watools_util')
   
    return outFC
   
#   
##Code for getting lc distances in script 2:
#    ##THIS DOESN'T WORK.  Any time a lcp goes through an intermediate airspace, it is shorter than the lcp going through the adjacency boundary.
#    #Therefore LCDists reported by this script will sometimes be longer than actual.
#    #try with euclidean?
#
#WILL WORK for grouping though- because always focusing on nearest first.  Nearest neighbor will always have correct distance.  
#
#    #use as first guess, and LIMIT cwd calcs to this value plus 100 km??
#
#
#
## Import system modules
#import sys, string, os, arcgisscripting
#from numpy import *
#from watools_util import *
#
## Create the Geoprocessor object
#gp = arcgisscripting.create(9.3)
#gp.CheckOutExtension("spatial")
##
##gp.AddToolbox("C:/Program Files/ArcGIS/ArcToolbox/Toolboxes/Data Management Tools.tbx")
##gp.AddToolbox("C:/Program Files/ArcGIS/ArcToolbox/Toolboxes/Analysis Tools.tbx")
#
#gp.OverwriteOutput = 1
#
#alloc= "C:\\watools\\Project1\\Adj\\cwd_alloc_ras"
#cwd = "C:\\watools\\Project1\\Adj\\cwd"
#res = "C:\\watools\\test_data\\resist_ras"
#
##alloc= "C:\\li\\Adj\\cwd_alloc_ras"
##cwd = "C:\\li\\Adj\\cwd"
##res = "C:\\li\\li_crs_tmp"
#
#cellSize = gp.Describe(res).MeanCellHeight
#gp.CellSize = cellSize
#
#print 'Cell height is ',str(gp.Cellsize)
#posShift = gp.CellSize
#negShift = -1*float(gp.CellSize)
#scratchDir = "c:\\testoutput\\adj3\\"
#
#gp.workspace = scratchDir
#
#gp.addmessage('Calculating least-cost distances crossing horizontal allocation boundaries...')
#startTime=time.clock()
#gp.Shift_management(alloc, "alloc_r", posShift, "0")
#gp.Shift_management(cwd, "cwd_r", posShift, "0")
#gp.Shift_management(res, "res_r", posShift, "0")
#
#expression1 = "(con(" + alloc + " - alloc_r <> 0, 1))"
#expression = expression1 + " * ( " + cwd + " + cwd_r + (( " + res + " + res_r ) * " + str(cellSize) + " / 2 ))"
#gp.SingleOutputMapAlgebra_sa(expression, "lcdist_r")
#
#alloc_r = "alloc_r"
#lcDistTable_r,adjTable_r = getLeastCostDistancesFromShift(gp.workspace,alloc,alloc_r,"lcdist_r")
#startTime,hours,mins,secs = elapsedTime(startTime)
#
#
#gp.addmessage('Calculating least-cost distances crossing upper-left diagonal allocation boundaries...')
#gp.Shift_management(alloc, "alloc_ul", negShift, posShift)
#gp.Shift_management(cwd, "cwd_ul", negShift, posShift)
#gp.Shift_management(res, "res_ul", negShift, posShift)
#
#expression1 = "(con(" + alloc + " - alloc_ul <> 0, 1))"
#expression = expression1 + " * ( " + cwd + " + cwd_ul + (( " + res + " + res_ul ) * " + str(cellSize) + " * Sqrt(2) / 2 ))"
#gp.SingleOutputMapAlgebra_sa(expression, "lcdist_ul")
#
#alloc_ul = "alloc_ul"
#lcDistTable_ul,adjTable_ul = getLeastCostDistancesFromShift(gp.workspace,alloc,alloc_ul,"lcdist_ul")
#startTime,hours,mins,secs = elapsedTime(startTime)
#
#
#gp.addmessage('Calculating least-cost distances crossing upper-right diagonal allocation boundaries...')
#gp.Shift_management(alloc, "alloc_ur", posShift, posShift)
#gp.Shift_management(cwd, "cwd_ur", posShift, posShift)
#gp.Shift_management(res, "res_ur", posShift, posShift)
#
#expression1 = "(con(" + alloc + " - alloc_ur <> 0, 1))"
#expression = expression1 + " * ( " + cwd + " + cwd_ur + (( " + res + " + res_ur ) * " + str(cellSize) + " * Sqrt(2) / 2 ))"
#gp.SingleOutputMapAlgebra_sa(expression, "lcdist_ur")
#
#alloc_ur = "alloc_ur"
#lcDistTable_ur,adjTable_ur = getLeastCostDistancesFromShift(gp.workspace,alloc,alloc_ur,"lcdist_ur")
#startTime,hours,mins,secs = elapsedTime(startTime)
#
#
#gp.addmessage('Calculating least-cost distances crossing vertical allocation boundaries...')
#gp.Shift_management(alloc, "alloc_u", "0", posShift)
#gp.Shift_management(cwd, "cwd_u", "0", posShift)
#gp.Shift_management(res, "res_u", "0", posShift)
#
#expression1 = "(con(" + alloc + " - alloc_u <> 0, 1))"
#expression = expression1 + " * ( " + cwd + " + cwd_u + (( " + res + " + res_u ) * " + str(cellSize) + " / 2 ))"
#gp.SingleOutputMapAlgebra_sa(expression, "lcdist_u")
#
#alloc_u = "alloc_u"
#lcDistTable_u,adjTable_u = getLeastCostDistancesFromShift(gp.workspace,alloc,alloc_u,"lcdist_u")
#startTime,hours,mins,secs = elapsedTime(startTime)
#
#lcDistTable = combineLeastCostDistances(lcDistTable_r,lcDistTable_u,lcDistTable_ur,lcDistTable_ul)
#adjTable = combine_adjacency_tables(adjTable_r,adjTable_u,adjTable_ur,adjTable_ul)
#
#
#
#outcsvfile = scratchDir+'TEMP.CSV' #this will replace adj files...
#outfile = open(outcsvfile,"w")
#coreIds = 'HCAXXXX' #fixme temp
#outfile.write ("#Edge" + "," + str(coreIds) + "," + str(coreIds) + "_1" + "\n")
#for x in range(0,len(adjTable)):
#    outfile.write ( str(x) + "," + str(adjTable[x,0]) + "," + str(adjTable[x,1]) + "\n" )
#outfile.close()
#
#outcsvLogfile = scratchDir+'TEMP.CSV'
#outfile = open(outcsvLogfile,"w") #log file
#for x in range(0,len(adjTable)):
#    outfile.write ( str(x) + "," + str(adjTable[x,0]) + "," + str(adjTable[x,1]) + "\n" )
#outfile.close()
#
#
#print 'final table:'
#print lcDistTable
#    
#

#def getAdjList(workspace,alloc_ras,coreShapefile,coreIds): #Fixme: good method, but takes 200 secs instead of 20 with LI data.  Not implemented.
#    try:
#        gp.workspace=workspace
#        allocFC = "alloc.shp"
#        count,statement = 0, 'gp.rastertopolygon_conversion(alloc_ras,allocFC,"SIMPLIFY")'
#        while True:
#            try: exec statement
#            except:
#                count,tryAgain = hiccup_test(count,statement)
#                if not tryAgain: exec statement
#            else: break
#        
#        gp.workspace = 'in_memory'
#        
#        allocFL = "allocFeatureLayer"   
#        #desc = gp.Describe(inFC)    
#        gp.MakeFeatureLayer_management(workspace+ "\\"+allocFC, allocFL , "", "", "")
#        
#        coreList = getCoreList(coreShapefile,coreIds,coreIds)
#        cores = coreList[:,0]
#        
#        field = "GRIDCODE"
#           
#        adjList = zeros((0,2),dtype='int32')
#        
#        for core in cores:
#            #startTime = time.clock()
#    
#            #gp.addmessage('core '+ str(core))
#            adjacentCores = []
#            expression = field + " = " + str(core)
#    
#    #THIS IS WHAT TAKES SO LONG
#            gp.SelectLayerByAttribute_management(allocFL, "NEW_SELECTION", expression)
#            gp.SelectLayerByLocation_management(allocFL, "BOUNDARY_TOUCHES", allocFL, "", "NEW_SELECTION" )
#    ###########################
#    
#            rows = gp.SearchCursor(allocFL)
#            row = rows.Next()
#            while row:
#                value = row.GetValue(field)
#                if value > core:
#                    adjacentCores = append(adjacentCores,int(value))
#                row = rows.Next()
#           
#            uniqueAdjacentCores = unique(adjacentCores)
#            newAdjacentCores=zeros((len(uniqueAdjacentCores),2),dtype='int32')
#            newAdjacentCores[:,0]=core
#            newAdjacentCores[:,1]=uniqueAdjacentCores
#            adjList = append(adjList,newAdjacentCores,axis=0)
#            
#        return adjList
#    except arcgisscripting.ExecuteError:
#        raise_geoproc_error('watools_util')
#    except:
#        raise_python_error('watools_util')
#        
        
        
#def writeLinkMapsOLD(workspace,linkTableFile,coreShapefile,coreIds,coreAreaIds,step):  #Old code, not used, but save for future grouping script
#    try:
#        startTime = time.clock()
#        useEsriCentroids = True #Use new centroids that actually fall inside polygons.
#        gp.workspace=workspace # "in_memory"
#        gp.OverwriteOutput = 1
#        #########################################################
#        ##Relevant linkTable INDEXES
#        linkIdCol=0 #Link ID
#        HCA1Col=1 #1st 1st core area link connects
#        HCA2Col=2 #2nd core area link connects
#        core1Col=3 #core ID of1st core area link connects
#        core2Col=4 #core ID of2nd core area link connects
#        linkTypeCol=5 #0=no link. 2=corridor, 3=intermediate core area detected, 4=too long EucDist, 5=too long lcDist
#        eucDistCol=6
#        cwDistCol=7
#        #########################################################
#        
#        #Define Coordinate System for output shapefiles
#        desc = gp.Describe
#        SR = desc(coreShapefile).SpatialReference
#        ## ------------------------------------------
#        
#        linkTable = loadLinkTable(linkTableFile)
#        numLinks = linkTable.shape[0]
#        linkTable[:,linkTypeCol]=where(linkTable[:,linkTypeCol]==11,1,linkTable[:,linkTypeCol]) #Convert 11 grouplinks to 1
#        #def createLinkMap(linkTable,coreShapefile,coreAreaIds,coreIds):
#        gp.toolbox = "management"
#        
#        coresForLinework = "cores_for_linework.shp"
#        HCAsForLinework="cores_for_linework.shp"
#
#
##old method for mean_center_states to get geometric center
#        if useEsriCentroids == True: #New methods
#            pointArray = getCentroids(coreShapefile,coreAreaIds) 
#            makePoints(gp.workspace,pointArray,coresForLinework,coreAreaIds)
#            pointArray = getCentroids(coreShapefile,coreIds) 
#            makePoints(gp.workspace,pointArray,HCAsForLinework,coreIds)
#
#        else: #Use old method for mean_center_states to get geometric center
#            #Prepare for dissolve operations
#            groupDissolveShapefile = "tempgroups.shp" #old method for mean_center_states to get geometric center
#            coreDissolveShapefile = "tempcores.shp"#old method for mean_center_states to get geometric center
#            if gp.exists(coreDissolveShapefile):
#                gp.delete_management(coreDissolveShapefile)
#            if gp.exists(groupDissolveShapefile):#old method for mean_center_states to get geometric center
#                gp.delete_management(groupDissolveShapefile)#old method for mean_center_states to get geometric center
#            
#            ##Dissolve core polygons and get geographic centers
#            gp.Dissolve_management(coreShapefile, coreDissolveShapefile, coreIds)          
#            gp.MakeFeatureLayer(coreDissolveShapefile,"fcores")
#            gp.MeanCenter_stats("fcores", HCAsForLinework, "#", coreIds, "#") #old method, gave geometric center
#
#            ##Dissolve group polygons and get geographic centers
#            gp.Dissolve_management(coreShapefile, groupDissolveShapefile, coreAreaIds)#old method for mean_center_states to get geometric center           
#            gp.MakeFeatureLayer(groupDissolveShapefile,"fgroups")
#            gp.MeanCenter_stats("fgroups", coresForLinework, "#", coreAreaIds, "#") #old method, gave geometric center
#
#        #split linkTable into coreLinks, which are links BETWEEN GROUPS, and HCALinks, which are WITHIN GROUPS.
#        HCALinks=zeros((0,10))
#        coreLinks=zeros((0,10))
#        
#        #Pull out HCAlinks and grouplinks
#        numLinks=linkTable.shape[0]
#        rows,cols = where(linkTable[:,linkTypeCol:linkTypeCol+1]== 2)
#        coreLinks=linkTable[rows,:]
#        rows,cols=where(linkTable[:,linkTypeCol:linkTypeCol+1]== 10)
#        componentLinks=linkTable[rows,:]
#        coreLinks=append(coreLinks,componentLinks,axis=0)
#        
#        del componentLinks
#
#        rows,cols=where(linkTable[:,linkTypeCol:linkTypeCol+1]== 76)
#        componentLinks=linkTable[rows,:]
#        coreLinks=append(coreLinks,componentLinks,axis=0)
#
#        rows,cols=where(linkTable[:,linkTypeCol:linkTypeCol+1]== 1)
#        HCALinks=linkTable[rows,:]
#        
#        #create coreCoords array, with geographic centers of groups
#        cur = gp.SearchCursor(coresForLinework)
#        row = cur.Next()
#        
#        uniqueGroups = unique(coreLinks[:,core1Col:core2Col+1])
#        uniqueCores = unique(HCALinks[:,HCA1Col:HCA2Col+1])
#        
#        coreCoords = zeros((len(uniqueCores),6))
#        #coreId HCAx HCAy groupId corex corey
#        
#        coreCoords[:,0]=uniqueCores[:]
#        
#        coreList = getCoreList(coreShapefile,coreIds,coreAreaIds)
#        
#        #Get core coordinates into coreCoords
#        for i in range(0,len(coreCoords)):
#            
#            core=coreCoords[i,0]#ID of core we're operating on
#        
#            #get ID of group core belongs to
#            for j in range(0,len(coreList)):
#                if coreList[j,0]==core:
#                    group=coreList[j,1]
#                    coreCoords[i,3]=group
#                    break
#            
#            #get geographic center of core    
#            cur = gp.SearchCursor(HCAsForLinework)
#            row = cur.Next()
#            while row:
#                coreShape=row.GetValue(coreIds)
#                if coreShape==core:
#                    coreCoords[i,1]=row.GetValue("XCoord") 
#                    coreCoords[i,2]=row.GetValue("YCoord") 
#                    break
#                row = cur.Next()
#            del cur, row
#        
#            #get geographic center of group
#            cur = gp.SearchCursor(coresForLinework)
#            row = cur.Next()
#            while row:
#                groupShape=row.GetValue(coreAreaIds)
#                if groupShape==group:
#                    coreCoords[i,4]=row.GetValue("XCoord") 
#                    coreCoords[i,5]=row.GetValue("YCoord")
#                    break
#                row = cur.Next()
#            del cur, row
#               
#        #Create linkCoords array, which holds following info:
#        #linkIdCol core1Col core2Col eucdist lcdist core1x core1y core2x core2y 
#            
#        linkCoords=zeros((len(coreLinks),10))
#        linkCoords[:,0]=coreLinks[:,linkIdCol] #link IDs
#        linkCoords[:,1:3]=sort(coreLinks[:,core1Col:core2Col+1]) #core1Col and core2Col IDs, sorted left-to-right
#        linkCoords[:,3:5]=coreLinks[:,eucDistCol:cwDistCol+1] #euc and lc distances
#        linkCoords[:,9]=coreLinks[:,linkTypeCol] #2=minNN,10=component
#        
#        if len(coreLinks) > 0:   
#            ind=lexsort((linkCoords[:,2],linkCoords[:,1])) #sort by 1st groupId then by 2nd groupId
#            linkCoords = linkCoords[ind]    
#           
#        #Get rid of duplicate group pairs, retaining minimum distances between groups
#        numLinks=len(linkCoords)
#        x=1
#        while x<numLinks:
#            if linkCoords[x,1]==linkCoords[x-1,1] and linkCoords[x,2]==linkCoords[x-1,2]:
#                linkCoords[x-1,3]=min(linkCoords[x-1,3],linkCoords[x,3]) #Retain minimum euc distance between group pair
#                linkCoords[x-1,4]=min(linkCoords[x-1,4],linkCoords[x,4]) #Retain minimum lc distance between group pair
#                #linkCoords=deleterow(linkCoords,x)
#                #numLinks=numLinks-1
#                #x=x-1
#                linkCoords[x,0]=0 #mark for deletion  
#            x=x+1
#        if numLinks>0:
#            delRows=asarray(where(linkCoords[:,0]==0))
#            delRowsVector = zeros((delRows.shape[1]),dtype="int32")
#            delRowsVector[:] = delRows[0,:]
#            linkCoords=deleterow(linkCoords,delRowsVector)        
#            del delRows
#            del delRowsVector
#            
#        
#        #Get group coordinates into linkCoords
#        #linkCoords arranged by:
#        #linkIdCol core1Col core2Col eucdist lcdist core1x core1y core2x core2y 
#        for i in range(0,len(linkCoords)):
#            grp1=linkCoords[i,1]
#            grp2=linkCoords[i,2]
#            cur = gp.SearchCursor(coresForLinework)
#            row = cur.Next()
#            while row:
#                groupShape=row.GetValue(coreAreaIds)
#                if groupShape==grp1:
#                    linkCoords[i,5]=row.GetValue("XCoord") 
#                    linkCoords[i,6]=row.GetValue("YCoord") 
#                elif groupShape==grp2:
#                    linkCoords[i,7]=row.GetValue("XCoord") 
#                    linkCoords[i,8]=row.GetValue("YCoord") 
#                row = cur.Next()
#            del cur, row
#        
#        if len(coreLinks) > 0:  
#            ind=argsort((linkCoords[:,linkIdCol])) #sort by linkIdCol
#            linkCoords = linkCoords[ind]
#        
#       
#        linkTypeCol=5
#        linkTypes = linkTable[:,linkTypeCol]
#        numGroupLinks = sum(linkTypes==1)
#        numGroupLinks2 = sum(linkTypes==11)
#        numGroupLinks =numGroupLinks +numGroupLinks2
#        if numGroupLinks > 0:
#            #make coreGroupings.shp using coreCoords table
#            # will contain linework between each core and its group center
#            coreGroupingsShapefile='coreGroupings_step' + str(step) + '.shp'
#            
#            #fixme: may want to do calcs "in memory" and save to disk later for speed
#            gp.CreateFeatureclass(gp.workspace, coreGroupingsShapefile, "POLYLINE")
#            gp.defineprojection(coreGroupingsShapefile, SR)
#    
#            #Add fields to core groupings shapefile
#            gp.addfield(coreGroupingsShapefile,"Core_ID","SHORT")
#            gp.addfield(coreGroupingsShapefile,"Group_ID","SHORT")            
#            
#            
#            #Open a cursor to insert rows into the shapefile.
#            cur = gp.InsertCursor(coreGroupingsShapefile)
#            
#            #Create an Array and Point object.
#            lineArray = gp.CreateObject("Array")
#            pnt = gp.CreateObject("Point")
#            
#            #coreCoords indices:
#            #HCAId HCAx HCAy coreId corex corey
#            numLinks=len(coreCoords)
#            ##Loop through each record in coreCoords table
#            for i in range(0,numLinks):
#            
#                #Set the X and Y coordinates for origin vertex.
#                pnt.x = coreCoords[i,1]
#                pnt.y = coreCoords[i,2]
#                #Insert it into the line array    
#                lineArray.add(pnt)
#            
#                #Set the X and Y coordinates for destination vertex
#                pnt.x = coreCoords[i,4]
#                pnt.y = coreCoords[i,5]
#                #Insert it into the line array
#                lineArray.add(pnt)
#            
#                #Insert the new poly into the feature class.
#                feat = cur.NewRow()
#                feat.shape = lineArray
#                cur.InsertRow(feat)
#                lineArray.RemoveAll()
#            del cur
#            
#            
#            #Add attribute data to core groupings shapefile
#            rows=gp.UpdateCursor(coreGroupingsShapefile)
#            row = rows.Next()
#            line = 0
#            while row:
#                #coreCoords indices:
#                #HCAId HCAx HCAy coreId corex corey
#                row.SetValue("ID", line)
#                row.SetValue("Core_ID", coreCoords[line,0])
#                row.SetValue("Group_ID", coreCoords[line,3])
#                rows.UpdateRow(row)
#                row = rows.Next()
#                line = line + 1
#            # delete cursor and row points to remove locks on the data
#            del row, rows      
#        
#        #
#        coreLinksShapefile = 'sticks_step' + str(step) + '.shp'
#        
#        #make coreLinks.shp using linkCoords table
#        # will contain linework between each pair of connected groups
#        gp.CreateFeatureclass(gp.workspace, coreLinksShapefile, "POLYLINE")
#        
#        
#        #Define Coordinate System
#        desc = gp.Describe
#        SR = desc(coreShapefile).SpatialReference
#        gp.defineprojection(coreLinksShapefile, SR)
#        
#        #ADD ATTRIBUTES
#        gp.AddField_management(coreLinksShapefile, "Link_ID", "SHORT")
#        gp.AddField_management(coreLinksShapefile, "Link_Info", "TEXT")
#        gp.AddField_management(coreLinksShapefile, "From_Core", "SHORT")
#        gp.AddField_management(coreLinksShapefile, "To_Core", "SHORT")
#        gp.AddField_management(coreLinksShapefile, "Euc_Dist", "FLOAT")
#        gp.AddField_management(coreLinksShapefile, "CW_Dist", "FLOAT")
#        gp.AddField_management(coreLinksShapefile, "cwdToEucDistRatio", "FLOAT")                        
#        #
#        
#        #Open a cursor to insert rows into the shapefile.
#        cur = gp.InsertCursor(coreLinksShapefile)
#        
#        #Create an Array and Point object.
#        lineArray = gp.CreateObject("Array")
#        pnt = gp.CreateObject("Point")
#        
#        
#        #linkCoords indices:
#        #linkIdCol core1Col core2Col eucdist lcdist core1x core1y core2x core2y 
#        numLinks=len(linkCoords)
#        ##Loop through each record in linkCoords table
#        for i in range(0,numLinks):
#        
#            #Set the X and Y coordinates for origin vertex.
#            pnt.x = linkCoords[i,5]
#            pnt.y = linkCoords[i,6]
#            #Insert it into the line array    
#            lineArray.add(pnt)
#        
#            #Set the X and Y coordinates for destination vertex
#            pnt.x = linkCoords[i,7]
#            pnt.y = linkCoords[i,8]
#            #Insert it into the line array
#            lineArray.add(pnt)
#        
#            #Insert the new poly into the feature class.
#            feature = cur.NewRow()
#            feature.shape = lineArray
#            cur.InsertRow(feature)
#           
#            lineArray.RemoveAll()
#        del cur 
#            
#        #Add attribute data to link shapefile
#        rows=gp.UpdateCursor(coreLinksShapefile)
#        row = rows.Next()
#        line = 0
#        while row:
#            #linkCoords indices:
#            #linkIdCol core1Col core2Col eucdist lcdist core1x core1y core2x core2y
#            row.SetValue("Link_ID", linkCoords[line,0])
#            if linkCoords[line,9] == 2:
#                row.SetValue("Link_Info", "Connects_cores")
#            elif linkCoords[line,9] == 10:
#                row.SetValue("Link_Info", "Connects_constellations")
#            row.SetValue("From_Core", linkCoords[line,1])
#            row.SetValue("To_Core", linkCoords[line,2])
#            row.SetValue("Euc_Dist", linkCoords[line,3])
#            row.SetValue("CW_Dist", linkCoords[line,4])
#            distRatio1 = linkCoords[line,4]/linkCoords[line,3]
#            if linkCoords[line,4] <=0 or linkCoords[line,3] <= 0:
#                row.SetValue("cwdToEucDistRatio", -1)
#            else:
#                row.SetValue("cwdToEucDistRatio", linkCoords[line,4]/linkCoords[line,3])
#
#            rows.UpdateRow(row)
#            row = rows.Next()
#            line = line + 1
#        # delete cursor and row points to remove locks on the data
#        del row, rows
#        
#        #clean up temp files
#        if useEsriCentroids == False:
#            gp.delete_management(coreDissolveShapefile) #old method for mean_center_states to get geometric center
#            gp.delete_management(groupDissolveShapefile) #old method for mean_center_states to get geometric center
#        gp.delete_management(coresForLinework)
#        gp.delete_management(HCAsForLinework)
#        
##For in_memory operation... in progress
#        #coreLinksFINAL = outputWS + '\\'+'sticks_step' + str(step) + '.shp'
#        #coreGroupingsFINAL = outputWS + '\\' + 'coreGroupings_step' + str(step) + '.shp'  
#        #gp.copy_management(coreLinksShapefile,coreLinksFINAL)        
#        #gp.copy_management(coreGroupingsShapefile,coreGroupingsFINAL )
#        
#        startTime,hours,mins,secs = elapsedTime(startTime)
#        return
#    except arcgisscripting.ExecuteError:
#        raise_geoproc_error('watools_util')
#    except:
#        raise_python_error('watools_util')
#


From script s3_calcCwds.py:
                    #code to copy results to file geodatabase (TOO SLOW)
                #startTime = time.clock()
                #outDistanceRaster2 = cwdPath2 + "\cwd_" + str(int(sourceCore)) 
                #gp.CopyRaster_management (outDistanceRaster, outDistanceRaster2)#, "raster")
                #endTime = time.clock()
                #processTime = round((endTime - startTime), 2)
                #gp.addmessage('Done with copy. Elapsed time = '+str(processTime)+' seconds.')
    
#valuelist method
                    #gp.mask = "lcp"
                    #expression = core_ras + " * 1"
                    #gp.SingleOutputMapAlgebra_sa(expression, "testRas") #too slow- 13 sec vs faster zonalstats
                    #gp.mask = boundResis
                    #
                    #values = CreateValueList("testRas")
                    #if values == 'Failed':
                    #    pass
                    #else:
                    #    for value in range(len(values)):
                    #        if values[value] != int(sourceCore) and values[value] != int(targetCore):
                    #            break



From s5_calcLccs.py:
            #
        #gp.addmessage('\n*************************************************')
        #gp.addmessage('**********NOTE: Geodatabase writing disabled, writing old-fashioned grids for speed!')
        #gp.addmessage('**********ArcGIS cannot write quickly to directories with many grids.  ')
        #gp.addmessage('*************************************************\n')
        #
    #ANYWAY YOU SLICE IT GDBS SEEM SLOW- Regardless of where they are, what workspace set to, etc.
        #lccGDBname = "rawlccs"
        #gdbPath=projectDir+"\\"+ lccGDBname+".gdb"
        #
        #if gp.Exists(gdbPath):
        #    gp.addmessage('\nDeleting and re-creating ' + lccGDBname+".gdb")
        #    gp.delete_management(projectDir + "\\" + lccGDBname+".gdb")
        #try:
        #    gp.createfilegdb(projectDir, lccGDBname+".gdb")                   
        #except:
        #    pass
        #
    
                    
                    #it's not garbage collection that is slowing
                    #startTime = time.clock()
                    #gc.collect()
                    #endTime = time.clock()
                    #processTime = round((endTime - startTime), 2)        
                    #gp.addmessage('gc = ' + str(processTime) +' seconds.\n')
    
    
    #SLOW SLOW SLOW                            
                    #gdbraster = gdbPath + "\\TEST_" + str(corex) + "_" + str(corey)
                    #startTime=time.clock()
                    #gp.CopyRaster_management(lccRaster,gdbraster) 
                    #endTime = time.clock()
                    #processTime = round((endTime - startTime), 2)        
                    #gp.addmessage('COPY TIME = ' + str(processTime) +' seconds.\n')
                    
                    
                    

#cannot use SOMA for mosaic (unless jump thru hoops- search for basedemgrd and tempfillgrd to get thread).  


#This code is much slower than mosaic- maybe pyramids aren't that bad afterall?
#NOTE: code is not working, probably wrong syuntax for con statement
#See also code at: http://forums.arcgis.com/threads/8579-Adding-Together-Rasters

                #mosaicRaster2 = 'c:\\tempmos'
                #if numGridsWritten == 0 and dirCount == 0:
                #    gp.CopyRaster_management(lccNormRaster,mosaicRaster2)
                #else:                    
                #
                #    tempRas = 'c:\\tempx'
                #    expression = "con(isnull(" + lccNormRaster + ")," + mosaicRaster2 + "," + lccNormRaster + ")"
                #    gp.SingleOutputMapAlgebra_sa(expression, tempRas)
                #    
                #    tempRas2 = 'c:\\temp2'
                #    expression = "min(" + tempRas + "," + mosaicRaster2 + ")"
                #    gp.SingleOutputMapAlgebra_sa(expression, tempRas2)
                #    print 'test code done'
                #    startTime,hours,mins,secs = elapsedTime(startTime)
                #    gp.delete_management(mosaicRaster2) #May not be needed/
                #    gp.CopyRaster_management(tempRas2,mosaicRaster2)
                