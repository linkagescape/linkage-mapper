##*****************************************************************
## 2011_0128
## NAME: lm_util.py
##
## SUMMARY: Contains functions called by linkage mapper scripts
##
## SOFTWARE: ArcGIS 9.3 (requires Spatial Analyst extension)
##           Python 2.5
##
##*****************************************************************

import os
import string
import csv
from string import split
from time import localtime, strftime
import traceback

import arcgisscripting, sys, time, shutil
from numpy import *
from numpy import loadtxt, where, delete, arange

import lm_config

PROJECTDIR = lm_config.PROJECTDIR
OUTPUTDIR = lm_config.OUTPUTDIR
LOGDIR = lm_config.LOGDIR
SCRATCHDIR = lm_config.SCRATCHDIR
DATAPASSDIR = lm_config.DATAPASSDIR
DATAPASSARCHDIR = lm_config.DATAPASSARCHDIR
CWDBASEDIR = lm_config.CWDBASEDIR
CWDSUBDIR_MAME = lm_config.CWDSUBDIR_MAME
COREFC = lm_config.COREFC
COREFN = lm_config.COREFN
OUTPUTGDB = lm_config.OUTPUTGDB
STEP1 = lm_config.STEP1
STEP2 = lm_config.STEP1
STEP3 = lm_config.STEP3
STEP4 = lm_config.STEP4
STEP5 = lm_config.STEP5
S2EUCDISTFILE = lm_config.S2EUCDISTFILE
FCORES = lm_config.FCORES
GP = lm_config.GP

def get_linktable_row(linkId, linkTable):
    try:
        linkIdCol=0 #Link ID
        if linkTable[linkId-1,linkIdCol]==linkId: #Most likely.  linkTables tend to be in order with no skipped links.
            linkTableRow = linkId-1
            return linkTableRow
        else:
            numLinks = linkTable.shape[0]
            for linkTableRow in range (0,numLinks):
                if linkTable[linkTableRow,linkIdCol]==linkId:
                    return linkTableRow
        return -1 # Not found
    except:
        raise_python_error('lm_util')   

def get_link_type_desc(linkTypeCode):
    """For a linkType code in linkTable, returns description to attribute link in lcp and link maps
    NOTE: These are being overhauled to be more descriptive, particularly for nearest neighbor and 
    cluster (constellation) links identified in step 4.
    
    """
    activeLink = '0'
    if linkTypeCode == -2:
        linkTypeDesc = "Not_nearest_N_neighbors"
    elif linkTypeCode == -20:
        linkTypeDesc = "User_removed"        
    elif linkTypeCode == 1:
        linkTypeDesc = "Within-core"
    elif linkTypeCode == 2:
        linkTypeDesc = "Connects_cores"
        activeLink = '1'
    elif linkTypeCode == 3:
        linkTypeDesc = "Intermediate_core_detected"        
    elif linkTypeCode == 4:        
        linkTypeDesc = "Too_long_Euclidean_dist"
    elif linkTypeCode == 5:        
        linkTypeDesc = "Too_long_least_cost_dist"        
    elif linkTypeCode == 6:        
        linkTypeDesc = "Too_short_Euclidean_dist"
    elif linkTypeCode == 7:        
        linkTypeDesc = "Too_short_least_cost_dist"  
    elif linkTypeCode == 10:        
        linkTypeDesc = "Connects_constellations"
        activeLink = '1'        
    else:
        linkTypeDesc='Unknown'
    return activeLink, linkTypeDesc        


def get_links_from_core_pairs(linkTable, firstCore, secondCore):
    """Given two cores, finds the row in the link table corresponding to the link that connects them"""
    try:
        rows=zeros((0),dtype="int32")
        core1Col=1 #core ID of1st core area link connects
        core2Col=2 #core ID of2nd core area link connects        
        numLinks = linkTable.shape[0]
        for link in range (0,numLinks):
            corex=int(linkTable[link,core1Col])
            corey=int(linkTable[link,core2Col])
            if int(corex)==int(firstCore) and int(corey)==int(secondCore):
                rows = append(rows, link)
            elif int(corex)==int(secondCore) and int(corey)==int(firstCore):
                rows = append(rows, link)
        return rows
    except:
        raise_python_error('lm_util') 


def drop_links(linkTable, maxeud, mineud, maxcwd, mincwd, disableLeastCostNoVal):
    """Inactivates links that fail to meet min or max length criteria"""
    try:
        dashline(1)
        #set linkTable INDEXES
        #########################################################
        ##linkTable INDEXES
        linkIdCol = 0 #Link ID
        core1Col = 1 #core ID of1st core area link connects
        core2Col = 2 #core ID of2nd core area link connects
        linkTypeCol = 5 #0=no link. 2=corridor, 3=intermediate core area detected, 4=too long EucDist, 5=too long lcDist 6=too short EucDist, 7=too short cwDist, 10 = cluster link
        eucDistCol = 6
        cwDistCol = 7
        #########################################################

        numLinks = linkTable.shape[0]
        numDroppedLinks=0
        coreList=linkTable[:,core1Col:core2Col+1]
        coreList=sort(coreList)

        if disableLeastCostNoVal:
            for x in range(0,numLinks):
                linkId = str(int(linkTable[x,linkIdCol]))                            
                if linkTable[x,cwDistCol] == -1:
                    if linkTable[x,linkTypeCol]==2 or linkTable[x,linkTypeCol]==10: #Check only enabled corridor links
                        corex=str(int(linkTable[x,core1Col]))
                        corey=str(int(linkTable[x,core2Col]))
                        GP.addmessage("The least-cost corridor between " + str(corex) + " and " + str(corey)+" (link #"+linkId+ ") has an unknown length in cost distance units.  This means it is longer than the max cost-weighted distance specified in the 'Calc CWDs' script OR it passes through NODATA cells and will be dropped.\n")
                        
                        linkTable[x,linkTypeCol] = 5 #Disable link
                        numDroppedLinks = numDroppedLinks + 1

        #Check for corridors that are too long in Euclidean or cost-weighted distance
        if maxeud is not None or maxcwd is not None: 
            for x in range(0,numLinks):
                linkId = str(int(linkTable[x,linkIdCol]))
                if maxeud is not None:
                    if linkTable[x,eucDistCol] > maxeud:
                        if linkTable[x,linkTypeCol] == 2 or linkTable[x,linkTypeCol] == 10: #Check only enabled corridor links
                            corex=str(int(coreList[x,0]))
                            corey=str(int(coreList[x,1]))
                            GP.addmessage("Link #"+linkId+ " connecting cores " + str(corex) + " and " + str(corey) + " is  " + str(linkTable[x,eucDistCol]) + " units long- too long in Euclidean distance.")
                            linkTable[x,linkTypeCol]=4 #Disable link
                            numDroppedLinks=numDroppedLinks+1
                    
                    if maxcwd is not None:
                        if (linkTable[x,linkTypeCol] == 2
                            or linkTable[x,linkTypeCol] == 10): #Check only enabled corridor links
                            if (linkTable[x,cwDistCol] > maxcwd): #BHM 8/11/10 added check for -1 cwDistCol
                                corex=str(int(linkTable[x,core1Col]))
                                corey=str(int(linkTable[x,core2Col]))
                                GP.addmessage("Link #" + linkId + " connecting" 
                                              "cores " + str(corex) + " and " 
                                              + str(corey) + " is  " 
                                              + str(linkTable[x,cwDistCol]) 
                                              + " units long- too long in cost" 
                                              "distance units.")
                                linkTable[x,linkTypeCol]=5 #Disable link
                                numDroppedLinks=numDroppedLinks+1
        
        if mineud is not None or mincwd is not None: 
            for x in range(0,numLinks):
                linkId = str(int(linkTable[x,linkIdCol]))
                if linkTable[x,linkTypeCol]==2 or linkTable[x,linkTypeCol]==10: #Check only enabled corridor links
                    if mineud is not None:
                        if linkTable[x,eucDistCol] < mineud:
                            corex=str(int(coreList[x,0]))
                            corey=str(int(coreList[x,1]))
                            GP.addmessage("Link #"+linkId+ " connecting cores " + str(corex) + " and " + str(corey) + " is only " + str(linkTable[x,eucDistCol]) + " units long- too short in Euclidean distance.")
                            linkTable[x,linkTypeCol]=6 #Disable link
                            numDroppedLinks=numDroppedLinks+1
        
                    if mincwd is not None:
                        if (linkTable[x,cwDistCol] < mincwd) and (linkTable[x,cwDistCol]) !=-1:
                            if linkTable[x,linkTypeCol]==2 or linkTable[x,linkTypeCol]==10:
                                corex=str(int(linkTable[x,core1Col]))
                                corey=str(int(linkTable[x,core2Col]))
                                GP.addmessage("Link #"+linkId+ " connecting cores " + str(corex) + " and " + str(corey) + " is only " + str(linkTable[x,cwDistCol]) + " units long- too short in cost distance units.")
                                linkTable[x,linkTypeCol]=7 #Disable link
                                numDroppedLinks=numDroppedLinks+1
        return linkTable,numDroppedLinks
    except:
        raise_python_error('lm_util')   
        

def get_zonal_minimum(dbfFile):
    """Finds the minimum value in a table generated by zonal statistics"""
    try:
        rows = GP.searchcursor(dbfFile)
        row = rows.Next()
        try:
            coreMin = row.Min
        except:
            return 'Failed'
        while row:
            if coreMin > row.Min:
                coreMin = row.Min 
            row = rows.next()
        del row
        del rows
        return coreMin
    except arcgisscripting.ExecuteError:
        raise_geoproc_error('lm_util')
    except:
        raise_python_error('lm_util')   

  
def get_core_list():
    """Returns a list of core area IDs from polygon file"""
    try:
        #Get the number of cores
        coreCount = int(GP.GetCount_management(COREFC).GetOutput(0))#FIXME: I think this returns number of shapes, not number of unique cores.
        #Get core data into numpy array
        coreList = zeros((coreCount, 2))
        cur = GP.SearchCursor(COREFC)
        row = cur.Next()
        i = 0
        while row:
            coreList[i,0] = row.GetValue(COREFN)
            coreList[i,1] = row.GetValue(COREFN)  
            row = cur.Next()
            i = i + 1
        del cur, row
    except arcgisscripting.ExecuteError:
        raise_geoproc_error('lm_util')
    except:
        raise_python_error('lm_util')
    return coreList

    
def get_core_targets(core, linkTable):
    """For a given core area, returns a list of other core areas the core area is connected to"""
    try:
        #Relevant linkTable index values:
        core1Col=1 #core ID of 1st core area link connects
        core2Col=2 #core ID of 2nd core area link connects
        linkTypeCol=5 #0=no link. 2=corridor, 3=intermediate core area detected, 4=too long EucDist, 5=too long lcDist
        #########################################################
        
        targetList=zeros((len(linkTable),2),dtype="int32")
        targetList[:,0]=linkTable[:,core1Col]#possible targets 1st column=core1Col
        targetList[:,1]=linkTable[:,core2Col]#possible targets 2nd column=core2Col
        validPair=linkTable[:,linkTypeCol] #Copy of linkTypeCol column.  
        validPair=where(validPair==2,1,0) #2=map corridor.
        targetList[:,0]=multiply(targetList[:,0],validPair)
        targetList[:,1]=multiply(targetList[:,1],validPair)
        
        rows,cols=where(targetList==int(core))
        targetList=targetList[rows, 1-cols]
        targetList=unique(asarray(targetList))
    except arcgisscripting.ExecuteError:
        raise_geoproc_error('lm_util')
    except:
        raise_python_error('lm_util')
    return targetList
   

def elapsed_time(startTime):
    """Returns elapsed time given a start time"""
    try:
        now = time.clock()
        elapsed = now - startTime
        secs = int(elapsed)
        mins = int(elapsed / 60)
        hours = int(mins / 60)
        mins = mins-hours * 60
        secs = secs-mins * 60 - hours * 3600
        if mins == 0:
            GP.addmessage('That took ' + str(secs) + ' seconds.\n')
        elif hours == 0:
            GP.addmessage('That took ' + str(mins) + ' minutes and ' 
                           + str(secs) + ' seconds.\n')
        else:
            GP.addmessage('That took ' + str(hours) + ' hours ' + str(mins) 
                           + ' minutes and ' + str(secs) + ' seconds.\n')
        return now, hours, mins, secs
    except:
        raise_python_error('lm_util')


def report_pct_done(current, goal):
    """Reports percent done"""
    try:
        goal = int(ceil(goal / 10))
        goal = float(goal + 1) * 10
        pctdone = ((float(current) / goal) * 100)
        if pctdone/10 == int(pctdone/10):
            if pctdone > 0:
                GP.addmessage(str(int(pctdone))+ " percent done.")
    except:
        raise_python_error('lm_util')


def report_links(linkTable):
    """Prints number of links in a link table"""
    try:
        numLinks = linkTable.shape[0]
        GP.addmessage('There are ' + str(numLinks) + ' links in the table.')
        linkTypeCol=5
        linkTypes = linkTable[:,linkTypeCol]
        numCorridorLinks = sum(linkTypes==2)
        numGrpLinks = sum(linkTypes==1)+ sum(linkTypes==11)
        numComponentLinks = sum(linkTypes==10)
        if numComponentLinks  > 0:
            GP.addmessage('This includes ' + str(numCorridorLinks) + ' potential corridor links and '+ str(numComponentLinks) + ' component links.')
        elif numCorridorLinks > 0:
            GP.addmessage('This includes ' + str(numCorridorLinks) + ' potential corridor links.')
        else:
            GP.addmessage('\n***NOTE: There are NO corridors to map!')
        dashline(2)
        
    except arcgisscripting.ExecuteError:
        raise_geoproc_error('lm_util')
    except:
        raise_python_error('lm_util')
    return

      
############################################################################    
## Adjacency and allocation functions ##########################
############################################################################

def get_adj_using_shift_method(alloc):
    """Returns table listing adjacent core areas by shifting allocation grid one pixel
    and then looking for pixels with different allocations across shifted grids.
    """
    global GP
    cellSize = GP.Describe(alloc).MeanCellHeight
    GP.CellSize = cellSize
    
    posShift = GP.CellSize
    negShift = -1*float(GP.CellSize)
    
    GP.workspace = SCRATCHDIR
    
    GP.addmessage('Calculating adjacencies crossing horizontal allocation boundaries...')
    startTime=time.clock()
    GP.Shift_management(alloc, "alloc_r", posShift, "0")
       
    alloc_r = "alloc_r"
    adjTable_r = get_allocs_from_shift(GP.workspace,alloc,alloc_r)
    startTime, hours, mins, secs = elapsed_time(startTime)
    
    GP.addmessage('Calculating adjacencies crossing upper-left diagonal allocation boundaries...')
    GP.Shift_management(alloc, "alloc_ul", negShift, posShift)
    
    alloc_ul = "alloc_ul"
    adjTable_ul = get_allocs_from_shift(GP.workspace,alloc,alloc_ul)
    startTime, hours, mins, secs = elapsed_time(startTime)
      
    GP.addmessage('Calculating adjacencies crossing upper-right diagonal allocation boundaries...')
    GP.Shift_management(alloc, "alloc_ur", posShift, posShift)

    alloc_ur = "alloc_ur"
    adjTable_ur  = get_allocs_from_shift(GP.workspace,alloc,alloc_ur)
    startTime, hours, mins, secs = elapsed_time(startTime)
    
    GP.addmessage('Calculating adjacencies crossing vertical allocation boundaries...')
    GP.Shift_management(alloc, "alloc_u", "0", posShift)
    
    alloc_u = "alloc_u"
    adjTable_u  = get_allocs_from_shift(GP.workspace,alloc,alloc_u)
    startTime, hours, mins, secs = elapsed_time(startTime)

    adjTable = combine_adjacency_tables(adjTable_r,adjTable_u,adjTable_ur,adjTable_ul)
        
    return adjTable
    
            
def combine_adjacency_tables(adjTable_r, adjTable_u, adjTable_ur, adjTable_ul):
    """Combines tables describing whether core areas are adjacent based on
    allocation zones that touch on horizontal, vertical, and diagonal axes
    """
    try:
        adjTable = append(adjTable_r,adjTable_u,axis=0)
        adjTable = append(adjTable,adjTable_ur,axis=0)
        adjTable = append(adjTable,adjTable_ul,axis=0)
            
        pairs = sort(adjTable[:,0:2])
        adjTable[:,0:2] = pairs
       
        ind=lexsort((adjTable[:,1],adjTable[:,0])) #sort by 1st core Id then by 2nd core Id
        adjTable = adjTable[ind]    
        
        numDists=len(adjTable)
        x=1
        while x<numDists:
            if adjTable[x,0]==adjTable[x-1,0] and adjTable[x,1]==adjTable[x-1,1]:
                adjTable[x-1,0]=0 #mark for deletion
            x=x+1
        
        if numDists>0:
            delRows=asarray(where(adjTable[:,0]==0))
            delRowsVector = zeros((delRows.shape[1]),dtype="int32")
            delRowsVector[:] = delRows[0,:]
            adjTable=delete_row(adjTable,delRowsVector)        
            del delRows
            del delRowsVector
       
        return adjTable   
    except arcgisscripting.ExecuteError:
        raise_geoproc_error('lm_util')
    except:
        raise_python_error('lm_util')      
    
               
def get_allocs_from_shift(workspace, alloc, alloc_sh):
    """Returns a table of adjacent allocation zones using grid shift method"""
    try:
        combine_ras = os.path.join(GP.workspace, "combine")
        count,statement = 0, 'GP.SingleOutputMapAlgebra_sa("combine("+alloc+","+ alloc_sh+")", combine_ras, alloc, alloc_sh)'
        while True:
            try: exec statement
            except:
                count,tryAgain = hiccup_test(count,statement)
                if not tryAgain: exec statement
            else: break
        allocLookupTable = get_alloc_lookup_table(GP.workspace,combine_ras)        
        return allocLookupTable[:,1:3]

    except arcgisscripting.ExecuteError:
        raise_geoproc_error('lm_util')
    except:
        raise_python_error('lm_util')   


def get_alloc_lookup_table(workspace, combine_ras): # Try this in script 4 too... and anything else with minras...
    """Returns a table of adjacent allocation zones given a raster with allocation zone attributes"""    
    try:
        desc = GP.describe(combine_ras)
        fldlist = GP.listfields(combine_ras)
        
        valFld = fldlist[1].name
        allocFld = fldlist[3].name
        allocFld_sh= fldlist[4].name
        
        allocLookupTable = zeros((0,3),dtype="int32")
        appendRow = zeros((1,3),dtype="int32")
        
        rows = GP.searchcursor(combine_ras)
        row = rows.next()
        while row:
            alloc = row.getvalue(allocFld)
            alloc_sh = row.getvalue(allocFld_sh)
            if alloc != alloc_sh:
                appendRow[0,0] = row.getvalue(valFld)
                appendRow[0,1] = alloc
                appendRow[0,2] = alloc_sh
                allocLookupTable = append(allocLookupTable, appendRow, axis=0)
            row = rows.next()
        del row
        del rows
    
        return allocLookupTable
    except arcgisscripting.ExecuteError:
        raise_geoproc_error('lm_util')
    except:
        raise_python_error('lm_util')   

            
############################################################################    
## Bounding Circle Functions ##########################
############################################################################

def new_extent(fc, field, value):
    """Returns the maximum area extent of features where field == value"""
    try:
        shapeFieldName = GP.describe(fc).shapefieldname
        
        #searchRows = GP.searchcursor(fc, " "'"core_ID"'" = " + str(currentID))
        searchRows = GP.searchcursor(fc, field + ' = ' + str(value))
        searchRow = searchRows.next()
        #get the 1st features extent
        extentObj = searchRow.getvalue(shapeFieldName).extent
        xMin = extentObj.xmin
        yMin = extentObj.ymin
        xMax = extentObj.xmax
        yMax = extentObj.yMax
        searchRow = searchRows.next()#now move on to the other features
        while searchRow:
           extentObj = searchRow.getvalue(shapeFieldName).extent
           if extentObj.xmin < xMin:
              xMin = extentObj.xmin 
           if extentObj.ymin < yMin:
              yMin = extentObj.ymin
           if extentObj.xmax > xMax:
              xMax = extentObj.xmax
           if extentObj.ymax > yMax:
              yMax = extentObj.ymax
           searchRow = searchRows.next()
        del searchRow
        del searchRows
    except arcgisscripting.ExecuteError:
        raise_geoproc_error('lm_util')
    except:
        raise_python_error('lm_util')   
 
    #return [xMin,yMin,xMax,yMax]
    #strExtent = str(xMin)+' '+str(yMin)+' '+str(xMax)+' '+str(yMax)#"%s %s %s %s" % tuple(lstExt)
    #strExtent
    return  str(xMin),str(yMin),str(xMax),str(yMax)


def get_centroids(shapefile, field):
    """Returns centroids of features"""
    try:
        pointArray = zeros((0,3),dtype="float32")
        xyCumArray = zeros((0,3),dtype="float32")
        xyArray = zeros((1,3),dtype="float32")
        rows = GP.SearchCursor(shapefile)
        row = rows.Next()
        while row:
            #list1.append(row.field)
            feat = row.shape
            center = feat.Centroid
            center = str(center)
            xy = center.split(" ")
            xyArray[0,0] = float(xy[0])
            xyArray[0,1] = float(xy[1])
            value = row.GetValue(field)
            xyArray[0,2] = int(value)            
            xyCumArray = append(xyCumArray, xyArray, axis=0)
            row = rows.Next()
        pointArray = append(pointArray,xyCumArray,axis=0)
        
        return pointArray 

    except arcgisscripting.ExecuteError:
        raise_geoproc_error('lm_util')
    except:
        raise_python_error('lm_util')   


def get_bounding_circle_data(extentBoxList, corex, corey, distbuff):
    """Returns centroid and radius of circles bounding the extent boxes of two core areas"""
    try:
        circlePointData = zeros((1,5),dtype='float32')                        
        numBoxes = extentBoxList.shape[0]
        for line in range (numBoxes):  #doing loop because can't get where to work correctly on this list
            if extentBoxList[line,0] == corex:
                corexData = extentBoxList[line,:]
            if extentBoxList[line,0] == corey:
                coreyData = extentBoxList[line,:]
                break
        
        x_ulx = corexData[1] 
        x_lrx = corexData[2]
        x_uly = corexData[3]
        x_lry = corexData[4]
        y_ulx = coreyData[1] 
        y_lrx = coreyData[2]
        y_uly = coreyData[3]
        y_lry = coreyData[4]
        
        xmin = min(x_ulx,y_ulx)
        ymin = min(x_lry,y_lry)       
        xmax = max(x_lrx,y_lrx)
        ymax = max(x_uly,y_uly)

        centX = xmin + (xmax - xmin)/2
        centY = ymin + (ymax - ymin)/2
        radius = sqrt(pow((xmax - xmin)/2,2) + pow((ymax - ymin)/2,2))
        if distbuff != 0:
            radius = radius + int(distbuff)
        circlePointData[0,:] = [centX, centY, corex, corey, radius]
    except arcgisscripting.ExecuteError:
        raise_geoproc_error('lm_util')
    except:
        raise_python_error('lm_util')     
    
    return circlePointData


def get_extent_box_coords(fieldValue=None):
    """Get coordinates of bounding box that contains selected features"""
    try:        
        # get all features, not just where COREFN = fieldValue
        if fieldValue is None:
            fieldValue = 1
            desc = GP.Describe
            extent=desc(FCORES).extent
            lr = extent.lowerright
            ul = extent.upperleft
            ulx=ul.x
            uly=ul.y
            lrx=lr.x
            lry=lr.y
        else:        
            ulx,lry,lrx,uly =  new_extent(FCORES, COREFN, fieldValue)
        
        ulx=float(ulx)
        lrx=float(lrx)
        uly=float(uly)
        lry=float(lry)
        boxData = zeros((1,5),dtype='float32')
        boxData[0,:] = [fieldValue, ulx, lrx, uly, lry]
    except arcgisscripting.ExecuteError:
        raise_geoproc_error('lm_util')
    except:
        raise_python_error('lm_util')
        
    return boxData
    

def make_points(workspace, pointArray, outFC):
    """Creates a shapefile with points specified by coordinates in pointArray
    
       outFC is just the filename, not path. 
       pointArray is x,y,corex,corey,radius
       
    """
    try:
        global GP
        wkspbefore = GP.workspace
        GP.workspace = workspace
        if GP.exists(outFC):
            GP.delete_management(outFC)
        GP.CreateFeatureclass_management(workspace, outFC, "POINT")
        #for field in fieldArray:
        if pointArray.shape[1] > 3:
            GP.addfield(outFC, "corex", "SHORT")
            GP.addfield(outFC, "corey", "SHORT")
            GP.addfield(outFC, "radius", "DOUBLE")
            GP.addfield(outFC, "cores_x_y", "TEXT")
        else:
            GP.addfield(outFC, "XCoord", "DOUBLE")
            GP.addfield(outFC, "YCoord", "DOUBLE")                        
            GP.addfield(outFC, COREFN, "SHORT")
        rows = GP.InsertCursor(outFC)
        
        numPoints = pointArray.shape[0]
        for i in range (numPoints):
            point = GP.CreateObject("Point")
            point.ID = i
            point.X = float(pointArray[i,0])
            point.Y = float(pointArray[i,1])
            row = rows.NewRow()
            row.shape = point
            row.SetValue("ID", i)
            if pointArray.shape[1]>3:            
                row.SetValue("corex", int(pointArray[i,2]))
                row.SetValue("corey", int(pointArray[i,3]))
                row.SetValue("cores_x_y", str(int(pointArray[i,2])) + '_' + str(int(pointArray[i,3])))
                row.SetValue("radius", float(pointArray[i,4]))
            else:
                row.SetValue("XCoord", float(pointArray[i,0]))
                row.SetValue("YCoord", float(pointArray[i,1]))
                row.SetValue(COREFN, float(pointArray[i,2]))

            rows.InsertRow(row)
        del row
        del point
        del rows
        GP.workspace = wkspbefore
    
    except arcgisscripting.ExecuteError:
        raise_geoproc_error('lm_util')
    except:
        raise_python_error('lm_util')
    return 

   
############################################################################    
## LCP Shapefile Functions #################################################
############################################################################

def create_lcp_shapefile(linkTable, sourceCore, targetCore, lcpLoop, SR):
    """Creates lcp shapefile, which shows locations of least-cost path lines attributed with corridor info/status"""
    try:
        startTime = time.clock()                 
        linkIdCol = 0
        linkTypeCol = 5
        eucDistCol = 6
        cwDistCol = 7
        rows = get_links_from_core_pairs(linkTable,sourceCore,targetCore) 
        link = rows[0]
        
        # lcpline = os.path.join(SCRATCHDIR + "lcpline.shp" 
        #                     + str(int(sourceCore)) + "_" 
        #                     + str(int(targetCore)))
        lcpline = os.path.join(SCRATCHDIR, "lcpline.shp")
        GP.RasterToPolyline_conversion("lcp", lcpline, "NODATA", "", 
                                       "NO_SIMPLIFY") 

        lcplineDslv = os.path.join(SCRATCHDIR, "lcplineDslv.shp")
        GP.Dissolve_management(lcpline,lcplineDslv)

        GP.AddField_management(lcplineDslv,"Link_ID","SHORT","5")
        GP.CalculateField_management(lcplineDslv,"Link_ID",
                                     int(linkTable[link,linkIdCol]))

        linkTypeCode = linkTable[link,linkTypeCol]
        activeLink, linkTypeDesc = get_link_type_desc(linkTypeCode)
        
        GP.AddField_management(lcplineDslv,"Active","SHORT")
        GP.CalculateField_management(lcplineDslv,"Active", activeLink)

        GP.AddField_management(lcplineDslv,"Link_Info","TEXT")
        GP.CalculateField_management(lcplineDslv,"Link_Info",linkTypeDesc)

        GP.AddField_management(lcplineDslv,"From_Core","SHORT","5")
        GP.CalculateField_management(lcplineDslv,"From_Core",int(sourceCore))
        GP.AddField_management(lcplineDslv,"To_Core","SHORT","5")
        GP.CalculateField_management(lcplineDslv,"To_Core",int(targetCore))
        
        GP.AddField_management(lcplineDslv,"Euc_Dist","DOUBLE","10","2")
        GP.CalculateField_management(lcplineDslv,"Euc_Dist",linkTable[link,eucDistCol])
        
        GP.AddField_management(lcplineDslv,"CW_Dist","DOUBLE","10","2")
        GP.CalculateField_management(lcplineDslv,"CW_Dist",linkTable[link,cwDistCol])
        GP.AddField_management(lcplineDslv,"LCP_Length","DOUBLE","10","2")

        rows=GP.UpdateCursor(lcplineDslv)
        row = rows.Next()
        while row:
            feat=row.shape
            lcpLength=int(feat.length)
            row.SetValue("LCP_Length",lcpLength)
            rows.UpdateRow(row)
            row = rows.Next()
        del row, rows
        
        distRatio1 = (float(linkTable[link,cwDistCol]) 
                      / float(linkTable[link,eucDistCol]))
        GP.AddField_management(lcplineDslv,"cwd2Euc_R","DOUBLE","10","2")
        GP.CalculateField_management(lcplineDslv,"cwd2Euc_R",distRatio1)

        distRatio2 = float(linkTable[link,cwDistCol]) / float(lcpLength)
        GP.AddField_management(lcplineDslv,"cwd2Path_R","DOUBLE","10","2")
        GP.CalculateField_management(lcplineDslv,"cwd2Path_R",distRatio2)

        lcpLoop=lcpLoop+1
        lcpShapefile = os.path.join(DATAPASSDIR, "lcplines_step3.shp")
        if lcpLoop == 1:
            if GP.Exists(lcpShapefile):
                try:
                    GP.Delete(lcpShapefile)

                except:
                    dashline(1)
                    msg = ('ERROR: Could not remove LCP shapefile ' 
                           + lcpShapefile + '. Is it open in ArcMap?')
                    GP.AddError(msg)
                    exit(1)

            GP.copy_management(lcplineDslv,lcpShapefile)
        else:
            GP.Append_management(lcplineDslv,lcpShapefile,"TEST")
            
        GP.defineprojection(lcpShapefile, SR)

        ## remove below?  Is it worth the time?        
        # logLcpShapefile = os.path.join(LOGDIR, "lcplines_step3.shp")
        # if GP.exists(logLcpShapefile):
            # try:
                # GP.delete_management(logLcpShapefile)
            # except:
                # dashline(1)
                # msg = 'ERROR: Could not remove lcp shapefile from log directory: ' + logLcpShapefile + '. Is it open in ArcMap?'
                # GP.AddError(msg)
                # exit(1)
        # GP.copy_management(lcpShapefile,logLcpShapefile)               
        
        return lcpLoop

    except arcgisscripting.ExecuteError:
        raise_geoproc_error('lm_util')
    except:
        raise_python_error('lm_util')
        

def update_lcp_shapefile(linkTable, lastStep, thisStep):
    """Updates lcp shapefiles with new link information/status"""
    try:
        startTime = time.clock()        
        
        #set linkTable INDEXES
        linkTypeCol = 5
        eucDistCol = 6
        cwDistCol=7
        lcpLengthCol = 10
        cwdToEucRatioCol = 11
        cwdToPathRatioCol = 12
        
        numLinks = linkTable.shape[0]
        extraCols=zeros((numLinks,3),dtype="float64") #g1' g2' THEN c1 c2
        linkTableTemp = append(linkTable,extraCols,axis=1)
        del extraCols

        linkTableTemp[:,lcpLengthCol] = -1
        linkTableTemp[:,cwdToEucRatioCol] = -1
        linkTableTemp[:,cwdToPathRatioCol] = -1

        lcpShapefile = os.path.join(DATAPASSDIR, "lcplines_step" + str(thisStep) 
                        + ".shp")
        if lastStep != thisStep:
            if thisStep == 5:
                oldLcpShapefile = os.path.join(DATAPASSDIR, "lcplines_step" 
                                   + str(lastStep) + ".shp")
                # If last step wasn't step 4 then must be step 3
                if not GP.exists(oldLcpShapefile):
                    oldLcpShapefile = os.path.join(DATAPASSDIR, "lcplines_step" 
                                       + str(lastStep-1) + ".shp") # step 3
            else:         
                oldLcpShapefile = os.path.join(DATAPASSDIR, "lcplines_step" 
                                   + str(lastStep) + ".shp")
                
            if GP.exists(lcpShapefile):
                try:
                    GP.delete_management(lcpShapefile)
                except:
                    dashline(1)
                    msg = ('ERROR: Could not remove LCP shapefile ' 
                           + lcpShapefile + '. Is it open in ArcMap?')
                    GP.AddError(msg)
                    exit(1)   
            GP.copy_management(oldLcpShapefile,lcpShapefile)        


        rows=GP.UpdateCursor(lcpShapefile)
        row = rows.Next()
        line = 0
        while row:
            linkId = row.getvalue("Link_ID")
            linkTypeCode = linkTable[linkId-1,linkTypeCol]
            activeLink, linkTypeDesc = get_link_type_desc(linkTypeCode)
            row.SetValue("Link_Info", linkTypeDesc)
            row.SetValue("Active", activeLink)
            rows.UpdateRow(row)

            linkTableRow = get_linktable_row(linkId, linkTableTemp)
            linkTableTemp[linkTableRow,lcpLengthCol] = row.getvalue("LCP_Length")
            linkTableTemp[linkTableRow,cwdToEucRatioCol] = row.getvalue("cwd2Euc_R")
            linkTableTemp[linkTableRow,cwdToPathRatioCol] = row.getvalue("cwd2Path_R")
            
            row = rows.Next()
            line = line + 1
        # delete cursor and row points to remove locks on the data
        del row, rows
        
        # logLcpShapefile = os.path.join(LOGDIR, "lcplines_step" 
        #                             + str(thisStep) + ".shp"
        # if GP.exists(logLcpShapefile):
            # try:
                # GP.delete_management(logLcpShapefile)
            # except:
                # dashline(1)
                # msg = 'ERROR: Could not remove lcp shapefile from log directory: ' + logLcpShapefile + '. Is it open in ArcMap?'
                # GP.AddError(msg)
                # exit(1)
        # GP.copy_management(lcpShapefile,logLcpShapefile)

        outputLcpShapefile = os.path.join(OUTPUTDIR, "lcplines_step" + str(thisStep) + ".shp")
        if GP.exists(outputLcpShapefile):
            try:
                GP.delete_management(outputLcpShapefile)
            except:
                dashline(1)
                msg = 'ERROR: Could not remove lcp shapefile from output directory: ' + outputLcpShapefile + '. Is it open in ArcMap?'
                GP.AddError(msg)
                exit(1)
        GP.copy_management(lcpShapefile,outputLcpShapefile)
        
        # oldOutputLcpShapefile = os.path.join(OUTPUTDIR, "lcplines_step" + str(lastStep) + ".shp")
        # if GP.exists(oldOutputLcpShapefile):
            # try:
                # GP.delete_management(oldOutputLcpShapefile) #clean up to keep clutter down in output dir
            # except:
                # pass
        
        return linkTableTemp
    
    except arcgisscripting.ExecuteError:
        raise_geoproc_error('lm_util')
    except:
        raise_python_error('lm_util')   

############################################################################    
## Graph Functions ########################################################
############################################################################

def delete_row(A, delrow):
    """Deletes rows from a matrix"""
    try:
        m = A.shape[0]
        n = A.shape[1]
        keeprows = delete (arange(0, m), delrow)
        keepcols = arange(0, n)
    except:
        raise_python_error('lm_util')
    return A[keeprows][:,keepcols]


def delete_col(A, delcol):
    """Deletes columns from a matrix"""
    try:
        m = A.shape[0]
        n = A.shape[1]
        keeprows = arange(0, m)
        keepcols = delete (arange(0, n), delcol)
    except:
        raise_python_error('lm_util')
    return A[keeprows][:,keepcols]


def delete_row_col(A, delrow, delcol):
    """Deletes rows and columns from a matrix"""
    try:
        m = A.shape[0]
        n = A.shape[1]

        keeprows = delete (arange(0, m), delrow)
        keepcols = delete (arange(0, n), delcol)
    except:
        raise_python_error('lm_util')
    return A[keeprows][:,keepcols]
    

def components_no_sparse(G):
    """Returns components of a graph while avoiding use of sparse matrices"""
    #from gapdt.py by Viral Shah
    #        G = sparse.coo_matrix(G) ############
    U,V= where(G)
    
    n = G.shape[0]
    D = arange (0, n, dtype='int32')

    while True:
        D = conditional_hooking(D, U, V)
        star = check_stars (D)
        
        if (sum(star) == n):
            return relabel(D, 1)
            break

        D = pointer_jumping(D)


def relabel( oldlabel, offset=0):#from gapdt.py by Viral Shah
    """Relabels components"""
    newlabel = zeros(size(oldlabel), dtype='int32')
    s = sort(oldlabel)
    perm = argsort(oldlabel)
    f = where(diff(concatenate(([s[0]-1], s))))
    newlabel[f] = 1
    newlabel = cumsum(newlabel)
    newlabel[perm] = copy(newlabel)
    return newlabel-1+offset


def conditional_hooking (D, u, v):#from gapdt.py by Viral Shah
    """Utility for components code"""
    Du = D[u]
    Dv = D[v]
    
    hook = where ((Du == D[Du]) & (Dv < Du))
    Du = Du[hook]
    Dv = Dv[hook]
    
    D[Du] = Dv
    return D


def check_stars (D):#from gapdt.py by Viral Shah
    """Utility for components code"""    
    n = D.size
    star = ones (n, dtype='int32')

    notstars = where (D != D[D])
    star[notstars] = 0
    Dnotstars = D[notstars]
    star[Dnotstars] = 0
    star[D[Dnotstars]] = 0

    star = star[D]
    return star

def pointer_jumping (D):#from gapdt.py by Viral Shah
    """Utility for components code"""    
    n = D.size
    Dold = zeros(n, dtype='int32');

    while any(Dold != D):
        Dold = D
        D = D[D]
    return D
       
    
    
    
############################################################################    
## Input Functions ########################################################
############################################################################

def load_link_table(linkTableFile):
    """Reads link table created by previous step """    
    try:
        linkTable1 = loadtxt(linkTableFile, dtype = 'Float64', comments='#', delimiter=',')
        if len(linkTable1) == linkTable1.size: #Just one connection
            linkTable = zeros((1,len(linkTable1)),dtype='Float64')
            linkTable[:,0:len(linkTable1)] = linkTable1[0:len(linkTable1)]
        else:
            linkTable=linkTable1
        return linkTable
    except:
        raise_python_error('lm_util')   


        
############################################################################    
## Output Functions ########################################################
############################################################################
def write_link_table(linkTable, outlinkTableFile):
    """Writes link tables to pass link data between steps """
    try:
        
        numLinks = linkTable.shape[0]
        outFile = open(outlinkTableFile,"w")
        
        if linkTable.shape[1] == 10:
            outFile.write ("#link,coreId1,coreId2,cluster1,cluster2,linkType,eucDist,lcDist,eucAdj,cwdAdj\n")
        
            for x in range(0,numLinks):
                for y in range(0,9):
                    outFile.write (str(linkTable[x,y]) + "," )
                outFile.write (str(linkTable[x,9]))
                outFile.write ("\n")

        else:
            
            outFile.write ("#link,coreId1,coreId2,cluster1,cluster2,linkType,eucDist,lcDist,eucAdj,cwdAdj,lcpLength,cwdToEucRatio,cwdToPathRatio\n")
        
            for x in range(0,numLinks):
                for y in range(0,12):
                    outFile.write (str(linkTable[x,y]) + "," )
                outFile.write (str(linkTable[x,12]))
                outFile.write ("\n")

        outFile.close()
    except arcgisscripting.ExecuteError:
        raise_geoproc_error('lm_util')
    except:
        raise_python_error('lm_util')
    return

def write_adj_file(outcsvfile, adjTable):
    """Writes csv file listing adjacent core areas to pass adjacency info between steps"""
    outfile = open(outcsvfile,"w")
    outfile.write ("#Edge" + "," + str(COREFN) + "," + str(COREFN) + "_1" + "\n")
    for x in range(0,len(adjTable)):
        outfile.write ( str(x) + "," + str(adjTable[x,0]) + "," + str(adjTable[x,1]) + "\n" )
    outfile.close()         

def write_link_maps(linkTableFile, step):
    """Writes stick maps (aka link maps)
    
    These are vector line files showing links connecting core area pairs. Links
    contain attributes about corridor characteristics.
    
    """
    try:
        global GP
        GP.workspace = OUTPUTDIR

        #########################################################
        ##Relevant linkTable INDEXES
        linkIdCol=0
        core1Col=1 #core ID of 1st core area link connects
        core2Col=2 #core ID of 2nd core area link connects
        # 0=no link. 1=core, 2=corridor, 3=intermediate core area detected, 
        # 4=too long EucDist, 5=too long lcDist
        linkTypeCol=5 
        eucDistCol=6
        cwDistCol=7
        #########################################################
        
        linkTable = load_link_table(linkTableFile)
        #linkTable = loadtxt(linkTableFile, dtype = 'Float64', comments='#', 
        #                    delimiter=',')
        numLinks = linkTable.shape[0]
        GP.toolbox = "management"
        
        coresForLinework = "cores_for_linework.shp"

        #Preferred method to get geometric center
        pointArray = get_centroids(COREFC, COREFN)
                
        make_points(GP.workspace, pointArray, coresForLinework) 
        numLinks=linkTable.shape[0]
        #rows,cols = where(linkTable[:,linkTypeCol:linkTypeCol+1]== 2)

        coreLinks=linkTable
        
        #create coreCoords array, with geographic centers of cores
        coreCoords = zeros(pointArray.shape,dtype='float64')
        coreCoords[:,0]=pointArray[:,2]
        coreCoords[:,1]=pointArray[:,0]
        coreCoords[:,2]=pointArray[:,1]
        
        #Create linkCoords array, which holds following info:
        #linkIdCol core1Col core2Col eucdist lcdist core1x core1y core2x core2y 
            
        linkCoords=zeros((len(coreLinks),10))
        linkCoords[:,0]=coreLinks[:,linkIdCol] #link IDs
        linkCoords[:,1:3]=sort(coreLinks[:,core1Col:core2Col+1]) #core1Col and core2Col IDs, sorted left-to-right
        linkCoords[:,3:5]=coreLinks[:,eucDistCol:cwDistCol+1] #euc and lc distances
        linkCoords[:,9]=coreLinks[:,linkTypeCol] #2=minNN,10=component
        
        if len(coreLinks) > 0:   
            ind=lexsort((linkCoords[:,2],linkCoords[:,1])) #sort by 1st coreId then by 2nd coreId
            linkCoords = linkCoords[ind]    
           
        #Get core coordinates into linkCoords
        #linkCoords arranged by:
        #linkIdCol core1Col core2Col eucdist lcdist core1x core1y core2x core2y
        for i in range(0,len(linkCoords)):
            grp1=linkCoords[i,1]
            grp2=linkCoords[i,2]
            
            for core in range (0,len(coreCoords)):
                if coreCoords[core,0] == grp1:
                    linkCoords[i,5]=coreCoords[core,1]
                    linkCoords[i,6]=coreCoords[core,2]
                elif coreCoords[core,0] == grp2:                
                    linkCoords[i,7]=coreCoords[core,1]
                    linkCoords[i,8]=coreCoords[core,2]

        if len(coreLinks) > 0:  
            ind=argsort((linkCoords[:,0])) #sort by linkIdCol
            linkCoords = linkCoords[ind]
              
        linkTypes = linkTable[:,linkTypeCol]
        
        #
        coreLinksShapefile = 'sticks_step' + str(step) + '.shp'
        
        #make coreLinks.shp using linkCoords table
        # will contain linework between each pair of connected cores
        GP.CreateFeatureclass(GP.workspace, coreLinksShapefile, "POLYLINE")
        
        
        #Define Coordinate System for output shapefiles
        desc = GP.Describe
        SR = desc(COREFC).SpatialReference
        GP.defineprojection(coreLinksShapefile, SR)
        
        #ADD ATTRIBUTES
        GP.AddField_management(coreLinksShapefile, "Link_ID", "SHORT")
        GP.AddField_management(coreLinksShapefile, "Active", "SHORT")
        GP.AddField_management(coreLinksShapefile, "Link_Info", "TEXT")
        GP.AddField_management(coreLinksShapefile, "From_Core", "SHORT")
        GP.AddField_management(coreLinksShapefile, "To_Core", "SHORT")
        GP.AddField_management(coreLinksShapefile, "Euc_Dist", "FLOAT")
        GP.AddField_management(coreLinksShapefile, "CW_Dist", "FLOAT")
        GP.AddField_management(coreLinksShapefile, "cwd2Euc_R", "FLOAT")                        
        #
                
        #Create an Array and Point object.
        lineArray = GP.CreateObject("Array")
        pnt = GP.CreateObject("Point")
               
        #linkCoords indices:
        #linkIdCol core1Col core2Col eucdist lcdist core1x core1y core2x core2y 
        numLinks=len(linkCoords)

        #Open a cursor to insert rows into the shapefile.
        cur = GP.InsertCursor(coreLinksShapefile)

        ##Loop through each record in linkCoords table
        for i in range(0,numLinks):
        
            #Set the X and Y coordinates for origin vertex.
            pnt.x = linkCoords[i,5]
            pnt.y = linkCoords[i,6]
            #Insert it into the line array    
            lineArray.add(pnt)
        
            #Set the X and Y coordinates for destination vertex
            pnt.x = linkCoords[i,7]
            pnt.y = linkCoords[i,8]
            #Insert it into the line array
            lineArray.add(pnt)
        
            #Insert the new poly into the feature class.
            feature = cur.NewRow()
            feature.shape = lineArray
            cur.InsertRow(feature)
           
            lineArray.RemoveAll()
        del cur
                    
        #Add attribute data to link shapefile
        rows=GP.UpdateCursor(coreLinksShapefile)
        row = rows.Next()
        line = 0
        while row:
            #linkCoords indices:
            #linkIdCol core1Col core2Col eucdist lcdist core1x core1y core2x core2y
            row.SetValue("Link_ID", linkCoords[line,0])
            if linkCoords[line,9] == 2:
                row.SetValue("Link_Info", "Group_Pair")
            linkTypeCode = linkCoords[line,9]
            activeLink, linkTypeDesc = get_link_type_desc(linkTypeCode)
            row.SetValue("Active", activeLink)
            row.SetValue("Link_Info", linkTypeDesc)

            row.SetValue("From_Core", linkCoords[line,1])
            row.SetValue("To_Core", linkCoords[line,2])
            row.SetValue("Euc_Dist", linkCoords[line,3])
            row.SetValue("CW_Dist", linkCoords[line,4])
            distRatio1 = float(linkCoords[line,4])/float(linkCoords[line,3])
            if linkCoords[line,4] <=0 or linkCoords[line,3] <= 0:
                row.SetValue("cwd2Euc_R", -1)
            else:
                row.SetValue("cwd2Euc_R", linkCoords[line,4]/linkCoords[line,3])

            rows.UpdateRow(row)
            row = rows.Next()
            line = line + 1

        del row, rows
             
        #clean up temp files
        GP.delete_management(coresForLinework)     
        
        return
    
    except arcgisscripting.ExecuteError:
        raise_geoproc_error('lm_util')
    except:
        raise_python_error('lm_util')
      
        
############################################################################    
## File and Path Management Functions ######################################
############################################################################


def archive_datapass():    
    if os.path.exists(DATAPASSARCHDIR):
        shutil.rmtree(DATAPASSARCHDIR)
    if os.path.exists(DATAPASSDIR):
        shutil.copytree(DATAPASSDIR, DATAPASSARCHDIR)
    return
    
    
def get_cwd_path(core):
    """Returns the path for the cwd raster corresponding to a core area """          
    dirCount = int(core / 100)
    if dirCount > 0:    
        return os.path.join(CWDBASEDIR, CWDSUBDIR_MAME + str(dirCount), 
                         "cwd_" + str(core))
    else:
        return os.path.join(CWDBASEDIR, CWDSUBDIR_MAME, "cwd_" 
                         + str(core))


def check_project_dir():
    """Checks to make sure path name is not too long (can cause problems with ESRI grids) """
    if len(PROJECTDIR) > 100:
        msg = ('ERROR: Project directory "' + PROJECTDIR 
               + '" is too deep.  Please choose a shallow directory'
               '(something like "C:\ANBO").')
        GP.AddError(msg)
        GP.AddMessage(GP.GetMessages(2))
        exit(1) 
    return


def get_prev_step_link_table(step):
    """Returns the name of the link table created by the previous step"""
    try:
        prevStep = step - 1
        
        if step == 5:
            prevStepLinkTable = os.path.join(DATAPASSDIR, 
                                          'linkTable_step4.csv')
            GP.addmessage('\nLooking for '+ DATAPASSDIR
                          + 'linkTable_step4.csv')
            if os.path.exists(prevStepLinkTable) == True:
                return prevStepLinkTable
            else:
                prevStep = 3 #Can skip step 4
        
        prevStepLinkTable = os.path.join(DATAPASSDIR, 'linkTable_step' 
                                      + str(prevStep) + '.csv')
        GP.addmessage('\nLooking for '+ DATAPASSDIR + 'linkTable_step' 
                      + str(prevStep) + '.csv')
        if os.path.exists(prevStepLinkTable) == True:
            return prevStepLinkTable 
        else:
            msg = ('\nERROR: Could not find a linktable from step previous to' 
                   'step #' + str(step) + ' in datapass directory.  See above' 
                   'for valid linkTable files.')
            GP.AddError(msg)
            exit(1)
        
    except arcgisscripting.ExecuteError:
        raise_geoproc_error('lm_util')
    except:
        raise_python_error('lm_util') 


def get_this_step_link_table(step):
    """Returns name of link table to write for current step"""
    try:  
        filename = os.path.join(DATAPASSDIR, 'linkTable_step' + str(step) 
                             + '.csv')
        return filename
    
    except arcgisscripting.ExecuteError:
        raise_geoproc_error('lm_util')
    except:
        raise_python_error('lm_util') 

    
def clean_up_link_tables(step):
    """Remove link tables from previous runs."""
    try:
        for stepNum in range(step,7):
            filename = os.path.join(DATAPASSDIR, 'linkTable_step' 
                                 + str(stepNum) + '.csv')
            if os.path.isfile(filename):
                os.remove(filename)

        filename = os.path.join(OUTPUTDIR, 'linkTable_final.csv')
        if os.path.isfile(filename)==True:
            os.remove(filename)
        # newfilename = os.path.join(OUTPUTDIR, 'old_results', 
        #                         'linkTable_final.csv')
        # shutil.move(filename, newfilename)
            
    except arcgisscripting.ExecuteError:
        raise_geoproc_error('lm_util')
    except:
        raise_python_error('lm_util') 


def copy_final_link_maps():
    """Copies final link maps from datapass to the output directory """
    try:
        step = 5 #This is the number of the final step                      
        
        coreLinksShapefile = os.path.join(OUTPUTDIR, 'sticks_step' + str(step)
                                       + '.shp')
        lcpShapefile = os.path.join(DATAPASSDIR, "lcplines_step" + str(step) 
                                 + ".shp")
        
        if GP.exists(coreLinksShapefile):
            GP.MakeFeatureLayer(coreLinksShapefile,"flinks")
            field = "Active"
            expression = field + " = " + str(1)
            GP.selectlayerbyattribute("flinks", "NEW_SELECTION", expression)
                  
            activeLinksShapefile= os.path.join(OUTPUTGDB, 
                                            'Active_Sticks')
            GP.CopyFeatures_management("flinks", activeLinksShapefile)

            expression = field + " = " + str(0)
            GP.selectlayerbyattribute("flinks", "NEW_SELECTION", expression)
            inActiveLinksShapefile= os.path.join(OUTPUTGDB, 
                                              'Inactive_Sticks')
            GP.CopyFeatures_management("flinks",inActiveLinksShapefile)
            
            #move_map(coreLinksShapefile,activeLinksShapefile)
            #copy_map(activeLinksShapefile,inActiveLinksShapefile)


        if GP.exists(lcpShapefile):
            GP.MakeFeatureLayer(lcpShapefile,"flcp")
            field = "Active"
            expression = field + " = " + str(1)
            GP.selectlayerbyattribute("flcp", "NEW_SELECTION", expression)

            activeLcpShapefile= os.path.join(OUTPUTGDB, 'Active_LCPs')
            GP.CopyFeatures_management("flcp", activeLcpShapefile)

            expression = field + " = " + str(0)
            GP.selectlayerbyattribute("flcp", "NEW_SELECTION", expression)
            inActiveLcpShapefile= os.path.join(OUTPUTGDB, 
                                            'Inactive_LCPs')
            GP.CopyFeatures_management("flcp",inActiveLcpShapefile)

        # Move stick and lcp maps for each step to log directory to reduce clutter in output
        for i in range(2,6):
            oldLinkFile = os.path.join(OUTPUTDIR, 'sticks_step' + str(i) 
                                    + '.shp')
            logLinkFile = os.path.join(LOGDIR, 'sticks_step' + str(i) 
                                    + '.shp')
            if GP.exists(oldLinkFile):
                try:
                    move_map(oldLinkFile,logLinkFile)
                except:
                    pass
            oldLcpShapeFile = os.path.join(OUTPUTDIR, 'lcplines_step' + str(i)
                                        + '.shp')
            logLcpShapeFile = os.path.join(LOGDIR, 'lcplines_step' + str(i) 
                                        + '.shp')
            if GP.exists(oldLcpShapeFile):
                try:
                    move_map(oldLcpShapeFile,logLcpShapeFile)
                except:
                    pass       
        return
    except arcgisscripting.ExecuteError:
        raise_geoproc_error('lm_util')
    except:
        raise_python_error('lm_util') 
       
        
def move_map(oldMap, newMap):
    """Moves a map to a new location """
    if GP.exists(oldMap):
        if GP.exists(newMap):
            try:
                GP.delete_management(newMap)
            except:
                pass
        try:
            GP.CopyFeatures_management (oldMap, newMap)
            GP.delete_management(oldMap)
        except:
            pass
    return


############################################################################    
##Error Checking and Handling Functions ####################################
############################################################################

def print_conefor_warning():
    GP.addmessage('\nWARNING: At least one potential link was dropped because')
    GP.addmessage('there was no Euclidean distance value in the input Euclidean')
    GP.addmessage('distance file from Conefor extension.\n')
    GP.addmessage('This may just mean that there were core areas that were adjacent')
    GP.addmessage('but were farther apart than the optional maximum distance used ')
    GP.addmessage('when running Conefor.  But it can also mean that distances were')
    GP.addmessage('calculated using a different core area shapefile or the wrong field')
    GP.addmessage('in the same core area shapefile.\n')

    
def check_steps():
    """Check to make sure there are no skipped steps in a sequence of chosen 
    steps (except step 4 which is optional)
    
    """
    skipStep = False       
    if STEP1 and not STEP2 and STEP3: 
        skipStep = True    
    if STEP2 and not STEP3 and (STEP4 or STEP5):      
        skipStep = True           
    if skipStep:
        try:
            dashline(1)
            msg = ("Error: You can start or stop at different steps, but you "
                   "can't SKIP any except for step 4.\n")
            GP.AddError(msg)
            exit(0)                    
        except:
            raise_python_error('lm_util') 
    return

        
def check_dist_file():
    """Checks for Euclidean distance file from Conefor Inputs tool"""
    # Text file from conefor sensinode extension of edge-edge distances between core area pairs  
    if not os.path.exists(S2EUCDISTFILE):
        GP.AddMessage('\nERROR: Euclidean distance file not found: ' + S2EUCDISTFILE)
        exit(0)
    return


def hiccup_test(count, statement):
    """Re-tries some ArcGIS calls in case of server problems or 'other hiccups'"""
    try:
        if count < 10:
            sleepTime = 10 * count
            count = count + 1
            dashline(1)
            if count == 1:
                GP.addmessage('Failed to execute ' + statement + ' on try #' + str(count) + '.\n Could be an ArcGIS hiccup.')  
                dashline(2)
                GP.addmessage("Here's the error being reported: ")    
                import traceback
                tb = sys.exc_info()[2]  # get the traceback object
                tbinfo = traceback.format_tb(tb)[0]      # tbinfo contains the error's line number and the code
                line = tbinfo.split(", ")[1]

                err = traceback.format_exc().splitlines()[-1]

                for msg in range(0, GP.MessageCount):
                    if GP.GetSeverity(msg) == 2:
                        GP.AddReturnMessage(msg)
                    print GP.AddReturnMessage(msg)
                    dashline(2)
            else:
                GP.addmessage('Failed again executing ' + statement + ' on try #' + str(count) + '.\n Could be an ArcGIS hiccup- scroll up for error description.')  

                GP.addmessage('---------Trying again in ' + str(sleepTime) + ' seconds---------\n')
            time.sleep(sleepTime)
            return count,True
        else:
            sleepTime = 300
            count = count + 1
            dashline(1)
            GP.addmessage('Failed to execute ' + statement + ' on try #' + str(count) + '.\n Could be an ArcGIS hiccup.  Trying again in 5 minutes.\n')
            time.sleep(sleepTime)
            return count,True
    except:
        raise_python_error('lm_util') 

def raise_geoproc_error(filename):
    """Handle geoprocessor errors and provide details to user"""
    dashline(1)
    tb = sys.exc_info()[2]  # get the traceback object
    tbinfo = traceback.format_tb(tb)[0]      # tbinfo contains the error's line number and the code
    line = tbinfo.split(", ")[1]
    #filename = sys.argv[0]
    err = traceback.format_exc().splitlines()[-1]
    
    GP.AddError("Geoprocessing error on **" + line + "** of " + filename + " :")
    dashline(1)         
    for msg in range(0, GP.MessageCount):
        if GP.GetSeverity(msg) == 2:
            GP.AddReturnMessage(msg)
        dashline(2)
        print GP.AddReturnMessage(msg)
        dashline(2)
    exit(0)
        
def raise_python_error(filename):
    """Handle python errors and provide details to user"""
    dashline(1)  
    tb = sys.exc_info()[2]  # get the traceback object
    tbinfo = traceback.format_tb(tb)[0]      # tbinfo contains the error's line number and the code
    line = tbinfo.split(", ")[1]
    #filename = sys.argv[0]
    err = traceback.format_exc().splitlines()[-1]

    GP.AddError("Python error on **" + line + "** of " + filename)
    GP.AddError(err)
    dashline(2)
    exit(0)
    
def dashline(lspace=0):
    """Output dashed line in tool output dialog.
        
       0 = No empty line 
       1 = Empty line before 
       2 = Empty line after
        
    """
    if lspace == 1:
        GP.addmessage('\n')
    GP.addmessage('---------------------------------')
    if lspace == 2:
        GP.addmessage('\n')