import arcgisscripting, sys, time
from time import localtime, strftime
import os, string,csv    
from string import split
from watools_util import *

Version = '2010.10.18'

try:
    from numpy import *
except:    
    gp.AddError("Numpy is not installed.  Please get the Python 2.5-compatible version from http://sourceforge.net/projects/numpy/files/NumPy/1.4.1/numpy-1.4.1-win32-superpack-python2.5.exe/download ") 
    sys.exit(1) # end script process

from watools_util import *

gp = arcgisscripting.create(9.3)
gp.OverwriteOutput = 1

# need to have a Spatial Analyst extension
try:
    gp.CheckOutExtension("Spatial")
except RuntimeError:
    gp.AddError("You need to have a Spatial Analyst extension for this tool")
    sys.exit(1) # end script process

try:
    
    linkIdCol=0
    linkTypeCol = 5

    narg = len(sys.argv)
    if narg == 1: #if run from Python and not ArcGIS.
        projectDir = "C:\\WATOOLS\\Demo\\project"
    else:
        projectDir=gp.GetParameterAsText(0)

    gp.addmessage('\n---------------------------------')
    gp.addmessage('Running script updateLCPs.py')

    outputDir = projectDir + "\\output"
    linkTableFile=outputDir + '\\linkTable_Final.csv'

    linkTable = loadLinkTable(linkTableFile)
    
    linkTypeCol = 5
    eucDistCol = 6
    cwDistCol=7
    lcpLengthCol = 10
    cwdToEucRatioCol = 11
    cwdToPathRatioCol = 12
        
        
    if linkTable.shape[1] == 10:
        numLinks = linkTable.shape[0]
        extraCols=zeros((numLinks,3),dtype="float64") #g1' g2' THEN c1 c2
        linkTableTemp = append(linkTable,extraCols,axis=1)
        del extraCols
    
        linkTableTemp[:,lcpLengthCol] = -1
        linkTableTemp[:,cwdToEucRatioCol] = -1
        linkTableTemp[:,cwdToPathRatioCol] = -1
    else:
        linkTableTemp=linkTable
        

    for i in range(0,2):
        if i==0:
            lcpShapefile= outputDir + '\\Active_LCPs.shp'
        else:
            lcpShapefile= outputDir + '\\Inactive_LCPs.shp'
            
        #Add new metrics if lcp files are from old version
        field="LCP_Length"
        fields = gp.ListFields(lcpShapefile, field)
        if len(fields) == 0: #If created by older version
            gp.AddField_management(lcpShapefile,"LCP_Length","DOUBLE","10","2")
        field="Cwd2Euc_R"
        fields = gp.ListFields(lcpShapefile, field)
        if len(fields) == 0: 
            gp.AddField_management(lcpShapefile,"Cwd2Euc_R","DOUBLE","10","2")
        field="Cwd2Path_R"
        fields = gp.ListFields(lcpShapefile, field)
        if len(fields) == 0: 
            gp.AddField_management(lcpShapefile,"Cwd2Path_R","DOUBLE","10","2")
        
        field="Dist_Ratio"
        fields = gp.ListFields(lcpShapefile, field)
        if len(fields) != 0: 
            gp.deletefield (lcpShapefile, field)

        field="DistRatio1"
        fields = gp.ListFields(lcpShapefile, field)
        if len(fields) != 0: 
            gp.deletefield (lcpShapefile, field)

        field="DistRatio2"
        fields = gp.ListFields(lcpShapefile, field)
        if len(fields) != 0: 
            gp.deletefield (lcpShapefile, field)
    
        rows=gp.UpdateCursor(lcpShapefile)
        row = rows.Next()
        while row:
            feat=row.shape
            lcpLength=int(feat.length)
            row.SetValue("LCP_Length",lcpLength)
            eucDist = row.GetValue("Euc_Dist")               
            cwDist = row.GetValue("CW_Dist")
            cwdToEucRatio = float(cwDist) / float(eucDist)
            row.SetValue("Cwd2Euc_R",cwdToEucRatio)
            cwdToPathRatio = float(cwDist) / float(lcpLength)
            row.SetValue("Cwd2Path_R",cwdToPathRatio)
            
            rows.UpdateRow(row)
            row = rows.Next()
        del row, rows
    
        #Get info into linkTable
        rows=gp.UpdateCursor(lcpShapefile)
        row = rows.Next()
        line = 0
        while row:
            linkId = row.getvalue("Link_ID")
            linkTypeCode = linkTable[linkId-1,linkTypeCol]
            activeLink, linkTypeDesc = getLinkTypeDesc(linkTypeCode)
            row.SetValue("Link_Info", linkTypeDesc)
            row.SetValue("Active", activeLink)
            rows.UpdateRow(row)
    
            linkTableTemp[linkId-1,lcpLengthCol] = row.getvalue("LCP_Length")
            linkTableTemp[linkId-1,cwdToEucRatioCol] = row.getvalue("Cwd2Euc_R")
            linkTableTemp[linkId-1,cwdToPathRatioCol] = row.getvalue("Cwd2Path_R")
            
            row = rows.Next()
            line = line + 1
        # delete cursor and row points to remove locks on the data
        del row, rows
    
        gp.addmessage('\nFinal LCP shapfile updated '+lcpShapefile)
    
    
    coreIds = "core_ID"
    settingsHistory = 'dummy'
    writelinkTable(linkTableTemp,coreIds,linkTableFile,settingsHistory)
    gp.addmessage('\nUpdated link table written to '+linkTableFile)

    #Pull out corridor and constellation links
    numLinks=linkTableTemp.shape[0]
    rows,cols = where(linkTableTemp[:,linkTypeCol:linkTypeCol+1]== 2)
    coreLinks=linkTableTemp[rows,:]
    rows,cols=where(linkTableTemp[:,linkTypeCol:linkTypeCol+1]== 10)
    componentLinks=linkTableTemp[rows,:]
    activeLinkTable=append(coreLinks,componentLinks,axis=0)        
    del componentLinks
    
    ind=argsort((activeLinkTable[:,linkIdCol])) #sort by linkIdCol
    activeLinkTable = activeLinkTable[ind]
    
    activeLinkTableFile=projectDir+"\\output\\"+"linkTable_Final_Active_Links_Only.csv"   
    writelinkTable(activeLinkTable,coreIds,activeLinkTableFile,settingsHistory)   
    gp.addmessage('\nUpdated table of active links written to '+activeLinkTableFile)

    gp.addmessage('\nDONE!\n')
    
# Return GEOPROCESSING specific errors
except arcgisscripting.ExecuteError:
    filename =  'updateLCPs.py'
    raiseGeoprocError(filename)

# Return any PYTHON or system specific errors
except:
    filename =  'updateLCPs.py'
    raisePythonError(filename)
            
del gp

