##*****************************************************************
## 2011_0128
## NAME: s2_buildNetwork.py
##
## SUMMARY: Generates initial version of linkTable.csv based on euclidean distances and adjacencies of core areas
##
## SOFTWARE: ArcGIS 9.3 (requires Spatial Analyst extension)
##           Python 2.5
##
##*****************************************************************

# import required modules
import arcgisscripting, sys, time
from watools_util import *

from numpy import *        

def step2_build_network(gp,Version,options):
    """Generates initial version of linkTable.csv based on euclidean distances
    and adjacencies of core areas.
    
    """
    try:
        # ------------------------------------------------------------------            
        # Unpack options
        projectDir=options["projectDir"]
        coreShapefile=options["coreShapefile"]
        coreIds=options["coreIds"]
        useMinEucDist=options["useMinEucDist"]
        useMaxEucDist=options["useMaxEucDist"]
        maxEucDist=options["maxEucDist"]
        minEucDist=options["minEucDist"]
        keepLinkCriteria=options["keepLinkCriteria"]
        defaultNullValue = options["defaultNullValue"]		
        coreAreaIds = coreIds # We now set variable coreAreaIds field name to user input ID field name                
        # ------------------------------------------------------------------        
        
        outlinkTableFile = get_this_step_link_table(projectDir,step=2)
        
        # ------------------------------------------------------------------
        ## linkTable column numbers
        linkIdCol=0 # Link ID
        core1Col=1 # core ID of 1st core area link connects
        core2Col=2 # core ID of 2nd core area link connects
        cluster1Col=3
        cluster2Col=4
        linkTypeCol=5 # 0=no link. 2=corridor, 3=intermediate core area detected (so do not map)
        eucDistCol=6
        cwDistCol=7
        eucAdjCol=8
        cwdAdjCol=9
        # ------------------------------------------------------------------

        dropFlag=False  # This is a warning flag if distances are mising in conefor 

        # Folders and paths and workspace
        outputRootDir = projectDir + "\\output\\" 
        dataPassDir = projectDir + "\\dataPass\\"
        
        # ------------------------------------------------------------------        
        # Load text file from conefor sensinode extension of edge-edge distances between core area pairs
        dir, file= os.path.split(coreShapefile)
        base, extension= os.path.splitext(file)
        eucDistFile=dataPassDir + "distances_" + base + ".txt" # base is same base name as core area shapefile.
        dummyDistFile=dataPassDir + "distances_" + base + "_2.txt" # if Conefor tried to write to a directory that already had a distances file
        if os.path.exists(dummyDistFile): 
            gp.addmessage('WARNING: found file ' + dummyDistFile + ' which means Conefor extension tried to write to a directory')
            gp.addmessage('that already had a euclidean distance file.')
            gp.addmessage('If you are sure that ' + eucDistFile + 'is correct, then proceed.\n')
            raw_input('Hit any key to proceed, or click cancel to bail.\n') # Fixme: raw_input does not work when called from ArcMap   
        # ------------------------------------------------------------------        

        
        # ------------------------------------------------------------------                
        # Load euclidean adjacency file if needed
        if keepLinkCriteria == "Keep all adjacent links" or keepLinkCriteria == "Keep Euclidean adjacent links":
            eucAdjFile= dataPassDir + "eucAdj.csv" # adjacency file created from 2getAdjacencies.py
            if os.path.exists(eucAdjFile) == False: 
                gp.AddMessage('\nERROR: Euclidean adjacency file required to keep Euclidean adjacent links: ' + eucAdjFile)
                exit(0)        
        else:
            eucAdjFile = defaultNullValue
        # ------------------------------------------------------------------        

        
        # ------------------------------------------------------------------                
        # Load cwd adjacency file if needed
        if keepLinkCriteria == "Keep all adjacent links" or keepLinkCriteria == "Keep cost-weighted adjacent links":
            cwdAdjFile= dataPassDir + "cwdAdj.csv" # adjacency file created from 2getAdjacencies.py
            if os.path.exists(cwdAdjFile) == False:  
                gp.AddMessage('\nERROR: Cost-weighted adjacency file required to keep CWD adjacent links: ' + cwdAdjFile)
                exit(0)        
        else:
            cwdAdjFile = defaultNullValue
        # ------------------------------------------------------------------        
        
    
        #------------------------------------------------------------------------------
        # Load eucDists matrix from file and sort
        eucDists = loadtxt(eucDistFile, dtype = 'Float64', comments='#')
        numDists = eucDists.shape[0]
        gp.addmessage('\n---------------------')
        gp.addmessage('Core area distance list from Conefor Sensinode loaded.  \n')
        gp.addmessage('number of pairwise distances = ' + str(numDists))
        gp.addmessage('---------------------\n')
        eucDists[:,0:2]=sort(eucDists[:,0:2])
        
        # sort eucDists by 1st column then by 2nd then by 3rd    
        ind=lexsort((eucDists[:,2],eucDists[:,1],eucDists[:,0])) 
        eucDists = eucDists[ind]    
        # ------------------------------------------------------------------        

        
        #------------------------------------------------------------------------------
        # Get rid of duplicate pairs of cores, retaining MINIMUM distance between them
        numDistsOld=numDists
        for x in range(numDists-2,-1,-1):
            if eucDists[x,0]==eucDists[x+1,0] and (eucDists[x,1]==eucDists[x+1,1]):
                eucDists[x+1,0]=0
        delRows=asarray(where(eucDists[:,0]==0))
        delRowsVector = zeros((delRows.shape[1]),dtype="int32")
        delRowsVector[:] = delRows[0,:]
        eucDists=delete_row(eucDists,delRowsVector)        
        del delRows
        del delRowsVector    
        numDists = eucDists.shape[0]
        gp.addmessage('\n---------------------')
        gp.addmessage('Removed '+str(numDistsOld - numDists) + ' duplicate core pairs in Euclidean distance table.\n')                      
        maxEucDistId=max(eucDists[:,1])
        gp.addmessage('After removing duplicates and distances that exceed maximum, \nthere are ' + str(numDists) + ' pairwise distances.  Max core ID number is ' + str(int(maxEucDistId)) + '.')
        gp.addmessage('---------------------\n')
        #------------------------------------------------------------------------------

        
        # Begin creating and manipulating linktables
        distlinkTable = zeros((len(eucDists),10)) # zeros and many other array functions are imported from numpy 
        distlinkTable[:,1:3]=eucDists[:,0:2] # this kind of indexing is from numpy too
        distlinkTable[:,eucDistCol]=eucDists[:,2] # eucDistCol is just a number from the index table above.  It is just used to specify the column where euclidean distances are stored
        
        
        #------------------------------------------------------------------------------
        # Get adjacencies using adj files from step 1.
        if cwdAdjFile!= defaultNullValue or eucAdjFile!= defaultNullValue:  
            if cwdAdjFile!= defaultNullValue:
                adjList = loadtxt(cwdAdjFile, dtype = 'int32', comments='#', delimiter=',') # creates a numpy array
                
                if len(adjList) == adjList.size: # Just one connection
                    cwdAdjList = zeros((1,3),dtype='int32')
                    cwdAdjList[:,0:3] = adjList[0:3]                   
                else:
                    cwdAdjList=adjList
                cwdAdjList =cwdAdjList[:,1:3] # Drop first column
                cwdAdjList = sort(cwdAdjList)
                gp.addmessage('Cost-weighted adjacency file loaded.')
                maxCwdAdjCoreID=max(cwdAdjList[:,1])
            else:
                maxCwdAdjCoreID=0        
                
            if eucAdjFile!= defaultNullValue:
                adjList = loadtxt(eucAdjFile, dtype = 'int32', comments='#', delimiter=',') # creates a numpy array

                if len(adjList) == adjList.size: # Just one connection
                    eucAdjList = zeros((1,3),dtype='int32')
                    eucAdjList[:,0:3] = adjList[0:3]                   
                else:
                    eucAdjList=adjList
                eucAdjList=eucAdjList[:,1:3] # Drop first column
                eucAdjList = sort(eucAdjList)
                gp.addmessage('Euclidean adjacency file loaded')
                maxEucAdjCoreID=max(eucAdjList[:,1])
            else:
                maxEucAdjCoreID=0
            gp.addmessage('-------------------------\n')
    
            maxCoreId = max(maxEucAdjCoreID,maxCwdAdjCoreID,maxEucDistId)
            
            
            # FIXME: consider using a lookup table to reduce size of matrix when there are gaps in coreIDs.
            if cwdAdjFile!= defaultNullValue:
                cwdAdjMatrix=zeros((maxCoreId+1,maxCoreId+1),dtype='int32') 
                for x in range(0,len(cwdAdjList)):
                    cwdAdjMatrix[cwdAdjList[x,0],cwdAdjList[x,1]]=1
                cwdAdjMatrix[0,:]=0 # 0 values for core Ids are invalid
                cwdAdjMatrix[0,:]=0
    
            if eucAdjFile!= defaultNullValue:
                eucAdjMatrix=zeros((maxCoreId+1,maxCoreId+1),dtype='int32')
                for x in range(0,len(eucAdjList)):
                    eucAdjMatrix[eucAdjList[x,0],eucAdjList[x,1]]=1
                eucAdjMatrix[0,:]=0 # 0 values for core Ids are invalid
                eucAdjMatrix[0,:]=0
        
            distanceMatrix=zeros((maxCoreId+1,maxCoreId+1),dtype='int32')
            for x in range(0,len(eucDists)):
                distanceMatrix[eucDists[x,0],eucDists[x,1]]=int(eucDists[x,2])
############################# KNOCKING THIS OUT IN ORDER TO INCLUDE ALL LINKS, INCLUDING NON ADJACENT ONES ##############    
####            if keepLinkCriteria == "Keep cost-weighted adjacent links":
####                difference=where(distanceMatrix,1,0)
####                difference=difference-cwdAdjMatrix
####                if amin(difference) < 0:
####                    dropFlag=True
####                distanceMatrix = multiply(distanceMatrix,cwdAdjMatrix) # Drop anything not cwd adjacent
####                
####            elif keepLinkCriteria == "Keep Euclidean adjacent links":
####                difference=where(distanceMatrix,1,0)
####                difference=difference-eucAdjMatrix
####                if amin(difference) < 0:
####                    dropFlag=True
####                del difference           
####                distanceMatrix = multiply(distanceMatrix,eucAdjMatrix) # Drop anything not euc adjacent
####                
####            else: # "Keep all adjacent links"
####                adjMatrix = eucAdjMatrix
####                adjMatrix = adjMatrix + cwdAdjMatrix
####                adjMatrix =where(adjMatrix==2,1,adjMatrix)
####                difference=where(distanceMatrix,1,0)
####                difference=difference-adjMatrix
####                
####                if amin(difference) < 0:
####                    dropFlag=True
####       
####                del difference
####                distanceMatrix = multiply(distanceMatrix,adjMatrix) # Drop anything not adjacent                
        #------------------------------------------------------------------------------            

        
        #------------------------------------------------------------------------------    
        # OK, we have distance matrix (which now defines which pairs can be potential links).  Use it to create link table.   
        gp.addmessage('creating link table')
        distanceMatrix = delete_row_col(distanceMatrix,0,0) # Get rid of 0 index- we don't have any valid core ids with 0 values
        rows,cols=where(distanceMatrix)
        linkTable = zeros((len(rows),10),dtype='int32')
    
        for x in range(0,len(rows)):
            linkTable[x,core1Col]=rows[x]+1
            linkTable[x,core2Col]=cols[x]+1
            linkTable[x,eucDistCol]=distanceMatrix[rows[x],cols[x]]
        del distanceMatrix
        
        if eucAdjFile!= defaultNullValue:
            eucAdjMatrix = delete_row_col(eucAdjMatrix,0,0) # Get rid of 0 index- we don't have any valid core ids with 0 values
            for x in range(0,len(rows)):
                linkTable[x,eucAdjCol]=eucAdjMatrix[rows[x],cols[x]]
            del eucAdjMatrix
    
        if cwdAdjFile!= defaultNullValue:
            cwdAdjMatrix = delete_row_col(cwdAdjMatrix,0,0) # Get rid of 0 index- we don't have any valid core ids with 0 values
            for x in range(0,len(rows)):
                linkTable[x,cwdAdjCol]=cwdAdjMatrix[rows[x],cols[x]]
            del cwdAdjMatrix
    
        if dropFlag==True:
            gp.addmessage('\n---------------------')        
            gp.addmessage('NOTE: At least one adjacent link was dropped because there was no Euclidean ')
            gp.addmessage('distance value in the input distance file from Conefor extension.')
            gp.addmessage('---------------------\n')        

        # Get list of core IDs, based on core area shapefile.
        coreList=get_core_list(coreShapefile,coreIds,coreAreaIds)
        numCores = len(coreList)    
        
        linkTable[:,cluster1Col]=-1 # No clusters until later steps
        linkTable[:,cluster2Col]=-1
        
        if cwdAdjFile == defaultNullValue:
            linkTable[:,cwdAdjCol]=-1 #Euc adjacency not evaluated
        if eucAdjFile == defaultNullValue:
            linkTable[:,eucAdjCol]=-1 # Cost-weighted adjacency not evaluated
        linkTable[:,cwDistCol]=-1 # not evaluated yet. May eventually have ability to get lcdistances for adjacent cores from getAdjacencies.py
        
        # Update linkTable with new core IDs
        numCores = len(coreList)
        for core in range(numCores):
            if coreList[core,0]!=coreList[core,1]:
                linkTable[:,core1Col]=where(linkTable[:,core1Col] == coreList[core,0],coreList[core,1],linkTable[:,core1Col])
                linkTable[:,core2Col]=where(linkTable[:,core2Col] == coreList[core,0],coreList[core,1],linkTable[:,core2Col])
        
        # Set linkTypeCol to 2 for all valid corridors (those not yet grouped, i.e. not assigned value of 1)    
        linkTable[:,linkTypeCol]=where(linkTable[:,linkTypeCol]==0,2,1)
        
        # Make sure linkTable is sorted 
        ind=lexsort((linkTable[:,core2Col],linkTable[:,core1Col])) 
        linkTable = linkTable[ind]
        
        # Assign link IDs in order
        for x in range(len(linkTable)):
            linkTable[x,linkIdCol]=x+1
                   
        if len(unique(coreList[:,1]))<2:
            gp.addmessage('\n***WARNING: There are less than two core areas. \nThis means there is nothing to connect with linkages.')
        #------------------------------------------------------------------------------
            
            
        # Drop links that are too long 
        gp.addmessage('\nChecking for corridors that are too long to map.')
        disableLeastCostNoVal = False
        linkTable,numDroppedLinks = drop_links(linkTable,useMaxEucDist,maxEucDist,False,0,useMinEucDist,minEucDist,False,0,disableLeastCostNoVal)
        if numDroppedLinks > 0:
            gp.addmessage('\n-------------------------------')
            gp.addmessage('Removed '+str(numDroppedLinks) + ' links that were too long in Euclidean distance.')                                 
            gp.addmessage('-------------------------------\n')
        
        #Write linkTable to disk
        gp.addmessage('\nWriting ' + outlinkTableFile)
        write_link_table(linkTable,coreIds,outlinkTableFile)
        linkTableLogFile=projectDir+"\\log\\"+"linkTable_step2.csv"   
        write_link_table(linkTable,coreIds,linkTableLogFile)
        report_links(linkTable)

        gp.addmessage ('Creating shapefiles with linework for links.\n')
        write_link_maps(outputRootDir,outlinkTableFile,coreShapefile,coreIds,coreAreaIds,step=2)
        gp.addmessage('Linework shapefiles written.')
                   
        if dropFlag==True:
            print_conefor_warning
                
    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        gp.addmessage("\n--------------------")
        gp.addmessage('****Failed in step 2. Details follow.****')        
        filename =  __file__
        raise_geoproc_error(filename)
    
    # Return any PYTHON or system specific errors
    except:
        gp.addmessage("\n--------------------")
        gp.addmessage('****Failed in step 2. Details follow.****')        
        filename =  __file__
        raise_python_error(filename)
        
    del gp    
    return    
    
    

