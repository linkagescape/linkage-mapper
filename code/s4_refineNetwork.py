##*****************************************************************
## 2011_0128
## NAME: s4_refineNetwork.py
##
## SUMMARY: Allows user to only connect each core area to its N
## nearest neighbors, then connect any disjunct clusters ('constellations')
## of core areas to their nearest neighboring cluster
##
## SOFTWARE: ArcGIS 9.3 (requires Spatial Analyst extension)
##           Python 2.5
##
##*****************************************************************
#
# Note: because cwds calculated in step 3, constellation links just connect core pairs, 
# not all cores in 1 constellation to all cores in another.
# Could use previous (now discarded) combo code to mosaic CWDS if wanted.

# Import required modules
import arcgisscripting, sys, time
from numpy import *
from watools_util import *

def step4_refine_network(gp,Version,options):
    """Allows user to only connect each core area to its N
    nearest neighbors, then connect any disjunct clusters ('constellations')
    of core areas to their nearest neighboring cluster
    
    """
    try:
        gp.OverwriteOutput = 1
        
        # ------------------------------------------------------------------
        # Unpack options
        projectDir=options["projectDir"]
        coreShapefile=options["coreShapefile"]
        coreIds=options["coreIds"]
        useMinEucDist=options["useMinEucDist"]
        useMaxEucDist=options["useMaxEucDist"]
        maxEucDist=options["maxEucDist"]
        minEucDist=options["minEucDist"]
        useMinCostDist=options["useMinCostDist"]
        useMaxCostDist=options["useMaxCostDist"]
        maxCostDist=options["maxCostDist"]
        minCostDist=options["minCostDist"]
        connectNearestOnly=options["connectNearestOnly"]
        connectComponents=options["connectComponents"]
        maxNumNearestNeighbors=options["maxNumNearestNeighbors"]
        nearestNbrDistType=options["nearestNbrDistType"]
        step3=options["step3"]
        defaultNullValue = options["defaultNullValue"]
        coreAreaIds = coreIds # We now set variable coreAreaIds field name to user input ID field name        
        # ------------------------------------------------------------------        
        
        # ------------------------------------------------------------------
        # Folders, paths, filenames, and workspaces
        outputRootDir = projectDir + "\\output\\"
        dataPassDir = projectDir + "\\dataPass\\"
        gp.Workspace = outputRootDir        
        linkTableFile = get_prev_step_link_table(projectDir,step=4) 
        # ------------------------------------------------------------------    
                         
        # ------------------------------------------------------------------
        ## linkTable column numbers
        linkIdCol=0 # Link ID
        core1Col=1 # core ID of 1st core area link connects
        core2Col=2 # core ID of 2nd core area link connects
        cluster1Col=3 # component ID of 1st core area link connects
        cluster2Col=4 # component ID of 2nd core area link connects
        linkTypeCol=5 # 0=no link. 2=corridor, 3=intermediate core area detected, 4=too long EucDist, 5=too long lcDist
        eucDistCol=6
        cwDistCol=7
        #----------------------------------------------------------
        
        linkTable = load_link_table(linkTableFile)   
        numLinks = linkTable.shape[0]
        report_links(linkTable)

        if step3==False:
            # re-check for links that are too long in case script run out of sequence with more stringent settings
            gp.addmessage('Double-checking for corridors that are too long or too short to map.')
            disableLeastCostNoVal = True
            linkTable,numDroppedLinks = drop_links(linkTable,useMaxEucDist,maxEucDist,False,maxCostDist,useMinEucDist,minEucDist,False,minCostDist,disableLeastCostNoVal)

        rows,cols = where(linkTable[:,linkTypeCol:linkTypeCol+1] == 2) 
        corridorLinks=linkTable[rows,:]
        coresToProcess = unique(corridorLinks[:,core1Col:core2Col+1])
    
        if nearestNbrDistType == "Euclidean":
            distCol = eucDistCol
        else:
            distCol = cwDistCol

            
        if maxNumNearestNeighbors != defaultNullValue: 
            #----------------------------------------------------------
            # Flag links that do not connect any core areas to their nearest N neighbors. (N = maxNumNearestNeighbors)
            gp.addmessage('\n---------------------')
            gp.addmessage('Connecting each core area to its nearest '+str(maxNumNearestNeighbors)+' nearest neighbors.')
      
            for core in coresToProcess: # Code written assuming NO duplicate core pairs.
                rows,cols = where(corridorLinks[:,core1Col:core2Col+1] == core)
                distsFromCore=corridorLinks[rows,:]

                # Sort by distance from target core
                ind = argsort(distsFromCore[:,distCol])
                distsFromCore = distsFromCore[ind]
                
                # Set N nearest neighbor connections to 20
                maxRange = min(len(rows),maxNumNearestNeighbors)
                for link in range (0,maxRange):
                    linkId = distsFromCore[link,linkIdCol]
                    linkTable[linkId-1,linkTypeCol] = 20 # assumes linktable sequentially numbered with no gaps
            #----------------------------------------------------------
        

            #----------------------------------------------------------
            # Connect constellations (aka compoments or clusters)        
            # Fixme: needs testing.  Move to function.
            if connectComponents == True: 
                gp.addmessage('----------------------------------------------')
                gp.addmessage('Connecting constellations')
        
                # linkTableComp has 4 extra cols to track COMPONENTS
                numLinks = linkTable.shape[0]
                compCols=zeros((numLinks,4),dtype="int32") # g1' g2' THEN c1 c2
                linkTableComp = append(linkTable,compCols,axis=1)
                del compCols
                
                # renumber cores to save memory for this next step.  Place in columns 10 and 11
                for coreInd in range(0,len(coresToProcess)):
                    rows,cols=where(linkTableComp[:,core1Col:core2Col+1]==coresToProcess[coreInd]) # here, cols are 0 for core1Col and 1 for core2Col
                    linkTableComp[rows,cols+10]=coreInd # want results in cols 10 and 11- These are NEW core numbers (0-numcores)
                
                rows,cols = where(linkTableComp[:,linkTypeCol:linkTypeCol+1] == 20)
                corridorLinksComp=linkTableComp[rows,:] # The new, improved corridorLinks- only "20" links
          
                coresToProcess = unique(linkTableComp[:,10:12]) # These are NEW core numbers (range from 0 to numcores)
                
                #Create graph describing connected cores.  
                Graph = zeros((len(coresToProcess),len(coresToProcess)),dtype="int32")
                rows = corridorLinksComp[:,10].astype('int32')
                cols = corridorLinksComp[:,11].astype('int32')
                vals = where(corridorLinksComp[:,linkTypeCol] == 20, 1,0) # But aren't these all 1?
                Graph[rows,cols]=vals # why not just 1?
                Graph=Graph+Graph.T
        
                # Use graph to identify components (disconnected sub-groups) in core area network
                components = components_no_sparse(Graph)
                numComponents = len(unique(components))
        
                for coreInd in range(0,len(coresToProcess)):
                    rows,cols=where(linkTableComp[:,10:12]==coresToProcess[coreInd]) # In resulting cols, cols are 0 for core1Col and 1 for core2Col
                    linkTableComp[rows,cols+12]=components[coreInd] # want results in cols 12 and 13  Note: we've replaced new core numbers with COMPONENT numbers.

                # Additional column indexes for linkTableComp
                component1Col = 12
                component2Col = 13
                linkTableComp[:,cluster1Col] = linkTableComp[:,component1Col]
                linkTableComp[:,cluster2Col] = linkTableComp[:,component2Col]
                
                # Sort by distance 
                ind = argsort(linkTableComp[:,distCol])
                linkTableComp= linkTableComp[ind]
           
                # Connect constellations via shortest inter-constellation links, until all constellations connected.
                for row in range(0,numLinks):
                    if (linkTableComp[row,distCol] > 0) and (linkTableComp[row,linkTypeCol] == 2) and (linkTableComp[row,component1Col] != linkTableComp[row,component2Col]): 
                        linkTableComp[row,linkTypeCol] = 10 # Make this an inter-component link
                        newComp=min(linkTableComp[row,component1Col:component2Col+1])
                        oldComp=max(linkTableComp[row,component1Col:component2Col+1])
                        rows,cols=where(linkTableComp[:,component1Col:component2Col+1]== oldComp) # cols are 0 and 1
                        linkTableComp[rows,cols+12] = newComp # want results in cols 12 and 13
                        
                # Remove extra columns from link table
                linkTable = delete_col(linkTableComp,[10, 11, 12, 13])

                # Re-sort link table by link ID
                ind = argsort(linkTable[:,linkIdCol])
                linkTable= linkTable[ind]
                #----------------------------------------------------------    
            
            
            # At end, any non-comp links that are 2 get assigned -2 (too long to be in maxnn, not a component link)   
            rows=where(linkTable[:,linkTypeCol] == 2)
            linkTable[rows,linkTypeCol] = -2
            
            # set 20 links back to 2, get rid of extra columns, re-sort linktable
            rows=where(linkTable[:,linkTypeCol] == 20)
            linkTable[rows,linkTypeCol] = 2
       
        # Write linkTable to disk
        outlinkTableFile = get_this_step_link_table(projectDir,step=4)
        gp.addmessage('---------------------\n')
        gp.addmessage('\nWriting ' + outlinkTableFile)
        write_link_table(linkTable,coreIds,outlinkTableFile)
        linkTableLogFile=projectDir+"\\log\\"+"linkTable_step4.csv"   
        write_link_table(linkTable,coreIds,linkTableLogFile)
        
        startTime = time.clock()
        dummy = update_lcp_shapefile(projectDir,linkTable,lastStep=3,thisStep=4)
        startTime,hours,mins,secs = elapsed_time(startTime)
          
        gp.addmessage('---------------------')
        gp.addmessage('\nCreating shapefiles with linework for links.')
        write_link_maps(outputRootDir,linkTableFile,coreShapefile,coreIds,coreAreaIds,step=4)    
        
    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        gp.addmessage('****Failed in step 4. Details follow.****')        
        filename =  __file__
        raise_geoproc_error(filename)
    
    # Return any PYTHON or system specific errors
    except:
        gp.addmessage('****Failed in step 4. Details follow.****')                  
        filename =  __file__
        raise_python_error(filename)
        
    del gp
    return    
    

