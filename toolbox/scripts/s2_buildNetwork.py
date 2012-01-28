#!/usr/bin/env python2.5
# Authors: Brad McRae and Darren Kavanagh

"""Step 2: Build network.

Generates initial version of linkTable.csv based on euclidean distances and
adjacencies of core areas

"""

__filename__ = "s2_buildNetwork.py"
__version__ = "0.7.7beta"

import os.path as path

import arcgisscripting
import numpy as npy
import time
import gc

from lm_config import Config as Cfg
import lm_util as lu

###
SIMPLIFY_CORES = Cfg.SIMPLIFY_CORES
NEAR_TBL = path.join(Cfg.SCRATCHDIR, "neartbl.dbf")
DIST_FNAME = path.join(Cfg.PROJECTDIR, (Cfg.COREFC + "_dists.txt"))
NEAR_FN = "NEAR_DIST"
EUCADJMATRIXFILE = path.join(Cfg.SCRATCHDIR, "eucAdj.npy")
CWDADJMATRIXFILE = path.join(Cfg.SCRATCHDIR, "cwdcAdj.npy")
DISTMATRIXFILE = path.join(Cfg.SCRATCHDIR, "distances.npy")
EUCADJFILE = Cfg.EUCADJFILE
CWDADJFILE = Cfg.CWDADJFILE        
gp = Cfg.gp
if not Cfg.LOGMESSAGES:
    gprint = gp.addmessage
else:
    gprint = lu.gprint



def STEP2_build_network():
    """Generates initial version of linkTable.csv based on euclidean distances
    and adjacencies of core areas.

    """
    try:
        lu.dashline(1)
        gprint('Running script ' + __filename__)
        outlinkTableFile = lu.get_this_step_link_table(step=2)

        # Warning flag for missing distances in conefor file
        # dropFlag = False

        # ------------------------------------------------------------------
        # adjacency file created from s1_getAdjacencies.py
        if not path.exists(EUCADJFILE):
            msg = ('\nERROR: Euclidean adjacency file required from '
                  'Step 1: ' + EUCADJFILE)
            lu.raise_error(msg)
            
        # ------------------------------------------------------------------
        # adjacency file created from s1_getAdjacencies.py
        if not path.exists(CWDADJFILE):
            msg=('\nERROR: Cost-weighted adjacency file required from'
                              'Step 1: ' + CWDADJFILE)
            lu.raise_error(msg)
        #----------------------------------------------------------------------

        # Load eucDists matrix from file and npy.sort
        if Cfg.S2EUCDISTFILE is None:
            eucdist_file = generate_distance_file(CWDADJFILE,EUCADJFILE)
        else:
            eucdist_file = Cfg.S2EUCDISTFILE
        
        eucDists_in = npy.loadtxt(eucdist_file, dtype='Float64', comments='#')

        if eucDists_in.size == 3:  # If just one line in file
            eucDists = npy.zeros((1,3),dtype = 'Float64')
            eucDists[0,:] = eucDists_in
            numDists = 1

        else:
            eucDists = eucDists_in
            numDists = eucDists.shape[0]
        del eucDists_in
        eucDists[:, 0:2] = npy.sort(eucDists[:, 0:2])
        ind = npy.lexsort((eucDists[:, 2], eucDists[:, 1], eucDists[:, 0]))
        eucDists = eucDists[ind]
        gprint('Core area distance list loaded.')   
        gprint('number of pairwise distances = ' + str(numDists))            
        # sort eucDists by 1st column then by 2nd then by 3rd

        #----------------------------------------------------------------------
        # Get rid of duplicate pairs of cores, retaining MINIMUM distance
        # between them
        numDistsOld = numDists
        for x in range(numDists - 2, -1, -1):
            if (eucDists[x, 0] == eucDists[x + 1, 0]
                and (eucDists[x, 1] == eucDists[x + 1, 1])):
                eucDists[x + 1, 0] = 0
        delRows = npy.asarray(npy.where(eucDists[:, 0] == 0))
        delRowsVector = npy.zeros((delRows.shape[1]), dtype="int32")
        delRowsVector[:] = delRows[0, :]
        eucDists = lu.delete_row(eucDists, delRowsVector)
        del delRows
        del delRowsVector
        numDists = eucDists.shape[0]

        lu.dashline(1)
        gprint('Removed ' + str(numDistsOld - numDists) +
                          ' duplicate core pairs in Euclidean distance table.'
                          '\n')
        maxEucDistID = max(eucDists[:, 1])
        gprint('After removing duplicates and distances that exceed'
                          ' maximum, \nthere are ' + str(numDists) +
                          ' pairwise distances.  Max core ID number is ' +
                          str(int(maxEucDistID)) + '.')


        # Begin creating and manipulating linktables
        # zeros and many other array functions are imported from numpy
        linkTable = npy.zeros((len(eucDists), 10), dtype='int32')
        linkTable[:, 1:3] = eucDists[:, 0:2]
        linkTable[:, Cfg.LTB_EUCDIST] = eucDists[:, 2]

        #----------------------------------------------------------------------
        # Get adjacencies using adj files from step 1.
        cwdAdjTable = get_adj_list(CWDADJFILE)
        cwdAdjList=[]
        for i in range(0,len(cwdAdjTable)):
            listEntry=(str(cwdAdjTable[i,0])+'_'+str(cwdAdjTable[i,1]))
            cwdAdjList.append(listEntry)
        gprint('Cost-weighted adjacency file loaded.')
        maxCwdAdjCoreID = max(cwdAdjTable[:, 1])
        del cwdAdjTable
        
        eucAdjTable = get_adj_list(EUCADJFILE)
        eucAdjList=[]
        for i in range(0,len(eucAdjTable)):
            listEntry=(str(eucAdjTable[i,0])+'_'+str(eucAdjTable[i,1]))
            eucAdjList.append(listEntry)
        maxEucAdjCoreID = max(eucAdjTable[:, 1])      
        del eucAdjTable
        
        maxCoreId = max(maxEucAdjCoreID, maxCwdAdjCoreID, maxEucDistID)

        del eucDists
        
        gprint('Creating link table')        
        linkTable[:, Cfg.LTB_CWDADJ] = -1  # Euc adjacency not evaluated
        linkTable[:, Cfg.LTB_EUCADJ] = -1        
        for x in range(0, linkTable.shape[0]):
            listEntry=(str(linkTable[x, Cfg.LTB_CORE1])+'_'+str(linkTable[x, Cfg.LTB_CORE2]))
            if listEntry in cwdAdjList:
                linkTable[:, Cfg.LTB_CWDADJ] = 1
            else:
                linkTable[:, Cfg.LTB_CWDADJ] = 0                
            if listEntry in eucAdjList:
                linkTable[:, Cfg.LTB_EUCADJ] = 1
            else:
                linkTable[:, Cfg.LTB_EUCADJ] = 0
        
        if Cfg.S2ADJMETH_CW and Cfg.S2ADJMETH_EU:  # "Keep all adjacent links"
            gprint("\nKeeping all adjacent links\n")
            rows = []
            for row in range(0, linkTable.shape[0]):
                if linkTable[row, Cfg.LTB_EUCADJ] == 0 and linkTable[row, Cfg.LTB_CWDADJ] == 0:
                    rows.append(row)
            linkTable = lu.delete_row(linkTable, rows)    

        elif Cfg.S2ADJMETH_CW:
            gprint("\nKeeping cost-weighted adjacent links\n")
            delRows = npy.asarray(npy.where(linkTable[:, Cfg.LTB_CWDADJ] == 0))
            delRowsVector = npy.zeros((delRows.shape[1]), dtype="int32")
            delRowsVector[:] = delRows[0, :]            
            linkTable = lu.delete_row(linkTable, delRowsVector)                
            
        else:
            gprint("\nKeeping Euclidean adjacent links\n")
            delRows = npy.asarray(npy.where(linkTable[:, Cfg.LTB_EUCADJ] == 0))
            delRowsVector = npy.zeros((delRows.shape[1]), dtype="int32")
            delRowsVector[:] = delRows[0, :]            
            linkTable = lu.delete_row(linkTable, delRowsVector)                

        # if dropFlag:
            # lu.dashline(1)
            # gprint('NOTE: At least one adjacent link was dropped '
                          # 'because there was no Euclidean ')
            # gprint('distance value in the input distance file from '
                          # 'Conefor extension.')
            # lu.dashline(2)


        linkTable[:, Cfg.LTB_CLUST1] = -1  # No clusters until later steps
        linkTable[:, Cfg.LTB_CLUST2] = -1

        # not evaluated yet. May eventually have ability to get lcdistances
        # for adjacent cores from s1_getAdjacencies.py
        linkTable[:, Cfg.LTB_CWDIST] = -1

        
        # Get list of core IDs, based on core area shapefile.
        coreList = lu.get_core_list(Cfg.COREFC, Cfg.COREFN)
        if len(npy.unique(coreList[:, 1])) < 2:
            lu.dashline(1)
            msg =('\nERROR: There are less than two core '
                  'areas.\nThis means there is nothing to connect '
                  'with linkages. Bailing.')
            lu.raise_error(msg)
            
        # Set Cfg.LTB_LINKTYPE to valid corridor code
        linkTable[:, Cfg.LTB_LINKTYPE] = Cfg.LT_CORR
        # Make sure linkTable is sorted
        ind = npy.lexsort((linkTable[:, Cfg.LTB_CORE2],
              linkTable[:, Cfg.LTB_CORE1]))
        if len(linkTable) == 0:
            msg = ('\nERROR: There are no valid core area '
                            'pairs. This can happen when core area numbers in '
                            'your Conefor distances text file do not match '
                            'those in your core area feature class.')
            lu.raise_error(msg)

        linkTable = linkTable[ind]

        # Assign link IDs in order
        for x in range(len(linkTable)):
            linkTable[x, Cfg.LTB_LINKID] = x + 1

        #----------------------------------------------------------------------

        # Drop links that are too long
        gprint('\nChecking for corridors that are too long to map.')
        DISABLE_LEAST_COST_NO_VAL = False
        linkTable, numDroppedLinks = lu.drop_links(linkTable, Cfg.MAXEUCDIST,
                                                   0, Cfg.MINEUCDIST, 0,
                                                   DISABLE_LEAST_COST_NO_VAL)
        if numDroppedLinks > 0:
            lu.dashline(1)
            gprint('Removed ' + str(numDroppedLinks) +
                              ' links that were too long in Euclidean '
                              'distance.')
            # lu.dashline(2)

        # Write linkTable to disk
        gprint('Writing ' + outlinkTableFile)
        lu.write_link_table(linkTable, outlinkTableFile)
        linkTableLogFile = path.join(Cfg.LOGDIR, "linkTable_s2.csv")
        lu.write_link_table(linkTable, linkTableLogFile)
        lu.report_links(linkTable)

        gprint('Creating shapefiles with linework for links.\n')
        try:
            lu.write_link_maps(outlinkTableFile, step=2)
        except:
            lu.write_link_maps(outlinkTableFile, step=2)
        gprint('Linework shapefiles written.')
        
        # if dropFlag:
            # print_conefor_warning()

    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 2. Details follow.****')
        lu.raise_geoproc_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 2. Details follow.****')
        lu.raise_python_error(__filename__)

    return

# Fixme: routine below could be used for other operations in code above.
def get_adj_list(adjFile): 
    try:
        inAdjList = npy.loadtxt(adjFile, dtype='int32', comments='#',
                          delimiter=',')  # creates a numpy array
        if len(inAdjList) == inAdjList.size:  # Just one connection
            outAdjList = npy.zeros((1, 3), dtype='int32')
            outAdjList[:, 0:3] = inAdjList[0:3]
        else:
            outAdjList = inAdjList
        outAdjList = outAdjList[:, 1:3]  # Drop first column
        outAdjList = npy.sort(outAdjList) #sorts left-right
        return outAdjList

    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 2. Details follow.****')
        lu.raise_geoproc_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 2. Details follow.****')
        lu.raise_python_error(__filename__)
    
    
def generate_distance_file(CWDADJFILE,EUCADJFILE):
    """Use ArcGIS to create Conefor distance file

    Requires ArcInfo license.

    """
    try:
        #gp.Extent = gp.Describe(Cfg.COREFC).Extent
        gp.CellSize = gp.Describe(Cfg.RESRAST).MeanCellHeight    
                
        if SIMPLIFY_CORES == True:
            gprint('Simplifying polygons for core pair distance calculations')
            COREFC_SIMP = path.join(Cfg.SCRATCHDIR, "CoreFC_Simp.shp")
            tolerance = float(gp.CellSize) / 3
            
            arc10 = True
            try:
                import arcpy.cartography as CA
            except:
                arc10 = False
            if arc10 == True:
                CA.SimplifyPolygon(Cfg.COREFC, COREFC_SIMP,"POINT_REMOVE", 
                                    tolerance, "#", "NO_CHECK")
            else:
                gp.SimplifyPolygon(Cfg.COREFC, COREFC_SIMP,"POINT_REMOVE", 
                                    tolerance, "#", "NO_CHECK")
                                         
                                          
            S2COREFC = COREFC_SIMP
        else:
            S2COREFC = Cfg.COREFC
           
        gp.workspace = Cfg.SCRATCHDIR
        FS2COREFC = "fcores"
        FS2COREFC2 = "fcores2"
        gp.MakeFeatureLayer(S2COREFC, FS2COREFC)
        gp.MakeFeatureLayer(S2COREFC, FS2COREFC2)
        
        output = []
        csvseparator = "\t"

        adjList = get_full_adj_list(CWDADJFILE,EUCADJFILE)
        sourceCores = npy.unique(adjList[:,0])

        gprint('\nFinding distances between cores using Generate Near Table.')
        
        # gprint('old method')
        # start_time = time.clock()
        # gp.generateneartable(S2COREFC, S2COREFC, NEAR_TBL, "#",
                           # "NO_LOCATION", "NO_ANGLE", "ALL", "0")
        # start_time = lu.elapsed_time(start_time)
        
        gprint('There are '+str(len(adjList))+' core pairs to process.')
        pctDone = 0
        start_time = time.clock()
        for x in range(0,len(adjList)):
            
            pctDone = lu.report_pct_done(x, len(adjList), pctDone)
            sourceCore = adjList[x,0]
            targetCore = adjList[x,1]
            expression = Cfg.COREFN + " = " + str(sourceCore)
            gp.selectlayerbyattribute(FS2COREFC, "NEW_SELECTION", expression)
            expression = Cfg.COREFN + " = " + str(targetCore)
            gp.selectlayerbyattribute(FS2COREFC2, "NEW_SELECTION", expression)
    
            gp.generateneartable(FS2COREFC, FS2COREFC2, NEAR_TBL, "#",
                               "NO_LOCATION", "NO_ANGLE", "ALL", "0")        

            rows = gp.searchcursor(NEAR_TBL)
            row = rows.Next()
            if row: #May be running on selected core areas in step 2
                dist = row.getvalue(NEAR_FN)
                if dist <= 0: #In case simplified polygons abut one another
                    dist = gp.CellSize
                outputrow = []           
                outputrow.append(str(sourceCore))
                outputrow.append(str(targetCore))                    
                outputrow.append(str(dist))
                output.append(csvseparator.join(outputrow))
                del row
            del rows   
        start_time = lu.elapsed_time(start_time)           
 
        dist_file = open(DIST_FNAME, 'w')
        dist_file.write('\n'.join(output))
        dist_file.close()
        gprint('Distance file ' + DIST_FNAME + ' generated.\n')
        
        return DIST_FNAME

    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 2. Details follow.****')
        lu.raise_geoproc_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 2. Details follow.****')
        lu.raise_python_error(__filename__)

        
def print_conefor_warning():
    """Warns that some links have no euclidean distances in conefor file."""
    gprint('\nWARNING: At least one potential link was dropped because\n'
        'there was no Euclidean distance value in the input Euclidean\n' 
        'distance file from Conefor extension.\n'
        '   This may just mean that there were core areas that were adjacent\n'
        'but were farther apart than the optional maximum distance used\n'
        'when running Conefor.  But it can also mean that distances  were\n'
        'calculated using a different core area shapefile or the wrong field\n'
        'in the same core area shapefile.\n')

def get_full_adj_list(CWDADJFILE,EUCADJFILE):       
    try:    
        cwdAdjList = get_adj_list(CWDADJFILE)
        eucAdjList = get_adj_list(EUCADJFILE)
        adjList = npy.append(eucAdjList, cwdAdjList, axis=0)
        adjList = npy.sort(adjList)
        
                # sort by 1st core Id then by 2nd core Id
        ind = npy.lexsort((adjList[:, 1], adjList[:, 0]))
        adjList = adjList[ind]
        
        numDists = len(adjList)
        x = 1
        while x < numDists:
            if (adjList[x, 0] == adjList[x - 1, 0] and
                adjList[x, 1] == adjList[x - 1, 1]):
                adjList[x - 1, 0] = 0  # mark for deletion
            x = x + 1

        if numDists > 0:
            delRows = npy.asarray(npy.where(adjList[:, 0] == 0))
            delRowsVector = npy.zeros((delRows.shape[1]), dtype="int32")
            delRowsVector[:] = delRows[0, :]
            adjList = lu.delete_row(adjList, delRowsVector)
            del delRows
            del delRowsVector

        return adjList

    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 2. Details follow.****')
        lu.raise_geoproc_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 2. Details follow.****')
        lu.raise_python_error(__filename__)




        
