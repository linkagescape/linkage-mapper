#!/usr/bin/env python2.5
# Authors: Brad McRae and Darren Kavanagh

"""Step 2: Build network.

Generates initial version of linkTable.csv based on euclidean distances and
adjacencies of core areas

"""

__filename__ = "s2_buildNetwork.py"
__version__ = "0.6.6"

import os.path as path

import arcgisscripting
import numpy as npy
import time

from lm_config import Config as Cfg
import lm_util as lu

SIMPLIFY_CORES = True #fixme: move to config

NEAR_TBL = path.join(Cfg.SCRATCHDIR, "neartbl.dbf")
DIST_FNAME = path.join(Cfg.PROJECTDIR, (Cfg.COREFC + "_dist.txt"))
FID_FN = "FID"
INFID_FN = "IN_FID"
NEARID_FN = "NEAR_FID"
NEAR_FN = "NEAR_DIST"
S2CORE_RAS = "s2core_ras"
S2COREFN = Cfg.COREFN
NEAR_COREFN = S2COREFN + "_1"


gp = Cfg.gp
gprint = gp.addmessage


def STEP2_build_network():
    """Generates initial version of linkTable.csv based on euclidean distances
    and adjacencies of core areas.

    """
    try:
        lu.dashline(1)
        gprint('Running script' + __filename__)
        outlinkTableFile = lu.get_this_step_link_table(step=2)

        # Warning flag for missing distances in conefor file
        dropFlag = False

        # ------------------------------------------------------------------
        # adjacency file created from s1_getAdjacencies.py
        eucAdjFile = path.join(Cfg.DATAPASSDIR, "eucAdj.csv")
        if not path.exists(eucAdjFile):
            gprint('\nERROR: Euclidean adjacency file required from '
                  'Step 1: ' + eucAdjFile)
            exit(0)

        # ------------------------------------------------------------------
        # adjacency file created from s1_getAdjacencies.py
        cwdAdjFile = path.join(Cfg.DATAPASSDIR, "cwdAdj.csv")
        if not path.exists(cwdAdjFile):
            gprint('\nERROR: Cost-weighted adjacency file required from'
                              'Step 1: ' + cwdAdjFile)
            exit(0)
        #----------------------------------------------------------------------

        # Load eucDists matrix from file and npy.sort
        if Cfg.S2EUCDISTFILE is None:
            eucdist_file = generate_distance_file(cwdAdjFile,eucAdjFile)
        else:
            eucdist_file = Cfg.S2EUCDISTFILE

        eucDists = npy.loadtxt(eucdist_file, dtype='Float64', comments='#')
        numDists = eucDists.shape[0]
        # lu.dashline()
        gprint('Core area distance list '
                          'loaded. \n')
        gprint('number of pairwise distances = ' + str(numDists))
        # lu.dashline(2)
        eucDists[:, 0:2] = npy.sort(eucDists[:, 0:2])

        # sort eucDists by 1st column then by 2nd then by 3rd
        ind = npy.lexsort((eucDists[:, 2], eucDists[:, 1], eucDists[:, 0]))
        eucDists = eucDists[ind]

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
        maxeudistid = max(eucDists[:, 1])
        gprint('After removing duplicates and distances that exceed'
                          ' maximum, \nthere are ' + str(numDists) +
                          ' pairwise distances.  Max core ID number is ' +
                          str(int(maxeudistid)) + '.')
        # lu.dashline(2)

        # Begin creating and manipulating linktables
        # zeros and many other array functions are imported from numpy
        distlinkTable = npy.zeros((len(eucDists), 10))
        # this kind of indexing is from numpy too
        distlinkTable[:, 1:3] = eucDists[:, 0:2]
        # Cfg.LTB_EUCDIST is just a number from the index table above.
        # It is just used to specify the column where euclidean distances are
        # stored
        distlinkTable[:, Cfg.LTB_EUCDIST] = eucDists[:, 2]

        #----------------------------------------------------------------------
        # Get adjacencies using adj files from step 1.
        adjList = npy.loadtxt(cwdAdjFile, dtype='int32', comments='#',
                          delimiter=',')  # creates a numpy array

        if len(adjList) == adjList.size:  # Just one connection
            cwdAdjList = npy.zeros((1, 3), dtype='int32')
            cwdAdjList[:, 0:3] = adjList[0:3]
        else:
            cwdAdjList = adjList
        cwdAdjList = cwdAdjList[:, 1:3]  # Drop first column
        cwdAdjList = npy.sort(cwdAdjList)
        gprint('Cost-weighted adjacency file loaded.')
        maxCwdAdjCoreID = max(cwdAdjList[:, 1])

        # Create Numpy array
        adjList = npy.loadtxt(eucAdjFile, dtype='int32', comments='#',
                              delimiter=',')

        if len(adjList) == adjList.size:  # Just one connection
            eucAdjList = npy.zeros((1, 3), dtype='int32')
            eucAdjList[:, 0:3] = adjList[0:3]
        else:
            eucAdjList = adjList
        eucAdjList = eucAdjList[:, 1:3]  # Drop first column
        eucAdjList = npy.sort(eucAdjList)
        gprint('Euclidean adjacency file loaded')
        maxEucAdjCoreID = max(eucAdjList[:, 1])

        maxCoreId = max(maxEucAdjCoreID, maxCwdAdjCoreID, maxeudistid)

        # FIXME: consider using a lookup table to reduce size of matrix
        # when there are gaps in core areas
        cwdAdjMatrix = npy.zeros((maxCoreId + 1, maxCoreId + 1),
                             dtype='int32')
        for x in range(0, len(cwdAdjList)):
            cwdAdjMatrix[cwdAdjList[x, 0], cwdAdjList[x, 1]] = 1
        cwdAdjMatrix[0, :] = 0  # 0 values for core Ids are invalid
        cwdAdjMatrix[0, :] = 0

        eucAdjMatrix = npy.zeros((maxCoreId + 1, maxCoreId + 1),
                             dtype='int32')
        for x in range(0, len(eucAdjList)):
            eucAdjMatrix[eucAdjList[x, 0], eucAdjList[x, 1]] = 1
        eucAdjMatrix[0, :] = 0  # 0 values for core Ids are invalid
        eucAdjMatrix[0, :] = 0

        distanceMatrix = npy.zeros((maxCoreId + 1, maxCoreId + 1),
                                   dtype='int32')

        for x in range(0, len(eucDists)):
            distanceMatrix[eucDists[x, 0], eucDists[x, 1]] = (
                int(eucDists[x, 2]))

        if Cfg.S2ADJMETH_CW:
            difference = npy.where(distanceMatrix, 1, 0)
            difference = difference - cwdAdjMatrix
            if npy.amin(difference) < 0:
                dropFlag = True
            # Drop anything not cwd adjacent
            distanceMatrix = npy.multiply(distanceMatrix, cwdAdjMatrix)

        elif Cfg.S2ADJMETH_EU:
            difference = npy.where(distanceMatrix, 1, 0)
            difference = difference - eucAdjMatrix
            if npy.amin(difference) < 0:
                dropFlag = True
            del difference
            # Drop anything not euc adjacent
            distanceMatrix = npy.multiply(distanceMatrix, eucAdjMatrix)

        else:  # "Keep all adjacent links"
            adjMatrix = eucAdjMatrix
            adjMatrix = adjMatrix + cwdAdjMatrix
            adjMatrix = npy.where(adjMatrix == 2, 1, adjMatrix)
            difference = npy.where(distanceMatrix, 1, 0)
            difference = difference - adjMatrix

            if npy.amin(difference) < 0:
                dropFlag = True

            del difference
            # Drop anything not adjacent
            distanceMatrix = npy.multiply(distanceMatrix, adjMatrix)

        #----------------------------------------------------------------------
        # OK, we have distance matrix (which now defines which pairs can be
        # potential links).  Use it to create link table.
        gprint('creating link table')
        # Get rid of 0 index- we don't have any valid core ids with 0 values
        distanceMatrix = lu.delete_row_col(distanceMatrix, 0, 0)
        rows, cols = npy.where(distanceMatrix)
        linkTable = npy.zeros((len(rows), 10), dtype='int32')

        for x in range(0, len(rows)):
            linkTable[x, Cfg.LTB_CORE1] = rows[x] + 1
            linkTable[x, Cfg.LTB_CORE2] = cols[x] + 1
            linkTable[x, Cfg.LTB_EUCDIST] = distanceMatrix[rows[x], cols[x]]
        del distanceMatrix

        # Get rid of 0 index- we don't have any valid core ids with 0
        # values
        eucAdjMatrix = lu.delete_row_col(eucAdjMatrix, 0, 0)
        for x in range(0, len(rows)):
            linkTable[x, Cfg.LTB_EUCADJ] = eucAdjMatrix[rows[x], cols[x]]
        del eucAdjMatrix

        # Get rid of 0 index- we don't have any valid core ids with 0
        # values
        cwdAdjMatrix = lu.delete_row_col(cwdAdjMatrix, 0, 0)
        for x in range(0, len(rows)):
            linkTable[x, Cfg.LTB_CWDADJ] = cwdAdjMatrix[rows[x], cols[x]]
        del cwdAdjMatrix

        if dropFlag:
            lu.dashline(1)
            gprint('NOTE: At least one adjacent link was dropped '
                          'because there was no Euclidean ')
            gprint('distance value in the input distance file from '
                          'Conefor extension.')
            # lu.dashline(2)

        # Get list of core IDs, based on core area shapefile.
        coreList = lu.get_core_list()
        numCores = len(coreList)

        linkTable[:, Cfg.LTB_CLUST1] = -1  # No clusters until later steps
        linkTable[:, Cfg.LTB_CLUST2] = -1

        linkTable[:, Cfg.LTB_CWDADJ] = -1  # Euc adjacency not evaluated
        linkTable[:, Cfg.LTB_EUCADJ] = -1

        # not evaluated yet. May eventually have ability to get lcdistances
        # for adjacent cores from s1_getAdjacencies.py
        linkTable[:, Cfg.LTB_CWDIST] = -1

        # Update linkTable with new core IDs
        numCores = len(coreList)
        for core in range(numCores):
            if coreList[core, 0] != coreList[core, 1]:
                linkTable[:, Cfg.LTB_CORE1] = npy.where(
                    linkTable[:, Cfg.LTB_CORE1] == coreList[core, 0],
                    coreList[core, 1], linkTable[:, Cfg.LTB_CORE1])
                linkTable[:, Cfg.LTB_CORE2] = npy.where(
                    linkTable[:, Cfg.LTB_CORE2] == coreList[core, 0],
                    coreList[core, 1], linkTable[:, Cfg.LTB_CORE2])

        # Set Cfg.LTB_LINKTYPE to valid corridor code
        linkTable[:, Cfg.LTB_LINKTYPE] = Cfg.LT_CORR
        # Make sure linkTable is sorted
        ind = npy.lexsort((linkTable[:, Cfg.LTB_CORE2],
              linkTable[:, Cfg.LTB_CORE1]))
        if len(linkTable) == 0:
            gp.Adderror('\nERROR: There are no valid core area '
                            'pairs. This can happen when core area numbers in '
                            'your Conefor distances text file do not match '
                            'those in your core area feature class.')
            exit(0)

        linkTable = linkTable[ind]

        # Assign link IDs in order
        for x in range(len(linkTable)):
            linkTable[x, Cfg.LTB_LINKID] = x + 1

        if len(npy.unique(coreList[:, 1])) < 2:
            gprint('\nERROR: There are less than two core '
                              'areas.\nThis means there is nothing to connect '
                              'with linkages. Bailing.')
            exit(0)
        #----------------------------------------------------------------------

        # Drop links that are too long
        gprint('\nChecking for corridors that are too long to map.')
        disableLeastCostNoVal = False
        linkTable, numDroppedLinks = lu.drop_links(linkTable, Cfg.MAXEUCDIST,
                                                   0, Cfg.MINEUCDIST, 0,
                                                   disableLeastCostNoVal)
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
        lu.write_link_maps(outlinkTableFile, step=2)
        gprint('Linework shapefiles written.')

        if dropFlag:
            lu.print_conefor_warning()

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
        outAdjList = npy.sort(outAdjList)
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
    
    
def generate_distance_file(cwdAdjFile,eucAdjFile):
    """Use ArcGIS to create Conefor distance file

    Requires ArcInfo license.

    """
    try:
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

        adjList = get_full_adj_list(cwdAdjFile,eucAdjFile)
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
            expression = S2COREFN + " = " + str(sourceCore)
            gp.selectlayerbyattribute(FS2COREFC, "NEW_SELECTION", expression)
            expression = S2COREFN + " = " + str(targetCore)
            gp.selectlayerbyattribute(FS2COREFC2, "NEW_SELECTION", expression)
    
            gp.generateneartable(FS2COREFC, FS2COREFC2, NEAR_TBL, "#",
                               "NO_LOCATION", "NO_ANGLE", "ALL", "0")        

            rows = gp.searchcursor(NEAR_TBL)
            row = rows.Next()
            if row: #May be running on selected core areas in step 2
                dist = row.getvalue(NEAR_FN)
                if dist == 0: #In case simplified polygons abut one another
                    dist = gp.CellSize
                outputrow = []           
                outputrow.append(str(sourceCore))
                outputrow.append(str(targetCore))                    
                outputrow.append(str(dist))
                output.append(csvseparator.join(outputrow))
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

        

def get_full_adj_list(cwdAdjFile,eucAdjFile):       
    try:    
        cwdAdjList = get_adj_list(cwdAdjFile)
        eucAdjList = get_adj_list(eucAdjFile)
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




        
