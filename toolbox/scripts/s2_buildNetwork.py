#!/usr/bin/env python2.5
# Authors: Brad McRae and Darren Kavanagh

"""Step 2: Build network.

Generates initial version of linkTable.csv based on euclidean distances and
adjacencies of core areas

"""

import os.path as path

import arcgisscripting
import numpy as npy
import time

from lm_config import tool_env as cfg
import lm_util as lu
from lm_retry_decorator import retry

_SCRIPT_NAME = "s2_buildNetwork.py"

gp = cfg.gp
gprint = lu.gprint


def STEP2_build_network():
    """Generates initial version of linkTable.csv based on euclidean distances
    and adjacencies of core areas.

    """
    try:
        lu.dashline(1)
        gprint('Running script ' + _SCRIPT_NAME)
        outlinkTableFile = lu.get_this_step_link_table(step=2)

        # Warning flag for missing distances in conefor file
        # dropFlag = False

        # ------------------------------------------------------------------
        # adjacency file created from s1_getAdjacencies.py
        if cfg.S2ADJMETH_EU and not path.exists(cfg.EUCADJFILE):
            msg = ('\nERROR: Euclidean adjacency file required from '
                  'Step 1: ' + cfg.EUCADJFILE)
            lu.raise_error(msg)

        # ------------------------------------------------------------------
        # adjacency file created from s1_getAdjacencies.py
        if cfg.S2ADJMETH_CW and not path.exists(cfg.CWDADJFILE):
            msg = ('\nERROR: Cost-weighted adjacency file required from'
                              'Step 1: ' + cfg.CWDADJFILE)
            lu.raise_error(msg)
        #----------------------------------------------------------------------

        # Load eucDists matrix from file and npy.sort
        if cfg.S2EUCDISTFILE is None:
            eucdist_file = generate_distance_file()
        else:
            eucdist_file = cfg.S2EUCDISTFILE

        eucDists_in = npy.loadtxt(eucdist_file, dtype='Float64', comments='#')

        if eucDists_in.size == 3:  # If just one line in file
            eucDists = npy.zeros((1, 3), dtype='Float64')
            eucDists[0, :] = eucDists_in
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
        linkTable[:, cfg.LTB_EUCDIST] = eucDists[:, 2]

        #----------------------------------------------------------------------
        # Get adjacencies using adj files from step 1.
        if cfg.S2ADJMETH_CW or cfg.S2ADJMETH_EU:  # Keep ALL links
            cwdAdjList = []
            eucAdjList = []
            if cfg.S2ADJMETH_CW:
                cwdAdjTable = get_adj_list(cfg.CWDADJFILE)
                cwdAdjList = []
                for i in range(0, len(cwdAdjTable)):
                    listEntry = (str(cwdAdjTable[i, 0]) + '_' + str(cwdAdjTable[i, 1]))
                    cwdAdjList.append(listEntry)
                gprint('Cost-weighted adjacency file loaded.')
                maxCwdAdjCoreID = max(cwdAdjTable[:, 1])
                del cwdAdjTable

            if cfg.S2ADJMETH_EU:
                eucAdjTable = get_adj_list(cfg.EUCADJFILE)
                eucAdjList = []
                for i in range(0, len(eucAdjTable)):
                    listEntry = (str(eucAdjTable[i, 0]) + '_' + str(eucAdjTable[i, 1]))
                    eucAdjList.append(listEntry)
                maxEucAdjCoreID = max(eucAdjTable[:, 1])
                del eucAdjTable

        # maxCoreId = max(maxEucAdjCoreID, maxCwdAdjCoreID, maxEucDistID)

        del eucDists

        gprint('Creating link table')
        linkTable[:, cfg.LTB_CWDADJ] = -1  # Euc adjacency not evaluated
        linkTable[:, cfg.LTB_EUCADJ] = -1
        if cfg.S2ADJMETH_CW or cfg.S2ADJMETH_EU:  
            for x in range(0, linkTable.shape[0]):
                listEntry = (str(linkTable[x, cfg.LTB_CORE1]) + '_' +
                             str(linkTable[x, cfg.LTB_CORE2]))
                if listEntry in cwdAdjList:
                    linkTable[x, cfg.LTB_CWDADJ] = 1
                else:
                    linkTable[x, cfg.LTB_CWDADJ] = 0
                if listEntry in eucAdjList:
                    linkTable[x, cfg.LTB_EUCADJ] = 1
                else:
                    linkTable[x, cfg.LTB_EUCADJ] = 0

        if cfg.S2ADJMETH_CW and cfg.S2ADJMETH_EU:  # "Keep all adjacent links"
            gprint("\nKeeping all adjacent links\n")
            rows = []
            for row in range(0, linkTable.shape[0]):
                if (linkTable[row, cfg.LTB_EUCADJ] == 0
                    and linkTable[row, cfg.LTB_CWDADJ] == 0):
                    rows.append(row)
            linkTable = lu.delete_row(linkTable, rows)

        elif cfg.S2ADJMETH_CW:
            gprint("\nKeeping cost-weighted adjacent links\n")
            delRows = npy.asarray(npy.where(linkTable[:, cfg.LTB_CWDADJ] == 0))
            delRowsVector = npy.zeros((delRows.shape[1]), dtype="int32")
            delRowsVector[:] = delRows[0, :]
            linkTable = lu.delete_row(linkTable, delRowsVector)

        elif cfg.S2ADJMETH_EU:
            gprint("\nKeeping Euclidean adjacent links\n")
            delRows = npy.asarray(npy.where(linkTable[:, cfg.LTB_EUCADJ] == 0))
            delRowsVector = npy.zeros((delRows.shape[1]), dtype="int32")
            delRowsVector[:] = delRows[0, :]
            linkTable = lu.delete_row(linkTable, delRowsVector)

        else:  # For Climate Corridor tool
            gprint("\nIgnoring adjacency and keeping all links\n")

        # if dropFlag:
            # lu.dashline(1)
            # gprint('NOTE: At least one adjacent link was dropped '
                          # 'because there was no Euclidean ')
            # gprint('distance value in the input distance file from '
                          # 'Conefor extension.')
            # lu.dashline(2)

        linkTable[:, cfg.LTB_CLUST1] = -1  # No clusters until later steps
        linkTable[:, cfg.LTB_CLUST2] = -1

        # not evaluated yet. May eventually have ability to get lcdistances
        # for adjacent cores from s1_getAdjacencies.py
        linkTable[:, cfg.LTB_CWDIST] = -1

        # Get list of core IDs, based on core area shapefile.
        coreList = lu.get_core_list(cfg.COREFC, cfg.COREFN)
        if len(npy.unique(coreList[:, 1])) < 2:
            lu.dashline(1)
            msg = ('\nERROR: There are less than two core '
                  'areas.\nThis means there is nothing to connect '
                  'with linkages. Bailing.')
            lu.raise_error(msg)

        # Set cfg.LTB_LINKTYPE to valid corridor code
        linkTable[:, cfg.LTB_LINKTYPE] = cfg.LT_CORR
        # Make sure linkTable is sorted
        ind = npy.lexsort((linkTable[:, cfg.LTB_CORE2],
              linkTable[:, cfg.LTB_CORE1]))
        if len(linkTable) == 0:
            msg = ('\nERROR: There are no valid core area '
                            'pairs. This can happen when core area numbers in '
                            'your Conefor distances text file do not match '
                            'those in your core area feature class.')
            lu.raise_error(msg)

        linkTable = linkTable[ind]

        # Assign link IDs in order
        for x in range(len(linkTable)):
            linkTable[x, cfg.LTB_LINKID] = x + 1

        #----------------------------------------------------------------------

        if cfg.CONNECTFRAGS:               
            connect_clusters(linkTable)     
        else:
            # Drop links that are too long
            gprint('\nChecking for corridors that are too long to map.')
            DISABLE_LEAST_COST_NO_VAL = False
            linkTable, numDroppedLinks = lu.drop_links(linkTable, cfg.MAXEUCDIST,
                                                       0, cfg.MINEUCDIST, 0,
                                                       DISABLE_LEAST_COST_NO_VAL)
            if numDroppedLinks > 0:
                lu.dashline(1)
                gprint('Removed ' + str(numDroppedLinks) +
                                  ' links that were too long in Euclidean '
                                  'distance.')

            # Write linkTable to disk
            gprint('Writing ' + outlinkTableFile)
            lu.write_link_table(linkTable, outlinkTableFile)
            linkTableLogFile = path.join(cfg.LOGDIR, "linkTable_s2.csv")
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
        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 2. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)

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
        outAdjList = npy.sort(outAdjList)  # sorts left-right
        return outAdjList

    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 2. Details follow.****')
        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 2. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)


def generate_distance_file():
    """Use ArcGIS to create Conefor distance file

    Requires ArcInfo license.

    """
    try:
        #gp.Extent = gp.Describe(cfg.COREFC).Extent
        gp.CellSize = gp.Describe(cfg.RESRAST).MeanCellHeight
        S2COREFC = cfg.COREFC
        if cfg.SIMPLIFY_CORES:
            try:
                gprint('Simplifying polygons for core pair distance calculations')
                COREFC_SIMP = path.join(cfg.SCRATCHDIR, "CoreFC_Simp.shp")
                tolerance = float(gp.CellSize) / 3

                try:
                    import arcpy
                    import arcpy.cartography as CA
                except:
                    arcpy = False
                if arcpy:
                    CA.SimplifyPolygon(cfg.COREFC, COREFC_SIMP, "POINT_REMOVE",
                                        tolerance, "#", "NO_CHECK")

                else:
                    gp.SimplifyPolygon(cfg.COREFC, COREFC_SIMP, "POINT_REMOVE",
                                        tolerance, "#", "NO_CHECK")

                S2COREFC = COREFC_SIMP
            except:
                pass # In case point geometry is entered for core area FC


        gp.workspace = cfg.SCRATCHDIR
        FS2COREFC = "fcores"
        FS2COREFC2 = "fcores2"
        gp.MakeFeatureLayer(S2COREFC, FS2COREFC)
        gp.MakeFeatureLayer(S2COREFC, FS2COREFC2)

        output = []
        csvseparator = "\t"
        

        adjList = get_full_adj_list()
        # sourceCores = npy.unique(adjList[:, 0])

        gprint('\nFinding distances between cores using Generate Near Table.')
#        gp.OutputCoordinateSystem = gp.describe(cfg.COREFC).SpatialReference
        near_tbl = path.join(cfg.SCRATCHDIR, "neartbl.dbf")
        # gprint('old method')
        # start_time = time.clock()
        # gp.generateneartable(S2COREFC, S2COREFC, near_tbl, "#",
                           # "NO_LOCATION", "NO_ANGLE", "ALL", "0")
        # start_time = lu.elapsed_time(start_time)

        gprint('There are ' + str(len(adjList)) + ' adjacent core pairs to '
               'process.')
        pctDone = 0
        start_time = time.clock()
        for x in range(0, len(adjList)):

            pctDone = lu.report_pct_done(x, len(adjList), pctDone)
            sourceCore = adjList[x, 0]
            targetCore = adjList[x, 1]
            expression = cfg.COREFN + " = " + str(sourceCore)
            gp.selectlayerbyattribute(FS2COREFC, "NEW_SELECTION", expression)
            expression = cfg.COREFN + " = " + str(targetCore)
            gp.selectlayerbyattribute(FS2COREFC2, "NEW_SELECTION", expression)

            gp.generateneartable(FS2COREFC, FS2COREFC2, near_tbl, "#",
                               "NO_LOCATION", "NO_ANGLE", "ALL", "0")

            rows = gp.searchcursor(near_tbl)
            row = rows.Next()
            minDist = 1e20
            if row:  # May be running on selected core areas in step 2
                while row:
                    dist = row.getvalue("NEAR_DIST")
                    if dist <= 0:  # In case simplified polygons abut one another
                        dist = float(gp.CellSize)
                    if dist < minDist:
                        minDist = dist
                        outputrow = []
                        outputrow.append(str(sourceCore))
                        outputrow.append(str(targetCore))
                        outputrow.append(str(dist))
                    del row
                    row = rows.Next()              
            del rows
            output.append(csvseparator.join(outputrow))  
              
        start_time = lu.elapsed_time(start_time)

        # In case coreFC is grouped in TOC, get coreFN for non-Arc statement
        group,coreFN = path.split(cfg.COREFC)

        dist_fname = path.join(cfg.PROJECTDIR, (coreFN + "_dists.txt"))
        dist_file = open(dist_fname, 'w')
        dist_file.write('\n'.join(output))
        dist_file.close()
        gprint('Distance file ' + dist_fname + ' generated.\n')

        return dist_fname

    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 2. Details follow.****')
        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 2. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)


def print_conefor_warning():
    """Warns that some links have no euclidean distances in conefor file."""
    gprint('\nWARNING:')
    gprint('At least one potential link was dropped because\n'
        'there was no Euclidean distance value in the input Euclidean\n'
        'distance file from Conefor extension.\n'
        '   This may just mean that there were core areas that were adjacent\n'
        'but were farther apart than the optional maximum distance used\n'
        'when running Conefor.  But it can also mean that distances  were\n'
        'calculated using a different core area shapefile or the wrong field\n'
        'in the same core area shapefile.\n')


def get_full_adj_list():
    try:
        if not cfg.S2ADJMETH_CW and not cfg.S2ADJMETH_EU:  # Keep ALL links
            coreList = lu.get_core_list(cfg.COREFC, cfg.COREFN)
            coreList = coreList[:,0]
            numCores = len(coreList)
            adjList = npy.zeros((numCores*(numCores-1)/2,2), dtype="int32")
            pairIndex = 0
            for sourceIndex in range(0,numCores-1):
                for targetIndex in range(sourceIndex + 1, numCores):
                    adjList[pairIndex,0]=coreList[sourceIndex]
                    adjList[pairIndex,1]=coreList[targetIndex]
                    pairIndex = pairIndex + 1
            return adjList
        eucAdjList = get_adj_list(cfg.EUCADJFILE)
        if cfg.S2ADJMETH_CW:
            cwdAdjList = get_adj_list(cfg.CWDADJFILE)
            adjList = npy.append(eucAdjList, cwdAdjList, axis=0)
        else:
            adjList = eucAdjList
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
        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 2. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)

def connect_clusters(linkTable):
        # CUSTOM Fragment connecting code
    try:
        clusterFC = path.join(cfg.SCRATCHDIR,"Cores_Grouped_dist"+str(int(cfg.MAXEUCDIST))+".shp")
        gp.CopyFeatures_management(cfg.COREFC,clusterFC)
#        coreFC_orig = cfg.COREFC 
#        clusterFC = tempShapefile
        
        gprint('Running custom fragment connecting code.')
        numLinks = linkTable.shape[0]

        fieldList = gp.ListFields(clusterFC)
        cluster_ID = 'clus' + str(int(cfg.MAXEUCDIST))
        for field in fieldList:
            if str(field.Name) == cluster_ID:
                gp.DeleteField_management(clusterFC, cluster_ID)
        gp.AddField_management(clusterFC, cluster_ID, "LONG")
        
        rows = gp.UpdateCursor(clusterFC)
        row = rows.Next()
        while row:
            # linkCoords indices
            fragID = row.GetValue(cfg.COREFN) 
            row.SetValue(cluster_ID, fragID)
            rows.UpdateRow(row)
            row = rows.Next()
        del row, rows


        linkTable[:, cfg.LTB_CLUST1] = linkTable[:, cfg.LTB_CORE1]
        linkTable[:, cfg.LTB_CLUST2] = linkTable[:, cfg.LTB_CORE2]
        
        #if frags less than cutoff set cluster_ID equal.        
        for x in range(0,numLinks):
            # if linkTable[x, cfg.LTB_LINKTYPE] == cfg.LT_CORR:            
                gprint("link #"+str(x+1))
                # Set newfragmentID of 2nd fragment to that of frag 1
                frag1ID = linkTable[x, cfg.LTB_CLUST1]
                frag2ID = linkTable[x, cfg.LTB_CLUST2]
                if frag1ID == frag2ID:
                    continue
                eucDist = linkTable[x, cfg.LTB_EUCDIST]

                # cur = gp.SearchCursor(clusterFC)
                # row = cur.Next()
                # while row:
                    # if row.GetValue(cluster_ID) == frag1ID:
                        # frag1CoreID =  row.GetValue(cluster_ID)
                    # elif row.GetValue(cluster_ID) == frag2ID:
                        # frag2CoreID =  row.GetValue(cluster_ID)
                    # row = cur.Next()
                    # i = i + 1
                # del cur, row

                if eucDist < cfg.MAXEUCDIST:
#                    gprint("Joining fragments "+str(frag1ID)+" and "+str(frag2ID)+" in Core #"+str(frag2CoreID)+" separated by distance "+str(eucDist))
                    gprint("Joining fragments "+str(frag1ID)+" and "+str(frag2ID)+" separated by distance "+str(eucDist))
                    # update linktable to new fragment ID in cluster field
                    rows = npy.where(linkTable[:,cfg.LTB_CLUST1] == frag2ID)
                    linkTable[rows, cfg.LTB_CLUST1] = frag1ID
                    rows = npy.where(linkTable[:,cfg.LTB_CLUST2] == frag2ID)
                    linkTable[rows, cfg.LTB_CLUST2] = frag1ID
                    del rows

                    # update shapefile to new fragment ID in cluster_ID field
                    rows = gp.UpdateCursor(clusterFC)
                    row = rows.Next()
                    while row:
                        if row.GetValue(cluster_ID) == frag2ID:
                            row.SetValue(cluster_ID, frag1ID)
                        rows.UpdateRow(row)
                        row = rows.Next()
                    del row, rows
                # else:   
                    # gprint("Adjacent fragments not close enough ("+str(frag1CoreID)+" and "+str(frag2CoreID)+")")
        
        gprint('Done Joining.  Creating output shapefiles.')
        
        coreBaseName = path.splitext(path.basename(cfg.COREFC))[0]

        # outputFN = coreBaseName + "Cluster_Dist"+ str(int(cfg.MAXEUCDIST)) + ".shp"
        # outputShapefile = path.join(cfg.PROJECTDIR, outputFN)
        # gp.CopyFeatures_management(clusterFC,outputShapefile)
        outputFN = coreBaseName + "_Cluster"+str(int(cfg.MAXEUCDIST))+"_dissolve.shp"
        outputShapefile = path.join(cfg.SCRATCHDIR,outputFN)
        gp.Dissolve_management(clusterFC, outputShapefile, cluster_ID)
        outputFN = coreBaseName + "_Cluster"+str(int(cfg.MAXEUCDIST))+"_dissolve_area.shp"
        coreFCWithArea = path.join(cfg.SCRATCHDIR,outputFN)
        gp.CalculateAreas_stats(outputShapefile, coreFCWithArea)

        outputFN = coreBaseName + "_Cluster"+str(int(cfg.MAXEUCDIST))+".shp"
        clusterFCFinal = path.join(cfg.PROJECTDIR,outputFN)
        gp.CopyFeatures_management(clusterFC,clusterFCFinal)
        
        # Update final core featureclass with cluster ID and area
        gp.AddField_management(clusterFCFinal, cluster_ID, "LONG")
        gp.AddField_management(clusterFCFinal, "clust_area", "DOUBLE")

        #run through rows- get cluster id, then get area from coreFCwitharea using searchcursor

        rows = gp.UpdateCursor(clusterFCFinal)
        row = rows.Next()
        while row:
            # linkCoords indices
            clustID = row.GetValue(cluster_ID)
            rows2 = gp.searchcursor(coreFCWithArea)
            row2 = rows2.Next()
            while row2:
                if row2.GetValue(cluster_ID) == clustID:
                    fArea = "F_AREA"
                    clustArea = row2.GetValue(fArea)
                    break
                row2 = rows2.Next()
            row.SetValue("clust_area", clustArea)
            rows.UpdateRow(row)
            row = rows.Next()
        del row, rows, row2, rows2
        gprint('Cores with cluster ID and cluster area written to: ' 
                + clusterFCFinal)

        outlinkTableFile = path.join(cfg.DATAPASSDIR,'linktable_clusters.csv')
        gprint('Writing ' + outlinkTableFile)
        lu.write_link_table(linkTable, outlinkTableFile)            
                
        
        
        
        ##########################################################    

    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 2. Details follow.****')
        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 2. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)
        