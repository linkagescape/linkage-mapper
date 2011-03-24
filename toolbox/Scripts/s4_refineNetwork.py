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
# Note: because cwds calculated in step 3, constellation links just connect
# core pairs, not all cores in 1 constellation to all cores in another.
# Could use previous (now discarded) combo code to mosaic CWDS if wanted.

# Import required modules
import sys
import os.path as path
import time

import arcgisscripting
from numpy import *

import lm_config
import lm_util as lu

OUTPUTDIR = lm_config.OUTPUTDIR
LOGDIR = lm_config.LOGDIR
STEP3 = lm_config.STEP3
S4CONNECT = lm_config.S4CONNECT
S4DISTTYPE_EU = lm_config.S4DISTTYPE_EU
S4MAXNN = lm_config.S4MAXNN
MAXEUCDIST = lm_config.MAXEUCDIST
MINEUCDIST = lm_config.MINEUCDIST
GP = lm_config.GP

def step4_refine_network():
    """Allows user to only connect each core area to its N
    nearest neighbors, then connect any disjunct clusters ('constellations')
    of core areas to their nearest neighboring cluster

    """
    try:
        global GP
        GP.Workspace = OUTPUTDIR

        linkTableFile = lu.get_prev_step_link_table(step=4)

        ## linkTable column numbers
        linkIdCol = 0 # Link ID
        core1Col = 1 # core ID of 1st core area link connects
        core2Col = 2 # core ID of 2nd core area link connects
        cluster1Col = 3 # component ID of 1st core area link connects
        cluster2Col = 4 # component ID of 2nd core area link connects
        # 0 = no link. 2 = corridor, 3 = intermediate core area detected,
        # 4 = too long EucDist, 5 = too long lcDist
        linkTypeCol = 5
        eucDistCol = 6
        cwDistCol = 7

        linkTable = lu.load_link_table(linkTableFile)
        numLinks = linkTable.shape[0]
        lu.report_links(linkTable)

        if not STEP3:
            # re-check for links that are too long in case script run out of
            # sequence with more stringent settings
            GP.addmessage('Double-checking for corridors that are too long or '
                          'too short to map.')
            disableLeastCostNoVal = True
            linkTable,numDroppedLinks = lu.drop_links(linkTable, MAXEUCDIST,
                                                   0, MINEUCDIST, 0,
                                                   disableLeastCostNoVal)

        rows,cols = where(linkTable[:,linkTypeCol:linkTypeCol+1] == 2)
        corridorLinks=linkTable[rows,:]
        coresToProcess = unique(corridorLinks[:,core1Col:core2Col+1])

        if S4DISTTYPE_EU:
            distCol = eucDistCol
        else:
            distCol = cwDistCol

        # Flag links that do not connect any core areas to their nearest
        # N neighbors. (N = S4MAXNN)
        lu.dashline(1)
        GP.addmessage('Connecting each core area to its nearest ' +
                      str(S4MAXNN) + ' nearest neighbors.')

        # Code written assuming NO duplicate core pairs
        for core in coresToProcess:
            rows,cols = where(corridorLinks[:,core1Col:core2Col+1] == core)
            distsFromCore=corridorLinks[rows,:]

            # Sort by distance from target core
            ind = argsort(distsFromCore[:,distCol])
            distsFromCore = distsFromCore[ind]

            # Set N nearest neighbor connections to 20
            maxRange = min(len(rows), S4MAXNN)
            for link in range (0,maxRange):
                linkId = distsFromCore[link,linkIdCol]
                # assumes linktable sequentially numbered with no gaps
                linkTable[linkId-1,linkTypeCol] = 20

        # Connect constellations (aka compoments or clusters)
        # Fixme: needs testing.  Move to function.
        if S4CONNECT:
            lu.dashline(1)
            GP.addmessage('Connecting constellations')

            # linkTableComp has 4 extra cols to track COMPONENTS
            numLinks = linkTable.shape[0]
            compCols=zeros((numLinks,4),dtype="int32") # g1' g2' THEN c1 c2
            linkTableComp = append(linkTable,compCols,axis=1)
            del compCols

            # renumber cores to save memory for this next step.  Place in
            # columns 10 and 11
            for coreInd in range(0, len(coresToProcess)):
                # here, cols are 0 for core1Col and 1 for core2Col
                rows,cols=where(linkTableComp[:,core1Col:core2Col+1] ==
                                coresToProcess[coreInd])
                # want results in cols 10 and 11- These are NEW core numbers
                # (0-numcores)
                linkTableComp[rows,cols+10] = coreInd

            rows,cols = where(linkTableComp[:,linkTypeCol:linkTypeCol+1] == 20)
            # The new, improved corridorLinks- only "20" links
            corridorLinksComp=linkTableComp[rows,:]
            # These are NEW core numbers (range from 0 to numcores)
            coresToProcess = unique(linkTableComp[:,10:12])

            #Create graph describing connected cores.
            Graph = zeros((len(coresToProcess), len(coresToProcess)),
                          dtype="int32")
            rows = corridorLinksComp[:,10].astype('int32')
            cols = corridorLinksComp[:,11].astype('int32')
            # But aren't these all 1?
            vals = where(corridorLinksComp[:,linkTypeCol] == 20, 1,0)
            Graph[rows,cols]=vals # why not just 1?
            Graph=Graph+Graph.T

            # Use graph to identify components (disconnected sub-groups) in
            # core area network
            components = lu.components_no_sparse(Graph)
            numComponents = len(unique(components))

            for coreInd in range(0,len(coresToProcess)):
                # In resulting cols, cols are 0 for core1Col and 1 for core2Col
                rows, cols = (where(linkTableComp[:,10:12] ==
                              coresToProcess[coreInd]))
                # want results in cols 12 and 13  Note: we've replaced new core
                # numbers with COMPONENT numbers.
                linkTableComp[rows,cols+12]=components[coreInd]

            # Additional column indexes for linkTableComp
            component1Col = 12
            component2Col = 13
            linkTableComp[:,cluster1Col] = linkTableComp[:,component1Col]
            linkTableComp[:,cluster2Col] = linkTableComp[:,component2Col]

            # Sort by distance
            ind = argsort(linkTableComp[:,distCol])
            linkTableComp= linkTableComp[ind]

            # Connect constellations via shortest inter-constellation links,
            # until all constellations connected.
            for row in range(0,numLinks):
                if ((linkTableComp[row,distCol] > 0) and
                    (linkTableComp[row,linkTypeCol] == 2) and
                    (linkTableComp[row,component1Col] !=
                     linkTableComp[row,component2Col])):
                    # Make this an inter-component link
                    linkTableComp[row,linkTypeCol] = 10
                    newComp = min(linkTableComp
                                  [row,component1Col:component2Col+1])
                    oldComp = max(linkTableComp
                                  [row,component1Col:component2Col+1])
                    # cols are 0 and 1
                    rows, cols = where(linkTableComp
                                       [:,component1Col:component2Col+1]
                                       == oldComp)
                    # want results in cols 12 and 13
                    linkTableComp[rows,cols+12] = newComp

            # Remove extra columns from link table
            linkTable = lu.delete_col(linkTableComp,[10, 11, 12, 13])

            # Re-sort link table by link ID
            ind = argsort(linkTable[:,linkIdCol])
            linkTable= linkTable[ind]

        # At end, any non-comp links that are 2 get assigned -2
        # (too long to be in S4MAXNN, not a component link)
        rows = where(linkTable[:,linkTypeCol] == 2)
        linkTable[rows,linkTypeCol] = -2

        # set 20 links back to 2, get rid of extra columns, re-sort linktable
        rows = where(linkTable[:,linkTypeCol] == 20)
        linkTable[rows,linkTypeCol] = 2

        # Write linkTable to disk
        outlinkTableFile = lu.get_this_step_link_table(step=4)
        lu.dashline(2)
        GP.addmessage('\nWriting ' + outlinkTableFile)
        lu.write_link_table(linkTable, outlinkTableFile)
        linkTableLogFile = path.join(LOGDIR, "linkTable_step4.csv")
        lu.write_link_table(linkTable, linkTableLogFile)

        startTime = time.clock()
        dummy = lu.update_lcp_shapefile(linkTable, lastStep=3, thisStep=4)
        startTime, hours, mins, secs = lu.elapsed_time(startTime)

        lu.dashline()
        GP.addmessage('\nCreating shapefiles with linework for links.')
        lu.write_link_maps(linkTableFile, step=4)

    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        GP.addmessage('****Failed in step 4. Details follow.****')
        filename =  __file__
        lu.raise_geoproc_error(filename)

    # Return any PYTHON or system specific errors
    except:
        GP.addmessage('****Failed in step 4. Details follow.****')
        filename =  __file__
        lu.raise_python_error(filename)

    return