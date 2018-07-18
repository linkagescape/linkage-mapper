#!/usr/bin/env python2.5
# Authors: Brad McRae and Darren Kavanagh

"""Step 4: Refine network.

Allows user to only connect each core area to its N nearest neighbors, then
connect any disjunct clusters ('constellations') of core areas to their
nearest neighboring cluster

"""

# Note: because cwds calculated in step 3, constellation links just connect
# core pairs, not all cores in 1 constellation to all cores in another.
# Could use previous (now discarded) combo code to mosaic CWDS if wanted.

import os.path as path
import time

import arcgisscripting
import numpy as npy

from lm_config import tool_env as cfg
import lm_util as lu

_SCRIPT_NAME = "s4_refineNetwork.py"

gprint = lu.gprint

def STEP4_refine_network():
    """Allows user to only connect each core area to its N
    nearest neighbors, then connect any disjunct clusters ('constellations')
    of core areas to their nearest neighboring cluster

    """
    try:

        lu.dashline(1)
        gprint('Running script ' + _SCRIPT_NAME)
        cfg.gp.Workspace = cfg.OUTPUTDIR

        linkTableFile = lu.get_prev_step_link_table(step=4)

        linkTable = lu.load_link_table(linkTableFile)
        numLinks = linkTable.shape[0]
        lu.report_links(linkTable)

        if not cfg.STEP3:
            # re-check for links that are too long in case script run out of
            # sequence with more stringent settings
            gprint('Double-checking for corridors that are too long'
                              ' or too short to map.')
            DISABLE_LEAST_COST_NO_VAL = True
            linkTable,numDroppedLinks = lu.drop_links(
                linkTable, cfg.MAXEUCDIST, 0, cfg.MAXCOSTDIST, 0,
                DISABLE_LEAST_COST_NO_VAL)

        rows, cols = npy.where(
                     linkTable[:,cfg.LTB_LINKTYPE:cfg.LTB_LINKTYPE + 1] > 0)
        # == cfg.LT_CORR
            # or
            # linkTable[:,cfg.LTB_LINKTYPE:cfg.LTB_LINKTYPE + 1] == cfg.LT_KEEP)
        corridorLinks = linkTable[rows,:]
        coresToProcess = npy.unique(
            corridorLinks[:, cfg.LTB_CORE1:cfg.LTB_CORE2 + 1])

        if cfg.S4DISTTYPE_EU:
            distCol = cfg.LTB_EUCDIST
        else:
            distCol = cfg.LTB_CWDIST

        # Flag links that do not connect any core areas to their nearest
        # N neighbors. (N = cfg.S4MAXNN)
        lu.dashline(1)
        # optionally ignore max nearest neighbor setting
        if cfg.IGNORES4MAXNN:
            gprint('Connecting each core area to all its neighbors.')
        else:
            gprint('Connecting each core area to its nearest ' +
                          str(cfg.S4MAXNN) + ' nearest neighbors.')

        # Code written assuming NO duplicate core pairs
        for core in coresToProcess:
            rows,cols = npy.where(
                corridorLinks[:,cfg.LTB_CORE1:cfg.LTB_CORE2+1] == core)
            distsFromCore = corridorLinks[rows,:]

            # Sort by distance from target core
            ind = npy.argsort(distsFromCore[:,distCol])
            distsFromCore = distsFromCore[ind]

            # Set N nearest neighbor connections to Nearest Neighbor (NNCT)
            # optionally ignore max nearest neighbor setting
            if cfg.IGNORES4MAXNN:
                maxRange = len(rows)
            else:
                maxRange = min(len(rows), cfg.S4MAXNN)
            for link in range (0, maxRange):
                linkId = distsFromCore[link, cfg.LTB_LINKID]
                # assumes linktable sequentially numbered with no gaps
                linkTable[linkId - 1, cfg.LTB_LINKTYPE] = cfg.LT_NNCT

        # Connect constellations (aka components or clusters)
        # Fixme: needs testing.  Move to function.
        if cfg.S4CONNECT:
            lu.dashline(1)
            gprint('Connecting constellations')

            # linkTableComp has 4 extra cols to track COMPONENTS
            numLinks = linkTable.shape[0]
            # g1' g2' THEN c1 c2
            compCols = npy.zeros((numLinks, 4), dtype="int32")
            linkTableComp = npy.append(linkTable, compCols, axis=1)
            del compCols

            # renumber cores to save memory for this next step.  Place in
            # columns 10 and 11
            for coreInd in range(0, len(coresToProcess)):
                # here, cols are 0 for cfg.LTB_CORE1 and 1 for cfg.LTB_CORE2
                rows, cols = npy.where(
                    linkTableComp[:,cfg.LTB_CORE1:cfg.LTB_CORE2+1] ==
                    coresToProcess[coreInd])
                # want results in cols 10 and 11- These are NEW core numbers
                # (0 - numcores)
                linkTableComp[rows, cols + 10] = coreInd

            rows, cols = npy.where(
                linkTableComp[:, cfg.LTB_LINKTYPE:cfg.LTB_LINKTYPE + 1] ==
                cfg.LT_NNCT)
            # The new, improved corridorLinks- only NN links
            corridorLinksComp = linkTableComp[rows, :]
            # These are NEW core numbers (range from 0 to numcores)
            coresToProcess = npy.unique(linkTableComp[:, 10:12])

            #Create graph describing connected cores.
            Graph = npy.zeros((len(coresToProcess), len(coresToProcess)),
                          dtype="int32")
            rows = corridorLinksComp[:,10].astype('int32')
            cols = corridorLinksComp[:,11].astype('int32')
            vals = npy.where(corridorLinksComp[:,cfg.LTB_LINKTYPE] ==
                         cfg.LT_NNCT, cfg.LT_CORR, 0)

            Graph[rows,cols] = vals
            Graph = Graph + Graph.T

            # Use graph to identify components (disconnected sub-groups) in
            # core area network
            components = lu.components_no_sparse(Graph)

            for coreInd in range(0,len(coresToProcess)):
                # In resulting cols, cols are 0 for LTB_CORE1 and 1 for
                # LTB_CORE2
                rows, cols = (npy.where(linkTableComp[:,10:12] ==
                              coresToProcess[coreInd]))
                # want results in cols 12 and 13  Note: we've replaced new core
                # numbers with COMPONENT numbers.
                linkTableComp[rows,cols+12] = components[coreInd]
            # Additional column indexes for linkTableComp
            component1Col = 12
            component2Col = 13
            linkTableComp[:,cfg.LTB_CLUST1] = linkTableComp[:,component1Col]
            linkTableComp[:,cfg.LTB_CLUST2] = linkTableComp[:,component2Col]

            # Sort by distance
            ind = npy.argsort(linkTableComp[:,distCol])
            linkTableComp = linkTableComp[ind]

            # Connect constellations via shortest inter-constellation links,
            # until all constellations connected.
            for row in range(0,numLinks):
                if ((linkTableComp[row,distCol] > 0) and
                    ((linkTableComp[row,cfg.LTB_LINKTYPE] == cfg.LT_CORR) or
                    (linkTableComp[row,cfg.LTB_LINKTYPE] == cfg.LT_KEEP)) and
                    (linkTableComp[row,component1Col] !=
                     linkTableComp[row,component2Col])):
                    # Make this an inter-component link
                    linkTableComp[row,cfg.LTB_LINKTYPE] = cfg.LT_CLU
                    newComp = min(linkTableComp
                                  [row,component1Col:component2Col + 1])
                    oldComp = max(linkTableComp
                                  [row,component1Col:component2Col + 1])
                    # cols are 0 and 1
                    rows, cols = npy.where(linkTableComp
                                       [:,component1Col:component2Col + 1]
                                       == oldComp)
                    # want results in cols 12 and 13
                    linkTableComp[rows,cols + 12] = newComp

            # Remove extra columns from link table
            linkTable = lu.delete_col(linkTableComp,[10, 11, 12, 13])

            # Re-sort link table by link ID
            ind = npy.argsort(linkTable[:,cfg.LTB_LINKID])
            linkTable = linkTable[ind]

        # At end, any non-constellation links that are not NN's get dropped
        # (too long to be in cfg.S4MAXNN, not a component link)
        rows = npy.where(linkTable[:,cfg.LTB_LINKTYPE] == cfg.LT_CORR)
        linkTable[rows,cfg.LTB_LINKTYPE] = cfg.LT_CPLK

        # set NNCT links to NN corridor links (NNC), get rid
        # of extra columns, re-sort linktable
        rows = npy.where(linkTable[:,cfg.LTB_LINKTYPE] == cfg.LT_NNCT)
        linkTable[rows,cfg.LTB_LINKTYPE] = cfg.LT_NNC

        # Write linkTable to disk
        outlinkTableFile = lu.get_this_step_link_table(step=4)
        # lu.dashline(1)
        gprint('\nWriting ' + outlinkTableFile)
        lu.write_link_table(linkTable, outlinkTableFile)
        linkTableLogFile = path.join(cfg.LOGDIR, "linkTable_s4.csv")
        lu.write_link_table(linkTable, linkTableLogFile)

        start_time = time.clock()
        lu.update_lcp_shapefile(linkTable, lastStep=3, thisStep=4)
        start_time = lu.elapsed_time(start_time)

        # lu.dashline()
        gprint('Creating shapefiles with linework for links.')
        try:
            lu.write_link_maps(outlinkTableFile, step=4)
        except:
            lu.write_link_maps(outlinkTableFile, step=4)

    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        gprint('****Failed in step 4. Details follow.****')
        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except:
        gprint('****Failed in step 4. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)

    return