# Need to wipe results of barrier and pinch/network if lm is run
#UPDATE CIRCUITSCAPE WITH VIRAL'S NEW GAPDT!!!!!!!!!!
# for component in components
# indices = where component == component
# delrows
# write graph and solve

# append prefixes to links shapefiles etc.
# move all raster into 1st script.  careful about circuit gdb- may need 2...
#capitalize constants like CONFIGDIR
# solidify output naming.  Consolidate active/inactive links?
# put option check in master.
#separate into 2 scripts.
#check in code this week
# node currents- write copy of core shapefile
# ratio of eff_r to r?
#break s7 into two?
#network mode requires one component!

# ini file?  Can set 
# calculate ratios for all linktables

# remove stick shapefile at start of script?
# output shapefiles to circuitscape gdb
#clean up network current maps
# eventually have nework connectedness measures.  these would have all pairs as sources and targets, 2 measures using eff resis with r's parm'd using eff resis and cwd vals

#TODO:
#circuitscape- don't import wx.  
# put all outputgdb raster names in cfg
# results tab in ug


#!/usr/bin/env python2.5

"""Maps pinch points using Circuitscape given CWD calculations from
       s3_calcCwds.py.
"""

__filename__ = "s8_centrality.py"
__version__ = "CENTRALITY TEST"

import os.path as path
import os
import time
import shutil
import arcgisscripting
import numpy as npy

from lm_config import Config as Cfg
import lm_util as lu

# Set local references to objects and constants from lm_config
gp = Cfg.gp
gprint = gp.addmessage

LTB_CORE1 = Cfg.LTB_CORE1
LTB_CORE2 = Cfg.LTB_CORE2
LTB_LINKID = Cfg.LTB_LINKID
LTB_LINKTYPE = Cfg.LTB_LINKTYPE
LTB_CWDIST = Cfg.LTB_CWDIST

    
def STEP8_calc_centrality():
    """Maps network using Circuitscape given Linkage Mapper outputs

    """
    try:        
        lu.dashline(0)
        gprint('Running script ' + __filename__)
            
        gp.workspace = Cfg.SCRATCHDIR

        
        # set the analysis extent and cell size to that of the resistance
        # surface


        linkTableFile = lu.get_prev_step_link_table(step=7)
        linkTable = lu.load_link_table(linkTableFile)
        numLinks = linkTable.shape[0]
        numCorridorLinks = lu.report_links(linkTable)
        if numCorridorLinks == 0:
            lu.dashline()
            gprint('\nThere are no linkages. Bailing.')
            time.sleep(5)
            return
            
        if linkTable.shape[1] < 11: # If linktable has no entries from prior
                                    # centrality or pinchpint analyses
            extraCols = npy.zeros((numLinks, 5), dtype="float64")
            linkTable = npy.append(linkTable, extraCols, axis=1)
            linkTable[:, Cfg.LTB_LCPLEN] = -1
            linkTable[:, Cfg.LTB_CWDEUCR] = -1
            linkTable[:, Cfg.LTB_CWDPATHR] = -1
            linkTable[:, Cfg.LTB_EFFRESIST] = -1
        
        linkTable[:, Cfg.LTB_CURRENT] = -1
        
        # set up directory for centrality output

        PREFIX = path.basename(Cfg.PROJECTDIR)
        

        coresToProcess = npy.unique(linkTable[:, LTB_CORE1:LTB_CORE2 + 1])
        maxCoreNum = max(coresToProcess)
        del coresToProcess

        lu.dashline(0)

        coreList = linkTable[:,LTB_CORE1:LTB_CORE2+1]
        coreList = npy.sort(coreList)
        #gprint('There are ' + str(len(npy.unique(coreList))) ' core areas.')
            
        INCENTRALITYDIR = Cfg.CENTRALITYBASEDIR
        OUTCENTRALITYDIR = path.join(Cfg.CENTRALITYBASEDIR, 
                                     Cfg.CIRCUITOUTPUTDIR_NM)
        CONFIGDIR = path.join(INCENTRALITYDIR, Cfg.CIRCUITCONFIGDIR_NM)              
        
        # Create output geodatabase
        gp.createfilegdb(Cfg.OUTPUTDIR, path.basename(Cfg.CENTRALITYGDB))


        # Set Circuitscape options and write config file
        options = lu.setCircuitscapeOptions() 
        options['data_type']='network'         
        options['habitat_file'] = path.join(INCENTRALITYDIR,
                                            'Circuitscape_graph.txt')
        # Setting point file equal to graph to do all pairs in Circuitscape
        options['point_file'] = path.join(INCENTRALITYDIR,
                                          'Circuitscape_graph.txt') 
        outputFN = 'Circuitscape_network.out'
        options['output_file'] = path.join(OUTCENTRALITYDIR, outputFN)
        configFN = 'Circuitscape_network.ini'
        outConfigFile = path.join(CONFIGDIR, configFN)
        lu.writeCircuitscapeConfigFile(outConfigFile, options)

        
        delRows = npy.asarray(npy.where(linkTable[:,LTB_LINKTYPE] < 1))
        delRowsVector = npy.zeros((delRows.shape[1]), dtype="int32")
        delRowsVector[:] = delRows[0, :]
        LT = lu.delete_row(linkTable, delRowsVector)
        del delRows
        del delRowsVector
        graphList = npy.zeros((LT.shape[0],3), dtype="float64")
        graphList[:,0] = LT[:,LTB_CORE1]
        graphList[:,1] = LT[:,LTB_CORE2]
        graphList[:,2] = LT[:,LTB_CWDIST]

        # Construct graph from list of nodes and resistances        
        # NOTE graph is in NODE NAME ORDER
        G, nodeNames = make_graph_from_list(graphList)  
        
        # Find connected components in graph
        components = lu.components_no_sparse(G) 
        # Currently, need to solve one component at a time with Circuitscape
        for component in range (1, npy.max(components) + 1):
            inComp = npy.where(components[:] == component)
            notInComp = npy.where(components[:] != component)
            componentGraph = lu.delete_row_col(G,notInComp,notInComp)
            componentNodeNames = nodeNames[inComp]
            rows,cols, = npy.where(componentGraph)
            data = componentGraph[rows,cols]
            numEntries = len(data)
            
            graphListComponent = npy.zeros((numEntries,3),dtype = 'Float64')
            graphListComponent[:,0] = componentNodeNames[rows]
            graphListComponent[:,1] = componentNodeNames[cols]
            graphListComponent[:,2] = data


            write_graph(options['habitat_file'] ,graphListComponent)                
            systemCall = path.join(os.path.dirname(
                os.path.realpath( __file__ )), 'circuitscape\\cs_run.exe')
            systemCall = systemCall + ' ' + outConfigFile  
            import subprocess
            gprint('Calculating current flow centrality...')
            subprocess.check_call(systemCall, shell=True)
            
            outputFN = 'Circuitscape_network_branch_currents_cum.txt'
            currentList = path.join(OUTCENTRALITYDIR, outputFN)
            currents = load_graph(currentList,graphType='graph/network',
                                  datatype='float64')
                       
            numLinks = currents.shape[0]
            for x in range(0,numLinks):
                corex = currents[x,0] 
                corey = currents[x,1] 
                
                #linkId = LT[x,LTB_LINKID]
                row = lu.get_links_from_core_pairs(linkTable, corex, corey)
                #row = lu.get_linktable_row(linkId, linkTable)
                linkTable[row,Cfg.LTB_CURRENT] = currents[x,2] 
        
        finalLinkTable = lu.update_lcp_shapefile(linkTable, lastStep=5,
                                                  thisStep=8)
        linkTableFile = path.join(Cfg.DATAPASSDIR, "linkTable_s7_s8.csv")
        lu.write_link_table(finalLinkTable, linkTableFile)
        linkTableFinalFile = path.join(Cfg.OUTPUTDIR, 
                                       PREFIX + "_linkTable_s7_s8.csv")
        lu.write_link_table(finalLinkTable, linkTableFinalFile)
        gprint('Copy of final linkTable written to '+
                          linkTableFinalFile)
        
        gprint('Creating shapefiles with linework for links.')
        lu.write_link_maps(linkTableFinalFile, step=8)
        
        #FIXME- copy this to centralitygdb.  use copy_final_ code...

        # Clean up temporary files
        lu.delete_dir(Cfg.SCRATCHDIR)
        
        
    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 8. Details follow.****')
        lu.raise_geoproc_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 8. Details follow.****')
        lu.raise_python_error(__filename__)

    return
    
    
    
def write_graph(filename,graphList):
    npy.savetxt(filename,graphList)
    return

def load_graph(filename,graphType,datatype):
    if os.path.isfile(filename)==False:
        raise RuntimeError('File "'  + filename + '" does not exist')
    f = open(filename, 'r')
    try:
        graphObject = npy.loadtxt(filename, dtype = 'Float64', comments='#') 
    except:
        try:
            graphObject = npy.loadtxt(filename, dtype = 'Float64', 
                                      comments='#', delimiter=',')
        except:
            raise RuntimeError('Error reading',type,
                               'file.  Please check file format')
    return graphObject
    

def make_graph_from_list(graphList):
    try:
        nodes = lu.delete_col(graphList,2) 
        nodeNames = (npy.unique(npy.asarray(nodes))).astype('int32')
        nodes[npy.where(nodes >= 0)] = (lu.relabel(nodes[npy.where(nodes >= 0)], 0))
        node1 = (nodes[:,0]).astype('int32')
        node2 = (nodes[:,1]).astype('int32')
        data = graphList[:,2]
        numnodes = nodeNames.shape[0]
        #FIXME!!! CHECK GRAPH ORDER!!!!
        g_graph = npy.zeros((numnodes,numnodes),dtype = 'float64')
        g_graph[node1,node2] = data
        g_graph = g_graph + g_graph.T
        #G = sparse.csr_matrix((conductances, (node1, node2)), shape = (numnodes, numnodes))
        #g_graph = npy.asarray(g_graph)
        return g_graph, nodeNames    

        # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 8. Details follow.****')
        lu.raise_python_error(__filename__)
            
    
    