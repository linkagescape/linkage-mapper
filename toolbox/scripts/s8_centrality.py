#!/usr/bin/env python2.5

"""Maps pinch points using Circuitscape given CWD calculations from
       s3_calcCwds.py.
"""

__filename__ = "s8_centrality.py"
__version__ = "0.7.5"

import os.path as path
import os
import time
import shutil
import arcgisscripting
import numpy as npy
import subprocess
            
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
PREFIX = Cfg.PREFIX
    
def STEP8_calc_centrality():
    """ Experimental code to analyze network centrality using Circuitscape 
        given Linkage Mapper outputs

    """
    try:        
        lu.dashline(0)
        gprint('Running script ' + __filename__)
            
        gp.workspace = Cfg.SCRATCHDIR

        csPath = lu.getCsPath()
        if csPath == None:
            msg = ('Cannot find an installation of Circuitscape 3.5.5' 
                    '\nor greater in your Program Files directory.') 
            gp.AddError(msg)
            gp.AddMessage(gp.GetMessages(2))
            exit(1)
        
        invalidFNs = ['fid','id','oid','shape']
        if Cfg.COREFN.lower() in invalidFNs:
        #if Cfg.COREFN == 'FID' or Cfg.COREFN == 'ID':
            lu.dashline(1)
            msg = ('ERROR: Core area field names ID, FID, SHAPE, and OID are '
                    'reserved for ArcGIS. Please choose another field- must be '
                    'a positive integer.')
            Cfg.gp.AddError(msg)
            exit(1)

        # set the analysis extent and cell size to that of the resistance
        # surface
        coreCopy =  path.join(Cfg.DATAPASSDIR, 
                                     'cores.shp')

        if gp.Exists(coreCopy):
            try:
                gp.Delete(coreCopy)

            except:
                dashline(1)
                msg = ('ERROR: Could not remove shapefile ' +
                       coreCopy + '. Was it open in ArcMap?\n You may '
                       'need to re-start ArcMap to release the file lock.')
                gp.AddError(msg)
                exit(1)


        gp.CopyFeatures_management(Cfg.COREFC, coreCopy)
        gp.AddField_management(coreCopy, "CF_Central", "DOUBLE", "10",
                                   "2")

        inLinkTableFile = lu.get_prev_step_link_table(step=8)
        linkTable = lu.load_link_table(inLinkTableFile)
        numLinks = linkTable.shape[0]
        numCorridorLinks = lu.report_links(linkTable)
        if numCorridorLinks == 0:
            lu.dashline()
            gprint('\nThere are no linkages. Bailing.')
            time.sleep(5)
            return
            
        if linkTable.shape[1] < 16: # If linktable has no entries from prior
                                    # centrality or pinchpint analyses
            extraCols = npy.zeros((numLinks, 6), dtype="float64")
            linkTable = linkTable[:,0:10]
            linkTable = npy.append(linkTable, extraCols, axis=1)
            linkTable[:, Cfg.LTB_LCPLEN] = -1
            linkTable[:, Cfg.LTB_CWDEUCR] = -1
            linkTable[:, Cfg.LTB_CWDPATHR] = -1
            linkTable[:, Cfg.LTB_EFFRESIST] = -1
            linkTable[:, Cfg.LTB_CWDTORR] = -1
            del extraCols

        linkTable[:, Cfg.LTB_CURRENT] = -1
        
        # set up directory for centrality output

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
            if npy.max(components) > 1:
                inComp = npy.where(components[:] == component)
                notInComp = npy.where(components[:] != component)
                componentGraph = lu.delete_row_col(G,notInComp,notInComp)
                componentNodeNames = nodeNames[inComp]
            else:
                componentGraph = G
                componentNodeNames = nodeNames

            rows,cols, = npy.where(componentGraph)
            data = componentGraph[rows,cols]
            numEntries = len(data)
            
            graphListComponent = npy.zeros((numEntries,3),dtype = 'Float64')
            graphListComponent[:,0] = componentNodeNames[rows]
            graphListComponent[:,1] = componentNodeNames[cols]
            graphListComponent[:,2] = data

            write_graph(options['habitat_file'] ,graphListComponent)                

            gprint('\nCalculating current flow centrality for component ' + str(component))
#            subprocess.check_call(systemCall, shell=True)
            subprocess.call([csPath, outConfigFile], shell=True)                     
            
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
        
            coreCurrentFN = 'Circuitscape_network_node_currents_cum.txt'
            nodeCurrentList = path.join(OUTCENTRALITYDIR, coreCurrentFN)
            nodeCurrents = load_graph(nodeCurrentList,graphType='graph/network',
                                  datatype='float64')        

            numNodeCurrents = nodeCurrents.shape[0]
            rows = gp.UpdateCursor(coreCopy)
            row = rows.Next()
            while row:
                coreID = row.getvalue(Cfg.COREFN)
                for i in range (0, numNodeCurrents):
                    if coreID == nodeCurrents[i,0]:
                        row.SetValue("CF_Central", nodeCurrents[i,1])
                        break
                rows.UpdateRow(row)
                row = rows.Next()
            del row, rows          
        
        finalLinkTable = lu.update_lcp_shapefile(linkTable, lastStep=5,
                                                  thisStep=8)
        linkTableFile = path.join(Cfg.DATAPASSDIR, "linkTable_s8.csv")
        lu.write_link_table(finalLinkTable, linkTableFile, inLinkTableFile)
        linkTableFinalFile = path.join(Cfg.OUTPUTDIR, 
                                       PREFIX + "_linkTable_s8.csv")
        lu.write_link_table(finalLinkTable, 
                            linkTableFinalFile, inLinkTableFile)
        gprint('Copy of final linkTable written to '+
                          linkTableFinalFile)
        
        

        # ListFeatureClasses (wildCard, geometryType)	Returns the feature classes in the current workspace
        # ListRasters (wildCard, rasterType)
        finalCoreFile = os.path.join(Cfg.CORECENTRALITYGDB,
                                     PREFIX + '_Cores')       
        #copy core area map to gdb.
        if not gp.exists(Cfg.CORECENTRALITYGDB):
            gp.createfilegdb(Cfg.OUTPUTDIR, path.basename(Cfg.CORECENTRALITYGDB))             
        gp.CopyFeatures_management(coreCopy, finalCoreFile)
        
        gprint('Creating shapefiles with linework for links.')
        lu.write_link_maps(linkTableFinalFile, step=8)
        
        # Copy final link maps to gdb.  
        lu.copy_final_link_maps(step=8)
        
        

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
        nodes[npy.where(nodes >= 0)] = (
            lu.relabel(nodes[npy.where(nodes >= 0)], 0))
        node1 = (nodes[:,0]).astype('int32')
        node2 = (nodes[:,1]).astype('int32')
        data = graphList[:,2]
        numnodes = nodeNames.shape[0]
        #FIXME!!! CHECK GRAPH ORDER!!!!
        g_graph = npy.zeros((numnodes,numnodes),dtype = 'float64')
        g_graph[node1,node2] = data
        g_graph = g_graph + g_graph.T

        return g_graph, nodeNames    

        # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 8. Details follow.****')
        lu.raise_python_error(__filename__)
            
    
    

    
    
