#!/usr/bin/env python2.5
# Author: Brad McRae

"""Maps pinch points using Circuitscape given CWD calculations from
       s3_calcCwds.py.
Reguired Software:
ArcGIS 10 with Spatial Analyst extension
Python 2.6
Numpy
"""

import os.path as path

import arcpy
import numpy as npy

from lm_config import tool_env as cfg
import lm_util as lu


_SCRIPT_NAME = "s7_centrality.py"

arcpy.CheckOutExtension("spatial")

gprint = lu.gprint


def STEP7_calc_centrality():
    """ Analyze network centrality using Circuitscape
        given Linkage Mapper outputs

    """
    try:
        lu.dashline(0)
        gprint('Running script ' + _SCRIPT_NAME)

        arcpy.env.workspace = cfg.SCRATCHDIR

        # Check for valid LCP shapefile
        prevLcpShapefile = lu.get_lcp_shapefile(None, thisStep = 7)
        if not arcpy.Exists(prevLcpShapefile):
            msg = ('Cannot find an LCP shapefile from step 5.  Please '
                    'rerun that step and any previous ones if necessary.')
            lu.raise_error(msg)

        # Remove lcp shapefile from this step if run previously
        lcpShapefile = path.join(cfg.DATAPASSDIR, "lcpLines_s7.shp")
        lu.delete_data(lcpShapefile)

        invalidFNs = ['fid','id','oid','shape']
        if cfg.COREFN.lower() in invalidFNs:
        #if cfg.COREFN == 'FID' or cfg.COREFN == 'ID':
            lu.dashline(1)
            msg = ('ERROR: Core area field names ID, FID, SHAPE, and OID are'
                    ' reserved for ArcGIS. \nPlease choose another field- must'
                    ' be a positive integer.')
            lu.raise_error(msg)

        lu.dashline(1)
        gprint('Mapping centrality of network cores and links'
                '\nusing Circuitscape....')
        lu.dashline(0)

        # set the analysis extent and cell size to that of the resistance
        # surface
        coreCopy =  path.join(cfg.SCRATCHDIR, 'cores.shp')

        arcpy.CopyFeatures_management(cfg.COREFC, coreCopy)
        exists = False
        field_names = [field.name for field in arcpy.ListFields(coreCopy)]
        if "CF_Central" in field_names:
            exists = True
        if not exists:
            # arcpy.AddField_management(coreCopy, "CF_Central", "DOUBLE", "10", "2")
            arcpy.AddField_management(coreCopy, "CF_Central", "DOUBLE")

        inLinkTableFile = lu.get_prev_step_link_table(step=7)
        linkTable = lu.load_link_table(inLinkTableFile)
        numLinks = linkTable.shape[0]
        numCorridorLinks = lu.report_links(linkTable)
        if numCorridorLinks == 0:
            lu.dashline(1)
            msg =('\nThere are no linkages. Bailing.')
            lu.raise_error(msg)

        if linkTable.shape[1] < 16: # If linktable has no entries from prior
                                    # centrality or pinchpint analyses
            extraCols = npy.zeros((numLinks, 6), dtype="float64")
            linkTable = linkTable[:,0:10]
            linkTable = npy.append(linkTable, extraCols, axis=1)
            linkTable[:, cfg.LTB_LCPLEN] = -1
            linkTable[:, cfg.LTB_CWDEUCR] = -1
            linkTable[:, cfg.LTB_CWDPATHR] = -1
            linkTable[:, cfg.LTB_EFFRESIST] = -1
            linkTable[:, cfg.LTB_CWDTORR] = -1
            del extraCols

        linkTable[:, cfg.LTB_CURRENT] = -1

        coresToProcess = npy.unique(linkTable[:, cfg.LTB_CORE1:cfg.LTB_CORE2 + 1])
        maxCoreNum = max(coresToProcess)
        del coresToProcess

        lu.dashline(0)

        coreList = linkTable[:,cfg.LTB_CORE1:cfg.LTB_CORE2+1]
        coreList = npy.sort(coreList)
        #gprint('There are ' + str(len(npy.unique(coreList))) ' core areas.')

        # set up directory for centrality
        INCENTRALITYDIR = cfg.CENTRALITYBASEDIR
        OUTCENTRALITYDIR = path.join(cfg.CENTRALITYBASEDIR,
                                     cfg.CIRCUITOUTPUTDIR_NM)
        CONFIGDIR = path.join(INCENTRALITYDIR, cfg.CIRCUITCONFIGDIR_NM)

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

        delRows = npy.asarray(npy.where(linkTable[:,cfg.LTB_LINKTYPE] < 1))
        delRowsVector = npy.zeros((delRows.shape[1]), dtype="int32")
        delRowsVector[:] = delRows[0, :]
        LT = lu.delete_row(linkTable, delRowsVector)
        del delRows
        del delRowsVector
        graphList = npy.zeros((LT.shape[0],3), dtype="float64")
        graphList[:,0] = LT[:,cfg.LTB_CORE1]
        graphList[:,1] = LT[:,cfg.LTB_CORE2]
        graphList[:,2] = LT[:,cfg.LTB_CWDIST]

        write_graph(options['habitat_file'] ,graphList)
        gprint('\nCalculating current flow centrality using Circuitscape...')
        
        memFlag = lu.call_circuitscape(cfg.CSPATH, outConfigFile)        
        
        outputFN = 'Circuitscape_network_branch_currents_cum.txt'
        currentList = path.join(OUTCENTRALITYDIR, outputFN)

        if not arcpy.Exists(currentList):
            write_graph(options['habitat_file'] ,graphList)
            gprint('\nCalculating current flow centrality using Circuitscape '
                   '(2nd try)...')
            memFlag = lu.call_circuitscape(cfg.CSPATH, outConfigFile)                    
            if not arcpy.Exists(currentList):
                lu.dashline(1)
                msg = ('ERROR: No Circuitscape output found.\n'
                       'It looks like Circuitscape failed.')
                arcpy.AddError(msg)
                lu.write_log(msg)
                exit(1)

                
        currents = load_graph(currentList,graphType='graph/network',
                              datatype='float64')

        numLinks = currents.shape[0]
        for x in range(0,numLinks):
            corex = currents[x,0]
            corey = currents[x,1]

            #linkId = LT[x,cfg.LTB_LINKID]
            row = lu.get_links_from_core_pairs(linkTable, corex, corey)
            #row = lu.get_linktable_row(linkId, linkTable)
            linkTable[row,cfg.LTB_CURRENT] = currents[x,2]

        coreCurrentFN = 'Circuitscape_network_node_currents_cum.txt'
        nodeCurrentList = path.join(OUTCENTRALITYDIR, coreCurrentFN)
        nodeCurrents = load_graph(nodeCurrentList,graphType='graph/network',
                              datatype='float64')

        numNodeCurrents = nodeCurrents.shape[0]
        rows = arcpy.UpdateCursor(coreCopy)
        row = rows.newRow()
        for row in rows:
            coreID = row.getValue(cfg.COREFN)
            for i in range (0, numNodeCurrents):
                if coreID == nodeCurrents[i,0]:
                    row.setValue("CF_Central", nodeCurrents[i,1])
                    break
            rows.updateRow(row)
            #row = rows.newRow()
        del row, rows
        gprint('Done with centrality calculations.')

        finalLinkTable = lu.update_lcp_shapefile(linkTable, lastStep=5,
                                                  thisStep=7)
        linkTableFile = path.join(cfg.DATAPASSDIR, "linkTable_s5_plus.csv")
        lu.write_link_table(finalLinkTable, linkTableFile, inLinkTableFile)
        linkTableFinalFile = path.join(cfg.OUTPUTDIR,
                                       cfg.PREFIX + "_linkTable_s5_plus.csv")
        lu.write_link_table(finalLinkTable,
                            linkTableFinalFile, inLinkTableFile)
        gprint('Copy of final linkTable written to '+
                          linkTableFinalFile)

        finalCoreFile = path.join(cfg.CORECENTRALITYGDB,
                                     cfg.PREFIX + '_Cores')
        #copy core area map to gdb.
        if not arcpy.Exists(cfg.CORECENTRALITYGDB):
            arcpy.CreateFileGDB_management(cfg.OUTPUTDIR,
                             path.basename(cfg.CORECENTRALITYGDB))
        arcpy.CopyFeatures_management(coreCopy, finalCoreFile)

        gprint('Creating shapefiles with linework for links.')
        lu.write_link_maps(linkTableFinalFile, step=7)

        # Copy final link maps to gdb and clean up.
        lu.copy_final_link_maps(step=7)

    # Return GEOPROCESSING specific errors
    except arcpy.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 7. Details follow.****')
        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except:
        lu.dashline(1)
        gprint('****Failed in step 7. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)

    return


def write_graph(filename,graphList):
    npy.savetxt(filename,graphList)
    return

def load_graph(filename,graphType,datatype):
    if not path.isfile(filename):
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
        gprint('****Failed in step 7. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)      
