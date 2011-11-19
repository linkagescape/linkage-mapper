to simplify core areas...
        S2CORE_RAS = "s2core_ras"
        gp.workspace = Cfg.SCRATCHDIR
        gp.CellSize = gp.Describe(Cfg.RESRAST).MeanCellHeight
        gp.extent = "MAXOF"
        gp.FeatureToRaster_conversion(Cfg.COREFC, Cfg.COREFN, S2CORE_RAS, gp.Cellsize)
        S2CORE_FC = "s2corefc" 
        gp.RasterToPolygon_conversion(S2CORE_RAS, S2CORE_FC, "SIMPLIFY")

then... create new field name same as core field name
        copy gridcode values to this field 
        
        
        
                                           
set options- write to project dir, let user know code will look there.

OPTIONS FILE FOUND
Options set to:

NO OPTIONS FILE FOUND
Default options will be used

trim down options- combine steps 1 and 2?


display cwd or ratio on top of corridors using sgticks?


        elif linktable.shape[1] == 13:
            outFile.write("#link,coreId1,coreId2,cluster1,cluster2,linkType,"
                           "eucDist,lcDist,eucAdj,cwdAdj,lcpLength,"
                           "cwdToEucRatio,cwdToPathRatio\n")
right now this is just for final linktable.  but calculated for all lcplines.

add lcplength lcp shapefile? right now is just ratio

# Min1=Rank(1,[acost1],[acost2],[acost3]....) where acostx is the
# accumulative cost to HCA x.
# Likewise, do:

# Min2=Rank(2,[acost1],[acost2],[acost3]....)
#can we do a rank with 100 acosts?  what about nodata? NO.
#So, just don't normalize and do regular min mosaic of cost distances.
Correct for euc dist?  cwd - pathlength or cwd/pathlength (minus doesn't give you wider corridors..)



                      

#For recursive deletes:
# from os.path import join, getsize
# for root, dirs, files in os.walk('python/Lib/email'):
    # print root, "consumes",
    # print sum(getsize(join(root, name)) for name in files),
    # print "bytes in", len(files), "non-directory files"
    # if 'CVS' in dirs:
        # dirs.remove('CVS')  # don't visit CVS directories

# In the next example, walking the tree bottom-up is essential: rmdir() doesn’t allow deleting a directory before the directory is empty:

# # Delete everything reachable from the directory named in "top",
# # assuming there are no symbolic links.
# # CAUTION:  This is dangerous!  For example, if top == '/', it
# # could delete all your disk files.
# import os
# for root, dirs, files in os.walk(top, topdown=False):
    # for name in files:
        # os.remove(os.path.join(root, name))
    # for name in dirs:
        # os.rmdir(os.path.join(root, name))


        
                                             
# Add defaults to dialog?  Or at least hints (bc's speed calcs).

#don't delete scratchdir if may be used by a next step

#set pinchpoint 0 current to nodata when input nodata
   #eventually do in cs
   #also in cs- show average current density in polygons
   #next release??
   
# add check for FID or ID as field name in toolbox.  Disallow.
# make sure core area field name is short integer

#Factorial least cost paths? (Consider all cores adjacent, add pixels)

#centrality- copy datapass/cores to outputgdb


#option 1 = centrality?
# Comment code
# Update demo

#barriers- allow map units or pixels?
# allow arbitrary resistance of restored habitat
# is it open in arcmap message
#Divide centrality by # cores?
# separate tools, and then a single all-in-one?
# capitalize constants like CONFIGDIR
# calculate ratios for all linktables
# remove stick shapefile at start of script?
# clean up network current maps
# eventually have nework connectedness measures.  these would have all pairs as sources and targets, 2 measures using eff resis with r's parm'd using eff resis and cwd vals

# put all outputgdb raster names in cfg
# results tab in ug

#UPDATE CIRCUITSCAPE WITH VIRAL'S NEW GAPDT
# Update circuitscape with multi component code

# global current flow? i.e. all pairs of cores with raster?

#!/usr/bin/env python2.5

"""Maps pinch points using Circuitscape given CWD calculations from
       s3_calcCwds.py.
"""

__filename__ = "s8_centrality.py"
__version__ = "0.6.6"

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
PREFIX = Cfg.PREFIX
    
def STEP8_calc_centrality():
    """ Experimental code to analyze network centrality using Circuitscape 
        given Linkage Mapper outputs

    """
    try:        
        lu.dashline(0)
        gprint('Running script ' + __filename__)
            
        gp.workspace = Cfg.SCRATCHDIR

        invalidFNs = {'fid','id','oid','shape'}
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
        
        
        
        
        #copy core area map to gdb.
        gp.createfilegdb(Cfg.OUTPUTDIR, path.basename(Cfg.CORECENTRALITYGDB))
        finalCoreFile = os.path.join(Cfg.CORECENTRALITYGDB,
                                     PREFIX + '_Cores')
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
            
    
    