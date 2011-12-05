##
## Circuitscape (C) 2008, Brad McRae and Viral B. Shah. 
##
## $Id: cs_compute.py 733 2010-12-04 16:11:42Z mcrae $
##

############## Enter True below to run in Stress Test mode
stressTest = False
stress_ncols = 4000
stress_nrows = 4000
#########################################################

import sys, time, string, os, math
import ConfigParser
import gc
import traceback
import wx

import numpy, scipy, pyamg

from sys import path
from string import split
from math import sqrt
from numpy import *
from scipy import sparse
from pyamg import *

from cs_util import *
from gapdt import *
import copy

pylab_available = False
try:
    import pylab
    import matplotlib
except ImportError:
    pylab_available = False

umfpack_available = False #Fixme- causes error on mac OS
try:
    import scipy.sparse.linalg.dsolve.umfpack as um
except ImportError:
    umfpack_available = False


# Set this to True to use -1 and +1 for currents
# Set this to False for the original
using_G_no_deleterow = True


print_timings_spaces = 0
print_timings = False

class cs_compute:
        
    def __init__(self, configFile, logger_func):
        gc.enable()
        self.state = {}
        self.state['amg_hierarchy'] = None
        self.options = readConfigFile(configFile)

###################################################################################
##############  OPTIONS FOR BETA CODE --EXPERIMENTAL AND NOT FOR USE ##############
        #self.options['data_type']='raster' # use 'network' for new network code
        
        self.options['write_baseline_results']=False
        if self.options['data_type']=='network': 
            self.options['write_baseline_results']=False
            self.options['graph_file'] = self.options['habitat_file']
            self.options['focal_node_file'] = self.options['point_file']

        self.options['write_graph'] = False #debug code 
        self.options['use_reclass_table'] = False
        self.options['reclass_file'] = './reclass.txt'        
        self.options['write_volt_drop_maps']=False
        self.options['set_null_voltages_to_nodata']=False
        self.options['write_max_cur_maps']=False #fixme: need to implement for network, also check memory error restart code 
        self.options['normalize_vdrop_maps']=False
#         self.options['use_included_pairs'] = True #debug code
#         self.options['included_pairs_file'] = './verify/1/include_matrix.txt' 
#         print'DEBUG CODE ACTIVATED. EXCLUDE SET TO',self.options['included_pairs_file']
#         
#         self.options['use_variable_source_strengths'] = True #debug code
#         self.options['variable_source_file'] = './verify/1/variable_source_list.txt' 
#         print'DEBUG CODE ACTIVATED. VARIABLE SOURCE FILE SET TO',self.options['variable_source_file']

#         self.options['use_mask'] = True #debug code
#         self.options['mask_file'] = './verify/1/mask.asc' 
#         print'DEBUG CODE ACTIVATED. MASK FILE SET TO',self.options['mask_file']
###################################################################################
        
        if umfpack_available:
            self.options['solver'] = 'umfpack'
        
        self.gapdt = gapdt()
        if pylab_available:
            self.fig = pylab.figure()

        global logger
        logger = logger_func

        global print_timings
        print_timings = self.options['print_timings']
        
        
    def print_timing(func):
        def wrapper(*arg):
            global print_timings_spaces
            print_timings_spaces +=  2
            t1 = time.time()
            res = func(*arg)
            t2 = time.time()
            print_timings_spaces -=  2
            if print_timings:
                print'%10d sec: ' % (t2-t1),
                for i in range(0,print_timings_spaces):
                    print" ",
                print'%s' % func.func_name
            return res
        return wrapper


    def cs_log(self, text,col):
        return
        (hours,mins,secs) = elapsed_time(self.state['lastUpdateTime'])
        if secs > 10: #FIXME- what is best time?
            self.state['lastUpdateTime'] = time.time()
            wx.SafeYield(None, True)  
            wx.GetApp().Yield(True)
        if logger:           
            logger(text,col)
        return


    @print_timing
    def compute(self):
        self.state['startTime'] = time.time()
        self.state['lastUpdateTime'] = time.time()        
        self.cs_log('',1)
        self.cs_log('',2)
        #Test write privileges by writing config file to output directory
        fileName = self.options['output_file']
        outputDir, outputFile = os.path.split(fileName)
        outputBase, outputExtension = os.path.splitext(outputFile)
        configFile = outputDir + '//' + outputBase + '.ini' 
        if os.path.isdir(outputDir):
            writeConfigFile(configFile, self.options)
        else:
            raise RuntimeError('Directory ' + outputDir + ' does not exist!')    
        
        if stressTest==True:
            print 'Calling stress test code'
            self.run_stress_test(stress_ncols,stress_nrows)
        elif self.options['data_type']=='network':
            result,solver_failed = self.network_module()
            self.logCompleteJob()
            return result,solver_failed #Fixme: add in solver failed check
        else:
            self.load_maps()
        
        if self.options['scenario']=='pairwise':
            resistances,solver_failed = self.pairwise_module(self.state['g_map'], self.state['poly_map'], self.state['points_rc'])
            self.logCompleteJob()
            return resistances,solver_failed     

        elif self.options['scenario'] == 'advanced':
            self.cs_log ('Calling solver module.',1)
            voltages,current_map,solver_failed = self.advanced_module(self.state['g_map'], self.state['poly_map'], self.state['source_map'], self.state['ground_map'],None,None,None,None,None)
            self.logCompleteJob()
            if solver_failed == True:
                print'Solver failed.\n'
            return voltages, solver_failed

        else:
            resistance_vector,solver_failed = self.one_to_all_module(self.state['g_map'], self.state['poly_map'], self.state['points_rc'])
            self.logCompleteJob()
            return resistance_vector,solver_failed    
    
    def grid_to_graph (self, x, y, node_map):
        return node_map[x, y] - 1
        
    def logCompleteJob(self):
        (hours,mins,secs) = elapsed_time(self.state['startTime'])
        if hours>0:
            self.cs_log('Job took ' + str(hours) +' hours ' + str(mins) + ' minutes to complete.',2)
        else:
            self.cs_log('Job took ' + str(mins) +' minutes ' + str(secs) + ' seconds to complete.',2)
 
  
    def get_overlap_polymap(self,point,point_map,poly_map_temp,new_poly_num): 
        point_poly = where(point_map == point, 1, 0) 
        poly_point_overlap = multiply(point_poly,poly_map_temp)
        overlap_vals = unique(asarray(poly_point_overlap))
        (rows) = where(overlap_vals>0)
        overlap_vals = overlap_vals[rows] #LIST OF EXISTING POLYGONS THAT OVERLAP WITH POINT
        for a in range (0, overlap_vals.size):
            poly_map_temp = where(poly_map_temp==overlap_vals[a],new_poly_num,poly_map_temp)
        poly_map_temp = where(point_map == point, new_poly_num, poly_map_temp) 
        return poly_map_temp


    @print_timing
    def network_module(self): 
        print 'Network mode- this is BETA CODE with no guarantees.'
        solver_failed = False
        (g_graph,nodeNames) = self.read_graph(self.options['graph_file'])
        C = self.gapdt.components(g_graph) 
        print max(C),'COMPONENTS:'
        print C
        if max(C)>1:
            raise RuntimeError('Input graph has more than one component.  Multiple components are not yet implemented in graph/network mode.')
        del C    
        G = self.laplacian(g_graph)
        del g_graph
        if self.options['scenario']=='advanced':
            (sources,grounds)= self.readSourcesGroundsNetwork(G, nodeNames, self.options['source_file'],self.options['ground_file'])
            #FIXME: retool readSourceGroundNodes to read both files and return complete sources and grounds
            result,solver_failed = self.advanced_module_network(G,sources,grounds,nodeNames)        

        else:
            if self.options['use_included_pairs']==True:
                self.state['includedPairs'] = self.readincludedPairs(self.options['included_pairs_file'])
            focalNodes = self.readFocalNodes(self.options['focal_node_file'])

            if self.options['scenario']=='pairwise':
                result,solver_failed = self.pairwise_module_network(G,focalNodes,nodeNames)
            else:
                result,solver_failed = self.one_to_all_module_network(G,focalNodes,nodeNames)
        return result,solver_failed #Fixme: need to check solver failed.
        
        
    def advanced_module_network(self,G,sources,grounds,nodeNames):
        (sources, grounds, finitegrounds) = self.resolve_conflicts(sources, grounds)
        try:
            voltages = self.multiple_solver(G, sources, grounds, finitegrounds)
            solver_failed = False
        except:
            solver_failed = True
        if self.options['write_cur_maps'] == True:
            (node_currents,branch_currents) = self.getCurrentsNetwork(G,voltages,finitegrounds)
            fileadd=''
            self.writeCurrentsNetwork(branch_currents, node_currents, nodeNames, fileadd)
        if self.options['write_volt_maps'] == True:
            self.writeVoltagesNetwork(voltages, nodeNames, fileadd='')
        
        return voltages, solver_failed
        
        
    def one_to_all_module_network(self,G,focalNodes,nodeNames):
        solver_failed = False
        if self.options['use_included_pairs']==True: #Prune points
            focalNodes = self.pruneIncludedPairsNetwork(focalNodes)          
            includedPairs = self.state['includedPairs'] 

        numpoints = focalNodes.size
        numnodes = G.shape[0]
        resistances = zeros((focalNodes.size,2),dtype = 'float64')
        resistances[:,0] = focalNodes
        if self.options['write_cur_maps'] == True:
            cum_node_currents = zeros((nodeNames.size,1),dtype = 'float64')
            cum_branch_currents = sparse.csr_matrix((G.shape))
        
        finitegrounds = [-9999]
        focalNodeLocs = self.namesToNodes(nodeNames,focalNodes) 
        
        sources = zeros(numnodes,dtype = 'float64')
        grounds = zeros(numnodes,dtype = 'float64')
        
        variableSources = self.getVariableSources(numnodes,focalNodes)
        
        if self.options['scenario'] == 'one-to-all':
            grounds[focalNodeLocs] = Inf
        else:
            sources = variableSources

        x = 0
        lastWriteTime = time.time()
        for i in range(0, numpoints):
            if self.options['use_included_pairs']==True: #Prune points
                groundsTemp=grounds
                sourcesTemp=sources
                for pair in range(0, focalNodes.size): #loop thru exclude[point,:], delete included pairs of focal point from point_map and points_rc_unique_temp 
                    if includedPairs[i+1,pair+1]==0 and i !=  pair:
                        dropNode = focalNodeLocs[pair]
                        groundsTemp[dropNode] = 0
                        sourcesTemp[dropNode] = 0
                        
            x = x+1
            (hours,mins,secs) = elapsed_time(self.state['startTime'])
            self.cs_log ('At ' + str(hours) +' hr ' + str(mins) + ' min solving focal node ' + str(x) + ' of '+ str(numpoints) + '.',1)

            node = focalNodeLocs[i]
            if self.options['scenario'] == 'one-to-all':
                grounds[node] = 0
                sources[node] = variableSources[node]
                try:
                    voltages = self.multiple_solver(G, sources, grounds, finitegrounds)
                    resistances[i,1] = voltages[node]
                except:
                    solver_failed = True
                    resistances[i,1] = -777

                grounds[node] = Inf
                sources[node] = 0
            else:
                grounds[node] = Inf
                sources[node] = 0
                try:
                    voltages = self.multiple_solver(G, sources, grounds, finitegrounds)
                except:
                    solver_failed = True
                    resistances[i,1] = -777
                grounds[node] = 0
                sources[node] = variableSources[node]               

            (hours,mins,secs) = elapsed_time(lastWriteTime)
            if secs > 120: 
                lastWriteTime = time.time()
                self.writeResistancesOneToAll(resistances,'_incomplete')       
            
            if self.options['write_cur_maps'] == True:
                (node_currents,branch_currents) = self.getCurrentsNetwork(G,voltages,finitegrounds)
                cum_node_currents=cum_node_currents+node_currents
                cum_branch_currents=cum_branch_currents+branch_currents
                if self.options['write_cum_cur_map_only']==False:     
                    self.writeCurrentsNetwork(branch_currents, node_currents, nodeNames, fileadd=str(int(focalNodes[i])))
            if self.options['write_volt_maps'] == True:
                self.writeVoltagesNetwork(voltages, nodeNames, fileadd=str(int(focalNodes[i])))
                
        if self.options['write_cur_maps'] == True:
            self.writeCurrentsNetwork(cum_branch_currents, cum_node_currents, nodeNames, 'cum')
            
        self.writeResistancesOneToAll(resistances,'')
        
        #Need to add row and column headers and write currents and voltagesto disk

        print 'focal nodes'
        print focalNodes
        return resistances,solver_failed


    def getVariableSources(self,numnodes,focalNodes):       
        variableSources = ones(numnodes,dtype = 'float64')
        if self.options['use_variable_source_strengths']==True:
            pointStrengths = self.readPointStrengths(self.options['variable_source_file']) 
            variableSourceNames = pointStrengths[:,0]
            variableSourceNodes = self.namesToNodes(focalNodes,variableSourceNames)
            print 'strengths'
            print pointStrengths
            print 'nodes'
            print variableSourceNodes

            try:
                for i in range (0,len(variableSourceNames)):
                    variableSources[variableSourceNodes[i]] = pointStrengths[variableSourceNodes[i],1]
            except IndexError:
                raise RuntimeError('Error assinging variable source strengths. Please make sure focal node names match.')                
        return variableSources

    
    def writeResistancesOneToAll(self,resistances,string):
        fileName = self.options['output_file']
        outputDir, outputFile = os.path.split(fileName)
        outputBase, outputExtension = os.path.splitext(outputFile)
        if self.options['scenario'] == 'one-to-all':
            outputFile = outputDir + '//' + outputBase + '_resistances' + string + outputExtension
        else:
            outputFile = outputDir + '//' + outputBase + '_results' + string + outputExtension     
        savetxt (outputFile, resistances) 
        
        #remove partial result file        
        if string=='':        
            if self.options['scenario'] == 'one-to-all':
                oldFile = outputDir + '//' + outputBase + '_resistances_incomplete' + string + outputExtension
            else:
                oldFile = outputDir + '//' + outputBase + '_results_incomplete' + string + outputExtension
            try:
                os.remove(oldFile)
            except:
                pass 
        return
    
   
    def pairwise_module_network(self,G,focalNodes,nodeNames):
        Gtemp=G.todense()
        solver_failed = False
        if self.options['use_included_pairs']==True: #Prune points
            focalNodes = self.pruneIncludedPairsNetwork(focalNodes)          
            includedPairs = self.state['includedPairs'] 
        else:
            includedPairs = ones((focalNodes.size+1,focalNodes.size+1),dtype = 'int32')

        numpoints = focalNodes.size
        resistances = -1*ones((focalNodes.size,focalNodes.size),dtype = 'float64')
        if self.options['write_cur_maps'] == True:
            cum_node_currents = zeros((G.shape[0],1),dtype = 'float64')
            cum_branch_currents = sparse.csr_matrix((G.shape))
        if (self.options['write_cur_maps'] == True) or (self.options['write_volt_maps'] == True) or (self.options['use_included_pairs']==True):        
            useResistanceCalcShortcut = False
        else:
            useResistanceCalcShortcut = True
            print 'Using shortcut'
            shortcutResistances = -1 * ones((numpoints, numpoints), dtype = 'float64') #temp bhm

        if (useResistanceCalcShortcut==True):
            voltmatrix = zeros((numpoints,numpoints),dtype = 'float64')     
        dstPoint = 0
        anchorPoint = 0            
        
        x = 0
        for i in range(0, numpoints):

#moved below to j loop           
#             x = x+1
#             (hours,mins,secs) = elapsed_time(self.state['startTime'])
#             if useResistanceCalcShortcut==True:
#                 y = numpoints
#                 self.cs_log ('At ' + str(hours) +' hr ' + str(mins) + ' min solving focal node ' + str(x) + ' of '+ str(y) + '.',1)
#             else:
#                 y = numpoints*(numpoints-1)/2
#                 self.cs_log ('At ' + str(hours) +' hr ' + str(mins) + ' min solving focal pair ' + str(x) + ' of '+ str(y) + '.',1)
           
            if range(i, numpoints) == []:
                break
            if (useResistanceCalcShortcut==True) and (dstPoint>0):
                break #No need to continue, we've got what we need to calculate resistances
         
            dstPoint = dstPoint+1
            
            if using_G_no_deleterow: #Fixme: right now this is assumed in network mode
                dst = self.nameToNode(nodeNames,focalNodes[i])
                G_dst_dst = G[dst, dst] 
                G[dst,dst] = 0
                self.state['amg_hierarchy'] = None
                gc.collect()
                self.create_amg_hierarchy(G)         
            
            for j in range(i+1, numpoints):
                x = x+1
                (hours,mins,secs) = elapsed_time(self.state['startTime'])
                if useResistanceCalcShortcut==True:
                    y = numpoints
                    self.cs_log ('At ' + str(hours) +' hr ' + str(mins) + ' min solving focal node ' + str(x) + ' of '+ str(y) + '.',1)
                else:
                    y = numpoints*(numpoints-1)/2
                    self.cs_log ('At ' + str(hours) +' hr ' + str(mins) + ' min solving focal pair ' + str(x) + ' of '+ str(y) + '.',1)            
                if includedPairs[i+1,j+1]==1:
                    src = self.nameToNode(nodeNames,focalNodes[j])
                    try:
                        voltages = self.single_ground_solver(G, src, dst)

                        resistances[i, j] = voltages[src] - voltages[dst]
                        resistances[j, i] = voltages[src] - voltages[dst]
                    except:
                        solver_failed = True
                        resistances[i, j] = -777
                        resistances[j, i] = -777
    
                    if (useResistanceCalcShortcut==True) and (dstPoint==1) and (solver_failed==False): #this occurs for first i that is in component
                        anchorPoint = i #for use later in shortcult resistance calc
                        voltmatrix = self.getVoltmatrixNetwork(i,j,numpoints,voltages,resistances,voltmatrix,focalNodes,nodeNames) 
    
                    if (self.options['write_cur_maps'] == True) and (solver_failed==False):
                        finitegrounds = [-9999]
                        (node_currents,branch_currents) = self.getCurrentsNetwork(G,voltages,finitegrounds)
                        cum_node_currents=cum_node_currents+node_currents
                        cum_branch_currents=cum_branch_currents+branch_currents
                        if self.options['write_cum_cur_map_only']==False:                
                            fileadd=str(int(focalNodes[i])) + '_' + str(int(focalNodes[j]))
                            self.writeCurrentsNetwork(branch_currents, node_currents, nodeNames, fileadd)
    
                    if (self.options['write_volt_maps'] == True) and (solver_failed==False):
                        fileadd=str(int(focalNodes[i])) + '_' + str(int(focalNodes[j]))

                        self.writeVoltagesNetwork(voltages, nodeNames, fileadd)
    
                    if solver_failed==True:
                        solver_failed = False
                        solver_failed_somehwere = True                 
            if (useResistanceCalcShortcut==True) and (i==anchorPoint): #this happens once per component. Anchorpoint is the first i in component
                shortcutResistances = self.getShortcutResistances(anchorPoint,voltmatrix,numpoints,resistances,shortcutResistances)
             
            if using_G_no_deleterow:
                G[dst,dst] = G_dst_dst
             
        if self.options['write_cur_maps'] == True:
          
            fileadd='cum'
            self.writeCurrentsNetwork(cum_branch_currents, cum_node_currents, nodeNames, fileadd)
            
        if useResistanceCalcShortcut==True:
            resistances = shortcutResistances
        for i in range(0,numpoints):
            resistances[i, i] = 0
        resistances = self.writeResistances(focalNodes, resistances)
        Gtemp=G.todense()
        return resistances,solver_failed


    def writeCurrentsNetwork(self, branch_currents, node_currents, nodeNames, fileadd):
        outputDir, outputFile = os.path.split(self.options['output_file'])
        outputBase, outputExtension = os.path.splitext(outputFile)
        if branch_currents!=None:
            file = outputDir + '//' + outputBase + '_branch_currents_' + fileadd + '.txt'
            self.writeGraph(file,branch_currents,nodeNames)
        outputNodeCurrents=zeros((len(node_currents),2),dtype='float64')
        outputNodeCurrents[:,0]=nodeNames[:]
        try:
            outputNodeCurrents[:,1]=node_currents[:,0]
        except:
            outputNodeCurrents[:,1]=node_currents[:]
        file = outputDir + '//' + outputBase + '_node_currents_' + fileadd + '.txt'
        savetxt(file,outputNodeCurrents)      
        return


    def writeVoltagesNetwork(self, voltages, nodeNames, fileadd):
        outputVoltages=zeros((len(voltages),2),dtype='float64')
        outputVoltages[:,0]=nodeNames[:]
        outputVoltages[:,1]=voltages[:]
        outputDir, outputFile = os.path.split(self.options['output_file'])
        outputBase, outputExtension = os.path.splitext(outputFile)
        file = outputDir + '//' + outputBase + '_voltages_' + fileadd + '.txt'
        savetxt(file,outputVoltages)      
        return        

        
    def getCurrentsNetwork(self,G,voltages,finitegrounds):
        G =  G.tocoo()
        node_currents = self.get_node_currents(voltages, G, finitegrounds)
        node_currents_col = zeros((node_currents.shape[0],1),dtype = 'float64')
        node_currents_col[:,0] = node_currents[:]
        branch_currents = self.get_branch_currents(G,voltages,True) 


        branch_currents = absolute(branch_currents) 
        G = G.tocsr()
        node_currents = node_currents_col
        return node_currents,branch_currents
    
    
    def getVoltmatrixNetwork(self,i,j,numpoints,voltages,resistances,voltmatrix,focalNodes,nodeNames):                                            
       voltvector = zeros((numpoints,1),dtype = 'float64')  
       for point in range(1,numpoints):
           node=self.nameToNode(nodeNames,focalNodes[point])
           voltageAtPoint = voltages[node] 
           voltageAtPoint = 1-(voltageAtPoint/resistances[i, j])
           voltvector[point] = voltageAtPoint
       voltmatrix[:,j] = voltvector[:,0] 
       return voltmatrix             
             
             
    def nameToNode(self, nodeNames,name):
        nodeNames = nodeNames.tolist()
        node = nodeNames.index(name)
        return node


    def namesToNodes(self, nodeNames,names):
        nodeNames = nodeNames.tolist()
        nodes = zeros(len(names),dtype = 'int32')
#         print 'nodes'
#         print nodes
#         print 'nodenames',nodeNames
        for i in range (0,len(names)):
           nodes[i] = nodeNames.index(names[i])
        return nodes

        
    def read_graph(self, filename):
        graphList = self.load_graph(filename,graphType='graph/network',datatype='float64')

        try:
            zeros_in_resistance_graph = False
            print graphList
            nodes = self.deletecol(graphList,2) 
            nodeNames = unique(asarray(nodes))
            nodes[where(nodes>= 0)] = self.gapdt.relabel(nodes[where(nodes>= 0)], 0)
            node1 = nodes[:,0]
            node2 = nodes[:,1]
            data = graphList[:,2]
            ######################## Reclassification code
            if self.options['use_reclass_table']==True:
                try:
                    reclassTable = self.readPointStrengths(self.options['reclass_file'])    
                except:
                    raise RuntimeError('Error reading reclass table.  Please check file format.')
                for i in range (0,reclassTable.shape[0]):
                    data = where(data==reclassTable[i,0],reclassTable[i,1],data)
                print'\n***** Reclassified habitat graph using', self.options['reclass_file'],'*****'
            ########################            
            if self.options['habitat_map_is_resistances'] == True:
                zeros_in_resistance_graph = (where(data==0, 1, 0)).sum() > 0
                print 'zeros', zeros_in_resistance_graph
                conductances = 1/data
            else:
                conductances = data
            numnodes = nodeNames.shape[0]
            G = sparse.csr_matrix((conductances, (node1, node2)), shape = (numnodes, numnodes))
            g_graph = G + G.T       
            print 'Graph'
            print g_graph
            print 'node names'
            print nodeNames
    
        except:
            raise RuntimeError('Error processing graph/network file.  Please check file format')
        if zeros_in_resistance_graph == True:
            raise RuntimeError('Error: zero resistance values are not currently allowed in habitat network/graph input file.')
        return g_graph, nodeNames


    def readFocalNodes(self, filename):
        focalNodes = self.load_graph(filename,graphType='focal node',datatype='int32')
        try:    
            if filename==self.options['graph_file']:#If graph was used as focal node file, then just want first two columns for focalNodes.
                focalNodes = self.deletecol(focalNodes, 2)
            focalNodes = unique(asarray(focalNodes))
            print 'focal nodes'
            print focalNodes
        except:
            raise RuntimeError('Error processing focal node file.  Please check file format')
        return focalNodes
        
        
    def load_graph(self,filename,graphType,datatype):
        if os.path.isfile(filename)==False:
            raise RuntimeError('File "'  + filename + '" does not exist')
        f = open(filename, 'r')
#        firstLine = f.readline()  #Can't use this- I think it stores numeric data as a string
#        if type(firstLine).__name__=='str':  #This fails sometimes with numeric data.

#         if firstLineType == 'str':
#             try:
#                 graphObject = loadtxt(filename, dtype = datatype, comments='#', skiprows=1) 
#                 print 'str!'
#             except:
#                 try:
#                     graphObject = loadtxt(filename, dtype = datatype, comments='#', skiprows=1, delimiter=',')
#                 except:
#                     raise RuntimeError('Error reading',type,'file.  Please check file format')
#         else:
        try:
            graphObject = loadtxt(filename, dtype = 'Float64', comments='#') 
        except:
            try:
                graphObject = loadtxt(filename, dtype = 'Float64', comments='#', delimiter=',')
            except:
                raise RuntimeError('Error reading',type,'file.  Please check file format')
        return graphObject
        

    def readSourcesGroundsNetwork(self, G, nodeNames, sourceFile,groundFile):
        if os.path.isfile(sourceFile)==False:
            raise RuntimeError('File "'  + sourceFile + '" does not exist')   
        if os.path.isfile(groundFile)==False:
            raise RuntimeError('File "'  + groundFile + '" does not exist')               
        rawSources = loadtxt(sourceFile, dtype = 'float64')
        rawGrounds = loadtxt(groundFile, dtype = 'float64')

        if self.options['ground_file_is_resistances']==True:
            rawGrounds[:,1] = 1 / rawGrounds[:,1]
        if self.options['use_direct_grounds']==True:
            rawGrounds[:,1] = where(rawGrounds[:,1],Inf,0)
        
        numnodes = G.shape[0]
        sourceNodes = self.namesToNodes(nodeNames,rawSources[:,0]) 
        sources = zeros((numnodes),dtype = 'float64')
        sources[sourceNodes[:]] = rawSources[:,1] 
        
        groundNodes = self.namesToNodes(nodeNames,rawGrounds[:,0]) 
        grounds = zeros((numnodes),dtype = 'float64')
        grounds[groundNodes[:]] = rawGrounds[:,1] 

        conflicts = logical_and(sources,grounds)
        if self.options['remove_src_or_gnd']=='rmvsrc':
            sources = where(conflicts,0,sources)
        elif self.options['remove_src_or_gnd']=='rmvgnd':
            grounds = where(conflicts,0,grounds)
        elif self.options['remove_src_or_gnd']=='rmvall':
            sources = where(conflicts,0,sources)
            grounds = where(conflicts,0,grounds)
        if size(where(sources)) == 0:
            raise RuntimeError('No valid sources detected. Please check source file') 
        if size(where(grounds)) == 0:
            raise RuntimeError('No valid grounds detected. Please check ground file') 
        return sources, grounds
        return sources,grounds


    @print_timing
    def one_to_all_module(self, g_map, poly_map, points_rc):
        lastWriteTime = time.time()

        if self.options['use_included_pairs']==True: #Prune points
            points_rc = self.pruneIncludedPairs(points_rc)          
            includedPairs = self.state['includedPairs'] 
        point_ids = unique(asarray(points_rc[:,0]))
        points_rc_unique = self.get_points_rc_unique(point_ids,points_rc)      
        
        resistance_vector = zeros((point_ids.size,2),float)
        solver_failed_somewhere = False
        if self.options['write_cur_maps'] == True:
            cum_current_map = zeros((self.state['nrows'],self.state['ncols']),dtype = 'float64')         
        if self.options['write_max_cur_maps']==True:
            max_current_map=cum_current_map
        oneToAllStreamline = False        
        if self.options['use_included_pairs']==False: #Will do this each time later if using included pairs
            point_map = numpy.zeros((self.state['nrows'],self.state['ncols']),int)
            point_map[points_rc[:,1],points_rc[:,2]] = points_rc[:,0]
       
            #combine point map and poly map
            poly_map_temp = self.get_poly_map_temp(poly_map,point_map,point_ids,None,None)
            unique_point_map = numpy.zeros((self.state['nrows'],self.state['ncols']),int)
            unique_point_map[points_rc_unique[:,1],points_rc_unique[:,2]] = points_rc_unique[:,0]

            (strengthMap,strengths_rc) = self.getstrengthMap(points_rc_unique,self.state['pointStrengths'])            

            #Are all points in same component?  If so, can streamline.  Create G here, just once
#             FIXME: place code below into module 
            node_map = self.construct_node_map(g_map, poly_map_temp) #******WANT POLYS BURNED IN 
            (component_map, components) = self.construct_component_map(g_map, node_map)
            pointComponents = where(unique_point_map, component_map, 0)#TODO use this to choose component to operate on below, unless all points in same component
            self.cs_log('Graph has ' + str(node_map.max()) + ' nodes and '+ str(components.max()) + ' components.',2)
            del node_map
            uniquePointComponentList = unique(pointComponents)
            numComponentsWithPoints = uniquePointComponentList.size
            if uniquePointComponentList[0]==0:
                numComponentsWithPoints = numComponentsWithPoints-1
            if numComponentsWithPoints==1:
                componentWithPoints = max(uniquePointComponentList)
                (G, node_map) = self.node_pruner(g_map, poly_map_temp, component_map, componentWithPoints) #FIXME drops polymaptemp node 1... #We'll keep node_map around in this case.  Only need to create it and G once.
                oneToAllStreamline = True
            else:
                del component_map, components  
        else:
            node_map = self.construct_node_map(g_map, poly_map) #******WANT POLYS BURNED IN 
            (component_map, components) = self.construct_component_map(g_map, node_map)
            self.cs_log('Graph has ' + str(node_map.max()) + ' nodes and '+ str(components.max()) + ' components.',2)
            del node_map, component_map, components

        for i in range(0, point_ids.size): #These are the 'src' nodes, i.e. the 'one' in all-to-one and one-to-all
            (hours,mins,secs) = elapsed_time(self.state['startTime'])
            self.cs_log ('At ' + str(hours) +' hr ' + str(mins) + ' min solving focal node ' + str(i+1) + ' of ' + str(point_ids.size) + '.',1)

            if self.options['use_included_pairs']==True: # Done above otherwise    
                #######################   
                points_rc_unique_temp = copy.copy(points_rc_unique)
                point_map = numpy.zeros((self.state['nrows'],self.state['ncols']),int)
                point_map[points_rc[:,1],points_rc[:,2]] = points_rc[:,0]       

                for pair in range(0, point_ids.size): #loop thru exclude[point,:], delete included pairs of focal point from point_map and points_rc_unique_temp 
                    if includedPairs[i+1,pair+1]==0 and i !=  pair:
                            id = point_ids[pair]
                            point_map = where(point_map==id,0,point_map)
                            points_rc_unique_temp[pair,0] = 0 #point will not be burned in to unique_point_map

                poly_map_temp = self.get_poly_map_temp2(poly_map,point_map,points_rc_unique_temp,includedPairs,i)

                unique_point_map = numpy.zeros((self.state['nrows'],self.state['ncols']),int)
                unique_point_map[points_rc_unique_temp[:,1],points_rc_unique_temp[:,2]] = points_rc_unique_temp[:,0]        

                (strengthMap,strengths_rc) = self.getstrengthMap(points_rc_unique_temp,self.state['pointStrengths'])
                ###########################################
                
            src = point_ids[i]
            sumConnections = unique_point_map.sum()
            sumConnections = sumConnections.sum()
            if sumConnections!= point_ids[i]: #If there are points to connect with src point
                if self.options['scenario'] == 'one-to-all':
                    if self.options['use_variable_source_strengths']==True:
                        strength = strengths_rc[i,0]
                    else:
                        strength = 1
                    source_map = where(unique_point_map == src,strength,0)
                    ground_map = point_map
                    ground_map =  where(unique_point_map == src,0,unique_point_map)
                    ground_map =  where(ground_map,Inf,0) 
                    self.options['remove_src_or_gnd'] = 'rmvgnd'
                else:
                    if self.options['use_variable_source_strengths']==True:
                        source_map = where(unique_point_map == src,0,strengthMap)
                    else:
                        source_map = where(unique_point_map,1,0)
                        source_map = where(unique_point_map == src,0,source_map)
                        
                    ground_map = point_map
                    ground_map =  where(unique_point_map == src,1,0)
                    ground_map =  where(ground_map,Inf,0)   

                    self.options['remove_src_or_gnd'] = 'rmvsrc'
                #FIXME: right now one-to-all *might* fail if there is just one node that is not grounded (I haven't encountered this problem lately BHM Nov 2009).

                if oneToAllStreamline==False:
                    G = None
                    node_map = None
                    component_map = None
                    componentWithPoints = None
                    
                resistance,current_map,solver_failed = self.advanced_module(self.state['g_map'], poly_map_temp, source_map, ground_map, src, G, node_map, component_map, componentWithPoints)
                if solver_failed == False:
                    if self.options['write_cur_maps'] == True:
                        cum_current_map = cum_current_map+current_map
                        if self.options['write_max_cur_maps']==True:
                            max_current_map=maximum(max_current_map,current_map)
                else:
                    print'Solver failed for at least one focal node.  \nFocal nodes with failed solves will be marked with value of -777 \nin output resistance list.\n'
    
                resistance_vector[i,0] = src
                resistance_vector[i,1] = resistance
                    
                if solver_failed==True:
                    solver_failed_somewhere = True
            else:
                resistance_vector[i,0] = src
                resistance_vector[i,1] = -1            

            (hours,mins,secs) = elapsed_time(lastWriteTime)
            if secs > 120: 
                lastWriteTime = time.time()
                self.writeResistancesOneToAll(resistance_vector,'_incomplete')
      
        if solver_failed_somewhere==False:
            if self.options['write_cur_maps'] == True:
                self.write_aaigrid('cum_curmap', '', cum_current_map)
                if self.options['write_max_cur_maps']==True:
                    self.write_aaigrid('max_curmap', '', max_current_map)

        #remove partial result file        
         
        self.writeResistancesOneToAll(resistance_vector,'')
       
        return resistance_vector,solver_failed_somewhere 

    def get_poly_map_temp(self,poly_map,point_map,point_ids,includedPairs,point1):
        if poly_map == []:
            poly_map_temp = point_map
        else:
            poly_map_temp = poly_map
            new_poly_num = numpy.max(poly_map)
            for point2 in range(0, point_ids.shape[0]):
                if (self.options['use_included_pairs']==True) and (includedPairs[point1+1,point2+1] == 1):
                   new_poly_num = new_poly_num+1
                   poly_map_temp = self.get_overlap_polymap(point_ids[point2],point_map,poly_map_temp, new_poly_num) 
                else: 
                    new_poly_num = new_poly_num+1
                    poly_map_temp = self.get_overlap_polymap(point_ids[point2],point_map,poly_map_temp, new_poly_num)                 
                
        return poly_map_temp

    #FIXME: Reconcile this module with the one above, consolidate.
    def get_poly_map_temp2(self,poly_map,point_map,points_rc_unique_temp,includedPairs,i):
        if poly_map == []:
            poly_map_temp = point_map
        else:
            poly_map_temp = poly_map
            new_poly_num = numpy.max(poly_map)
            #burn in src point to polygon map
            poly_map_temp = self.get_overlap_polymap(points_rc_unique_temp[i,0],point_map,poly_map_temp, new_poly_num)                     
            for point in range(0, points_rc_unique_temp.shape[0]): #burn in dst points to polygon map
                if includedPairs[i+1,point+1] == 1:  
                   new_poly_num = new_poly_num+1
                   poly_map_temp = self.get_overlap_polymap(points_rc_unique_temp[point,0],point_map,poly_map_temp, new_poly_num) 
        return poly_map_temp


    def pairwise_module(self, g_map, poly_map, points_rc):
        if self.options['write_cur_maps'] == True:
            cum_current_map = zeros((self.state['nrows'],self.state['ncols']),dtype = 'float64') 
            if self.options['write_max_cur_maps']==True:
                max_current_map=cum_current_map
            else:
                max_current_map=[]
        else:
            cum_current_map = []
            max_current_map=[]
        if self.options['write_volt_drop_maps'] == True:
            cum_vdrop_map = zeros((self.state['nrows'],self.state['ncols']),dtype = 'float64') 
        else:
            cum_vdrop_map = []
            
            

        # If there are no focal regions, pass all points to single_ground_all_pair_resistances,
        # otherwise, pass one point at a time.
        if self.options['point_file_contains_polygons']==False:
            if points_rc.shape[0]!= (unique(asarray(points_rc[:,0]))).shape[0]:
                raise RuntimeError('At least one focal node contains multiple cells.  If this is what you really want, then choose focal REGIONS in the pull-down menu') 
            
            else:
                if self.options['use_included_pairs']==True:
                    points_rc = self.pruneIncludedPairs(points_rc)
#                     includedPairs = self.state['includedPairs']
#                 else:
#                     numpoints = points_rc.shape[0]
#                     includedPairs = ones((numpoints+1,numpoints+1),dtype = 'int32')
                
                reportStatus = True
                try:
                    (resistances,cum_current_map,max_current_map,cum_vdrop_map,solver_failed) = self.single_ground_all_pair_resistances(g_map, poly_map, points_rc,cum_current_map,max_current_map,cum_vdrop_map,reportStatus)
                except MemoryError: #Give it a try, but starting again never seems to helps even from GUI.
                    self.enable_low_memory(True) #This doesn't seem to really clear out memory or truly restart.
                    if self.options['write_cur_maps'] == True:
                        cum_current_map = zeros((self.state['nrows'],self.state['ncols']),dtype = 'float64') 
                        if self.options['write_max_cur_maps']==True:
                            max_current_map=cum_current_map
                        else:
                            max_current_map=[]
                    else:
                        cum_current_map = []
                        max_current_map=[]
                        
                    #Note: This does not go through when it should.
                    (resistances,cum_current_map,max_current_map,cum_vdrop_map,solver_failed) = self.single_ground_all_pair_resistances(g_map, poly_map, points_rc,cum_current_map,max_current_map,cum_vdrop_map,reportStatus)
                if solver_failed == True:
                    print'Solver failed for at least one focal node pair.  \nPairs with failed solves will be marked with value of -777 \nin output resistance matrix.\n'

                point_ids = points_rc[:,0]

        else:
            if self.options['use_included_pairs']==True:
                points_rc = self.pruneIncludedPairs(points_rc)
                includedPairs = self.state['includedPairs']
            else:
                numpoints = points_rc.shape[0]
#Nov13_2010                includedPairs = ones((numpoints+1,numpoints+1),dtype = 'int32')
#Nov13_2010                numpoints = points_rc.shape[0]
#Nov13_2010                includedPairs = ones((numpoints+1,numpoints+1),dtype = 'int32')
                point_ids = unique(asarray(points_rc[:,0]))
                points_rc_unique = self.get_points_rc_unique(point_ids,points_rc)#Nov13_2010   #Fixme: can just use point ids for index size
                numUniquepoints= points_rc_unique.shape[0]#Nov13_2010
                includedPairs = ones((numUniquepoints+1,numUniquepoints+1),dtype = 'int32')#Nov13_2010

            point_map = numpy.zeros((self.state['nrows'],self.state['ncols']),int)
            point_map[points_rc[:,1],points_rc[:,2]] = points_rc[:,0]

            point_ids = unique(asarray(points_rc[:,0]))
            points_rc_unique = self.get_points_rc_unique(point_ids,points_rc)

            #numUniquepoints= points_rc_unique.shape[0]#Nov13_2010
            #includedPairs = ones((numUniquepoints+1,numUniquepoints+1),dtype = 'int32')#Nov13_2010      
      

            resistances = -1*ones((point_ids.size,point_ids.size),dtype = 'float64')
            x = 0
            for i in range(0, point_ids.size-1):
                for j in range(i+1, point_ids.size):
                    if includedPairs[i+1,j+1]==1:                
                        if poly_map == []:
                            poly_map_temp = zeros((self.state['nrows'],self.state['ncols']),int)
                            new_poly_num = 1
                        else:
                            poly_map_temp = poly_map
                            new_poly_num = numpy.max(poly_map)+1
                        point = point_ids[i]
                        poly_map_temp = self.get_overlap_polymap(point_ids[i],point_map,poly_map_temp, new_poly_num) 
                        poly_map_temp = self.get_overlap_polymap(point_ids[j],point_map,poly_map_temp, new_poly_num+1) 
    
                        #Get first instance of each point in points_rc
                        points_rc_temp = zeros((2,3),int)
                        points_rc_temp[0,:] = points_rc_unique[i,:]
                        points_rc_temp[1,:] = points_rc_unique[j,:]
    
                        numpoints = point_ids.size
                        x = x+1
                        y = numpoints*(numpoints-1)/2
                        (hours,mins,secs) = elapsed_time(self.state['startTime'])
                        self.cs_log ('At ' + str(hours) +' hr ' + str(mins) + ' min solving focal pair ' + str(x) + ' of '+ str(y) + '.',1)
                        reportStatus = False
                        
                        (pairwise_resistance,cum_current_map,max_current_map,cum_vdrop_map,solver_failed) = self.single_ground_all_pair_resistances(g_map, poly_map_temp, points_rc_temp,cum_current_map,max_current_map,cum_vdrop_map,reportStatus)
    
                        del poly_map_temp
                        if solver_failed == True:
                            print'Solver failed for at least one focal node pair.  \nPairs with failed solves will be marked with value of -777 \nin output resistance matrix.\n'
    
                        resistances[i,j] = pairwise_resistance[0,1]
                        resistances[j,i] = pairwise_resistance[0,1]

        for i in range(0,resistances.shape[0]): #Set diagonal to zero
            resistances[i, i] = 0

        #Add row and column headers and write resistances to disk
        resistances = self.writeResistances(point_ids, resistances)

        #BASELINE RESULTS FOR NETWORK CODE
        if self.options['write_baseline_results']==True:
            cum_currents=cum_current_map.flatten()
            fileadd='cum_baseline'
            nodeNames=arange(1,11) #Change for each scenario
            self.writeCurrentsNetwork(None, cum_currents, nodeNames, fileadd)
        ############
        
        if self.options['write_cur_maps'] == True:
            if solver_failed==False:
                if self.options['log_transform_maps'] == True:
                    cum_current_map = where(cum_current_map>0,log10(cum_current_map),self.state['nodata']) 
                    if self.options['write_max_cur_maps']==True:
                        max_current_map = where(max_current_map>0,log10(max_current_map),self.state['nodata']) 
                self.write_aaigrid('cum_curmap', '', cum_current_map)
                if self.options['write_max_cur_maps']==True:      
                   self.write_aaigrid('max_curmap', '', max_current_map)

        if self.options['write_volt_drop_maps']==True:  #Fixme: need to add vdrop maps into one-to-all mode
            if solver_failed==False:
                self.write_aaigrid('cum_volt_gradient_map', '', cum_vdrop_map)

        return resistances,solver_failed
    
    
    @print_timing
    def single_ground_all_pair_resistances(self, g_map, poly_map, points_rc,cum_current_map,max_current_map,cum_vdrop_map,reportStatus):
        lastWriteTime = time.time()
        numpoints = points_rc.shape[0]
        if (self.options['use_included_pairs']==False) or (self.options['point_file_contains_polygons']==True):
            includedPairs = ones((numpoints+1,numpoints+1),dtype = 'int32')
        else:
            includedPairs = self.state['includedPairs']
        
        if (self.options['point_file_contains_polygons']==True) or (self.options['write_volt_drop_maps'] == True) or (self.options['write_cur_maps'] == True) or (self.options['write_volt_maps'] == True) or (self.options['use_included_pairs']==True): 
           useResistanceCalcShortcut = False
        else:     
           useResistanceCalcShortcut = True #BHM We use this when there are no focal regions.  It saves time when we are also not creating maps
           shortcutResistances = -1 * ones((numpoints, numpoints), dtype = 'float64') #temp bhm
           
        solver_failed_somewhere = False
        node_map = self.construct_node_map(g_map, poly_map) #******WANT POLYS BURNED IN 
        if self.options['write_graph']==True:
            focalNodes = zeros((points_rc.shape[0],1),dtype = 'int32')
            
            for i in range(points_rc.shape[0]-1,-1,-1):
                focalNodes[i] = self.grid_to_graph (points_rc[i,1], points_rc[i,2], node_map)            
                if focalNodes[i]==-1:
                    focalNodes = self.deleterow(focalNodes, i)
            savetxt(self.options['focal_node_file'],focalNodes)
            
        (component_map, components) = self.construct_component_map(g_map, node_map)
        
        self.cs_log('Graph has ' + str(node_map.max()) + ' nodes and '+ str(components.max())+ ' components.',2)
        resistances = -1 * ones((numpoints, numpoints), dtype = 'float64')         #Inf creates trouble in python 2.5 on Windows. Use -1 instead.
        
        x = 0
        for c in range(1,int(components.max()+1)):
            points_in_this_component = self.checkPointsInComponent(c,numpoints,components,points_rc,node_map)
                        
            if points_in_this_component:
                (G, local_node_map) = self.node_pruner(g_map, poly_map, component_map, c)
                if c==int(components.max()):
                    del component_map 
                
                if (useResistanceCalcShortcut==True):
                    voltmatrix = zeros((numpoints,numpoints),dtype = 'float64')     #For resistance calc shortcut
                
                dstPoint = 0
                anchorPoint = 0
                
                ##############
                for i in range(0, numpoints):
                    if range(i, numpoints) == []:
                        break

                    if (useResistanceCalcShortcut==True) and (dstPoint>0): 
                        break #No need to continue, we've got what we need to calculate resistances

                    dst = self.grid_to_graph (points_rc[i,1], points_rc[i,2], node_map)
                    local_dst = self.grid_to_graph (points_rc[i,1], points_rc[i,2], local_node_map)
                    if (dst >=  0 and components[dst] == c):
                        dstPoint = dstPoint+1
                        Gsolve = []
                    
                        if using_G_no_deleterow:
                            G_dst_dst = G[local_dst, local_dst] 
                            G[local_dst,local_dst] = 0
    
                        else:
                            Gsolve = self.gapdt.deleterowcol(G, delrow = local_dst, delcol = local_dst)
    
                        if using_G_no_deleterow:
                            self.state['amg_hierarchy'] = None
                            gc.collect()
                            self.create_amg_hierarchy(G)
                        else:
                            self.create_amg_hierarchy(Gsolve)

                        ################    
                        for j in range(i+1, numpoints):
                            if includedPairs[i+1,j+1]==1: #Test for pair in includedPairs
                                if self.state['amg_hierarchy']==None: #Called in case of memory error in current mapping
                                    if using_G_no_deleterow:
                                        self.create_amg_hierarchy(G)
                                    else:
                                        self.create_amg_hierarchy(Gsolve)                            
                               
                                if reportStatus==True:
                                    x = x+1
                                    (hours,mins,secs) = elapsed_time(self.state['startTime'])
                                    if useResistanceCalcShortcut==True:
                                        y = numpoints
                                        self.cs_log ('At ' + str(hours) +' hr ' + str(mins) + ' min solving focal node ' + str(x) + ' of '+ str(y) + '.',1)
                                    else:
                                        y = numpoints*(numpoints-1)/2
                                        self.cs_log ('At ' + str(hours) +' hr ' + str(mins) + ' min solving focal pair ' + str(x) + ' of '+ str(y) + '.',1)
                                src = self.grid_to_graph (points_rc[j,1], points_rc[j,2], node_map)
                                local_src = self.grid_to_graph (points_rc[j,1], points_rc[j,2], local_node_map)
                                
                                if (src >=  0 and components[src] == c):
                                    # Solve resistive network (Kirchoff's laws)
                                    solver_failed = False
                                                
                                    try:
                                        if using_G_no_deleterow:
                                            voltages = self.single_ground_solver(G, local_src, local_dst)
                                        else:
                                            voltages = self.single_ground_solver(Gsolve, local_src, local_dst)
                                    except:
                                        solver_failed = True
                                        solver_failed_somewhere = True
                                        resistances[i, j] = -777
                                        resistances[j, i] = -777
        
                                    if solver_failed == False:
                                        if self.options['low_memory_mode']==True or self.options['point_file_contains_polygons']==True:
                                            self.state['amg_hierarchy'] = None
                                            gc.collect()    
                                        
                                        resistances[i, j] = voltages[local_src] - voltages[local_dst]
                                        resistances[j, i] = voltages[local_src] - voltages[local_dst]
                                        # Write maps to files
                                        frompoint = str(points_rc[i,0])
                                        topoint = str(points_rc[j,0])
                                        
                                        if useResistanceCalcShortcut==True:
                                            if dstPoint==1: #this occurs for first i that is in component
                                                anchorPoint = i #for use later in shortcult resistance calc
                                                voltmatrix = self.getVoltmatrix(i,j,numpoints,local_node_map,voltages,points_rc,resistances,voltmatrix)                                          

                                        #BASELINE RESULTS FOR NETWORK CODE
                                        if self.options['write_baseline_results']==True:
                                            fileadd=frompoint + '_' + topoint +'_baseline'
                                            nodeNames=arange(1,node_map.max()+1)
                                            self.writeVoltagesNetwork(voltages, nodeNames, fileadd)
                                        ########################    


                                        if self.options['write_volt_maps'] == True:
                                            if reportStatus==True:
                                                (hours,mins,secs) = elapsed_time(self.state['startTime'])
                                                self.cs_log ('At ' + str(hours) +' hr ' + str(mins) + ' min writing voltage map ' + str(x) + ' of ' + str(y) + '.',1)
                                            voltage_map = self.create_voltage_map(local_node_map,voltages) 
                                            self.write_aaigrid('voltmap', '_' + frompoint + '_' + topoint, voltage_map)
                                            del voltage_map
                                        if self.options['write_cur_maps'] == True:
                                            if reportStatus==True:
                                                (hours,mins,secs) = elapsed_time(self.state['startTime'])
                                                self.cs_log ('At ' + str(hours) +' hr ' + str(mins) + ' min writing current map ' + str(x) + ' of ' + str(y) + '.',1)
                                            finitegrounds = [-9999] #create dummy value for pairwise case
                                            
                                            try:
                                                G = G.tocoo()#Can cause memory error
                                            except MemoryError:
                                                self.enable_low_memory(False)
                                                G = G.tocoo()
                                            try:
                                                current_map = self.create_current_map(voltages, G, local_node_map, finitegrounds)   
                                            except MemoryError:
                                                self.enable_low_memory(False)
                                                current_map = self.create_current_map(voltages, G, local_node_map, finitegrounds)                                         
                                            G = G.tocsr()


                                            #BASELINE RESULTS FOR NETWORK CODE            
                                            if self.options['write_baseline_results']==True:
                                                if i==0 and j==1:
                                                    cum_branch_currents = sparse.csr_matrix((G.shape))
                                                G =  G.tocoo()
                                                branch_currents = self.get_branch_currents(G,voltages,True) #Fixme: add up cumulative map too.
                                                branch_currents = absolute(branch_currents) 
                                                cum_branch_currents=cum_branch_currents+branch_currents
                                                G = G.tocsr()
                                                outputDir, outputFile = os.path.split(self.options['output_file'])
                                                outputBase, outputExtension = os.path.splitext(outputFile)
                                                file = outputDir + '//' + outputBase + '_branch_currents_' + frompoint + '_' + topoint +'_baseline.txt'
                                                nodeNames=arange(1,node_map.max()+1)
                                                self.writeGraph(file,branch_currents,nodeNames)
                                            #########################                 

                                            cum_current_map = cum_current_map + current_map 
                                            if self.options['write_max_cur_maps']==True:
                                                max_current_map = maximum(max_current_map, current_map) 
                                            if self.options['write_cum_cur_map_only']==False:
                                                if self.options['log_transform_maps']==True:
                                                    current_map = where(current_map>0,log10(current_map),self.state['nodata'])
                                                self.write_aaigrid('curmap', '_' + frompoint + '_' + topoint, current_map)
                                            del current_map    
                           


                                        if self.options['write_volt_drop_maps']==True:  #Fixme: need to add vdrop maps into one-to-all mode
                                            voltage_map = self.create_voltage_map(local_node_map,voltages) 
                                            vdrop_map=self.create_vdrop_map(voltage_map, node_map)
                                            self.write_aaigrid('volt_gradient_map', '_' + frompoint + '_' + topoint, vdrop_map)
                                            ind = node_map == 0
                                            vdrop_map[where(ind)] = 0                                            
                                            cum_vdrop_map = cum_vdrop_map + vdrop_map
                                            del vdrop_map
                                        (hours,mins,secs) = elapsed_time(lastWriteTime)
                                        if secs > 120: 
                                            lastWriteTime = time.time()
                                            self.saveIncompleteResistances(resistances)#save incomplete resistances
                        if (useResistanceCalcShortcut==True and i==anchorPoint): #this happens once per component. Anchorpoint is the first i in component
                            shortcutResistances = self.getShortcutResistances(anchorPoint,voltmatrix,numpoints,resistances,shortcutResistances)
                                                
                        if using_G_no_deleterow:#G here
                            G[local_dst, local_dst] = G_dst_dst

                    #End for
                    self.state['amg_hierarchy'] = None
                    gc.collect()

                #End if
                self.state['amg_hierarchy'] = None
                gc.collect()

        # Finally, resistance to self is 0.
        if useResistanceCalcShortcut==True: 
            resistances = shortcutResistances
        for i in range(0,numpoints):
            resistances[i, i] = 0


        if self.options['write_baseline_results']==True:
            file = outputDir + '//' + outputBase + '_branch_currents_cum_baseline.txt'
            nodeNames=arange(1,node_map.max()+1)
            self.writeGraph(file,cum_branch_currents,nodeNames)
        return resistances,cum_current_map,max_current_map,cum_vdrop_map,solver_failed_somewhere

    @print_timing
    def single_ground_solver(self, G, src, dst):#G here
        n = G.shape[0]
        rhs = zeros(n, dtype = 'float64')
        if using_G_no_deleterow:
            if src==dst:
                voltages = zeros(n, dtype = 'float64')
            else:
                rhs[dst] = -1
                rhs[src] = 1
                voltages = self.solve_linear_system (G, rhs)

        else:
            if src==dst:
                voltages = zeros(n+1, dtype = 'float64')
            else:
            
                if src <=  dst:
                    rhs[src] = 1
                else:
                    rhs[src-1] = 1
                
                x = self.solve_linear_system (G, rhs)
             
                voltages = zeros(n+1, dtype = 'float64')
                keep = delete (arange(0, n+1, dtype = 'int32'), dst)
                voltages[keep] = x

        return voltages

    @print_timing
    def advanced_module(self, g_map, poly_map, source_map, ground_map,source_id, G, node_map, component_map, componentWithPoints): 
        if node_map==None:
            oneToAllStreamline = False
        else:
            oneToAllStreamline = True
        solver_called = False
        solver_failed = False 

        if oneToAllStreamline==False:
            node_map = self.construct_node_map(g_map, poly_map)
            (component_map, components) = self.construct_component_map(g_map, node_map)
        if self.options['scenario']=='advanced':
            self.cs_log('Graph has ' + str(node_map.max()) + ' nodes and '+ str(components.max())+ ' components.',2)
            
        if self.options['write_cur_maps'] == True:
            cum_current_map = zeros((self.state['nrows'],self.state['ncols']),dtype = 'float64')         
        
        if self.options['write_volt_maps'] == True: 
            cum_voltage_map = zeros((self.state['nrows'],self.state['ncols']),dtype = 'float64') 
        elif self.options['scenario']=='one-to-all':
            cum_voltage_map = zeros((self.state['nrows'],self.state['ncols']),dtype = 'float64') 
        
        if oneToAllStreamline==False:        
            cmin = 1
            cmax = components.max()
        else:
            cmin = componentWithPoints
            cmax = componentWithPoints
            
        for c in range(cmin,int(cmax+1)):
            c_map = where(component_map == c, 1, 0)
            local_source_map = multiply(c_map, source_map)
            local_ground_map = where(c_map, ground_map, 0) 
            del c_map
            
            source_in_component = (where(local_source_map, 1, 0)).sum() > 0
            ground_in_component = (where(local_ground_map, 1, 0)).sum() > 0
            
            if (source_in_component) & (ground_in_component):
                (rows, cols) = where(local_source_map)
                values = local_source_map[rows,cols]
                local_sources_rc = c_[values,rows,cols]
                (rows, cols) = where(local_ground_map)
                values = local_ground_map[rows,cols]
                local_grounds_rc = c_[values,rows,cols]
                del rows, cols, values, local_source_map, local_ground_map 

                if oneToAllStreamline==False:
                    (G, node_map) = self.node_pruner(g_map, poly_map, component_map, c)#G here

                numnodes = node_map.max()
                sources = zeros(numnodes)
                grounds = zeros(numnodes)
                num_local_sources = local_sources_rc.shape[0]
                num_local_grounds = local_grounds_rc.shape[0]

                for source in range(0, num_local_sources):
                    src = self.grid_to_graph (local_sources_rc[source,1], local_sources_rc[source,2], node_map)
                    # Possible to have more than one source at a node when there are polygons
                    sources[src] = sources[src] + local_sources_rc[source,0] 

                for ground in range(0, num_local_grounds):
                    gnd = self.grid_to_graph (local_grounds_rc[ground,1], local_grounds_rc[ground,2], node_map)
                    # Possible to have more than one ground at a node when there are polygons
                    grounds[gnd] = grounds[gnd] + local_grounds_rc[ground,0] 

                (sources, grounds, finitegrounds) = self.resolve_conflicts(sources, grounds)

                solver_called = True
                try:
                    voltages = self.multiple_solver(G, sources, grounds, finitegrounds) #G here
                    del sources, grounds
                    if c==int(cmax):
                        del g_map, poly_map, ground_map
                        gc.collect()
                except MemoryError:
                    raise MemoryError
                except:
                    voltages = -777
                    solver_failed = True
                    
                if solver_failed==False:
                    ##Voltage and current mapping are cumulative, since there may be independent components.
                    if self.options['write_volt_maps'] == True: 
                        cum_voltage_map +=  self.create_voltage_map(node_map,voltages) 
                    elif self.options['scenario']=='one-to-all':
                        cum_voltage_map +=  self.create_voltage_map(node_map,voltages) 
                    if self.options['write_cur_maps'] == True:
                        G = G.tocoo()
                        cum_current_map +=  self.create_current_map(voltages, G, node_map, finitegrounds) #G here
                
        if self.options['write_volt_maps'] == True: 
            if solver_failed==False:
                filetext = 'voltmap'
                if self.options['scenario']=='advanced':
                    fileadd = ''
                else:
                    fileadd = '_'+str(source_id)
                self.write_aaigrid(filetext, fileadd, cum_voltage_map)  

        if self.options['scenario']=='advanced':
            if self.options['write_cur_maps'] == True: 
                if solver_failed==False:
                    if self.options['log_transform_maps']==True:
                        cum_current_map = where(cum_current_map>0,log10(cum_current_map),self.state['nodata']) 
                    filetext = 'curmap'
                    fileadd = ''
                    self.write_aaigrid(filetext, fileadd, cum_current_map) 
            else:
               cum_current_map = None

        else:
            if self.options['write_cur_maps'] == True: 
                if self.options['write_cum_cur_map_only']==False:
                    if solver_failed==False:
                        if self.options['log_transform_maps']==True:
                            cum_current_map = where(cum_current_map>0,log10(cum_current_map),self.state['nodata']) 
                        filetext = 'curmap'   
                        fileadd = '_'+str(source_id)   
                        self.write_aaigrid(filetext, fileadd, cum_current_map) 
            else:
                cum_current_map = None
            
        if self.options['scenario']=='one-to-all':
            if solver_failed==False:
                (row, col) = where(source_map>0)
                voltages = cum_voltage_map[row,col]/source_map[row,col] #allows for variable source strength
        elif self.options['scenario']=='all-to-one':
            if solver_failed==False:
                voltages = 0 #return 0 voltage/resistance for all-to-one mode           

        # Advanced mode will return voltages of the last component solved only for verification purposes.  
        if solver_called==False:
            voltages = -1

        return voltages,cum_current_map,solver_failed

    def resolve_conflicts(self, sources, grounds):

        finitegrounds = where(grounds<Inf,grounds,0)
        if (where(finitegrounds==0, 0, 1)).sum()==0:
            finitegrounds = [-9999]
        infgrounds = where(grounds==Inf,1,0)
        
        ##Resolve conflicts bewteen sources and grounds
        conflicts = logical_and(sources,grounds)
        if self.options['remove_src_or_gnd']=='rmvsrc':
            sources = where(conflicts,0,sources)
        elif self.options['remove_src_or_gnd']=='rmvgnd':
            grounds = where(conflicts,0,grounds)
        elif self.options['remove_src_or_gnd']=='rmvall':
            sources = where(conflicts,0,sources)
        infconflicts = logical_and(sources,infgrounds)
        grounds = where(infconflicts,0,grounds)
        if size(where(sources)) == 0:
            raise RuntimeError('All sources conflicted with grounds and were removed. There is nothing to solve.') 
        if size(where(grounds)) == 0:
            raise RuntimeError('All grounds conflicted with sources and were removed.  There is nothing to solve.') 

        return (sources, grounds, finitegrounds)

    @print_timing
    def multiple_solver(self, G, sources, grounds, finitegrounds):
        if finitegrounds[0]==-9999:#Fixme: no need to do this, right?
            finitegrounds = zeros(G.shape[0],dtype = 'int32') #create dummy vector for pairwise case
            Gsolve = G + sparse.spdiags(finitegrounds.T, 0, G.shape[0], G.shape[0]) 
            finitegrounds = [-9999]
        else:
            Gsolve = G + sparse.spdiags(finitegrounds.T, 0, G.shape[0], G.shape[0]) 
           
        ##remove infinite grounds from graph
        infgroundlist = where(grounds==Inf)
        infgroundlist = infgroundlist[0]
        numinfgrounds = infgroundlist.shape[0]
        
        dst_to_delete = []
        for ground in range(1, numinfgrounds+1):
            dst = infgroundlist[numinfgrounds-ground]
            dst_to_delete.append(dst)
            #Gsolve = self.gapdt.deleterowcol(Gsolve, delrow = dst, delcol = dst)
            keep = delete (arange(0, sources.shape[0]), dst)
            sources = sources[keep]            
        Gsolve = self.gapdt.deleterowcol(Gsolve, delrow = dst_to_delete, delcol = dst_to_delete)
        
        self.create_amg_hierarchy(Gsolve)
        voltages = self.solve_linear_system(Gsolve, sources)
        del Gsolve
        self.state['amg_hierarchy'] = None

        numinfgrounds = infgroundlist.shape[0]
        if numinfgrounds>0:
            #replace infinite grounds in voltage vector
            for ground in range(numinfgrounds,0, -1): 
                node = infgroundlist[numinfgrounds - ground] 
                voltages = asmatrix(insert(voltages,node,0)).T
        return asarray(voltages).reshape(voltages.size)
                    
    @print_timing
    def construct_node_map(self, g_map, poly_map):
        node_map = zeros(g_map.shape, dtype = 'int32')
        node_map[g_map.nonzero()] = arange(1, sum(g_map>0)+1, dtype = 'int32')

        if poly_map == []:
            return node_map

        # Remove extra points from poly_map that are not in g_map
        poly_map_pruned = zeros(g_map.shape, dtype = 'int32')
        poly_map_pruned[where(g_map)] = poly_map[where(g_map)]
        
        polynums = unique (poly_map)
   
        for i in range(0, polynums.size):
            polynum = polynums[i]
            if polynum !=  0:

                (pi, pj) = where (poly_map_pruned == polynum) #
                (pk, pl) = where (poly_map == polynum) #Added 040309 BHM                
                if len(pi)>0:  
                    node_map[pk, pl] = node_map[pi[0], pj[0]] #Modified 040309 BHM  
        node_map[where(node_map)] = self.gapdt.relabel(node_map[where(node_map)], 1) #BHM 072409

        return node_map

    @print_timing
    def construct_component_map(self, g_map, node_map):
        G = self.construct_g_graph(g_map, node_map) 
        C = self.gapdt.components(G) 

        (I, J) = where(node_map)
        nodes = node_map[I, J].flatten()

        component_map = zeros(node_map.shape, dtype = 'int32')
        component_map[I, J] = C[nodes-1]

        return (component_map, C)


    @print_timing
    def construct_g_graph(self, g_map, node_map):
        numnodes = node_map.max()
        (node1, node2, conductances) = self.get_conductances(g_map, node_map)
        G = sparse.csr_matrix((conductances, (node1, node2)), shape = (numnodes, numnodes))#G here, but not laplacian.  Memory hogging operation?
        if self.options['write_graph']==True:
            self.writeGraph('./graph_2x5.txt',G,None)
        g_graph = G + G.T
        return g_graph

    def writeGraph(self,filename,graph,nodeNames):
#         print 'nodenames size',nodeNames.shape
#         print 'graph shape',graph.shape
        Gcoo =  graph.tocoo()
        mask = Gcoo.data > 0
        
        #debug
#         test=Gcoo.row[mask]
#         test=Gcoo.col[mask]
        
        
#         graphNcol = zeros((mask.size,3),dtype = "float64")
        graphNcol = zeros((Gcoo.row[mask].size,3),dtype = "float64") #Fixme: this may result in zero-current nodes being left out.  Needed to make change so dimensions would match Gcoo.data[mask]
        
        if nodeNames==None:
#             print 'mask size:', mask.size
#             test= Gcoo.row[mask]
#             print 'Gcoo.row[mask] size:',test.shape
#             test= Gcoo.col[mask]
#             print 'Gcoo.col[mask] size:',test.shape            
#             
#             print 'graphNcol size:',graphNcol.shape
#             print 'mask:'
#             print mask
#             print 'Gcoo.row[mask]:'
#             print test
            graphNcol[:,0] = Gcoo.row[mask]
            graphNcol[:,1] = Gcoo.col[mask]
        else:
            graphNcol[:,0]=nodeNames[Gcoo.row[mask]]
            graphNcol[:,1]=nodeNames[Gcoo.col[mask]]
        graphNcol[:,2] = Gcoo.data[mask]
        savetxt(filename,graphNcol)
        return

    def get_horiz_neighbors(self, g_map):
        m = g_map.shape[0]
        n = g_map.shape[1]

        g_map_l = g_map[:, 0:(n-1)]
        g_map_r = g_map[:, 1:n]
        g_map_lr = double(logical_and(g_map_l, g_map_r))
        s_horiz = where(c_[g_map_lr, zeros((m,1), dtype = 'int32')].flatten())
        t_horiz = where(c_[zeros((m,1), dtype = 'int32'), g_map_lr].flatten())

        return (s_horiz, t_horiz)

    def get_vert_neighbors(self, g_map):
        m = g_map.shape[0]
        n = g_map.shape[1]

        g_map_u = g_map[0:(m-1), :]
        g_map_d = g_map[1:m    , :]
        g_map_ud = double(logical_and(g_map_u, g_map_d))
        s_vert = where(r_[g_map_ud, zeros((1,n), dtype = 'int32')].flatten())
        t_vert = where(r_[zeros((1,n), dtype = 'int32'), g_map_ud].flatten())
        
        return (s_vert, t_vert)

    def get_diag1_neighbors(self, g_map):
        m = g_map.shape[0]
        n = g_map.shape[1]

        z1 = zeros((m-1, 1), dtype = 'int32')
        z2 = zeros((1  , n), dtype = 'int32')
        
        g_map_ul  = g_map[0:m-1, 0:n-1]
        g_map_dr  = g_map[1:m  , 1:n  ]
        g_map_udr = double(logical_and(g_map_ul, g_map_dr)) 
        s_dr      = where(r_[c_[g_map_udr, z1], z2].flatten())
        t_dr      = where(r_[z2, c_[z1, g_map_udr]].flatten())
        
        return (s_dr, t_dr)

    def get_diag2_neighbors(self, g_map):
        m = g_map.shape[0]
        n = g_map.shape[1]

        z1 = zeros((m-1, 1), dtype = 'int32')
        z2 = zeros((1  , n), dtype = 'int32')

        g_map_ur  = g_map[0:m-1, 1:n  ]
        g_map_dl  = g_map[1:m  , 0:n-1]
        g_map_udl = double(logical_and(g_map_ur, g_map_dl)) 
        s_dl      = where(r_[c_[z1, g_map_udl], z2].flatten())
        t_dl      = where(r_[z2, c_[g_map_udl, z1]].flatten())
                        
        return (s_dl, t_dl)
        
    def get_conductances(self, g_map, node_map):
        (s_horiz, t_horiz) = self.get_horiz_neighbors(g_map)
        (s_vert,  t_vert)  = self.get_vert_neighbors(g_map)

        s = c_[s_horiz, s_vert].flatten()
        t = c_[t_horiz, t_vert].flatten()

        # Conductances
        g1 = g_map.flatten()[s]
        g2 = g_map.flatten()[t]

        if self.options['connect_using_avg_resistances'] == False:
            conductances = (g1+g2)/2
        else:
            conductances = 1 /((1/g1+1/g2)/2)

        if self.options['connect_four_neighbors_only'] == False:
            (s_dr, t_dr) = self.get_diag1_neighbors(g_map)
            (s_dl, t_dl) = self.get_diag2_neighbors(g_map)

            sd = c_[s_dr, s_dl].flatten()
            td = c_[t_dr, t_dl].flatten()

            # Conductances
            g1 = g_map.flatten()[sd]
            g2 = g_map.flatten()[td]

            if self.options['connect_using_avg_resistances'] == False:
                conductances_d = (g1+g2) / (2*sqrt(2))
            else:
                conductances_d =  1 / (sqrt(2)*(1/g1 + 1/g2) / 2)

            conductances = r_[conductances, conductances_d]

            s = r_[s, sd].flatten()
            t = r_[t, td].flatten()

        # Nodes in the g_graph. 
        # Subtract 1 for Python's 0-based indexing. Node numbers start from 1
        node1 = node_map.flatten()[s]-1
        node2 = node_map.flatten()[t]-1
        
        return (node1, node2, conductances)

    @print_timing
    def node_pruner(self, g_map, poly_map, component_map, keep_component):
        selector = component_map == keep_component
        
        g_map_pruned = selector * g_map
        poly_map_pruned = []
        if poly_map !=  []:
            poly_map_pruned = selector * poly_map

        node_map_pruned = self.construct_node_map (g_map_pruned, poly_map_pruned)
        G_pruned = self.construct_g_graph (g_map_pruned, node_map_pruned) #G here

        if pylab_available:
            pylab.spy(selector)
            pylab.draw()

        G = self.laplacian(G_pruned) #G here
        
        return (G, node_map_pruned)#G here

    @print_timing
    def laplacian(self, G): #G here
        n = G.shape[0]

        # FIXME: Potential for memory savings, if assignment is used
        G = G - sparse.spdiags(G.diagonal(), 0, n, n)
        G = -G + sparse.spdiags(G.sum(0), 0, n, n)

        return G

    @print_timing
    def create_amg_hierarchy(self, G): #G here
        if self.options['solver'] == 'amg' or self.options['solver'] == 'cg+amg':
            self.state['amg_hierarchy'] = None
            # construct the MG hierarchy
            ml = []
            if using_G_no_deleterow:
                ml = smoothed_aggregation_solver(G)
            else:
                ml = ruge_stuben_solver(G)    
            self.state['amg_hierarchy'] = ml
  
        return

    @print_timing
    def solve_linear_system(self, G, rhs): #G here
        gc.collect()
        # Solve G*x = rhs
        x = []
        if self.options['solver'] == 'cg+amg':
            ml = self.state['amg_hierarchy']
            G.psolve = ml.psolve
            (x, flag) = sparse.linalg.cg(G, rhs, tol = 1e-6, maxiter = 100000)
            if flag !=  0 or linalg.norm(G*x-rhs) > 1e-3:
                raise RuntimeError('CG did not converge. May need more iterations.') 

        if self.options['solver'] == 'umfpack':
            umfpack = um.UmfpackContext()
            x = umfpack( um.UMFPACK_A, G, rhs )
            
        if self.options['solver'] == 'amg':
            ml = self.state['amg_hierarchy']
            x = ml.solve(rhs, tol = 1e-6);
            
        if self.options['solver'] == 'superlu':
            x = sparse.linalg.spsolve(G, rhs)

        return x 

    @print_timing
    def create_voltage_map(self, node_map, voltages):
        voltage_map = numpy.zeros((self.state['nrows'], self.state['ncols']), dtype = 'float64')
        ind = node_map > 0
        voltage_map[where(ind)] = asarray(voltages[node_map[ind]-1]).flatten()
        if self.options['set_null_voltages_to_nodata']==True:
            ind = node_map == 0
            voltage_map[where(ind)] = self.state['nodata']
        return voltage_map


######################### BEGIN CURRENT MAPPING CODE ########################################
    @print_timing
    def create_current_map(self, voltages, G, node_map, finitegrounds):#G here
        gc.collect()
        node_currents = self.get_node_currents(voltages, G, finitegrounds)
        (rows, cols) = where(node_map)
        vals = node_map[rows, cols]-1
        current_map = zeros((self.state['nrows'],self.state['ncols']),dtype = 'float64')
        current_map[rows,cols] = node_currents[vals]

        return current_map

    def get_node_currents(self, voltages, G, finitegrounds):
        node_currents_pos = self.get_node_currents_posneg (G, voltages, finitegrounds, True) 
        node_currents_neg = self.get_node_currents_posneg (G, voltages, finitegrounds, False)
        node_currents = where(node_currents_neg > node_currents_pos, node_currents_neg, node_currents_pos)

        return asarray(node_currents)[0]

    def get_node_currents_posneg (self, G, voltages, finitegrounds, pos):
        branch_currents = self.get_branch_currents(G,voltages,pos)
        branch_currents = branch_currents-branch_currents.T #Can cause memory error
        
        branch_currents = branch_currents.tocoo() #Can cause memory error, but this and code below more memory efficient than previous version.
        mask = branch_currents.data > 0
        row  = branch_currents.row[mask]
        col  = branch_currents.col[mask]
        data = branch_currents.data[mask]
        del mask
        n = G.shape[0]
        branch_currents = sparse.csr_matrix((data, (row, col)), shape = (n,n))
           
        if finitegrounds[0]!= -9999:  
            finiteground_currents = multiply(finitegrounds, voltages)
            if pos == True:
                finiteground_currents = where(finiteground_currents < 0, -finiteground_currents, 0)
            else:
                finiteground_currents = where(finiteground_currents > 0, finiteground_currents, 0)  
            n = G.shape[0]
            branch_currents = branch_currents + sparse.spdiags(finiteground_currents.T, 0, n, n)        

        return branch_currents.sum(0)
    
    def get_branch_currents(self,G,voltages,pos):    
        branch_currents = self.get_branch_currents_posneg(G,voltages,pos)
        n = G.shape[0]
        mask = G.row < G.col
        branch_currents = sparse.csr_matrix((branch_currents, (G.row[mask], G.col[mask])), shape = (n,n)) #SQUARE MATRIX, SAME DIMENSIONS AS GRAPH
        return branch_currents

    def get_branch_currents_posneg(self,G,voltages,pos):
        mask = G.row < G.col
        if pos==True:
             vdiff = voltages[G.row[mask]]              
             vdiff -=  voltages[G.col[mask]]             

        else:
             vdiff = voltages[G.col[mask]]              
             vdiff -=  voltages[G.row[mask]]             

        conductances = where(G.data[mask] < 0, -G.data[mask], 0)
        del mask
        
        branch_currents = asarray(multiply(conductances,vdiff.T)).flatten()
        maxcur = max(branch_currents)
        branch_currents = where(absolute(branch_currents/maxcur) < 1e-8,0,branch_currents) #Delete very small branch currents to save memory
        return branch_currents
######################### END CURRENT MAPPING CODE ########################################        
        
        
    ### FILE I/O ###
    def write_aaigrid(self, type, fileadd, data):
        if type == 'voltmap':
            if self.options['write_volt_maps'] == False: 
                return
        elif type == 'curmap' or type == 'cum_curmap' or type == 'max_curmap':
            if self.options['write_cur_maps'] == False: 
                return
        elif type == 'volt_gradient_map' or type == 'cum_volt_gradient_map':
            if self.options['write_volt_drop_maps'] == False: 
                return
        else:
            return

        outputDir, outputFile = os.path.split(self.options['output_file'])
        outputBase, outputExtension = os.path.splitext(outputFile)
        file = outputDir + '//' + outputBase + '_' + type + fileadd +'.asc'

        writer(file, data, self.state, self.options['compress_grids'])
        
    def read_cell_map(self, filename):
        (ncols, nrows, xllcorner, yllcorner, cellsize, nodata) = read_header(filename)
        self.state['ncols'] = ncols
        self.state['nrows'] = nrows
        self.state['xllcorner'] = xllcorner
        self.state['yllcorner'] = yllcorner
        self.state['cellsize'] = cellsize
        if nodata==False:
            self.state['nodata'] = -9999
        else:
            self.state['nodata'] = nodata

        cell_map = reader(filename, 'float64')

        ######################## Reclassification code
        if self.options['use_reclass_table']==True:
            try:
                reclassTable = self.readPointStrengths(self.options['reclass_file'])    
            except:
                raise RuntimeError('Error reading reclass table')
            for i in range (0,reclassTable.shape[0]):
                cell_map = where(cell_map==reclassTable[i,0],reclassTable[i,1],cell_map)
            print'\n***** Reclassified habitat map using', self.options['reclass_file'],'*****'
        ########################
        
        if self.options['habitat_map_is_resistances'] == True:
            zeros_in_resistance_map = (where(cell_map==0, 1, 0)).sum() > 0
            if zeros_in_resistance_map == True: #FIXME: Should be easy to accomodate zeros in resistance map, just treat them like polygons.
                raise RuntimeError('Error: zero resistance values are not currently supported for habitat maps.  Use a short-circuit region file instead.')
            g_map = 1 / cell_map  
            g_map = where(cell_map == -9999,0,g_map)
        else:
            g_map = where(cell_map == -9999,0,cell_map)    
        g_map = where(g_map < 0,0,g_map)    
        return g_map


    def read_point_map(self, filename):
        if os.path.isfile(filename)==False:
            raise RuntimeError('File "'  + filename + '" does not exist')
        base, extension = os.path.splitext(filename)
        
        if extension == ".txt":
            try:
                points = loadtxt(filename)
            except ValueError:
                raise RuntimeError('File "'  + filename + '" is not in correct text list format. \n If it is an ASCII grid, please use .asc extension.')                

            
            points_rc = zeros(points.shape,dtype = 'int32')
            try:
                points_rc[:,0] = points[:,0]
                points_rc[:,1] = ceil((self.state['nrows']-(points[:,2]-self.state['yllcorner'])/self.state['cellsize']))-1
                points_rc[:,2] = ceil(((points[:,1]-self.state['xllcorner'])/self.state['cellsize']))-1
                i = argsort(points_rc[:,0])
                points_rc = points_rc[i]
            except IndexError:
                raise RuntimeError('Error extracting focal node locations. Please check file format.')                

        elif extension == ".asc":
            readingMask = False
            (ncols, nrows, xllcorner, yllcorner, cellsize, nodata) = read_header(filename)
            if cellsize!= self.state['cellsize']:
                print'\n********\nWarning: Focal node raster has different \ncell size than habitat raster. \nCircuitscape will try to crudely resample the focal node raster. \nWe recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.'
                point_map = self.resampleMap(filename,readingMask)
            elif ncols!= self.state['ncols']:
                print'\n********\nWarning: Focal node raster has different \nnumber of columns than habitat raster. \nCircuitscape will try to crudely resample the focal node raster. \nWe recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.'
                point_map = self.resampleMap(filename,readingMask)
            elif nrows!= self.state['nrows']:
                print'\n********\nWarning: Focal node raster has different \nnumber of rows than habitat raster. \nCircuitscape will try to crudely resample the focal node raster. \nWe recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.'
                point_map = self.resampleMap(filename,readingMask)
            elif xllcorner!= self.state['xllcorner']:
                print'\n********\nWarning: Focal node raster has different \nxllcorner than habitat raster. \nCircuitscape will try to crudely resample the focal node raster. \nWe recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.'
                point_map = self.resampleMap(filename,readingMask)
            elif yllcorner!= self.state['yllcorner']:
                print'\n********\nWarning: Focal node raster has different \nyllcorner than habitat raster. \nCircuitscape will try to crudely resample the focal node raster. \nWe recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.'
                point_map = self.resampleMap(filename,readingMask)
            else:
                point_map = reader(filename, 'int32')

            (rows, cols) = where(point_map > 0)

            values = zeros(rows.shape,dtype = 'int32') 
            for i in range(0, rows.size):
                values[i] = point_map[rows[i], cols[i]]
            points_rc = c_[values,rows,cols]
            try:            
                i = argsort(points_rc[:,0])
                points_rc = points_rc[i]
            except IndexError:
                raise RuntimeError('Error extracting focal node locations. Please check file format.')                
        else:
            raise RuntimeError('Focal node file must have a .txt or .asc extension')
        
        #Check to make sure points fall within cellmap
        if min(points_rc[:,1])<0 or min(points_rc[:,2])<0:
            raise RuntimeError('At least one focal node location falls outside of habitat map')
        elif max(points_rc[:,1])>self.state['nrows']-1:
            raise RuntimeError('At least one focal node location falls outside of habitat map')
        elif  max(points_rc[:,2])>self.state['ncols']-1:
            raise RuntimeError('At least one focal node location falls outside of habitat map')
        if (unique(asarray(points_rc[:,0]))).shape[0]<2:
            raise RuntimeError('Less than two valid focal nodes found. Please check focal node location file.')                    
        
        return points_rc
        
    def read_poly_map(self, filename,readingMask):  
        (ncols, nrows, xllcorner, yllcorner, cellsize, nodata) = read_header(filename)
        if cellsize!= self.state['cellsize']:
            print'\n********\nWarning: Short-circuit region or mask raster has different \ncell size than habitat raster. \nCircuitscape will try to crudely resample the raster. \nWe recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.'
            map = self.resampleMap(filename,readingMask)
        elif ncols!= self.state['ncols']:
            print'\n********\nWarning: Short-circuit region or mask raster has different \nnumber of columns than habitat raster. \nCircuitscape will try to crudely resample the raster. \nWe recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.'
            map = self.resampleMap(filename,readingMask)
        elif nrows!= self.state['nrows']:
            print'\n********\nWarning: Short-circuit region or mask raster has different \nnumber of rows than habitat raster. \nCircuitscape will try to crudely resample the raster. \nWe recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.'
            map = self.resampleMap(filename,readingMask)
        elif xllcorner!= self.state['xllcorner']:
            print'\n********\nWarning: Short-circuit region or mask raster has different \nxllcorner than habitat raster. \nCircuitscape will try to crudely resample the raster. \nWe recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.'
            map = self.resampleMap(filename,readingMask)
        elif yllcorner!= self.state['yllcorner']:
            print'\n********\nWarning: Short-circuit region or mask raster has different \nyllcorner than habitat raster. \nCircuitscape will try to crudely resample the raster. \nWe recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.'
            map = self.resampleMap(filename,readingMask)            
        else:
            map = reader(filename, 'int32')
            
        map = where(map == nodata,0,map)        

        if readingMask==True:
            map = where(map < 0, 0, map)        

        return map

    def read_source_and_ground_maps(self, source_filename, ground_filename): 
        #FIXME: reader does not currently handle infinite inputs for ground conductances.
        if os.path.isfile(source_filename)==False:
            raise RuntimeError('File "'  + source_filename + '" does not exist')
        base, extension = os.path.splitext(source_filename)
        if extension == ".txt":  #FIXME: probably want to roll code used for reading source, ground and point text files into single utility
            try:
                sources = loadtxt(source_filename)
            except ValueError:
                raise RuntimeError('File "'  + source_filename + '" is not in correct text list format. \n If it is an ASCII grid, please use .asc extension.')                

            sources_rc = zeros(sources.shape,dtype = 'int32')
            sources_rc[:,0] = sources[:,0]
            sources_rc[:,1] = ceil((self.state['nrows']-(sources[:,2]-self.state['yllcorner'])/self.state['cellsize']))-1
            sources_rc[:,2] = ceil(((sources[:,1]-self.state['xllcorner'])/self.state['cellsize']))-1
            source_map = zeros((self.state['nrows'],self.state['ncols']),dtype = 'float64')
            source_map[sources_rc[:,1],sources_rc[:,2]] = sources_rc[:,0]
        elif extension=='.asc':
            (ncols, nrows, xllcorner, yllcorner, cellsize, nodata) = read_header(source_filename)
            if cellsize!= self.state['cellsize']:
                raise RuntimeError('Current source raster must have same cell size and number of rows and columns as habitat raster') 
            if ncols!= self.state['ncols']:
                raise RuntimeError('Current source raster must have same cell size and number of rows and columns as habitat raster') 
            if nrows!= self.state['nrows']:
                raise RuntimeError('Current source raster must have same cell size and number of rows and columns as habitat raster') 
            if yllcorner!= self.state['yllcorner']:
                raise RuntimeError('Current source raster must have same xllcorner and yllcorner as habitat raster') 
            if xllcorner!= self.state['xllcorner']:
                raise RuntimeError('Current source raster must have same xllcorner and yllcorner as habitat raster') 
                
                
                
            source_map = reader(source_filename, 'float64')
            source_map = where(source_map == -9999,0,source_map)

        else:
            raise RuntimeError('Current source files must have a .txt or .asc extension')
        if self.options['use_unit_currents']==True:
            source_map = where(source_map,1,0)

        if os.path.isfile(ground_filename)==False:
            raise RuntimeError('File "'  + ground_filename + '" does not exist')
        base, extension = os.path.splitext(ground_filename)
        if extension == ".txt":
            try:
                grounds = loadtxt(ground_filename)
            except ValueError:
                raise RuntimeError('File "'  + ground_filename + '" is not in correct text list format. \n If it is an ASCII grid, please use .asc extension.')                

            grounds_rc = zeros(grounds.shape,dtype = 'int32')
            grounds_rc[:,0] = grounds[:,0]
            grounds_rc[:,1] = ceil((self.state['nrows']-(grounds[:,2]-self.state['yllcorner'])/self.state['cellsize']))-1
            grounds_rc[:,2] = ceil(((grounds[:,1]-self.state['xllcorner'])/self.state['cellsize']))-1
            ground_map_raw = -9999*ones((self.state['nrows'],self.state['ncols']),dtype = 'float64')
            ground_map_raw[grounds_rc[:,1],grounds_rc[:,2]] = grounds_rc[:,0]
        elif extension=='.asc':
            (ncols, nrows, xllcorner, yllcorner, cellsize, nodata) = read_header(ground_filename)
            if cellsize!= self.state['cellsize']:
                raise RuntimeError('Ground raster must have same cell size and number of rows and columns as habitat raster') 
            if ncols!= self.state['ncols']:
                raise RuntimeError('Ground raster must have same cell size and number of rows and columns as habitat raster') 
            if nrows!= self.state['nrows']:
                raise RuntimeError('Ground raster must have same cell size and number of rows and columns as habitat raster') 
            if yllcorner!= self.state['yllcorner']:
                raise RuntimeError('Ground raster must have same xllcorner and yllcorner as habitat raster') 
            if xllcorner!= self.state['xllcorner']:
                raise RuntimeError('Ground raster must have same xllcorner and yllcorner as habitat raster') 
                
            ground_map_raw = reader(ground_filename, 'float64')
        else:
            raise RuntimeError('Ground files must have a .txt or .asc extension')
        if self.options['ground_file_is_resistances']==True:
            ground_map = 1 / ground_map_raw
            ground_map = where(ground_map_raw == -9999,0,ground_map)
        else:
            ground_map = where(ground_map_raw == -9999,0,ground_map_raw)
        if self.options['use_direct_grounds']==True:
            ground_map = where(ground_map,Inf,0)

        conflicts = logical_and(source_map,ground_map)
        if self.options['remove_src_or_gnd']=='rmvsrc':
            source_map = where(conflicts,0,source_map)
        elif self.options['remove_src_or_gnd']=='rmvgnd':
            ground_map = where(conflicts,0,ground_map)
        elif self.options['remove_src_or_gnd']=='rmvall':
            source_map = where(conflicts,0,source_map)
            ground_map = where(conflicts,0,ground_map)
        if size(where(source_map)) == 0:
            raise RuntimeError('No valid sources detected. Please check source file') 
        if size(where(ground_map)) == 0:
            raise RuntimeError('No valid grounds detected. Please check ground file') 
        return source_map, ground_map


    def readincludedPairs(self, filename):
        if os.path.isfile(filename)==False:
            raise RuntimeError('File "'  + filename + '" does not exist')
        
        try:
            f = open(filename, 'r')
            [ign, minval] = string.split(f.readline())
            [ign, maxval] = string.split(f.readline())
            minval = float(minval)
            maxval = float(maxval)
            f.close()            
            
            includedPairs = loadtxt(filename, skiprows = 2, dtype = 'Float64')
            pointIds = includedPairs[:,0]
            includedPairs = where(includedPairs>maxval,0,includedPairs)
            includedPairs = where(includedPairs<minval,0,1)             
            includedPairs[:,0] = pointIds
            includedPairs[0,:] = pointIds
            includedPairs[0,0] = -1
            i = argsort(includedPairs[:,0])
            includedPairs = includedPairs[i]
            includedPairs = includedPairs.T            
            i = argsort(includedPairs[:,0])
            includedPairs = includedPairs[i] 
        
        except:
            raise RuntimeError('Error reading focal node include/exclude matrix. Please check file format.')                

        return includedPairs

    def readPointStrengths(self, filename):
        if os.path.isfile(filename)==False:
            raise RuntimeError('File "'  + filename + '" does not exist')
        
        try:
            pointStrengths = loadtxt(filename)
        except ValueError:
            raise RuntimeError('Error reading focal node source strength list. Please check file format.')                
           
        try:
            pointIds = pointStrengths[:,0]
            i = argsort(pointStrengths[:,0])
            pointStrengths = pointStrengths[i]

        except:
            raise RuntimeError('Error reading focal node source strength list. Please check file format.')                
        
        return pointStrengths


    def enable_low_memory(self, restart):
        self.state['amg_hierarchy'] = None
        gc.collect()
        if self.options['low_memory_mode']==True:
            if restart==False: #If this module has already been called
                raise MemoryError
        self.options['low_memory_mode'] = True
        print'\n**************\nMemory error reported.'        

        type, value, tb = sys.exc_info()
        info = traceback.extract_tb(tb)
        print'Full traceback:'
        print info
        print'***************'
        filename, lineno, function, text = info[-1] # last line only
        print"\n %s:%d: %s: %s (in %s)" %\
              (filename, lineno, type.__name__, str(value), function)

        type = value = tb = None # clean up
        print'\nWARNING: CIRCUITSCAPE IS RUNNING LOW ON MEMORY.'
        if restart==True:
            print'Restarting in low memory mode, which will take somewhat longer to complete.'
        else:
            print'Switching to low memory mode, which will take somewhat longer to complete.'            
        print'CLOSING OTHER PROGRAMS CAN HELP FREE MEMORY RESOURCES.'
        print'Please see the user guide for more information on memory requirements.\n'               
        if restart==True:
            print'***Restarting in low memory mode***\n'
        else:
            print'***Continuing in low memory mode***\n'
        return


    def resampleMap(self,filename,readingMask):
        try:
            (ncols, nrows, xllcorner, yllcorner, cellsize, nodata) = read_header(filename)
            map = reader(filename, 'int32')
            map = where(map == nodata,0,map)
            
            if readingMask==True:
                map = where(map>0, 1, 0)      
                map = 1-map #now zeros are areas to keep, ones are masked out
                
            (rows, cols) = where(map > 0)
            
            xcoords = (cols+0.5)*cellsize+xllcorner
            ycoords = (nrows-rows-0.5)*cellsize+yllcorner
            
            values = zeros(rows.shape,dtype = 'int32') 
            for i in range(0, rows.size):
                values[i] = map[rows[i], cols[i]]
            mapCoords = c_[values,xcoords,ycoords]
    
            i = argsort(mapCoords[:,0])
            mapCoords = mapCoords[i]
            
            #From mapCoords to mapRc:
            mapRc = zeros(mapCoords.shape,dtype = 'int32')
            mapRc[:,0] = mapCoords[:,0]
            mapRc[:,1] = ceil((self.state['nrows']-(mapCoords[:,2]-self.state['yllcorner'])/self.state['cellsize']))-1
            mapRc[:,2] = ceil(((mapCoords[:,1]-self.state['xllcorner'])/self.state['cellsize']))-1
            i = argsort(mapRc[:,0])
            mapRc = mapRc[i]
    
            rows = mapRc[:,1]
            (delrows) = asarray(where(rows<0))
            delrows2 = zeros(delrows.shape[1]) #turn into 1-d array
            delrows2[:] = delrows[:]
            if delrows2!= []:
                mapRc = self.deleterow(mapRc,delrows2)
            rows = mapRc[:,1] 
            (delrows) = asarray(where(rows>self.state['nrows']-1))
            delrows2 = zeros(delrows.shape[1]) #turn into 1-d array
            delrows2[:] = delrows[:]
            if delrows2!= []:
                mapRc = self.deleterow(mapRc,delrows2)
            cols = mapRc[:,2]
            (delrows) = asarray(where(cols<0))
            delrows2 = zeros(delrows.shape[1]) #turn into 1-d array
            delrows2[:] = delrows[:]
            if delrows2!= []:
                mapRc = self.deleterow(mapRc,delrows2)
            cols = mapRc[:,2]
            (delrows) = asarray(where(cols>self.state['ncols']-1))
            delrows2 = zeros(delrows.shape[1]) #turn into 1-d array
            delrows2[:] = delrows[:]
            if delrows2!= []:
                mapRc = self.deleterow(mapRc,delrows2)
            del delrows
            del delrows2
    
            #From mapRc to map:                
            map = numpy.zeros((self.state['nrows'],self.state['ncols']),int)
            map[mapRc[:,1],mapRc[:,2]] = mapRc[:,0] 
            if readingMask==True:
                map = 1-map #now zeros are areas to mask out, ones are kept
            
        except: 
            raise RuntimeError('Error resampling focal node, mask, or short-circuit region locations to match habitat map cell size and extent.  We recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.')   

        return map

    def getVoltmatrix(self,i,j,numpoints,local_node_map,voltages,points_rc,resistances,voltmatrix):                                            
       self.options['set_null_voltages_to_nodata']=False #Not mapping if using resistance shortcut, which calls this module
       voltvector = zeros((numpoints,1),dtype = 'float64')  
       voltage_map = self.create_voltage_map(local_node_map,voltages) 
       for point in range(1,numpoints):
           voltageAtPoint = voltage_map[points_rc[point,1], points_rc[point,2]]
           voltageAtPoint = 1-(voltageAtPoint/resistances[i, j])
           voltvector[point] = voltageAtPoint
       voltmatrix[:,j] = voltvector[:,0] 
       return voltmatrix


    def getShortcutResistances(self,anchorPoint,voltmatrix,numpoints,resistances,shortcutResistances): #FIXME: no solver failed capability
        point1 = anchorPoint
        for pointx in range(0, numpoints): #point1 is source node, i.e. the 1 in R12.  point 2 is the dst node. TEMP was range(1, numpoints-1)
            R1x = resistances[point1,pointx]
            if R1x!= -1:
                shortcutResistances[point1,pointx] = R1x
                shortcutResistances[pointx,point1] = R1x
                for point2 in range(pointx,numpoints):
                    R12 = resistances[point1,point2] 
                    if R12!= -1:
                        shortcutResistances[point2,point1] = R12
                        shortcutResistances[point1,point2] = R12
                        Vx = voltmatrix[pointx,point2]
                        R2x = 2*R12*Vx+R1x-R12
                        shortcutResistances[pointx,point2] = R2x
                        shortcutResistances[point2,pointx] = R2x   
        return shortcutResistances                        


    def deleterow(self, A, delrow):
        m = A.shape[0]
        n = A.shape[1]
        keeprows = delete (arange(0, m), delrow)
        keepcols = arange(0, n)
        return A[keeprows][:,keepcols]

    def deletecol(self, A, delcol):
        m = A.shape[0]
        n = A.shape[1]
        keeprows = arange(0, m)
        keepcols = delete (arange(0, n), delcol)
        return A[keeprows][:,keepcols]

    def writeResistances(self, point_ids, resistances):
        focal_labels = insert(point_ids, [0], 0, axis = 0)
        resistances = insert(resistances, [0], 0, axis = 0)
        resistances = insert(resistances, [0], 0, axis = 1)
        resistances[0,:] = focal_labels
        resistances[:,0] = focal_labels
                
        fileName = self.options['output_file']
        outputDir, outputFile = os.path.split(fileName)
        outputBase, outputExtension = os.path.splitext(outputFile)
        outputFile = outputDir + '//' + outputBase + '_resistances' + outputExtension 
        if self.options['write_baseline_results']==True:
            outputFile = outputDir + '//' + outputBase + '_resistances_baseline.txt' 
        savetxt (outputFile, resistances)

        numPoints = resistances.shape[0]-1
        numEntries = numPoints*(numPoints-1)/2
        resistances3columns = zeros((numEntries,3),dtype = 'float64') 
        x = 0
        for i in range(1,numPoints):
            for j in range(i+1,numPoints+1):
                resistances3columns[x,0] = resistances[i,0]    
                resistances3columns[x,1] = resistances[0,j]
                resistances3columns[x,2] = resistances[i,j]
                x = x+1
        outputFile = outputDir + '//' + outputBase + '_resistances_3columns' + outputExtension 
        if self.options['write_baseline_results']==True:
            outputFile = outputDir + '//' + outputBase + '_resistances_3columns_baseline.txt'
        savetxt (outputFile, resistances3columns)
        
        #remove partial result file        
        oldFile = outputDir + '//' + outputBase + '_resistances_incomplete' + outputExtension
        try:
            os.remove(oldFile)
        except:
            pass 
        return resistances
        
    def saveIncompleteResistances(self, resistances):
        fileName = self.options['output_file']
        outputDir, outputFile = os.path.split(fileName)
        outputBase, outputExtension = os.path.splitext(outputFile)
        outputFile = outputDir + '//' + outputBase + '_resistances_incomplete' + outputExtension
        savetxt (outputFile, resistances)
        return


    def pruneIncludedPairsNetwork(self,focalNodes):
        includedPairs = (self.state['includedPairs'])
        includeList = list(includedPairs[0,:])
        point = 0
        dropFlag = False
        while point < focalNodes.size: #Prune out any points not in includeList
            if focalNodes[point] in includeList: #match
                point = point+1
            else:
                dropFlag = True   
                focalNodes=delete(focalNodes,point)
         
        includeList = list(focalNodes[:])
        numConnectionRows = includedPairs.shape[0]
        row = 1
        while row <numConnectionRows: #Prune out any entries in includeList that are not in focalNodes
            if includedPairs [row,0] in includeList: #match
                row = row+1
            else:
                includedPairs = self.gapdt.deleterowcol(includedPairs,delrow = row,delcol = row)   
                dropFlag = True
                numConnectionRows = numConnectionRows-1

        self.state['includedPairs'] = includedPairs                     
#         if dropFlag==True:
#             print'\nNOTE: Code to exclude pairwise calculations is activated and \nsome entries did not match with focal node file.  \nSome focal nodes may have been dropped.'      
        return focalNodes


    def pruneIncludedPairs(self,points_rc):
        includedPairs = (self.state['includedPairs'])
        includeList = list(includedPairs[0,:])
        point = 0
        dropFlag = False
        while point < points_rc.shape[0]: #Prune out any points not in includeList
            if points_rc[point,0] in includeList: #match
                point = point+1
            else:
                dropFlag = True   
                points_rc = self.deleterow(points_rc,point)  
         
        includeList = list(points_rc[:,0])
        numConnectionRows = includedPairs.shape[0]
        row = 1
        while row <numConnectionRows: #Prune out any entries in includeList that are not in points_rc
            if includedPairs [row,0] in includeList: #match
                row = row+1
            else:
                includedPairs = self.gapdt.deleterowcol(includedPairs,delrow = row,delcol = row)   
                dropFlag = True
                numConnectionRows = numConnectionRows-1

        self.state['includedPairs'] = includedPairs                     
#         if dropFlag==True:
#             print'\nNOTE: Code to exclude pairwise calculations is activated and \nsome entries did not match with focal node file.  \nSome focal nodes may have been dropped.'      
        return points_rc
    
    
    def get_points_rc_unique(self,point_ids,points_rc):
        points_rc_unique = zeros((point_ids.size,3), int)
        for i in range(0, point_ids.size):
            for j in range(0, points_rc.shape[0]):
                if points_rc[j,0]==point_ids[i]:
                    points_rc_unique[i,:] = points_rc[j,:] 
                    break                    
        return points_rc_unique          
        
    def checkPointsInComponent(self,c,numpoints,components,points_rc,node_map):
        points_in_this_component = False            
        for pt1 in range(0, numpoints): 
            if points_in_this_component == False:
                src = self.grid_to_graph (points_rc[pt1,1], points_rc[pt1,2], node_map)
                for pt2 in range(pt1+1, numpoints):
                    dst = self.grid_to_graph (points_rc[pt2,1], points_rc[pt2,2], node_map)
                    if (src >=  0 and components[src] == c) and (dst >=  0 and components[dst] == c):
                        points_in_this_component = True
                        break        
        return points_in_this_component    

    def getstrengthMap(self,points_rc_unique,pointStrengths):
        if self.options['use_variable_source_strengths']==True:
            if self.options['scenario'] == 'one-to-all': 
                strengths_rc = self.get_strengths_rc(self.state['pointStrengths'],points_rc_unique)
                strengthMap = None
            else:
                strengths_rc = self.get_strengths_rc(pointStrengths,points_rc_unique)
                strengthMap = numpy.zeros((self.state['nrows'],self.state['ncols']),dtype = 'Float64')
                strengthMap[points_rc_unique[:,1],points_rc_unique[:,2]] = strengths_rc[:,0]     
            return strengthMap,strengths_rc
        else:
            return None,None
        
    def get_strengths_rc(self,pointStrengths,points_rc_unique):
        strengths_rc = zeros(points_rc_unique.shape,dtype = 'float64')
        strengths_rc[:,1] = points_rc_unique[:,1]
        strengths_rc[:,2] = points_rc_unique[:,2]
        for point in range(0,points_rc_unique.shape[0]):
            try:
                strengthIds = list(pointStrengths[:,0])
                strengthValues = list(pointStrengths[:,1])
                pointId = points_rc_unique[point,0]
                indx = strengthIds.index(pointId)
                strengths_rc[point,0] = strengthValues[indx]
            except ValueError:
                strengths_rc[point,0] = 1
        return strengths_rc        

    @print_timing
    def load_maps(self):
        self.cs_log('Reading maps',1)
        self.cs_log('',2)
        self.state['g_map'] = self.read_cell_map(self.options['habitat_file'])
        if self.options['use_polygons']:
            self.state['poly_map'] = self.read_poly_map(self.options['polygon_file'],readingMask = False)
        else:
            self.state['poly_map'] = []
 
        if self.options['use_mask']==True:
            mask = self.read_poly_map(self.options['mask_file'],readingMask = True)
            mask = where(mask !=  0, 1, 0) 
            self.state['g_map'] = multiply(self.state['g_map'],mask)
            
            sumGmap = (self.state['g_map']).sum()
            sumGmap = sumGmap.sum()
            if sumGmap==0:
                raise RuntimeError('All entries in habitat map have been dropped after masking with the mask file.  There is nothing to solve.')             
            del mask
        else:
            self.state['mask'] = []
            

        if self.options['scenario']=='advanced':
            self.state['points_rc'] = []
            (self.state['source_map'], self.state['ground_map']) = self.read_source_and_ground_maps(self.options['source_file'], self.options['ground_file'])

        else:        
            self.state['points_rc'] = self.read_point_map(self.options['point_file'])
            self.state['source_map'] = []
            self.state['ground_map'] = []


        if self.options['use_included_pairs']==True:
            self.state['includedPairs'] = self.readincludedPairs(self.options['included_pairs_file'])
        
        self.state['pointStrengths'] = None
        if self.options['use_variable_source_strengths']==True:
            self.state['pointStrengths'] = self.readPointStrengths(self.options['variable_source_file']) 
        
        if pylab_available:
            pylab.spy(self.state['g_map'], hold = True, cmap = matplotlib.colors.ListedColormap(['c','r']))
            pylab.draw()

            pylab.spy(self.state['poly_map'], hold = True,  cmap = matplotlib.colors.ListedColormap(['c','g']))
            pylab.draw()

        self.cs_log('Processing maps',1)
        return 
        
        
#################### BEGIN STRESS TESTING CODE ############################################
    def run_stress_test(self,stress_ncols,stress_nrows):
        ##### 
        pointfile_contains_polys = False
        #############################
        print'\n****************STRESS TEST CODE ENABLED *************************'
        self.state['ncols'] = stress_ncols
        self.state['nrows'] = stress_nrows
        print'nrows = ',self.state['ncols']
        print'ncols = ',self.state['nrows'],'\n'
        self.state['xllcorner'] = 0
        self.state['yllcorner'] = 0
        self.state['cellsize'] = 1
        self.state['nodata'] = -9999            
        self.state['g_map'] = ones((self.state['nrows'], self.state['ncols']), dtype = 'float64')
        g_map = self.state['g_map'] #TEMP
#         g_map[500,:] = 0.002
#         g_map[1001,:] = 0.002        
        self.state['g_map'] = g_map
        self.state['poly_map'] = []
        self.options['write_cur_maps'] = False
#         self.state['points_rc'] =  array([[1,2,2],[2,40,40],[3,60,60],[4,80,80],[5,100,100],[6,120,120],[7,140,140],[8,160,160],[9,180,180],[10,200,200],[11,1600,400],[12,1640,440],[13,1660,460],[14,1680,480],[15,1700,500],[16,1900,500],[17,1900,500],[18,1900,700],[19,1900,900],[20,980,4],[21,1840,900]])        
        self.state['points_rc'] =  array([[1,1,1],[2,20,20],[3,30,30],[4,40,40],[5,50,50],[6,60,60],[7,70,70],[8,80,80],[9,90,90],[10,100,100],[11,800,200],[12,820,220],[13,830,230],[14,840,240],[15,850,250],[16,950,250],[17,950,250],[18,950,350],[19,950,450],[20,490,2],[21,920,450]]) 
        if pointfile_contains_polys==True:            
            self.state['points_rc'] =  array([[1,1,1],[1,10,110],[2,100,100],[3,200,200]])             
            self.options['point_file_contains_polygons'] = True
            print'Point file contains polygons.'
        if self.options['scenario']=='advanced':    
            self.state['source_map'] = zeros((self.state['nrows'], self.state['ncols']), dtype = 'float64')
            self.state['ground_map'] = zeros((self.state['nrows'], self.state['ncols']), dtype = 'float64') 
            self.state['source_map'][10,10] = 1
            self.state['source_map'][100,100] = 2
            self.state['ground_map'][20,20] = 1
            self.state['ground_map'][200,200] = 2 
        return 
#################### END STRESS TESTING CODE ############################################

#    raw_input('Hit any key to continue')    #debug code
# Fixme: place -2 in resistances for included pairs

        # 
# #        	Graph.write_adjacency('filename')
#       	Graph.Read_adjacency('filename')
#         test =   Graph.Read_ncol('graph_ncol.txt')
# #     def writeGraph(self,G):
#         numnodes = G.shape[0]
#         graphList = []
#         for i in range(0,numnodes-1)
#             for j in range(i+1,numnodes)
#                 if G[i,j]>0:
#                  resistances = insert(resistances, [0], 0, axis = 0)
#                     x = x+1
#                     graphList[x,1] = i+1
#                     graphList[x,2] = j+1
#                 
#         node1 = graphList[:,0]
#         node2 = graphList[:,1]
#         data = graphList[:,2]
# 
#         savetxt ('graph_written.txt', graphList)



# todo: warning if header of mask file is different specifies mask file
# advanced- test with source ground files
#tests for all applicable cases?

# for now, assume one component
# include/exclude NO
# strengths YES
# reclass? Yes
# write currents and voltages
# write version#

#VERIFY
# single
# single, exclude
# one
# all, exclude, varsrc
# adv
# adv, rmvsrc
# adv rmvgrnd
# adv rmvall
#reclass for regular verify

#test shortcut resistances
#auto detect ascii grids

#test header check code

#drop current in polygons?
#get unique node numbers, (node_map == value).sum()
#create a reclass table to create a node map that has number of times a node is in the map.  Divide current by that. Zero out voltage dropes when n>1.
#percentile maps?
#voltage drop


#VOLTAGE DROP- right now it's the literal voltage drop.  Should amend to max i*r, and only use when average resistances?  That way drop happens in cell responsible for it, not spread across barrier and adhacent cells.

#NODATA for currents when input hab map=nodata NO NEED?  0 is valid
#NODATA for voltages when g=0  DO NEED.  Voltage is undefined

# incomplete resistances can slow down!
#add incomplete r's to pairwise network

    def create_vdrop_map(self, volt_map, node_map):
#         print 'nodes'
#         print node_map
        #FIXME: Not sure yet how polymap affects this... 
        if self.options['normalize_vdrop_maps']==True:
            volt_map=volt_map / min(volt_map)
        up=logical_and(node_map,self.shift_down(node_map))
        vdiff_up = where(up,volt_map-self.shift_down(volt_map),0)
#         r_ratio_up=where(up,((1/g_map)/((1/g_map)+(1/(self.shift_down(g_map))))),0) FIXME: need something like this, but problematic with polygons?  Maybe not....
        vdiff_up[0,:]=0

        left=logical_and(node_map,self.shift_right(node_map))
        vdiff_left = where(left,volt_map-self.shift_right(volt_map),0)
        vdiff_left[:,0]=0

        right=logical_and(node_map,self.shift_left(node_map))               
        vdiff_right = where(right,volt_map-self.shift_left(volt_map),0)
        vdiff_right[:,self.state['ncols']-1]=0

        down=logical_and(node_map,self.shift_up(node_map))
        vdiff_down = where(down,volt_map-self.shift_up(volt_map),0)
        vdiff_down[self.state['nrows']-1,:]=0
        
        if self.options['connect_four_neighbors_only'] == False:
            ul=logical_and(node_map,self.shift_right(self.shift_down(node_map)))
            vdiff_ul = where(ul,volt_map-self.shift_right(self.shift_down(volt_map)),0) #FIXME: doesn't consider that shifted cell may not be a node (or does it?)
            vdiff_ul[0,:]=0
            vdiff_ul[:,0]=0
    
            ur=logical_and(node_map,self.shift_left(self.shift_down(node_map)))
            vdiff_ur = where(ur,volt_map-self.shift_left(self.shift_down(volt_map)),0)
            vdiff_ur[0,:]=0
            vdiff_ur[:,self.state['ncols']-1]=0
    
            dl=logical_and(node_map,self.shift_right(self.shift_up(node_map)))
            vdiff_dl = where(dl,volt_map - self.shift_right(self.shift_up(volt_map)),0)
            vdiff_dl[self.state['nrows']-1,:]=0
            vdiff_dl[:,0]=0
    
            dr=logical_and(node_map,self.shift_left(self.shift_up(node_map)))
            vdiff_dr = where(dr,volt_map-self.shift_left(self.shift_up(volt_map)),0)
            vdiff_dr[self.state['nrows']-1,:]=0
            vdiff_dr[:,self.state['ncols']-1]=0

        else:
            vdiff_ul=zeros(volt_map.shape,dtype='int32')
            vdiff_ur=vdiff_ul
            vdiff_dl=vdiff_ul
            vdiff_dr=vdiff_ul
            
#         print 'vdiff_ul'
#         print vdiff_ul
#         print 'vdiff_up'
#         print vdiff_up
#         print 'vdiff_ur'
#         print vdiff_ur
#         print 'vdiff_left'
#         print vdiff_left
#         print 'vdiff_right'
#         print vdiff_right
#         print 'vdiff_dl'
#         print vdiff_dl
#         print 'vdiff_down'
#         print vdiff_down
#         print 'vdiff_dr'
#         print vdiff_dr
        
        vdrop_3d=numpy.zeros((8,self.state['nrows'],self.state['ncols']),dtype='float64')
        vdrop_3d[0,:,:]=vdiff_ul
        vdrop_3d[1,:,:]=vdiff_up
        vdrop_3d[2,:,:]=vdiff_ur
        
        vdrop_3d[3,:,:]=vdiff_left
        vdrop_3d[4,:,:]=vdiff_right
        vdrop_3d[5,:,:]=vdiff_dl
        vdrop_3d[6,:,:]=vdiff_down
        vdrop_3d[7,:,:]=vdiff_dr
        

        
        vdrop_3d=numpy.sort(vdrop_3d, axis=0)        
        
        vdrop_min=where(vdrop_3d[0,:,:]<0,vdrop_3d[0,:,:],0)
        vdrop_max=where(vdrop_3d[7,:,:]>0,vdrop_3d[7,:,:],0)
        vdrop_map=(vdrop_max-vdrop_min)/2 #FIXME: scale by per-cell resistances!

        ind = node_map == 0
        vdrop_map[where(ind)] = self.state['nodata']

        return vdrop_map


    def shift_up(self,cells):
        return concatenate((cells[1:], cells[:1]))

    def shift_up(self,cells):
        return concatenate((cells[1:], cells[:1]))
    
    def shift_down(self,cells):
        return concatenate((cells[-1:], cells[:-1]))
    
    def shift_left(self,cells):
        return transpose(self.shift_up(transpose(cells)))
    
    def shift_right(self,cells):
        return transpose(self.shift_down(transpose(cells)))
        














