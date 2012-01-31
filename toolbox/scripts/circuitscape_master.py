#!/usr/bin/env python2.5

"""Master script for circuitscape analysis in linkage mapper.

Reguired Software:
ArcGIS 9.3 with Spatial Analyst extension
Python 2.5
Numpy

"""

__filename__ = "Circuitscape_master.py"
__version__ = "0.7.7beta-a"

import os.path as path
import arcgisscripting
from lm_config import Config as Cfg
import lm_util as lu
import shutil

import os

import s8_pinchpoints as s8 
import s7_centrality as s7

gp = Cfg.gp
if not Cfg.LOGMESSAGES:
    gprint = gp.addmessage
else:
    gprint = lu.gprint

def circuitscape_master():
    """
    
    """
    try:
        lu.createfolder(Cfg.LOGDIR)
        lu.createfolder(Cfg.MESSAGEDIR)
        Cfg.logFile=lu.create_log_file(Cfg.MESSAGEDIR, Cfg.TOOL, Cfg.PARAMS)
        
        # Check core ID field.
        lu.check_cores(Cfg.COREFC, Cfg.COREFN) 
        
        gp.OutputCoordinateSystem = gp.describe(Cfg.COREFC).SpatialReference
        gp.pyramid = "NONE"
        gp.rasterstatistics = "NONE"              
                                               
        # Move adj and cwd results from earlier versions to datapass directory
        lu.move_old_results()
              
        lu.delete_dir(Cfg.SCRATCHDIR)
                     
        if Cfg.DOPINCH == False and Cfg.DOCENTRALITY == False:            
            msg = ('ERROR: Please choose at least one option: pinch point or\n'
                    'network centrality analysis.')
            lu.raise_error(msg)

        lu.createfolder(Cfg.SCRATCHDIR) 
        lu.createfolder(Cfg.ARCSCRATCHDIR)
        
        if Cfg.DO_ALLPAIRS == True:
            #  Fixme: move raster path to config
            S5CORRIDORRAS = os.path.join(Cfg.OUTPUTGDB,Cfg.PREFIX + 
                                         "_lcc_mosaic_int") 
            if not gp.Exists(S5CORRIDORRAS):
                msg = ('ERROR: Corridor raster created in step 5 is required'
                        '\nfor all-pair analyses, but was not found.')
                lu.raise_error(msg)
        if Cfg.DOPINCH == True:
            # Make a local grid copy of resistance raster-
            # will run faster than gdb.
            lu.delete_data(Cfg.RESRAST)
            if not gp.Exists(Cfg.RESRAST_IN):
                msg = ('ERROR: Resistance raster is required for pinch point'
                        ' analyses, but was not found.')
                lu.raise_error(msg)
                  
            arcpy.env.extent = Cfg.RESRAST_IN
            arcpy.env.snapRaster = Cfg.RESRAST_IN
            gprint('\nMaking local copy of resistance raster.')
            gp.CopyRaster_management(Cfg.RESRAST_IN, Cfg.RESRAST)          
                    
        
        if Cfg.DOCENTRALITY == True:             
            gprint("Creating output folder: " + Cfg.CENTRALITYBASEDIR)
            if path.exists(Cfg.CENTRALITYBASEDIR):
                shutil.rmtree(Cfg.CENTRALITYBASEDIR)
            lu.createfolder(Cfg.CENTRALITYBASEDIR)
            gp.CreateFolder_management(Cfg.CENTRALITYBASEDIR, 
                                        Cfg.CIRCUITOUTPUTDIR_NM)
            gp.CreateFolder_management(Cfg.CENTRALITYBASEDIR, 
                                        Cfg.CIRCUITCONFIGDIR_NM)    
            lu.clean_out_workspace(Cfg.CORECENTRALITYGDB)
            
            s7.STEP7_calc_centrality()
            if Cfg.SAVECENTRALITYDIR == False:
                lu.delete_dir(Cfg.CENTRALITYBASEDIR)
    
        if Cfg.DOPINCH == True:     
            gprint("Creating output folder: " + Cfg.CIRCUITBASEDIR)
            lu.delete_dir(Cfg.CIRCUITBASEDIR)
            lu.createfolder(Cfg.CIRCUITBASEDIR)
            gp.CreateFolder_management(Cfg.CIRCUITBASEDIR, 
                                        Cfg.CIRCUITOUTPUTDIR_NM)
            gp.CreateFolder_management(Cfg.CIRCUITBASEDIR, 
                                        Cfg.CIRCUITCONFIGDIR_NM)                 
            
            lu.clean_out_workspace(Cfg.PINCHGDB)
            lu.delete_data(Cfg.PINCHGDB)
            
            s8.STEP8_calc_pinchpoints()            

            lu.delete_dir(Cfg.SCRATCHDIR)
            if Cfg.SAVECIRCUITDIR == False:
                lu.delete_dir(Cfg.CIRCUITBASEDIR)
            
        gprint('\nDONE!\n')

    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.raise_geoproc_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        lu.raise_python_error(__filename__)
       
        
if __name__ == "__main__":
    circuitscape_master()
