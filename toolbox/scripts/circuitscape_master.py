#!/usr/bin/env python2.5

"""Master script for circuitscape analysis in linkage mapper.

Reguired Software:
ArcGIS 9.3 with Spatial Analyst extension
Python 2.5
Numpy

"""

__filename__ = "Circuitscape_master.py"
__version__ = "CIRCUITSCAPE TEST"

import os.path as path
import arcgisscripting
from lm_config import Config as Cfg
import lm_util as lu
import shutil

import os

gp = Cfg.gp
gprint = gp.addmessage

def pinch_master():
    """
    
    """
    try:
        def createfolder(lmfolder):
            """Creates folder if it doesn't exist."""
            if not path.exists(lmfolder):
                gp.CreateFolder_management(path.dirname(lmfolder),
                                               path.basename(lmfolder))               


        gp.RefreshCatalog(Cfg.OUTPUTDIR)
        
        if Cfg.DOPINCH == False and Cfg.DOCENTRALITY == False:            
            msg = ('ERROR: Please choose at least one option: pinch point or\n'
                    'network centrality analysis.')
            gp.AddError(msg)
            exit(1)    

            
        if Cfg.DOPINCH == True:     
            gprint("Creating output folder: " + Cfg.CIRCUITBASEDIR)
            if path.exists(Cfg.CIRCUITBASEDIR):
                shutil.rmtree(Cfg.CIRCUITBASEDIR)
            createfolder(Cfg.CIRCUITBASEDIR)
            gp.CreateFolder_management(Cfg.CIRCUITBASEDIR, 
                                        Cfg.CIRCUITOUTPUTDIR_NM)
            gp.CreateFolder_management(Cfg.CIRCUITBASEDIR, 
                                        Cfg.CIRCUITCONFIGDIR_NM)                 
            createfolder(Cfg.SCRATCHDIR)    
            
            import s7_pinchpoints as s7            
            s7.STEP7_calc_pinchpoints()            

        if Cfg.DOCENTRALITY == True:            
            gprint("Creating output folder: " + Cfg.CENTRALITYBASEDIR)
            if path.exists(Cfg.CENTRALITYBASEDIR):
                shutil.rmtree(Cfg.CENTRALITYBASEDIR)
            createfolder(Cfg.CENTRALITYBASEDIR)
            gp.CreateFolder_management(Cfg.CENTRALITYBASEDIR, 
                                        Cfg.CIRCUITOUTPUTDIR_NM)
            gp.CreateFolder_management(Cfg.CENTRALITYBASEDIR, 
                                        Cfg.CIRCUITCONFIGDIR_NM)    
            import s8_centrality as s8
            s8.STEP8_calc_centrality()

        gprint('\nDONE!\n')

        

    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.raise_geoproc_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        lu.raise_python_error(__filename__)
       
        
if __name__ == "__main__":
    pinch_master()
