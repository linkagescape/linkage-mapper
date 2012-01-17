#!/usr/bin/env python2.5

"""Master script for barrier analysis in linkage mapper.

Reguired Software:
ArcGIS 9.3 with Spatial Analyst extension
Python 2.5
Numpy

"""

__filename__ = "barrier_master.py"
__version__ = "0.7.6"

import os.path as path
import arcgisscripting
from lm_config import Config as Cfg
import lm_util as lu
import s6_barriers as s6 

gp = Cfg.gp
gprint = gp.addmessage

def bar_master():
    """ Experimental code to detect barriers using cost-weighted distance
    outputs from Linkage Mapper tool.
    
    """
    try:
        gp.RefreshCatalog(Cfg.OUTPUTDIR)
        
        # Move adj and cwd results from earlier versions to datapass directory
        lu.move_old_results()

        # Delete final ouptut geodatabase
        lu.delete_dir(Cfg.BARRIERGDB)
        if not gp.exists(Cfg.BARRIERGDB):
            # Create output geodatabase
            Cfg.gp.createfilegdb(Cfg.OUTPUTDIR, path.basename(Cfg.BARRIERGDB))        
                
        lu.createfolder(Cfg.OUTPUTDIR)
        lu.createfolder(Cfg.SCRATCHDIR) 

        gprint('\nMaking local copy of resistance raster.')
        try:
            gp.CopyRaster_management(Cfg.RESRAST_IN, Cfg.RESRAST)          
        except: # This sometimes fails due to bad file locks
            Cfg.RESRAST = Cfg.RESRAST_IN        
     
        s6.STEP6_calc_barriers()
        
        #clean up
        lu.delete_dir(Cfg.SCRATCHDIR)
        if Cfg.SAVEBARRIERDIR ==  False:
            lu.delete_dir(Cfg.BARRIERBASEDIR)
        gp.addmessage('\nDONE!\n')


    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.raise_geoproc_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        lu.raise_python_error(__filename__)

if __name__ == "__main__":
    bar_master()
