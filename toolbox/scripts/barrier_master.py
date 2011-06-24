#!/usr/bin/env python2.5

"""Master script for barrier analysis in linkage mapper.

Reguired Software:
ArcGIS 9.3 with Spatial Analyst extension
Python 2.5
Numpy

"""

__filename__ = "barrier_master.py"
__version__ = "BARRIER TEST"

import os.path as path
import arcgisscripting
from lm_config import Config as Cfg
import lm_util as lu
import s6_barriers as s6 

def bar_master():
    """
    
    """
    try:
        Cfg.gp.RefreshCatalog(Cfg.OUTPUTDIR)
        
        # Delete final ouptut geodatabase
        if Cfg.gp.Exists(Cfg.BARRIERGDB): 
            Cfg.gp.addmessage('Deleting geodatabase ' + Cfg.BARRIERGDB)
            try:
                Cfg.gp.delete_management(Cfg.BARRIERGDB)
            except:
                lu.dashline(1)
                msg = ('ERROR: Could not remove geodatabase ' +
                       Cfg.BARRIERGDB + '. Is it open in ArcMap?\n You may '
                       'need to re-start ArcMap to release the file lock.')
                Cfg.gp.AddError(msg)
                exit(1)
                
        def createfolder(lmfolder):
            """Creates folder if it doesn't exist."""
            if not path.exists(lmfolder):
                Cfg.gp.CreateFolder_management(path.dirname(lmfolder),
                                               path.basename(lmfolder))
        createfolder(Cfg.OUTPUTDIR)
        createfolder(Cfg.SCRATCHDIR)    
     
        s6.STEP6_calc_barriers()
        Cfg.gp.addmessage('\nDONE!\n')


    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.raise_geoproc_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        lu.raise_python_error(__filename__)

if __name__ == "__main__":
    bar_master()
