#!/usr/bin/env python2.5

"""Master script for linkage mapper.

Reguired Software:
ArcGIS 9.3 with Spatial Analyst extension
Python 2.5
Numpy

"""

__filename__ = "lm_master.py"
__version__ = "0.6.4"

import os.path as path

import arcgisscripting

from lm_config import Config as Cfg
import lm_util as lu
import s1_getAdjacencies as s1
import s2_buildNetwork as s2
import s3_calcCwds as s3
import s4_refineNetwork as s4
import s5_calcLccs as s5


def lm_master():
    """Main function for linkage mapper.

    Called by ArcMap with parameters or run from command line with parameters
    entered in script below.  Calls functions in dedicated scripts for each of
    5 processing steps.

    """
    try:
        if Cfg.gp.Exists(Cfg.OUTPUTDIR):
            Cfg.gp.RefreshCatalog(Cfg.OUTPUTDIR)
        
        # Delete final ouptut geodatabase
        if Cfg.gp.Exists(Cfg.OUTPUTGDB_OLD) and Cfg.STEP5:
            try:
                Cfg.gp.delete_management(Cfg.OUTPUTGDB_OLD)
            except:
                pass
        if Cfg.gp.Exists(Cfg.OUTPUTGDB) and Cfg.STEP5:
            Cfg.gp.addmessage('Deleting geodatabase ' + Cfg.OUTPUTGDB)
            try:
                Cfg.gp.delete_management(Cfg.OUTPUTGDB)
            except:
                lu.dashline(1)
                msg = ('ERROR: Could not remove geodatabase ' +
                       Cfg.OUTPUTGDB + '. Is it open in ArcMap?\n You may '
                       'need to re-start ArcMap to release the file lock.')
                Cfg.gp.AddError(msg)
                exit(1)

        # Delete final link map geodatabase
        if Cfg.gp.Exists(Cfg.LINKMAPGDB) and Cfg.STEP5:
            Cfg.gp.addmessage('Deleting geodatabase ' + Cfg.LINKMAPGDB)
            try:
                Cfg.gp.delete_management(Cfg.LINKMAPGDB)
            except:
                lu.dashline(1)
                msg = ('ERROR: Could not remove geodatabase ' +
                       Cfg.LINKMAPGDB + '. Is it open in ArcMap?\n You may '
                       'need to re-start ArcMap to release the file lock.')
                Cfg.gp.AddError(msg)
                exit(1)
                
                
                
                
        def createfolder(lmfolder):
            """Creates folder if it doesn't exist."""
            if not path.exists(lmfolder):
                Cfg.gp.CreateFolder_management(path.dirname(lmfolder),
                                               path.basename(lmfolder))
        createfolder(Cfg.OUTPUTDIR)
        createfolder(Cfg.LOGDIR)
        createfolder(Cfg.SCRATCHDIR)
        createfolder(Cfg.DATAPASSDIR)

        # Identify first step cleanup link tables from that point
        if Cfg.STEP1:
            firststep = 1
        elif Cfg.STEP2:
            firststep = 2
        elif Cfg.STEP3:
            firststep = 3
        elif Cfg.STEP4:
            firststep = 4
        elif Cfg.STEP5:
            firststep = 5
        lu.clean_up_link_tables(firststep)

        # Move adj and cwd results from earlier versions to datapass directory
        lu.move_old_results()

        # Run linkage mapper processing steps
        if Cfg.STEP1:
            s1.STEP1_get_adjacencies()
        if Cfg.STEP2:
            s2.STEP2_build_network()
        if Cfg.STEP3:
            s3.STEP3_calc_cwds()
        if Cfg.STEP4:
            s4.STEP4_refine_network()
        if Cfg.STEP5:
            s5.STEP5_calc_lccs()

        Cfg.gp.addmessage('\nDONE!\n')

    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.raise_geoproc_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        lu.raise_python_error(__filename__)

if __name__ == "__main__":
    lm_master()
