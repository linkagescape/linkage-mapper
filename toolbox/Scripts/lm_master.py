#!/usr/bin/env python2.5

"""Master script for linkage mapper

Reguired Software:
ArcGIS 9.3 with Spatial Analyst extension
Python 2.5
Numpy

"""

__version__ = "0.6.0"

import sys
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

        lu.check_project_dir()

        if set("# ").intersection(Cfg.PROJECTDIR):
            Cfg.gp.AddError("\nSpaces are not allowed in project directory "
                            "path. Please avoid using any special characters "
                            "in directory names.")
            sys.exit(1)  # end script process

        def createfolder(lmfolder):
            """Creates folder if it doesn't exist."""
            if not path.exists(lmfolder):
                Cfg.gp.CreateFolder_management(path.dirname(lmfolder),
                                               path.basename(lmfolder))
        createfolder(Cfg.OUTPUTDIR)
        createfolder(Cfg.LOGDIR)
        createfolder(Cfg.SCRATCHDIR)
        createfolder(Cfg.DATAPASSDIR)

        # Check to make sure no missing steps, and identify first step checked
        lu.check_steps()        
        if Cfg.STEP1:        
            firststep = 1
            lu.check_dist_file()
        elif Cfg.STEP2:            
            firststep = 2
            lu.check_dist_file()
        elif Cfg.STEP3:
            firststep = 3
        elif Cfg.STEP4:
            firststep = 4
        elif Cfg.STEP5:
            firststep = 5
        else:
            lu.dashline(1)
            msg = 'ERROR: Please check at least one step.'
            Cfg.gp.AddError(msg)
            exit(1)

        # Make a backup copy of datapass directory.
        # This will have lcp maps and link tables from previous run.
        lu.archive_datapass()

        lu.clean_up_link_tables(firststep)

        if Cfg.STEP5:
            if Cfg.gp.Exists(Cfg.OUTPUTGDB):
                Cfg.gp.RefreshCatalog(Cfg.OUTPUTGDB)
                Cfg.gp.addmessage('Deleting geodatabase ' + Cfg.OUTPUTGDB)
                try:
                    Cfg.gp.delete_management(Cfg.OUTPUTGDB)
                except:
                    lu.dashline(1)
                    msg = ('ERROR: Could not remove geodatabase ' +
                           Cfg.OUTPUTGDB + '. Is it open in ArcMap?')
                    Cfg.gp.AddError(msg)
                    exit(1)

        # Step 1
        if Cfg.STEP1:
            lu.dashline(1)
            Cfg.gp.addmessage('Running script s1_getAdjacencies.py')
            s1.STEP1_get_adjacencies()

        # Step 2
        if Cfg.STEP2:
            lu.dashline(1)
            Cfg.gp.addmessage('Running script s2_buildNetwork.py')
            s2.STEP2_build_network()

        # Step 3
        if Cfg.STEP3:
            lu.dashline(1)
            Cfg.gp.addmessage('Running script s3_calcCwds.py')
            s3.STEP3_calc_cwds()  # Run step 3

        # Step 4
        if Cfg.STEP4:
            lu.dashline(1)
            Cfg.gp.addmessage('Running script s4_refineNetwork.py')
            s4.STEP4_refine_network()  # Run step 4

        # Step 5
        if Cfg.STEP5:
            lu.dashline(1)
            Cfg.gp.addmessage('Running script s5_calcLccs.py')
            s5.STEP5_calc_lccs()  # Run step 5

        Cfg.gp.addmessage('\nDONE!\n')
        del Cfg.gp
        
    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        filename = 'lm_master.py'
        lu.raise_geoproc_error(filename)

    # Return any PYTHON or system specific errors
    except:
        filename = 'lm_master.py'
        lu.raise_python_error(filename)

if __name__ == "__main__":
    lm_master()
