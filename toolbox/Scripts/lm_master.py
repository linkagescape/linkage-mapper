"""Master script for linkage mapper

Reguired Software:
ArcGIS 9.3 with Spatial Analyst extension
Python 2.5
Numpy

"""

__version__ = "0.6.0"

import sys
import os.path as path
import time
import string
import csv

from time import localtime, strftime # remove DMK
from string import split # remove DMK

import arcgisscripting
from numpy import *

import lm_config
import lm_util as lu
import s1_getAdjacencies as s1
import s2_buildNetwork as s2
import s3_calcCwds as s3
import s4_refineNetwork as s4
import s5_calcLccs as s5

PROJECTDIR = lm_config.PROJECTDIR
OUTPUTDIR = lm_config.OUTPUTDIR
LOGDIR = lm_config.LOGDIR
SCRATCHDIR = lm_config.SCRATCHDIR
DATAPASSDIR = lm_config.DATAPASSDIR
STEP1 = lm_config.STEP1
STEP2 = lm_config.STEP2
STEP3 = lm_config.STEP3
STEP4 = lm_config.STEP4
STEP5 = lm_config.STEP5
OUTPUTGDB = lm_config.OUTPUTGDB
GP = lm_config.GP

def lm_master():
    """Main function for linkage mapper

    Called by ArcMap with parameters or run from command line with parameters
    entered in script below.  Calls functions in dedicated scripts for each of
    5 processing steps.

    """
    try:
        lu.check_project_dir()

        if set("# ").intersection(PROJECTDIR):
            GP.AddError("\nSpaces are not allowed in project directory path."
                        "Please avoid using any special characters in "
                        "directory names.")
            sys.exit(1) # end script process

        def createfolder(lmfolder):
            if not path.exists(lmfolder):
                GP.CreateFolder_management(path.dirname(lmfolder),
                                       path.basename(lmfolder))
        createfolder(OUTPUTDIR)
        createfolder(LOGDIR)
        createfolder(SCRATCHDIR)
        createfolder(DATAPASSDIR)

        # Check to make sure no missing steps, and identify first step checked
        lu.check_steps()
        if STEP5: firstStep = 5
        if STEP4:
            firstStep = 4
            #move_stick_maps(firstStep)
        if STEP3:
            firstStep = 3
            #move_stick_maps(firstStep)
        if STEP2:
            lu.check_dist_file()
            firstStep = 2
            #move_stick_maps(firstStep)
        if STEP1:
            firstStep = 1
            #move_stick_maps(firstStep)
        try:
            test=firstStep
        except:
            lu.dashline(1)
            msg = 'ERROR: Please check at least one step.'
            GP.AddError(msg)
            exit(1)

        # Make a backup copy of datapass directory.
        # This will have lcp maps and link tables from previous run.
        lu.archive_datapass()

        lu.clean_up_link_tables(firstStep)

        if STEP5:
            if GP.Exists(OUTPUTGDB):
                GP.addmessage('Deleting geodatabase ' + OUTPUTGDB)
                try:
                    GP.delete_management(OUTPUTGDB)
                except:
                    lu.dashline(1)
                    msg = ('ERROR: Could not remove geodatabase ' + OUTPUTGDB +
                          '. Is it open in ArcMap?')
                    GP.AddError(msg)
                    exit(1)
            #move_stick_maps(firstStep)

        # Step 1
        if STEP1:
            lu.dashline(1)
            GP.addmessage('Running script s1_getAdjacencies.py')
            s1.step1_get_adjacencies()

        # Step 2
        if STEP2:
            lu.dashline(1)
            GP.addmessage('Running script s2_buildNetwork.py')
            s2.step2_build_network()

        # Step 3
        if STEP3:
            lu.dashline(1)
            GP.addmessage('Running script s3_calcCwds.py')
            s3.step3_calc_cwds()  # Run step 3

        # Step 4
        if STEP4 :
            lu.dashline(1)
            GP.addmessage('Running script s4_refineNetwork.py')
            s4.step4_refine_network()  # Run step 4

        # Step 5
        if STEP5 :
            lu.dashline(1)
            GP.addmessage('Running script s5_calcLccs.py')
            s5.step5_calc_lccs() # Run step 5

        GP.addmessage('\nDONE!\n')

    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        filename =  'lm_master.py'
        lu.raise_geoproc_error(filename)

    # Return any PYTHON or system specific errors
    except:
        filename =  'lm_master.py'
        lu.raise_python_error(filename)

if __name__ == "__main__":
    lm_master()