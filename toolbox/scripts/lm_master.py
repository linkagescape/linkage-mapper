#!/usr/bin/env python2.5
# Authors: Brad McRae and Darren Kavanagh

"""Master script for Linkage Lapper.

Reguired Software:
ArcGIS 9.3 with Spatial Analyst extension
Python 2.5
Numpy

"""

__filename__ = "lm_master.py"
__version__ = "0.6.6"

import os.path as path

import arcgisscripting

from lm_config import Config as Cfg
import lm_util as lu
import s1_getAdjacencies as s1
import s2_buildNetwork as s2
import s3_calcCwds as s3
import s4_refineNetwork as s4
import s5_calcLccs as s5

gp = Cfg.gp
gprint = gp.addmessage

def lm_master():
    """Main function for linkage mapper.

    Called by ArcMap with parameters or run from command line with parameters
    entered in script below.  Calls functions in dedicated scripts for each of
    5 processing steps.

    """
    try:
        gprint('\nLinkage Mapper Version ' + str(__version__))
        
        # Check core ID field.
        lu.check_cores()
        
       
        if gp.Exists(Cfg.OUTPUTDIR):
            gp.RefreshCatalog(Cfg.OUTPUTDIR)
        
      
        def delete_final_gdb(finalgdb):
            if gp.Exists(finalgdb) and Cfg.STEP5:
                gp.addmessage('Deleting geodatabase ' + finalgdb)
                try:
                    gp.delete_management(finalgdb)
                except:
                    lu.dashline(1)
                    msg = ('ERROR: Could not remove geodatabase ' +
                           finalgdb + '. Is it open in ArcMap?\n You may '
                           'need to re-start ArcMap to release the file lock.')
                    gp.AddError(msg)
                    exit(1)       
                    
        # Delete final output geodatabase
        delete_final_gdb(Cfg.OUTPUTGDB_OLD)
        delete_final_gdb(Cfg.OUTPUTGDB)
        delete_final_gdb(Cfg.EXTRAGDB)
        delete_final_gdb(Cfg.LINKMAPGDB)        

  
        def createfolder(lmfolder):
            """Creates folder if it doesn't exist."""
            if not path.exists(lmfolder):
                gp.CreateFolder_management(path.dirname(lmfolder),
                                               path.basename(lmfolder))
        createfolder(Cfg.OUTPUTDIR)
        createfolder(Cfg.LOGDIR)
        createfolder(Cfg.DATAPASSDIR)
        # Create fresh scratch directory
        lu.delete_dir(Cfg.SCRATCHDIR)
        createfolder(Cfg.SCRATCHDIR)
        
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

        # Make a local grid copy of resistance raster- will run faster than gdb
        # Don't know if raster is in a gdb if entered from TOC
        lu.delete_data(Cfg.RESRAST)
        gprint('\nMaking local copy of resistance raster.')
        gp.CopyRaster_management(Cfg.RESRAST_IN, Cfg.RESRAST)          
        gp.SnapRaster = Cfg.RESRAST
        
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
        
        # Clean up
        lu.delete_dir(Cfg.SCRATCHDIR)

        gp.addmessage('\nDONE!\n')

    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.raise_geoproc_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        lu.raise_python_error(__filename__)

if __name__ == "__main__":
    lm_master()
