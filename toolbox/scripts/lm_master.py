#!/usr/bin/env python2.5
# Authors: Brad McRae and Darren Kavanagh

"""Master script for Linkage Lapper.

Reguired Software:
ArcGIS 9.3 with Spatial Analyst extension
Python 2.5
Numpy

"""

__filename__ = "lm_master.py"
__version__ = "0.7.6"

import os.path as path
import os

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
        installD = gp.GetInstallInfo("desktop")
        gprint('\nLinkage Mapper Version ' + str(__version__))
        try:
            gprint('on ArcGIS '+ installD['ProductName'] + ' ' + 
                installD['Version'] + ' Service Pack ' + installD['SPNumber'])
        except: pass   
    
        gp.OutputCoordinateSystem = gp.describe(Cfg.COREFC).SpatialReference

        # Check core ID field.
        lu.check_cores()        
       
        if gp.Exists(Cfg.OUTPUTDIR):
            gp.RefreshCatalog(Cfg.OUTPUTDIR)
        
      
        def delete_final_gdb(finalgdb):
            if gp.Exists(finalgdb) and Cfg.STEP5:
                try:
                    lu.clean_out_workspace(finalgdb)
#                    gp.addmessage('Deleting geodatabase ' + finalgdb)
 #                   gp.delete_management(finalgdb)

                except:
                    lu.dashline(1)
                    msg = ('ERROR: Could not remove contents of geodatabase ' +
                           finalgdb + '. Is it open in ArcMap?\n You may '
                           'need to re-start ArcMap to release the file lock.')
                    gp.AddError(msg)
                    exit(1)     
        # Delete final output geodatabase
        delete_final_gdb(Cfg.OUTPUTGDB_OLD)
        delete_final_gdb(Cfg.OUTPUTGDB)
        delete_final_gdb(Cfg.EXTRAGDB)
        delete_final_gdb(Cfg.LINKMAPGDB)        

  
        lu.createfolder(Cfg.OUTPUTDIR)
        lu.createfolder(Cfg.LOGDIR)
        lu.createfolder(Cfg.DATAPASSDIR)
        # Create fresh scratch directory
        lu.delete_dir(Cfg.SCRATCHDIR)
        lu.createfolder(Cfg.SCRATCHDIR)
        
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
        gp.OverwriteOutput = True
        
        if Cfg.STEP1 or Cfg.STEP3:
            # Make a local grid copy of resistance raster for cwd runs-
            # will run faster than gdb.
            # Don't know if raster is in a gdb if entered from TOC
            lu.delete_data(Cfg.RESRAST)
            gp.pyramid = "NONE"
            gp.rasterstatistics = "NONE"

            gprint('\nMaking local copy of resistance raster.')
            try:
                gp.CopyRaster_management(Cfg.RESRAST_IN, Cfg.RESRAST)          
            except: # This sometimes fails due to bad file locks
                Cfg.RESRAST = Cfg.RESRAST_IN
        else:
            Cfg.RESRAST = Cfg.RESRAST_IN
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
        lu.clean_out_workspace(Cfg.SCRATCHDIR)

        gp.addmessage('\nDONE!\n')

    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.raise_geoproc_error(__filename__)

    # Return any PYTHON or system specific errors
    except:
        lu.raise_python_error(__filename__)

if __name__ == "__main__":
    lm_master()
