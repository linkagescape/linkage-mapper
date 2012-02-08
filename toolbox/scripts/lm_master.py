#!/usr/bin/env python2.5
# Authors: Brad McRae and Darren Kavanagh

"""Master script for Linkage Lapper.

Reguired Software:
ArcGIS 9.3 with Spatial Analyst extension
Python 2.5
Numpy

"""

import os.path as path
import sys

import arcgisscripting

from lm_config import tool_env as cfg
import lm_util as lu
import s1_getAdjacencies as s1
import s2_buildNetwork as s2
import s3_calcCwds as s3
import s4_refineNetwork as s4
import s5_calcLccs as s5

_filename = path.basename(__file__)
#__version__ = "$Revision$"


def lm_master():
    """Main function for linkage mapper.

    Called by ArcMap with parameters or run from command line with parameters
    entered in script below.  Calls functions in dedicated scripts for each of
    5 processing steps.

    """
    # Setup global variables
    cfg.configure("linkage_mapper", sys.argv)
    gp = cfg.gp
    try:
        gprint = lu.gprint

        # Move results from earlier versions to new directory structure
        lu.move_old_results()
        gp.OverwriteOutput = True
        gp.pyramid = "NONE"
        gp.rasterstatistics = "NONE"

        # Create output directories if they don't exist
        if gp.Exists(cfg.OUTPUTDIR):
            gp.RefreshCatalog(cfg.OUTPUTDIR)
        lu.create_dir(cfg.OUTPUTDIR)
        lu.create_dir(cfg.LOGDIR)
        lu.create_dir(cfg.MESSAGEDIR)
        lu.create_dir(cfg.DATAPASSDIR)
        # Create fresh scratch directory if not restarting in midst of step 3
        # if cfg.S2EUCDISTFILE != None:
            # if cfg.S2EUCDISTFILE.lower() == "restart": pass
        # else:
        lu.delete_dir(cfg.SCRATCHDIR)
        lu.create_dir(cfg.SCRATCHDIR)
        lu.create_dir(cfg.ARCSCRATCHDIR)
        cfg.logFile = lu.create_log_file(cfg.MESSAGEDIR, cfg.TOOL, cfg.PARAMS)

        installD = gp.GetInstallInfo("desktop")
        gprint('\nLinkage Mapper Version ' + cfg.releaseNum)
        try:
            gprint('on ArcGIS ' + installD['ProductName'] + ' ' +
                installD['Version'] + ' Service Pack ' + installD['SPNumber'])
        except:
            pass

        # Check core ID field.
        lu.check_cores(cfg.COREFC, cfg.COREFN)

        # Identify first step cleanup link tables from that point
        lu.dashline(1)
        if cfg.STEP1:
            gprint('Starting at step 1.')
            firststep = 1
        elif cfg.STEP2:
            gprint('Starting at step 2.')
            firststep = 2
        elif cfg.STEP3:
            gprint('Starting at step 3.')
            firststep = 3
            linkTableFile = lu.get_prev_step_link_table(step=3)  # Check exists
        elif cfg.STEP4:
            gprint('Starting at step 4.')
            firststep = 4
            linkTableFile = lu.get_prev_step_link_table(step=4)  # Check exists
        elif cfg.STEP5:
            gprint('Starting at step 5.')
            firststep = 5
            linkTableFile = lu.get_prev_step_link_table(step=5)  # Check exists
        lu.clean_up_link_tables(firststep)

        # Make a local grid copy of resistance raster for cwd runs-
        # will run faster than gdb.
        # Don't know if raster is in a gdb if entered from TOC
        lu.delete_data(cfg.RESRAST)
        gprint('\nMaking temporary copy of resistance raster for this run.')
        gp.Extent = gp.Describe(cfg.RESRAST_IN).Extent
        gp.SnapRaster = cfg.RESRAST_IN
        gp.CopyRaster_management(cfg.RESRAST_IN, cfg.RESRAST)

        if (cfg.STEP1) or (cfg.STEP3):
            # Make core raster file
            gprint('\nMaking temporary raster of core file for this run.')
            lu.delete_data(cfg.CORERAS)
            gp.FeatureToRaster_conversion(cfg.COREFC, cfg.COREFN,
                          cfg.CORERAS, gp.Describe(cfg.RESRAST).MeanCellHeight)
         # #   gp.RasterToPolygon_conversion(cfg.CORERAS, cfg.COREFC,
                                              # "NO_SIMPLIFY")

        def delete_final_gdb(finalgdb):
            """Deletes final geodatabase"""
            if gp.Exists(finalgdb) and cfg.STEP5:
                lu.clean_out_workspace(finalgdb)
                lu.delete_data(finalgdb)
                if gp.Exists(finalgdb):
                    lu.dashline(1)
                    msg = ('ERROR: Could not remove contents of geodatabase ' +
                           finalgdb + '. \nIs it open in ArcMap? You may '
                           'need to re-start ArcMap to release the file lock.')
                    lu.raise_error(msg)

        # Delete final output geodatabase
        delete_final_gdb(cfg.OUTPUTGDB_OLD)
        delete_final_gdb(cfg.OUTPUTGDB)
        delete_final_gdb(cfg.EXTRAGDB)
        delete_final_gdb(cfg.LINKMAPGDB)

        gp.OutputCoordinateSystem = gp.describe(cfg.COREFC).SpatialReference
        # Run linkage mapper processing steps
        if cfg.STEP1:
            s1.STEP1_get_adjacencies()
        if cfg.STEP2:
            s2.STEP2_build_network()
        if cfg.STEP3:
            s3.STEP3_calc_cwds()
        if cfg.STEP4:
            s4.STEP4_refine_network()
        if cfg.STEP5:
            s5.STEP5_calc_lccs()

        # Clean up
        lu.delete_dir(cfg.SCRATCHDIR)
        lu.close_log_file()

        gp.addmessage('\nDONE!\n')

        # severity = gp.MaxSeverity
        # gprint(str(severity))
        # test=gp.GetMessages(2)
        # if severity > 1:
            # gprint('Linkage Mapper SUCCEEDED. You can ignore any failure '
                    # 'messages from ArcGIS below.')
         
        return
    # Return GEOPROCESSING specific errors
    except arcgisscripting.ExecuteError:
        lu.print_geoproc_error(_filename)

    # Return any PYTHON or system specific errors
    except:
        lu.print_python_error(_filename)

if __name__ == "__main__":
    lm_master()
