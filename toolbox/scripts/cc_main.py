#!/usr/bin/env python2.6
# Author: Darren Kavanagh

"""Climate Corridor

Reguired Software:
ArcGIS 10.x with Spatial Analyst extension
Python 2.6
Numpy 1.3

"""

# $Revision$

import os
import sys
import csv
import itertools
from datetime import datetime, timedelta
import traceback
import pdb

import arcpy
import arcpy.sa as sa

from cc_config import cc_env
import cc_util
import lm_master

_SCRIPT_NAME = "cc_main.py"
__version__ = "0.0.3"


FR_COL = "From_Core"
TO_COL = "To_Core"

def main(argv=None):
    """Main function for Climate Corridor tool"""
    tformat = "%m/%d/%y %H:%M:%S"
    stime = datetime.now()

    arcpy.AddMessage("CLIMATE CORRIDOR " + __version__)
    arcpy.AddMessage("Start time: " + stime.strftime(tformat))

    zonal_tbl = "zstats.dbf"

    if argv is None:
        argv = sys.argv
    try:
        cc_env.configure(argv)

        import cc_grass_cwd

        # Check out the ArcGIS Spatial Analyst extension license
        #arcpy.AddMessage(arcpy.ProductInfo())
        arcpy.CheckOutExtension("Spatial")

        # Setup workspace
        arcpy.AddMessage("Creating output folder - " + cc_env.out_dir)
        arcpy.env.overwriteOutput = True
        cc_util.mk_proj_dir(cc_env.out_dir)
        arcpy.env.workspace = cc_env.out_dir
        lm_outdir = cc_util.mk_proj_dir("lm_out")

        # Clip inputs and create project area raster
        cc_clip_inputs()

        # Get zonal statistics for cores and climate
        arcpy.AddMessage("\nCALCULATING ZONAL STATISTICS FROM CLIMATE RASTER")
        climate_stats = arcpy.sa.ZonalStatisticsAsTable(cc_env.prj_core_fc,
            cc_env.core_fld, cc_env.prj_climate_rast, zonal_tbl, "DATA", "ALL")

        # Create core parings table
        core_parings = create_pair_tbl()

        # Limit core pairs based upon cliamte threashold
        limit_cores(core_parings, climate_stats)

        # Calculate distances and generate link table for Linkage Mapper
        core_list = gen_link_table(core_parings, lm_outdir)

        # Create CWD using Grass
        cc_grass_cwd.main(core_list)

        # Run Linkage Mapper
        arcpy.AddMessage("\nRUNNING LINKAGE MAPPER TO CREATE CLIMATE CORRIDORS")
        arg = (_SCRIPT_NAME, lm_outdir, cc_env.prj_core_fc, cc_env.core_fld,
               cc_env.prj_resist_rast, "false", "false", "#", "#", "true", 
               "false", "false", "4", "Cost-Weighted", "true", "true", "#", 
               "#", "#")
        lm_master.lm_master(arg)

    except arcpy.ExecuteError:
        msgs = "ArcPy ERRORS\n" + arcpy.GetMessages(2) + "\n"
        arcpy.AddError(msgs)
        pymsg = ("PYTHON ERRORS\n" + traceback.format_exc())
        arcpy.AddError(pymsg)
    except Exception:
        # tb = sys.exc_info()[2]
        # tbinfo = traceback.format_tb(tb)[0]
        pymsg = ("PYTHON ERRORS\n" + traceback.format_exc())
        arcpy.AddError(pymsg)
    finally:
        arcpy.CheckInExtension("Spatial")
        etime = datetime.now()
        rtime = etime - stime
        hours, minutes = ((rtime.days * 24 + rtime.seconds//3600),
                       (rtime.seconds//60)%60)
        arcpy.AddMessage("End time: " + etime.strftime(tformat))
        arcpy.AddMessage('Elapsed time: %s hrs %s mins' % (hours, minutes))


def cc_clip_inputs():
    """Clip Climate Corridor inputs to smallest extent"""
    ext_poly = "ext_poly.shp"  # Extent polygon
    try:
        arcpy.AddMessage("\nCLIPPING TO SMALLEST EXTENT")

        # Set environment so that smallest extent is used
        arcpy.env.extent = "MINOF"

        # Check if resistance raster is an input and get full path of input
        # raster/s in case they are raster layers (required for conversion)
        if cc_env.resist_rast is None:
            rast_list = arcpy.Describe(cc_env.climate_rast).catalogPath
        else:
            rast_list = (arcpy.Describe(cc_env.climate_rast).catalogPath + ";"
                        + arcpy.Describe(cc_env.resist_rast).catalogPath)
            if arcpy.Exists(cc_env.prj_resist_rast):
                arcpy.Delete_management(cc_env.prj_resist_rast)
        if arcpy.Exists(cc_env.prj_climate_rast):
            arcpy.Delete_management(cc_env.prj_climate_rast)

        arcpy.RasterToOtherFormat_conversion(rast_list, cc_env.out_dir, "GRID")
        arcpy.env.extent = None

        # Create project area raster
        proj_area_rast = sa.Con(sa.IsNull(cc_env.prj_climate_rast),
                                sa.Int(cc_env.prj_climate_rast), 1)
        proj_area_rast.save(cc_env.prj_area_rast)

        # Clip core feature class
        arcpy.RasterToPolygon_conversion(proj_area_rast, ext_poly,
                                         "NO_SIMPLIFY", "VALUE")
        arcpy.Clip_analysis(cc_env.core_fc, ext_poly,
                            cc_env.prj_core_fc)
    except Exception:
        raise
    finally:
        if arcpy.Exists(ext_poly):
            arcpy.Delete_management(ext_poly)


def create_pair_tbl():
    """Create table with all possible core to core combinations"""
    cpair_tbl = "corepairs.dbf"

    try:
        arcpy.AddMessage("\nCREATING CORE PAIRINGS TABLE")
        arcpy.CreateTable_management(cc_env.out_dir, cpair_tbl, "", "")
        arcpy.AddField_management(cpair_tbl, FR_COL, "Short", "", "",
                                      "", "", "NON_NULLABLE")
        arcpy.AddField_management(cpair_tbl, TO_COL, "Short", "", "",
                                      "", "", "NON_NULLABLE")
        arcpy.DeleteField_management(cpair_tbl, "Field1")

        srow, srows, outputrow, irows = None, None, None, None

        srows = arcpy.SearchCursor(cc_env.prj_core_fc, "", "",
                                   cc_env.core_fld, cc_env.core_fld + " A")

        core_list = [srow.getValue(cc_env.core_fld) for srow in srows]
        cores_product = list(itertools.combinations(core_list, 2))

        arcpy.AddMessage("There are " + str(len(core_list)) + " unique "
            "cores and " + str(len(cores_product)) + " pairings")

        irows = arcpy.InsertCursor(cpair_tbl)
        for x in cores_product:
            outputrow = irows.newRow()
            outputrow.setValue(FR_COL, int(x[0]))
            outputrow.setValue(TO_COL, int(x[1]))
            irows.insertRow(outputrow)

        return cpair_tbl

    except:
            raise
    finally:
        if srow:
            del srow
        if srows:
            del srows
        if outputrow:
            del outputrow
        if irows:
            del irows


def limit_cores(pair_tbl, stats_tbl):
    """Limit core pairs based upon climate threashold"""
    pair_vw = "dist_tbvw"
    stats_vw = "stats_tbvw"
    core_id = cc_env.core_fld.upper()

    try:
        arcpy.AddMessage("\nLIMITING CORE PAIRS BASED UPON CLIMATE "
                         "THREASHOLD")

        arcpy.MakeTableView_management(pair_tbl, pair_vw)
        arcpy.MakeTableView_management(stats_tbl, stats_vw)

        def add_stats(fld_pre, table_vw, join_col):
            """Add zonal and calculated statistics to stick table"""
            tmp_mea = fld_pre + "_tmp_mea"
            tmp_std = fld_pre + "_tmp_std"
            umin2std = fld_pre + "umin2std"

            # Add fields to stick table - has to be done before join
            arcpy.AddField_management(table_vw, tmp_mea, "Float", "", "",
                                      "", "", "NULLABLE")
            arcpy.AddField_management(table_vw, tmp_std, "Float", "", "",
                                      "", "", "NULLABLE")
            arcpy.AddField_management(table_vw, umin2std, "Float", "", "",
                                      "", "", "NULLABLE")

            # Join distance table to zonal stats table
            arcpy.AddIndex_management(table_vw, FR_COL, "fridx", "NON_UNIQUE",
                                      "ASCENDING")
            arcpy.AddIndex_management(table_vw, TO_COL, "toidx", "NON_UNIQUE",
                                      "ASCENDING")
            arcpy.AddIndex_management(stats_vw, core_id, "coreidx", "UNIQUE",
                                      "ASCENDING")
            arcpy.AddJoin_management(table_vw, join_col, stats_vw, core_id)

            tbl_name = arcpy.Describe(table_vw).baseName
            stats_tbl_nm = arcpy.Describe(stats_vw).baseName

            # Insert values into fields
            mean_value = "!" + stats_tbl_nm + ".MEAN" + "!"
            std_value = "!" + stats_tbl_nm + ".STD" + "!"
            mea_fld = "!" + tbl_name + "." + tmp_mea + "!"
            std_fld = "!" + tbl_name + "." + tmp_std + "!"

            arcpy.CalculateField_management(table_vw, tmp_mea, mean_value,
                                            "PYTHON")
            arcpy.CalculateField_management(table_vw, tmp_std, std_value,
                                            "PYTHON")
            expression = mea_fld + " - " + std_fld + " - " + std_fld
            arcpy.CalculateField_management(table_vw, umin2std, expression,
                                            "PYTHON")

            # Remove join
            arcpy.RemoveJoin_management(table_vw, stats_tbl_nm)

        # Add basic stats to distance table
        arcpy.AddMessage("Joining to zonal statistics to parings table")
        add_stats("fr", pair_vw, TO_COL)
        add_stats("to", pair_vw, FR_COL)

        # Calculate difference of 2 std
        arcpy.AddMessage("Calculating difference of 2 std")
        diffu_2std = "diffu_2std"
        frumin2std_fld = "!frumin2std!"
        toumin2std_fld = "!toumin2std!"
        arcpy.AddField_management(pair_vw, diffu_2std, "Float", "", "",
                                   "", "", "NULLABLE")
        expression = "abs(" + frumin2std_fld + " - " + toumin2std_fld + ")"
        arcpy.CalculateField_management(pair_vw, diffu_2std, expression,
                                        "PYTHON")

        # Filter distance table based on inputed threashold and delete rows
        arcpy.AddMessage("Filtering table based on threashold")
        diffu2std_fld = arcpy.AddFieldDelimiters(pair_vw, diffu_2std)
        expression = diffu2std_fld + " <= " + str(cc_env.climate_threashold)
        arcpy.SelectLayerByAttribute_management(pair_vw, "NEW_SELECTION",
                                                expression)
        rows_del = int(arcpy.GetCount_management(pair_vw).getOutput(0))
        if rows_del > 0:
            arcpy.DeleteRows_management(pair_vw)
        arcpy.AddMessage(str(rows_del) + " rows deleted")

    except:
        raise
    finally:
        if arcpy.Exists(stats_vw):
            arcpy.Delete_management(stats_vw)
        if arcpy.Exists(pair_vw):
            arcpy.Delete_management(pair_vw)

def gen_link_table(parings, outdir):
    """Calculate core to core distances and create linkage table

    Requires ArcInfo license.

    """
    coresim = "coresim"
    fcore_vw = "fcore_vw"
    tcore_vw = "tcore_vw"
    ntblid_fn = "FID"
    fmid_fn = "IN_FID"
    toid_fn = "NEAR_FID"
    ndist_fn = "NEAR_DIST"
    jtocore_fn = cc_env.core_fld[:8] + "_1"  # dbf field length
    near_tbl = os.path.join(cc_env.out_dir, "neartbl.dbf")
    link_fld = os.path.join(outdir, "datapass")
    link_file = os.path.join(link_fld, "linkTable_s2.csv")

    link_tbl, srow, srows = None, None, None

    try:
        arcpy.AddMessage("\nGENERATING LINKAGE TABLE")
        if cc_env.simplyfy_cores:
            arcpy.AddMessage("Simplifying polygons to speed up core pair "
                             "distance calculations")
            corefc = os.path.join(cc_env.out_dir, coresim + ".shp")
            climate_rast = arcpy.Raster(cc_env.prj_climate_rast)
            tolerance = climate_rast.meanCellHeight / 3
            arcpy.cartography.SimplifyPolygon(cc_env.prj_core_fc, corefc,
                "POINT_REMOVE", tolerance, "#", "NO_CHECK")
        else:
            corefc = cc_env.prj_core_fc

        # Get list of  core pairing and 'from cores'
        frm_cores = set()
        core_pairs = []
        srows = arcpy.SearchCursor(parings, "", "", FR_COL + "; " + TO_COL)
        for srow in srows:
           from_core = srow.getValue(FR_COL)
           to_core = str(srow.getValue(TO_COL))
           frm_cores.add(from_core)
           core_pairs.append([str(from_core), to_core])
        # del srows, srow
        frm_cores = map(str, sorted(frm_cores))

        # Generate link table file from near table results
        if not os.path.exists(link_fld):
            os.mkdir(link_fld)

        link_tbl = open(link_file, 'wb')
        writer = csv.writer(link_tbl, delimiter=',')
        headings = ["# link", "coreId1", "coreId2", "cluster1", "cluster2",
                    "linkType", "eucDist", "lcDist", "eucAdj", "cwdAdj"]
        writer.writerow(headings)        
        
        core_list = set()
        i = 1

        coreid_fld = arcpy.AddFieldDelimiters(corefc, cc_env.core_fld)
        for frm_core in frm_cores:
            # From cores
            expression = coreid_fld + " = " + frm_core
            arcpy.MakeFeatureLayer_management(corefc, fcore_vw, expression)

            # To cores
            to_cores_lst = [x[1] for x in core_pairs if frm_core == x[0]]
            to_cores = ', '.join(to_cores_lst)
            expression = coreid_fld + " in (" +  to_cores + ")"
            arcpy.MakeFeatureLayer_management(corefc, tcore_vw, expression)
            arcpy.AddMessage("Calculating Euclidean distance/s from core "
                             + frm_core + " to " + str(len(to_cores_lst))
                             + " other cores")

            # Generate near table for these core pairings
            arcpy.GenerateNearTable_analysis(fcore_vw, tcore_vw, near_tbl,
                 cc_env.max_euc_dist, "NO_LOCATION", "NO_ANGLE", "ALL")

            # Join near table to core table
            arcpy.JoinField_management(near_tbl, fmid_fn, corefc,
                                       ntblid_fn, cc_env.core_fld)
            arcpy.JoinField_management(near_tbl, toid_fn, corefc,
                                       ntblid_fn, cc_env.core_fld)

            # Limit pairings based on inputed Euclidean distances
            srow, srows = None, None
            euc_dist_fld = arcpy.AddFieldDelimiters(near_tbl, ndist_fn)
            expression = (euc_dist_fld + " > " + str(cc_env.min_euc_dist))
            srows = arcpy.SearchCursor(near_tbl, expression, "",
                                       jtocore_fn + "; " + ndist_fn,
                                       jtocore_fn + " A; " + ndist_fn + " A")

            # Process near table and output into a link table
            srow = srows.next()
            
            if srow:
                core_list.add(int(frm_core))
                while srow:
                    to_coreid = srow.getValue(jtocore_fn)
                    dist_value = srow.getValue(ndist_fn)
                    writer.writerow([i, frm_core, to_coreid, -1, -1, 1,
                                     dist_value, -1, -1, -1])
                    core_list.add(to_coreid)
                    srow = srows.next()
                    i += 1
                # del srows, srow            
                
        return map(str, sorted(core_list))

    except:
        raise
    finally:
        simplify_points = coresim + "_Pnt.shp"
        if arcpy.Exists(simplify_points):
                arcpy.Delete_management(simplify_points)

        if link_tbl:
            link_tbl.close()
        if srow:
            del srow
        if srows:
            del srows


if __name__ == "__main__":
    sys.exit(main())
