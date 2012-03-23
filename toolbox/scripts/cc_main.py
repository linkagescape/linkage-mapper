#!/usr/bin/env python2.6
# Author: Darren Kavanagh

"""Climate Corridor

Reguired Software:
ArcGIS 10.x with Spatial Analyst extension
Python 2.6
Numpy 1.3

"""

# $Revision$

import csv
import os
import sys
import traceback
import pdb

import arcpy
import arcpy.sa as sa

from cc_config import cc_env
import cc_util
import lm_master
import cc_grass_cwd

_SCRIPT_NAME = "cc_main.py"
__version__ = "0.0.1"


def main(argv=None):
    """Main function for Climate Corridor tool"""

    arcpy.AddMessage("Climate Corridor " + __version__)

    zonal_tbl = "zonal_stats.dbf"

    if argv is None:
        argv = sys.argv
    try:
        cc_env.configure(argv)
        # Check out the ArcGIS Spatial Analyst extension license
        arcpy.CheckOutExtension("Spatial")

        # Setup workspace
        arcpy.AddMessage("Creating output folder - " + cc_env.out_dir)
        arcpy.env.overwriteOutput = True
        cc_util.mk_proj_dir(cc_env.out_dir)
        arcpy.env.workspace = cc_env.out_dir

        # Clip inputs and create project area raster
        cc_clip_inputs()

        # Create distance file in conefor format
        conefor_file = generate_distance_file(cc_env.search_radius)
        #conefor_file = os.path.join(cc_env.out_dir, "coredist.txt")

        # Get zonal statistics for cores and climate
        arcpy.AddMessage("Calculating zonal statistics from climate raster")
        climate_stats = arcpy.sa.ZonalStatisticsAsTable(cc_env.prj_core_fc,
            cc_env.core_fld, cc_env.prj_climate_rast, zonal_tbl, "DATA", "ALL")

        # Run Linkage Mapper to create core sticks file
        arcpy.AddMessage("Running Linkage Mapper to create core sticks file")
        lm_outdir = cc_util.mk_proj_dir("lm_tmp")
        arg = (_SCRIPT_NAME, lm_outdir, cc_env.prj_core_fc, cc_env.core_fld,
               cc_env.prj_resist_rast, "true", "true", "#", conefor_file,
               "false", "true", "false", "4", "Cost-Weighted", "false",
               "false", "#", "#", "#")
        lm_master.lm_master(arg)
        arcpy.env.workspace = cc_env.out_dir

        # Limit core pairs based upon input restrictions
        conefor_file, core_list = limit_cores(lm_outdir, climate_stats)

        # Create CWD using Grass
        cc_grass_cwd.main(core_list)

        # Run Linkage Mapper
        arcpy.AddMessage("Running Linkage Mapper to create climate corridors")
        lm_outdir = cc_util.mk_proj_dir("lm_out")
        arg = (_SCRIPT_NAME, lm_outdir, cc_env.prj_core_fc, cc_env.core_fld,
               cc_env.prj_resist_rast, "true", "true", "#", conefor_file,
               "true", "false", "false", "4", "Cost-Weighted", "true", "true",
               "#", "#", "#")
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


def cc_clip_inputs():
    """Clip Climate Corridor inputs to smallest extent"""
    # Clip rasters
    try:
        arcpy.AddMessage("Clipping to smallest extent")

        ext_poly = "ext_poly.shp"  # Extent polygon
        arcpy.env.extent = "MINOF"

        if cc_env.resist_rast is None:
            rast_list = cc_env.climate_rast
        else:
            rast_list = cc_env.climate_rast + ";" + cc_env.resist_rast
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


def generate_distance_file(search_radius='#'):
    """Use ArcGIS to create Conefor distance file

    Requires ArcInfo license.

    """
    arcpy.AddMessage("Generating Conefor distance file")

    near_tbl = os.path.join(cc_env.out_dir, "neartbl.dbf")
    dist_fname = os.path.join(cc_env.out_dir, "coredist.txt")

    FID_FN = "FID"
    INFID_FN = "IN_FID"
    NEARID_FN = "NEAR_FID"
    NEAR_FN = "NEAR_DIST"
    NEAR_COREFN = cc_env.core_fld + "_1"

    # Generate near table
    arcpy.GenerateNearTable_analysis(cc_env.prj_core_fc, cc_env.prj_core_fc,
         near_tbl, search_radius, "NO_LOCATION", "NO_ANGLE", "ALL", "0")

    # Join near table to core table
    arcpy.JoinField_management(near_tbl, INFID_FN, cc_env.prj_core_fc, FID_FN,
                               cc_env.core_fld)
    arcpy.JoinField_management(near_tbl, NEARID_FN, cc_env.prj_core_fc, FID_FN,
                               cc_env.core_fld)

    # Process near table and output in Conefor format
    row, rows = None, None
    processed_ids = []
    dist_file = open(dist_fname, 'wt')
    writer = csv.writer(dist_file, delimiter='\t')

    rows = arcpy.SearchCursor (near_tbl, "", "", cc_env.core_fld + "; " +
                           NEAR_COREFN + "; " + NEAR_FN, cc_env.core_fld +
                           " A; " + NEAR_COREFN + " A")

    for row in rows:
        core_id = str(row.getValue(cc_env.core_fld))
        near_id = str(row.getValue(NEAR_COREFN))
        current_id = near_id + " " + core_id
        if current_id not in processed_ids:
            outputrow = []
            outputrow.extend([core_id, near_id, row.getValue(NEAR_FN)])
            writer.writerow(outputrow)
            processed_ids.append(core_id + " " + near_id)

    dist_file.close()
    if row:
        del row
    if rows:
        del rows

    return dist_fname

def limit_cores(lm_folder, stats_tbl):
        """Limit core pairs based upon input restrictions"""
        arcpy.AddMessage("Limiting core pairs based upon input "
                         "restrictions")

        sticks_tbl = "sticks_tbl.dbf"
        stick_tview = "sticks_tbl"
        fr_col = "From_Core"
        to_col = "To_Core"
        core_id = cc_env.core_fld.upper()
        sticks_file = os.path.basename(lm_folder) + "_sticks_s2.shp"
        sticks_fc = os.path.join(lm_folder, "output", sticks_file)

        # Create table limiting cores based on inputed Euclidean distances
        euc_dist = "Euc_Dist"
        euc_dist_fld = arcpy.AddFieldDelimiters(sticks_fc, euc_dist)
        expression = (euc_dist_fld + " > " + str(cc_env.min_euc_dist)
            + " and " + euc_dist_fld + " <  " + str(cc_env.max_euc_dist))
        if arcpy.Exists(sticks_tbl):
            arcpy.Delete_management(sticks_tbl)
        arcpy.TableToTable_conversion(sticks_fc, cc_env.out_dir, sticks_tbl,
                                      expression)
        arcpy.MakeTableView_management(sticks_tbl, stick_tview)
        arcpy.AddIndex_management(stick_tview, fr_col, "fridx", "NON_UNIQUE",
                                  "ASCENDING")
        arcpy.AddIndex_management(stick_tview, to_col, "toidx", "NON_UNIQUE",
                                  "ASCENDING")
        arcpy.MakeTableView_management(stats_tbl, "stats_tview")
        arcpy.AddIndex_management("stats_tview", core_id, "coreidx",
                                   "UNIQUE", "ASCENDING")

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

            # Join sticks table to zonal stats table
            arcpy.AddJoin_management(table_vw, join_col, stats_tbl, core_id)

            tbl_name = arcpy.Describe(table_vw).baseName
            stats_tbl_nm = arcpy.Describe(stats_tbl).baseName

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

        # Add basic stats to stick table
        add_stats("fr", stick_tview, to_col)
        add_stats("to", stick_tview, fr_col)

        # Calculate difference of 2 std
        diffu_2std = "diffu_2std"
        frumin2std_fld = "!frumin2std!"
        toumin2std_fld = "!toumin2std!"
        arcpy.AddField_management(stick_tview, diffu_2std, "Float", "", "",
                                   "", "", "NULLABLE")
        expression = "abs(" + frumin2std_fld + " - " + toumin2std_fld + ")"
        arcpy.CalculateField_management(stick_tview, diffu_2std, expression,
                                        "PYTHON")

        # Filter sticks table based on inputed threashold
        row, rows = None, None
        diffu2std_fld = arcpy.AddFieldDelimiters(stick_tview, diffu_2std)
        expression = diffu2std_fld + " > " + str(cc_env.climate_threashold)
        flds = fr_col + "; " + to_col + "; " + euc_dist
        rows = arcpy.SearchCursor (stick_tview, expression, "", flds)
        ndist_fname = os.path.join(cc_env.out_dir, "min_diff_over"
                                   + str(cc_env.climate_threashold) + ".txt")
        ndist_file = open(ndist_fname, 'wt')
        writer = csv.writer(ndist_file, delimiter='\t')
        core_list = set()

        for row in rows:
            outputrow = []
            fr_value = row.getValue(fr_col)
            to_value = row.getValue(to_col)
            dist_value = row.getValue(euc_dist)
            core_list.update([str(fr_value), str(to_value)])
            outputrow.extend([fr_value, to_value, dist_value])
            writer.writerow(outputrow)

        ndist_file.close()
        if row:
            del row
        if rows:
            del rows

        return ndist_fname, sorted(core_list)


if __name__ == "__main__":
    sys.exit(main())
