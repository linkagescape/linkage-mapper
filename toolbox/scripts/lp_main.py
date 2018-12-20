#!/usr/bin/env python2
# Authors: John Gallo and Randal Greene 2017

"""Linkage Priority main module."""

import os
import sys
import glob
import traceback
from datetime import datetime as dt
from collections import namedtuple

import arcpy

from lm_config import tool_env as lm_env
import lm_util


_SCRIPT_NAME = "lp_main.py"
TFORMAT = "%m/%d/%y %H:%M:%S"

NM_SCORE = "SCORE_RANGE"  # Score range normalization
NM_MAX = "MAX_VALUE"  # Maximum value normalization

CoordPoint = namedtuple('Point', 'x y')


def print_runtime(stime):
    """Print process time when running from script."""
    etime = dt.now()
    rtime = etime - stime
    hours, minutes = ((rtime.days * 24 + rtime.seconds // 3600),
                      (rtime.seconds // 60) % 60)
    print "End time: %s" % etime.strftime(TFORMAT)
    print "Elapsed time: %s hrs %s mins" % (hours, minutes)


def delete_dataset(dataset):
    """Delete one dataset."""
    try:
        arcpy.Delete_management(dataset)
    except arcpy.ExecuteError:
        arcpy.AddWarning("Error deleting scratch/intermediate/temporary "
                         "dataset %s. Program will continue." % dataset)


def delete_datasets_in_workspace():
    """Delete all datasets in workspace."""
    datasets = arcpy.ListDatasets()
    for dataset in datasets:
        delete_dataset(dataset)


def delete_arcpy_temp_datasets():
    """Delete datasets left behind by arcpy."""
    if arcpy.Exists(lm_env.SCRATCHGDB):
        # datasets in scratch gdb (including core_resistance_stats and old
        # hangovers)
        arcpy.env.workspace = lm_env.SCRATCHGDB
        delete_datasets_in_workspace()
    # datasets in intermediate gdb (including old hangovers)
    if not lm_env.KEEPINTERMEDIATE:
        if arcpy.Exists(lm_env.INTERGDB):
            arcpy.env.workspace = lm_env.INTERGDB
            delete_datasets_in_workspace()
    # datasets in current OS directory
    arcpy.env.workspace = os.getcwd()
    delete_datasets_in_workspace()


def normalize_raster(in_raster, normalization_method=NM_MAX, invert=False):
    """Normalize values in in_raster.

    Normalize values in in_raster using score range or max score method,
    with optional inversion.
    """
    lm_util.build_stats(in_raster)
    result = arcpy.GetRasterProperties_management(in_raster, "MINIMUM")
    min_val = float(result.getOutput(0))
    result = arcpy.GetRasterProperties_management(in_raster, "MAXIMUM")
    max_val = float(result.getOutput(0))
    if max_val > 0:
        if normalization_method == NM_SCORE:
            if invert:
                return (max_val - in_raster) / (max_val - min_val)
            return (in_raster - min_val) / (max_val - min_val)
        else:  # Max score normalization
            if invert:
                return (max_val + min_val - in_raster) / max_val
            return in_raster / max_val
    else:
        return in_raster * 0


def blended_priority(norm_trunc_raster, lp_raster, bp_raster):
    """Calculate overall Blended Priority."""
    lm_util.gprint("Calculating overall Blended Priority")
    tmp_raster3 = ((lm_env.TRUNCWEIGHT * arcpy.sa.Raster(norm_trunc_raster)) +
                   (lm_env.LPWEIGHT * arcpy.sa.Raster(lp_raster)))
    if lm_env.NORMALIZEBP:
        bp_raster_tmp = normalize_raster(tmp_raster3, NM_SCORE)
    else:
        bp_raster_tmp = tmp_raster3
    arcpy.CopyRaster_management(bp_raster_tmp, bp_raster, None, None, None,
                                None, None, "32_BIT_FLOAT")


def norm_trunc(trunc_raster, norm_trunc_raster):
    """Invert and normalize truncated raster (NORMTRUNC)."""
    lm_util.gprint("Inverting and normalizing truncated raster (NORMTRUNC)")
    norm_trunc_raster_tmp = normalize_raster(arcpy.sa.Raster(trunc_raster),
                                             lm_env.TRUNCNORMETH, True)
    arcpy.CopyRaster_management(norm_trunc_raster_tmp, norm_trunc_raster, None,
                                None, None, None, None, "32_BIT_FLOAT")


def linkage_priority(rci_raster, trunc_raster, lp_raster):
    """Clip RCI to extent of truncated raster (LP)."""
    lm_util.gprint("Calculating overall Linkage Priority")
    # mask RCI based on extent of truncated raster
    # (which is based on CWDTHRESH)
    lp_raster_tmp = arcpy.sa.ExtractByMask(rci_raster, trunc_raster)
    if lm_env.NORMALIZELP:
        lp_raster_tmp = normalize_raster(lp_raster_tmp, NM_SCORE)
    arcpy.CopyRaster_management(lp_raster_tmp, lp_raster, None, None, None,
                                None, None, "32_BIT_FLOAT")


def rci(cpv_raster, rci_raster):
    """Clip CPV to the MINCPV and renormalize.

    Creates relative corridor "importance (RCI).

    """
    lm_util.gprint("Calculating overall Relative Corridor Importance (RCI)")
    tmp_raster = arcpy.sa.ExtractByAttributes(
        cpv_raster, " VALUE >= " + str(lm_env.MINCPV))
    if lm_env.NORMALIZERCI:
        rci_raster_tmp = normalize_raster(tmp_raster, NM_SCORE)
    else:
        rci_raster_tmp = tmp_raster
    arcpy.CopyRaster_management(rci_raster_tmp, rci_raster, None, None, None,
                                None, None, "32_BIT_FLOAT")


def cpv(sum_rasters, cnt_non_null_cells_rast, max_rasters, cpv_raster):
    """Combine CSPs using Max and Mean to create overall CPV.

    CPV - Corridor Priority Value
    """
    lm_util.gprint("Combining CSPs using Max and Mean to create overall "
                   "Corridor Priority Value (CPV)")
    sum_all_raster = arcpy.sa.CellStatistics(sum_rasters, "SUM", "DATA")
    cnt_nonnull_rast = (
        arcpy.sa.CellStatistics(cnt_non_null_cells_rast, "SUM", "DATA"))
    max_all_raster = arcpy.sa.CellStatistics(max_rasters, "MAXIMUM", "DATA")
    cpv_raster_tmp = ((max_all_raster * lm_env.MAXCSPWEIGHT) +
                      ((sum_all_raster / cnt_nonnull_rast)
                       * lm_env.MEANCSPWEIGHT))
    cpv_raster_tmp.save(cpv_raster)
    if not lm_env.KEEPINTERMEDIATE:
        # clean-up CSP input rasters
        # could be multiple nlc folders
        nlc_idx = 0
        prev_ws = arcpy.env.workspace
        while True:
            nlc_str = ""
            if nlc_idx > 0:
                nlc_str = str(nlc_idx)
            if not os.path.exists(os.path.join(lm_env.DATAPASSDIR, "nlcc",
                                               "nlc" + nlc_str)):
                break
            arcpy.env.workspace = os.path.join(lm_env.DATAPASSDIR, "nlcc",
                                               "nlc" + nlc_str, "inv_norm")
            for in_rast in arcpy.ListRasters():
                lm_util.delete_data(in_rast)
            nlc_idx += 1
        arcpy.env.workspace = prev_ws


def save_interm_rast(cp_rast, save_rast):
    """Save intermediate raster."""
    if lm_env.KEEPINTERMEDIATE:
        arcpy.CopyRaster_management(
            cp_rast,
            os.path.join(lm_env.INTERGDB,
                         '_'.join([lm_env.PREFIX, save_rast])),
            None, None, None, None, None, "32_BIT_FLOAT")


def check_add_field(feature_class, field_name, data_type):
    """Check if field exists, and if not then add."""
    field_names = [field.name for field in arcpy.ListFields(feature_class)]
    if field_name not in field_names:
        exists = False
        arcpy.AddField_management(feature_class, field_name, data_type)
    else:
        exists = True
    return exists


def clim_priority_combine(lcp_lines):
    """Combine climate priority (A & L) values.

    Combine climate priority values A & L in a weighted sum to yield Core
    Areas Climate Linkage Priority Value (O).
    """
    lm_util.gprint("Calculating Core Areas Climate Linkage Priority Value")
    check_add_field(lcp_lines, "Clim_Lnk_Priority", "float")
    lnk_rows = arcpy.UpdateCursor(
        lcp_lines,
        fields="Clim_Lnk_Priority; NCLPv_Analog; NCLPv_Prefer")
    for lnk_row in lnk_rows:
        lnk_row.setValue(
            "Clim_Lnk_Priority",
            (lnk_row.getValue("NCLPv_Analog") * lm_env.CANALOG_WEIGHT) +
            (lnk_row.getValue("NCLPv_Prefer") * lm_env.CPREF_WEIGHT))
        lnk_rows.updateRow(lnk_row)
    del lnk_rows


def normalize_field(in_table, in_field, out_field,
                    normalization_method=NM_MAX, invert=False):
    """Normalize values in in_field into out_field.

    Normalize values in in_field into out_field using score range or max
    score method, with optional inversion.
    """
    check_add_field(in_table, out_field, "DOUBLE")
    min_val, max_val = value_range(in_table, in_field)

    if max_val > 0:
        try:
            if normalization_method == NM_SCORE:
                if invert:
                    arcpy.CalculateField_management(
                        in_table, out_field,
                        "(" + str(max_val) + " - !" + in_field + "!) / " +
                        str(max_val - min_val), "PYTHON_9.3")
                else:
                    arcpy.CalculateField_management(
                        in_table, out_field,
                        "(!" + in_field + "! - " + str(min_val) + ") / " +
                        str(max_val - min_val), "PYTHON_9.3")
            else:  # Max score normalization
                if invert:
                    arcpy.CalculateField_management(
                        in_table, out_field,
                        "(" + str(max_val + min_val) + " - !" +
                        in_field + "!) / " + str(max_val), "PYTHON_9.3")
                else:
                    arcpy.CalculateField_management(
                        in_table, out_field,
                        "!" + in_field + "! / " + str(max_val), "PYTHON_9.3")
        except Exception:
            # other exception - assume it was caused by situation described in
            # exception message
            exc_value, exc_traceback = sys.exc_info()[1:]
            arcpy.AddError(exc_value)
            lm_util.gprint("Traceback (most recent call last):\n" +
                           "".join(traceback.format_tb(exc_traceback)))
            raise Exception(
                "ERROR! MOST LIKELY CAUSE: One or more core areas are smaller "
                "than a pixel in the Resistance layer and/or Raster Analysis "
                "Cell Size environment setting. Try enlarging small core "
                "areas, resampling the Resistance layer or adjusting the Cell"
                "Size environment setting.")
    else:
        arcpy.CalculateField_management(in_table, out_field, "0", "PYTHON_9.3")


def clim_priority_val_normal(lcp_lines):
    """Normalize climate priority values (A & L) for each core pair."""
    lm_util.gprint("Normalizing Climate Analog and Preference Linkage "
                   "Priority Values")
    normalize_field(lcp_lines, "CLPv_Analog", "NCLPv_Analog",
                    lm_env.CANALOGNORMETH)
    normalize_field(lcp_lines, "CLPv_Prefer", "NCLPv_Prefer",
                    lm_env.CPREFERNORMETH)


def sline_y_value(x_coord, slope_val, intercept_val):
    """For a x value on a straight line find its corresponding y value.

    The equation of a straight line is: y = mx + b
    where m is the slope of the line and b is the intercept.
    """
    return (slope_val * x_coord) + intercept_val


def intercept(point, slope_val):
    """Find the intercept of a straight line.

    Where (x,y) is a point on the line and b is the slope of the line.
    """
    return point.y - (point.x * slope_val)


def slope(point1, point2):
    """Calculate slope of a straight line.

    Where (x1,y1) and (x2,y2) are points on the line.
    """
    return (point2.y - point1.y) / (point2.x - point1.x)


def clim_lnk_value(xlnk, min_pnt, target_pnt, max_pnt):
    """Calculate climate linkage value using equation of straight line."""
    if target_pnt.y is None:
        line_slope = slope(min_pnt, max_pnt)
        line_intercept = intercept(max_pnt, line_slope)
    elif xlnk < target_pnt.y:
        line_slope = slope(min_pnt, target_pnt)
        line_intercept = intercept(target_pnt, line_slope)
    elif xlnk > target_pnt.y:
        line_slope = slope(target_pnt, max_pnt)
        line_intercept = intercept(max_pnt, line_slope)
    else:
        line_slope = line_intercept = 1
    return sline_y_value(xlnk, line_slope, line_intercept)


def value_range(layer, field):
    """Get value range of field in layer."""
    min_value = arcpy.SearchCursor(
        layer, fields=field, sort_fields=field + " A").next().getValue(field)
    max_value = arcpy.SearchCursor(
        layer, fields=field, sort_fields=field + " D").next().getValue(field)
    return min_value, max_value


def clim_priority_values(lcp_lines):
    """Save climate priority values for each core pair to LCP layer.

    Save Climate Analog Linkage Priority Value (A) and
    Climate Preference Linkage Priority Value (L) for each core pair
    to LCP Layer.
    """
    lm_util.gprint("Calculating Climate Analog and Preference Linkage "
                   "Priority Values")
    check_add_field(lcp_lines, "CLPv_Analog", "float")
    check_add_field(lcp_lines, "CLPv_Prefer", "float")

    canalog_r_min, canalog_r_max = value_range(lcp_lines, "CAnalog_Ratio")
    cprefer_r_min, cprefer_r_max = value_range(lcp_lines, "CPrefer_Ratio")

    canalog_min_pnt = CoordPoint(canalog_r_min, lm_env.CANALOG_MIN)
    canalog_target_pnt = CoordPoint(lm_env.CANALOG_PIORITY,
                                    lm_env.CANALOG_TARGET)
    canalog_max_pnt = CoordPoint(canalog_r_max, lm_env.CANALOG_MAX)

    cprefer_min_pnt = CoordPoint(cprefer_r_min, lm_env.CPREF_MIN)
    cprefer_target_pnt = CoordPoint(1, 1)
    cprefer_max_pnt = CoordPoint(cprefer_r_max, lm_env.CPREF_MAX)

    lnk_rows = arcpy.UpdateCursor(
        lcp_lines,
        fields="CAnalog_Ratio; CPrefer_Ratio; CLPv_Analog; CLPv_Prefer")
    for lnk_row in lnk_rows:
        lnk_row.setValue(
            "CLPv_Analog",
            clim_lnk_value(
                lnk_row.getValue("CAnalog_Ratio"),
                canalog_min_pnt, canalog_target_pnt, canalog_max_pnt))

        lnk_row.setValue(
            "CLPv_Prefer",
            clim_lnk_value(
                lnk_row.getValue("CPrefer_Ratio"),
                cprefer_min_pnt, cprefer_target_pnt, cprefer_max_pnt))
        lnk_rows.updateRow(lnk_row)
    del lnk_rows


def clim_env_read(core_lyr, core, field):
    """Read core climate envelope."""
    rows = arcpy.SearchCursor(
        core_lyr, fields=field,
        where_clause=lm_env.COREFN + " = " + str(core))
    return rows.next().getValue(field)


def clim_ratios(lcp_lines, core_lyr):
    """Calculate and save climate ratios to LCP layer.

    Calculate and save climate ratios to LCP layer and include start core and
    destination.
    """
    lm_util.gprint("Calculating climate ratios for each core paring")
    check_add_field(lcp_lines, "Core_Start", "integer")
    check_add_field(lcp_lines, "Core_End", "integer")
    check_add_field(lcp_lines, "CAnalog_Ratio", "float")
    check_add_field(lcp_lines, "CPrefer_Ratio", "float")

    lnk_rows = arcpy.UpdateCursor(
        lcp_lines,
        fields=("From_Core; To_Core; Core_Start; Core_End; CAnalog_Ratio;"
                "CPrefer_Ratio"))
    for lnk_row in lnk_rows:
        core_start = lnk_row.getValue("From_Core")
        core_dest = lnk_row.getValue("To_Core")
        clim_start = clim_env_read(core_lyr, core_start, "cclim_env")
        clim_dest = clim_env_read(core_lyr, core_dest, "cclim_env")

        if clim_start <= clim_dest and not lm_env.HIGHERCE_COOLER:
            core_start, core_dest = core_dest, core_start
            clim_start, clim_dest = clim_dest, clim_start

        if lm_env.FCERAST_IN:
            clim_dest = clim_env_read(core_lyr, core_dest, "fclim_env")

        lnk_row.setValue("Core_Start", core_start)
        lnk_row.setValue("Core_End", core_dest)
        lnk_row.setValue("CAnalog_Ratio", clim_dest/clim_start)
        lnk_row.setValue("CPrefer_Ratio", clim_dest/lm_env.CPREF_VALUE)
        lnk_rows.updateRow(lnk_row)
    del lnk_rows


def core_mean(in_rast, core_lyr, in_var):
    """Calculate the mean values of a raster within each core area."""
    tbl_name = "_".join(["core", in_var])
    mean_fld = ".".join([lm_env.CORENAME, in_var])
    mean_value = "".join(["[", tbl_name, ".MEAN]"])

    mean_tbl = arcpy.sa.ZonalStatisticsAsTable(
        lm_env.COREFC, lm_env.COREFN, in_rast,
        os.path.join(lm_env.SCRATCHGDB, tbl_name),
        statistics_type="MEAN")
    arcpy.AddJoin_management(core_lyr, lm_env.COREFN, mean_tbl,
                             lm_env.COREFN)
    arcpy.CalculateField_management(core_lyr, mean_fld, mean_value)
    arcpy.RemoveJoin_management(core_lyr)


def clim_envelope(core_lyr):
    """Determine Climate Envelope for each core."""
    lm_util.gprint("Calculating Climate Envelope for each core")
    core_mean(lm_env.CCERAST_IN, core_lyr, "cclim_env")
    if lm_env.FCERAST_IN:
        core_mean(lm_env.FCERAST_IN, core_lyr, "fclim_env")


def clim_linkage_priority(lcp_lines, core_lyr):
    """Calculate Core Areas Climate Linkage Priority Value (O)."""
    clim_envelope(core_lyr)
    clim_ratios(lcp_lines, core_lyr)
    clim_priority_values(lcp_lines)
    clim_priority_val_normal(lcp_lines)
    clim_priority_combine(lcp_lines)


def eciv():
    """Normalize Expert Corridor Importance Value (ECIV) for each corridor."""
    lm_util.gprint("Normalizing Expert Corridor Importance Value (ECIV) "
                   "for each corridor")
    normalize_field(lm_env.COREPAIRSTABLE_IN, lm_env.ECIVFIELD, "neciv",
                    NM_SCORE)


def inv_norm():
    """Invert and normalize each corridor."""
    lm_util.gprint("Inverting and normalizing each corridor")
    prev_ws = arcpy.env.workspace
    # could be multiple nlc folders
    nlc_idx = 0
    while True:
        nlc_str = ""
        if nlc_idx > 0:
            nlc_str = str(nlc_idx)
        if not os.path.exists(
                os.path.join(lm_env.DATAPASSDIR, "nlcc", "nlc" + nlc_str)):
            break
        arcpy.env.workspace = os.path.join(lm_env.DATAPASSDIR, "nlcc",
                                           "nlc" + nlc_str)
        # process each corridor raster in folder
        for input_raster in arcpy.ListRasters():
            # max score normalization with inversion
            inv_norm_raster = normalize_raster(arcpy.sa.Raster(input_raster),
                                               lm_env.NORMCORRNORMETH, True)
            if not os.path.exists(os.path.join(
                    lm_env.DATAPASSDIR, "nlcc", "nlc" + nlc_str, "inv_norm")):
                arcpy.CreateFolder_management(
                    os.path.join(lm_env.DATAPASSDIR, "nlcc", "nlc" + nlc_str),
                    "inv_norm")
            inv_norm_raster.save(
                os.path.join(lm_env.DATAPASSDIR, "nlcc", "nlc" + nlc_str,
                             "inv_norm", input_raster))
        nlc_idx += 1
    arcpy.env.workspace = prev_ws


def chk_weights():
    """Check weights."""
    if lm_env.CCERAST_IN:
        if (lm_env.CLOSEWEIGHT + lm_env.PERMWEIGHT + lm_env.CAVWEIGHT +
                lm_env.ECIVWEIGHT + lm_env.CEDWEIGHT != 1.0):
            lm_util.gprint("Warning: CLOSEWEIGHT + PERMWEIGHT + CAVWEIGHT + "
                           "ECIVWEIGHT + CEDWEIGHT <> 1.0")
    else:
        if (lm_env.CLOSEWEIGHT + lm_env.PERMWEIGHT +
                lm_env.CAVWEIGHT + lm_env.ECIVWEIGHT != 1.0):
            lm_util.gprint("Warning: CLOSEWEIGHT + PERMWEIGHT + CAVWEIGHT + "
                           "ECIVWEIGHT <> 1.0")
    if lm_env.CEDWEIGHT > 0 and not lm_env.CCERAST_IN:
        lm_util.gprint("Warning: CEDWEIGHT > 0, but no Current Climate "
                       "Envelope raster input provided")
    if lm_env.CEDWEIGHT == 0 and lm_env.CCERAST_IN:
        lm_util.gprint("Warning: Current Climate Envelope raster input "
                       "provided, but CEDWEIGHT = 0")
    if lm_env.ECIVWEIGHT > 0 and (
            (not lm_env.COREPAIRSTABLE_IN) or (not lm_env.ECIVFIELD)):
        lm_util.gprint("Warning: ECIVWEIGHT > 0, but no Expert Corridor "
                       "Importance Value field provided")
    if (lm_env.ECIVWEIGHT == 0 and lm_env.COREPAIRSTABLE_IN
            and lm_env.ECIVFIELD):
        lm_util.gprint("Warning: Expert Corridor Importance Value field "
                       "provided, but ECIVWEIGHT = 0")


def csp(sum_rasters, cnt_non_null_cells_rast, max_rasters, lcp_lines,
        core_lyr):
    """Calculate Corridor Specific Priority (CSP) for each corridor."""
    lm_util.gprint("Calculating Corridor Specific Priority (CSP) for each "
                   "corridor")
    chk_weights()

    # invert and normalize each corridor
    inv_norm()

    # normalize Expert Corridor Importance Value (ECIV)
    if lm_env.COREPAIRSTABLE_IN:
        eciv()

    # calc climate envelope and analog ratio
    if lm_env.CCERAST_IN:
        clim_linkage_priority(lcp_lines, core_lyr)

    prev_ws = arcpy.env.workspace

    # could be multiple folders
    nlc_idx = 0
    while True:
        nlc_str = ""
        if nlc_idx > 0:
            nlc_str = str(nlc_idx)
        if not os.path.exists(os.path.join(lm_env.DATAPASSDIR, "nlcc", "nlc" +
                                           nlc_str, "inv_norm")):
            break
        arcpy.env.workspace = (os.path.join(lm_env.DATAPASSDIR, "nlcc", "nlc" +
                                            nlc_str, "inv_norm"))

        csp_rasters, count_rasters = [], []
        # process each corridor raster in folder
        for in_rast in arcpy.ListRasters():
            # check for max 1
            result = arcpy.GetRasterProperties_management(in_rast, "MAXIMUM")
            max_in = float(result.getOutput(0))
            if max_in != 1.0:
                lm_util.gprint("Warning: maximum " + in_rast +
                               " value <> 1.0")

            # get cores from raster name
            from_core, to_core = in_rast.split("_")

            # check for corresponding link
            links = arcpy.SearchCursor(
                lcp_lines,
                where_clause="(From_Core = " + from_core + " AND To_Core = "
                + to_core + ") OR (From_Core = " + to_core +
                " AND To_Core = " + from_core + ")",
                fields="Rel_Close; Rel_Perm; clim_lnk_priority")

            if links:
                link = links.next()

                # get and avg CAVs for the core pair
                x_cav = arcpy.SearchCursor(
                    core_lyr,
                    where_clause=lm_env.COREFN + " = " + from_core,
                    fields="norm_cav").next().getValue("norm_cav")
                y_cav = arcpy.SearchCursor(
                    core_lyr,
                    where_clause=lm_env.COREFN + " = " + to_core,
                    fields="norm_cav").next().getValue("norm_cav")

                avg_cav = (x_cav + y_cav) / 2

                # calc weighted sum
                output_raster = (
                    (lm_env.CLOSEWEIGHT * link.getValue("Rel_Close")) +
                    (lm_env.PERMWEIGHT * link.getValue("Rel_Perm")) +
                    (0.0001 * arcpy.sa.Raster(in_rast)) +
                    (lm_env.CAVWEIGHT * avg_cav))

                # get ECIV for the core pair
                if lm_env.COREPAIRSTABLE_IN and lm_env.ECIVFIELD:
                    neciv = arcpy.SearchCursor(
                        lm_env.COREPAIRSTABLE_IN,
                        where_clause="(" + lm_env.FROMCOREFIELD + " = " +
                        from_core + " AND " + lm_env.TOCOREFIELD + " = " +
                        to_core + ") OR" + "(" + lm_env.TOCOREFIELD + " = " +
                        from_core + " AND " + lm_env.FROMCOREFIELD + " = " +
                        to_core + ")").next().getValue("neciv")
                    output_raster += lm_env.ECIVWEIGHT * neciv

                # Increment weighted sum with Climate Gradient
                if lm_env.CCERAST_IN:
                    output_raster += (link.getValue("Clim_Lnk_Priority") *
                                      lm_env.CEDWEIGHT)

                save_interm_rast(output_raster, '_'.join(["CSP", in_rast]))

                # get max and min
                lm_util.build_stats(output_raster)
                result = arcpy.GetRasterProperties_management(output_raster,
                                                              "MAXIMUM")
                max_csp = float(result.getOutput(0))
                result = arcpy.GetRasterProperties_management(output_raster,
                                                              "MINIMUM")
                min_csp = float(result.getOutput(0))

                # determine threshold value
                diff_csp = max_csp - min_csp
                thres_csp = max_csp - (diff_csp * lm_env.PROPCSPKEEP)

                # apply threshold
                con_raster = arcpy.sa.Con(output_raster, output_raster, "#",
                                          "VALUE > " + str(thres_csp))
                is_null_raster = arcpy.sa.IsNull(con_raster)
                count_raster = arcpy.sa.EqualTo(is_null_raster, 0)
                count_rasters.append(count_raster)
                csp_rasters.append(con_raster)

                save_interm_rast(con_raster, '_'.join(["CSP_TOP", in_rast]))

            del link, links

        # perform intermediate calculations on CSPs leading toward CPV
        sum_rasters.append(arcpy.sa.CellStatistics(csp_rasters, "SUM", "DATA"))
        cnt_non_null_cells_rast.append(
            arcpy.sa.CellStatistics(count_rasters, "SUM", "DATA"))
        max_rasters.append(
            arcpy.sa.CellStatistics(csp_rasters, "MAXIMUM", "DATA"))
        nlc_idx += 1

    arcpy.env.workspace = prev_ws


def cav(core_lyr):
    """Calculate Core Area Value (CAV) and its components for each core."""
    lm_util.gprint("Calculating Core Area Value (CAV) and its components for "
                   "each core")
    # check weights and warn if issues
    if lm_env.OCAVRAST_IN:
        if (lm_env.RESWEIGHT + lm_env.SIZEWEIGHT + lm_env.APWEIGHT +
                lm_env.ECAVWEIGHT + lm_env.CFCWEIGHT + lm_env.OCAVWEIGHT
                != 1.0):
            lm_util.gprint("Warning: RESWEIGHT + SIZEWEIGHT + APWEIGHT + "
                           "ECAVWEIGHT + CFCWEIGHT + OCAVWEIGHT <> 1.0")
    else:
        if (lm_env.RESWEIGHT + lm_env.SIZEWEIGHT + lm_env.APWEIGHT +
                lm_env.ECAVWEIGHT + lm_env.CFCWEIGHT != 1.0):
            lm_util.gprint("Warning: RESWEIGHT + SIZEWEIGHT + APWEIGHT + "
                           "ECAVWEIGHT + CFCWEIGHT <> 1.0")
    if lm_env.OCAVWEIGHT > 0 and not lm_env.OCAVRAST_IN:
        lm_util.gprint("Warning: OCAVWEIGHT > 0 but no OCAV raster input "
                       "provided")
    if lm_env.OCAVWEIGHT == 0 and lm_env.OCAVRAST_IN:
        lm_util.gprint("Warning: OCAV raster input provided, but "
                       "OCAVWEIGHT = 0")

    # check/add fields
    for field in ("mean_res", "norm_res", "area", "norm_size", "perimeter",
                  "ap_ratio", "norm_ratio", "cav", "norm_cav", "cclim_env",
                  "fclim_env", "ocav", "nocav"):
        check_add_field(lm_env.COREFC, field, "DOUBLE")

    if not check_add_field(lm_env.COREFC, "ecav", "DOUBLE"):
        if lm_env.ECAVWEIGHT > 0:
            lm_util.gprint("Warning: ECAVWEIGHT > 0 but no ecav field in  "
                           "Cores feature class")
        arcpy.CalculateField_management(lm_env.COREFC, "ecav", "0")
    check_add_field(lm_env.COREFC, "necav", "DOUBLE")

    # current flow centrality (CFC, CF_Central) is copied from
    # Centrality Mapper
    if not check_add_field(lm_env.COREFC, "CF_Central", "DOUBLE"):
        # default to 0s
        arcpy.CalculateField_management(lm_env.COREFC, "CF_Central", "0")
    if lm_env.CFCWEIGHT > 0:
        # copy values from Centrality Mapper output
        # (core_centrality.gdb.project_Cores) if available
        centrality_cores = os.path.join(lm_env.CORECENTRALITYGDB,
                                        lm_env.PREFIX + "_Cores")
        if arcpy.Exists(centrality_cores):
            arcpy.AddJoin_management(core_lyr, lm_env.COREFN,
                                     centrality_cores, lm_env.COREFN)
            arcpy.CalculateField_management(
                core_lyr, lm_env.CORENAME + ".CF_Central", "[" +
                lm_env.PREFIX + "_Cores.CF_Central]")
            arcpy.RemoveJoin_management(core_lyr)
        # ensure cores have at least one non-0 value for CFC (could have been
        # copied above or set earlier)
        max_val = value_range(lm_env.COREFC, "CF_Central")[1]
        if max_val is None or max_val == 0:
            raise Exception(
                "ERROR: A Current Flow Centrality Weight (CFCWEIGHT) was "
                "provided but no Current Flow Centrality (CF_Central) "
                "values are available. Please run Centrality Mapper on "
                "this project, then run Linkage Priority.")

    check_add_field(lm_env.COREFC, "ncfc", "DOUBLE")

    # calc mean resistance
    core_mean(lm_env.RESRAST_IN, core_lyr, "mean_res")

    # calc area, perimeter and ratio
    arcpy.CalculateField_management(core_lyr, "area", "!SHAPE.AREA!",
                                    "PYTHON_9.3")
    arcpy.CalculateField_management(core_lyr, "perimeter", "!SHAPE.LENGTH!",
                                    "PYTHON_9.3")
    arcpy.CalculateField_management(core_lyr, "ap_ratio",
                                    "!area! / !perimeter!", "PYTHON_9.3")

    # normalize CAV inputs
    # resistance - invert
    normalize_field(core_lyr, "mean_res", "norm_res", lm_env.RESNORMETH,
                    True)

    # size
    normalize_field(core_lyr, "area", "norm_size", lm_env.SIZENORMETH)

    # area/perimeter ratio
    normalize_field(core_lyr, "ap_ratio", "norm_ratio", lm_env.APNORMETH)

    # ecav
    normalize_field(core_lyr, "ecav", "necav", lm_env.ECAVNORMETH)

    # cfc
    normalize_field(core_lyr, "CF_Central", "ncfc", lm_env.CFCNORMETH)

    # calc OCAV
    if lm_env.OCAVRAST_IN:
        # get max and min
        lm_util.build_stats(lm_env.OCAVRAST_IN)
        result = arcpy.GetRasterProperties_management(lm_env.OCAVRAST_IN,
                                                      "MAXIMUM")
        max_ocav = float(result.getOutput(0))
        result = arcpy.GetRasterProperties_management(lm_env.OCAVRAST_IN,
                                                      "MINIMUM")
        min_ocav = float(result.getOutput(0))
        # calc score range normalization on input
        ocav_raster = ((arcpy.sa.Raster(lm_env.OCAVRAST_IN) - min_ocav)
                       / (max_ocav - min_ocav))
        # calc aerial mean ocav for each core
        core_mean(ocav_raster, core_lyr, "ocav")
        normalize_field(core_lyr, "ocav", "nocav", NM_SCORE)

        # calc CAV
        arcpy.CalculateField_management(
            core_lyr, "cav",
            "(!norm_res! * " + str(lm_env.RESWEIGHT) + ") + (!norm_size! * " +
            str(lm_env.SIZEWEIGHT) + ") + (!norm_ratio! * " +
            str(lm_env.APWEIGHT) + ") + (!necav! * " + str(lm_env.ECAVWEIGHT)
            + ") + (!ncfc! * " + str(lm_env.CFCWEIGHT) + ") + (!nocav! * " +
            str(lm_env.OCAVWEIGHT) + ")", "PYTHON_9.3")

    else:
        # calc CAV
        arcpy.CalculateField_management(
            core_lyr, "cav", "(!norm_res! * " + str(lm_env.RESWEIGHT) +
            ") + (!norm_size! * " + str(lm_env.SIZEWEIGHT) +
            ") + (!norm_ratio! * " + str(lm_env.APWEIGHT) +
            ") + (!necav! * " + str(lm_env.ECAVWEIGHT) + ") + (!ncfc! * " +
            str(lm_env.CFCWEIGHT) + ")", "PYTHON_9.3")

    normalize_field(core_lyr, "cav", "norm_cav", NM_SCORE)


def calc_closeness(lcp_lines):
    """Calculate relative closeness for each Least Cost Path."""
    lm_util.gprint("Calculating relative closeness for each LCP line")
    normalize_field(lcp_lines, "LCP_Length", "Rel_Close",
                    lm_env.RELCLOSENORMETH, True)


def calc_permeability(lcp_lines):
    """Calculate raw and relative permeability for each Least Cost Path."""
    # raw
    lm_util.gprint("Calculating raw permeability for each LCP line")
    check_add_field(lcp_lines, "Raw_Perm", "DOUBLE")
    arcpy.CalculateField_management(lcp_lines, "Raw_Perm",
                                    "!LCP_Length! / !CW_Dist!", "PYTHON_9.3")

    # relative
    lm_util.gprint("Calculating relative permeability for each LCP line")
    normalize_field(lcp_lines, "Raw_Perm", "Rel_Perm", lm_env.RELPERMNORMETH)


def add_output_path(in_str):
    """Append LinkMap GDB path to inputted value."""
    return os.path.join(lm_env.OUTPUTGDB,
                        '_'.join([lm_env.PREFIX, in_str]))


def create_run_gdbs():
    """Create scratch and if necessary intermediate GDB."""
    if not os.path.isdir(lm_env.SCRATCHDIR):
        os.makedirs(lm_env.SCRATCHDIR)

    if not arcpy.Exists(lm_env.SCRATCHGDB):
        arcpy.CreateFileGDB_management(
            lm_env.SCRATCHDIR, os.path.basename(lm_env.SCRATCHGDB))
    arcpy.env.scratchWorkspace = lm_env.SCRATCHGDB

    if lm_env.KEEPINTERMEDIATE:
        if not arcpy.Exists(lm_env.INTERGDB):
            arcpy.CreateFileGDB_management(
                lm_env.SCRATCHDIR, os.path.basename(lm_env.SCRATCHDIR))


def chk_lnk_tbls():
    """Check that LM finished with steps 3 and 5."""
    if (not os.path.isfile(os.path.join(lm_env.DATAPASSDIR,
                                        "linkTable_s3.csv"))
            or not os.path.isfile(os.path.join(lm_env.DATAPASSDIR,
                                               "linkTable_s5.csv"))):
        raise Exception("ERROR: Project directory must contain a successful "
                        "Linkage Mapper run with Steps 3 and 5.")


def run_analysis():
    """Run main Linkage Priority analysis."""
    lm_util.gprint("Checking inputs")

    chk_lnk_tbls()
    create_run_gdbs()

    # calc permeability
    lcp_lines = os.path.join(lm_env.LINKMAPGDB, lm_env.PREFIX + "_LCPs")
    calc_permeability(lcp_lines)

    # calc relative closeness
    calc_closeness(lcp_lines)

    # calc Core Area Value (CAV) and its components for each core
    core_lyr = arcpy.MakeFeatureLayer_management(lm_env.COREFC, "core_lyr")
    cav(core_lyr)

    # calc Corridor Specific Priority (CSP)
    prev_ws = arcpy.env.workspace
    sum_rasters, cnt_non_null_cells_rast, max_rasters = [], [], []
    csp(sum_rasters, cnt_non_null_cells_rast, max_rasters, lcp_lines, core_lyr)

    # calc Corridor Priority Value (CPV)
    cpv_raster = add_output_path("CPV")
    cpv(sum_rasters, cnt_non_null_cells_rast, max_rasters, cpv_raster)
    arcpy.env.workspace = prev_ws

    # calc Relative Corridor Importance (RCI)
    rci_raster = add_output_path("RCI")
    rci(cpv_raster, rci_raster)

    if lm_env.CALCLP:
        trunc_raster = add_output_path(
            '_'.join(["corridors_truncated_at", lm_util.cwd_threash_str()]))
        if not arcpy.Exists(trunc_raster):
            lm_util.raise_error(
                "Truncated corridors raster not found.\n"
                "When CALCLP = True, set WRITETRUNCRASTER = True in "
                "Linkage Mapper")

        # calc Linkage Priority (LP)
        lp_raster = add_output_path("linkage_priority")
        linkage_priority(rci_raster, trunc_raster, lp_raster)

        # calc Blended Priority (BP)
        if lm_env.CALCBP:
            norm_trunc_raster = add_output_path("NORMTRUNC")
            bp_raster = add_output_path("blended_priority")
            norm_trunc(trunc_raster, norm_trunc_raster)
            blended_priority(norm_trunc_raster, lp_raster, bp_raster)
    else:
        if lm_env.CALCBP:
            lm_util.raise_error("When CALCBP = True, set CALCLP = True")

    # save a copy of Cores as the "Output for ModelBuilder Precondition"
    if lm_env.OUTPUTFORMODELBUILDER:
        arcpy.CopyFeatures_management(lm_env.COREFC,
                                      lm_env.OUTPUTFORMODELBUILDER)


def read_lm_params(proj_dir):
    """Read Linkage Pathways input parameters from log file."""
    # get log file for last LM run
    spath = os.path.join(proj_dir, "run_history", "log",
                         "*_Linkage Mapper.txt")
    entries = sorted(glob.glob(spath), key=os.path.getctime, reverse=True)
    last_lm_log = next(iter(entries or []), None)
    if not last_lm_log:
        raise Exception("ERROR: Log file for last Linkage Mapper run not "
                        "found. Please ensure Linkage Mapper is run "
                        "for this project before running Linkage Priority.")

    # read parameters section from file and turn into tuple for passing
    parms = ""
    with open(last_lm_log) as log_file:
        for line in log_file:
            if line[0:13] == "Parameters:\t[":
                parms = line[13:len(line) - 3].replace("\\\\", "\\")
                break
    if parms == "":
        raise Exception("ERROR: Log file for last Linkage Mapper run does not "
                        "contain a Parameters line")

    return tuple(parms.replace("'", "").split(", "))


def get_lm_params(argv):
    """Get settings from Linkage Pathways inputs."""
    lm_params = read_lm_params(argv[1])  # Pass in project dir
    argv.append(lm_params[17])  # Get CWDTHRESH


def log_setup():
    """Set up Linkage Mapper logging."""
    lm_env.logFilePath = lm_util.create_log_file(lm_env.MESSAGEDIR,
                                                 lm_env.TOOL, lm_env.PARAMS)
    lm_util.write_custom_to_log(lm_env.LPCUSTSETTINGS_IN)


def main(argv=None):
    """Run Linkage Priority tool."""
    start_time = dt.now()
    print "Start time: %s" % start_time.strftime(TFORMAT)

    if argv is None:
        # get parameters passed from tool dialog, if any
        argv = sys.argv
    try:
        # preparation steps
        get_lm_params(argv)
        lm_env.configure(lm_env.TOOL_LP, argv)
        lm_util.gprint("\nLinkage Priority Version " + lm_env.releaseNum)
        lm_util.check_project_dir()
        log_setup()
        # primary analysis
        run_analysis()
    except arcpy.ExecuteError:
        # arcpy error
        msg = arcpy.GetMessages(2)
        arcpy.AddError(msg)
        lm_util.write_log(msg)
        exc_traceback = sys.exc_info()[2]
        lm_util.gprint(
            "Traceback (most recent call last):\n" +
            "".join(traceback.format_tb(exc_traceback)[:-1]))
    except Exception:
        # other exceptions
        exc_value, exc_traceback = sys.exc_info()[1:]
        arcpy.AddError(exc_value)
        lm_util.gprint(
            "Traceback (most recent call last):\n" +
            "".join(traceback.format_tb(exc_traceback)))
    finally:
        delete_arcpy_temp_datasets()
        arcpy.CheckInExtension("Spatial")
        print_runtime(start_time)
        lm_util.close_log_file()


if __name__ == "__main__":
    main()
