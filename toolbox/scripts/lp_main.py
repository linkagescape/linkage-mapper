# Authors: John Gallo and Randal Greene 2017

"""Linkage Priority main module."""

import os
import sys
import glob
import traceback
from collections import namedtuple

import arcpy

from lm_config import tool_env as lm_env
import lm_util


_SCRIPT_NAME = "lp_main.py"

NM_SCORE = "SCORE_RANGE"  # Score range normalization
NM_MAX = "MAX_VALUE"  # Maximum value normalization

CoordPoint = namedtuple('Point', 'x y')


class AppError(Exception):
    """Custom error class."""

    pass


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


def save_interm_rast(rast_list, base_name):
    """Save intermediate rasters if user chooses."""
    if lm_env.KEEPINTERMEDIATE:
        gdb_name = ".".join([base_name, "gdb"])
        arcpy.CreateFileGDB_management(lm_env.SCRATCHDIR, gdb_name)
        for rast_name, rast in rast_list:
            rast.save(os.path.join(lm_env.SCRATCHDIR, gdb_name, rast_name))


def blended_priority(rast_list, lcp_ncsp):
    """Calculate overall Blended Priority."""
    lm_util.gprint("-Calculating overall BP")

    blend_rast = []
    for fname, rast in rast_list:
        blend_rast.append([fname, (lm_env.TRUNCWEIGHT * rast) +
                           (lm_env.LPWEIGHT * lcp_ncsp[fname])])

    save_interm_rast(blend_rast, "bp_step4")

    overall_bp = arcpy.sa.CellStatistics(
        [rast[1] for rast in blend_rast],
        statistics_type="MAXIMUM", ignore_nodata="DATA")
    overall_bp.save(os.path.join(lm_env.OUTPUTGDB, "blended_priority"))


def inv_norm(rast_list):
    """Invert and normalize each corridor."""
    lm_util.gprint("-Inverting and normalizing each corridor")
    norm_rast_list = []

    for fname, rast in rast_list:
        norm_rast = normalize_raster(rast, lm_env.NORMCORRNORMETH, True)
        norm_rast_list.append([fname, norm_rast])

    save_interm_rast(norm_rast_list, "bp_step3")
    return norm_rast_list


def clip_nlcc_to_threashold(lcp_list):
    """Clip NLCC_A_B rasters to CWD threshold.

    Clip the normalized least cost corridors using the specified CWD
    Threshold.
    """
    lm_util.gprint("-Clipping NLCC rasters to CWD threshold")

    prev_workspace = arcpy.env.workspace
    nlcc_top_list = []

    for dirpath, _, filenames in arcpy.da.Walk(
            lm_env.LCCBASEDIR, topdown=True, datatype="RasterDataset"):
        arcpy.env.workspace = dirpath
        for filename in filenames:
            if lm_env.LCCNLCDIR_NM in dirpath and filename in lcp_list:
                rast = arcpy.sa.ExtractByAttributes(
                    filename,
                    "VALUE <= {}".format(lm_env.CWDTHRESH))

                # Raster names cannot begin with a number in a GDB
                nfilename = '_'.join(["nlc", filename])

                nlcc_top_list.append([nfilename, rast])

    if not nlcc_top_list:
        raise AppError("No normalized least cost corridor rasters found")

    arcpy.env.workspace = prev_workspace
    save_interm_rast(nlcc_top_list, "bp_step1")

    return nlcc_top_list


def lcp_csp_for_bp(lcp_lines):
    """Get list of LCPs and their normalized CSP values for BP raster creation.

    Filter LCPs based on CPS cutoff if applicable. Return list with LCP
    filenames and dictionary with LCP GDB filename as key and normalized CSP
    as value.
    """
    lcp_list = []
    lcp_ncsp = {}

    if lm_env.CPSNORM_CUTOFF is not None:
        cursor_filter = "CSP_Norm_Trim=1"
    else:
        cursor_filter = None
    lcp_rows = arcpy.SearchCursor(
        lcp_lines,
        where_clause=cursor_filter,
        fields="From_Core; To_Core; CSP_Norm")

    for lcp_row in lcp_rows:
        lcp_name = ("{}_{}".format(lcp_row.getValue("From_Core"),
                                   lcp_row.getValue("To_Core")))
        lcp_list.append(lcp_name)

        # Use raster GDB name for key value
        lcp_ncsp['_'.join(["nlc", lcp_name])] = lcp_row.getValue("CSP_Norm")

    return lcp_list, lcp_ncsp


def calc_blended_priority(lcp_lines):
    """Generate Blended Priority raster from NLCC rasters."""
    lm_util.gprint("Calculating Blended Priority (BP):")

    lcp_list, lcp_ncsp = lcp_csp_for_bp(lcp_lines)
    nlcc_rast = clip_nlcc_to_threashold(lcp_list)
    nlcc_rast = inv_norm(nlcc_rast)

    blended_priority(nlcc_rast, lcp_ncsp)


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
        except arcpy.ExecuteError:
            lm_util.gprint(
                "ERROR! MOST LIKELY CAUSE: One or more core areas are smaller "
                "than a pixel in the Resistance layer and/or Raster Analysis "
                "Cell Size environment setting. Try enlarging small core "
                "areas, resampling the Resistance layer or adjusting the Cell"
                "Size environment setting.")
            raise
    else:
        arcpy.CalculateField_management(in_table, out_field, "0", "PYTHON_9.3")


def clim_priority_val_normal(lcp_lines):
    """Normalize climate priority values (A & L) for each core pair."""
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
    check_add_field(lcp_lines, "CLPv_Analog", "float")
    check_add_field(lcp_lines, "CLPv_Prefer", "float")

    canalog_r_min, canalog_r_max = value_range(lcp_lines, "CAnalog_Ratio")
    cprefer_r_min, cprefer_r_max = value_range(lcp_lines, "CPrefer_Ratio")

    if lm_env.CANALOG_MINRMAX > canalog_r_max:
        canalog_r_max = lm_env.CANALOG_MINRMAX

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
    mean_value = "".join(["!", tbl_name, ".MEAN!"])

    mean_tbl = arcpy.sa.ZonalStatisticsAsTable(
        lm_env.COREFC, lm_env.COREFN, in_rast,
        os.path.join(lm_env.SCRATCHGDB, tbl_name),
        statistics_type="MEAN")
    arcpy.AddJoin_management(core_lyr, lm_env.COREFN, mean_tbl,
                             lm_env.COREFN)
    arcpy.CalculateField_management(core_lyr, mean_fld, mean_value,
                                    "PYTHON_9.3")
    arcpy.RemoveJoin_management(core_lyr)
    lm_util.delete_data(mean_tbl)


def clim_envelope(core_lyr):
    """Determine Climate Envelope for each core."""
    core_mean(lm_env.CCERAST_IN, core_lyr, "cclim_env")
    if lm_env.FCERAST_IN:
        core_mean(lm_env.FCERAST_IN, core_lyr, "fclim_env")


def clim_linkage_priority(lcp_lines, core_lyr):
    """Calculate Core Areas Climate Linkage Priority Value."""
    lm_util.gprint("-Calculating Core Areas Climate Linkage Priority Value")
    clim_envelope(core_lyr)
    clim_ratios(lcp_lines, core_lyr)
    clim_priority_values(lcp_lines)
    clim_priority_val_normal(lcp_lines)
    clim_priority_combine(lcp_lines)


def eciv():
    """Normalize Expert Corridor Importance Value (ECIV) for each corridor."""
    lm_util.gprint("-Normalizing Expert Corridor Importance Value (ECIV) "
                   "for each corridor")
    normalize_field(lm_env.COREPAIRSTABLE_IN, lm_env.ECIVFIELD, "neciv",
                    NM_SCORE)


def chk_csp_wts():
    """Check weights used in CSP calculation."""
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


def calc_csp(lcp_lines, core_lyr):
    """Calculate Corridor Specific Priority (CSP) for each linkage."""
    lm_util.gprint("Calculating Corridor Specific Priority (CSP) for each "
                   "linkage:")
    chk_csp_wts()

    # Normalize Expert Corridor Importance Value (ECIV)
    if lm_env.COREPAIRSTABLE_IN:
        eciv()

    # Calc climate envelope and analog ratio
    if lm_env.CCERAST_IN:
        clim_linkage_priority(lcp_lines, core_lyr)
        lnk_fields = ("From_Core; To_Core; Rel_Close; Rel_Perm; "
                      "clim_lnk_priority; CSP")
    else:
        lnk_fields = "From_Core; To_Core; Rel_Close; Rel_Perm; CSP"

    lm_util.gprint("-Calculating CSP")
    check_add_field(lcp_lines, "CSP", "float")

    lnk_rows = arcpy.UpdateCursor(
        lcp_lines,
        fields=lnk_fields)

    for lnk_row in lnk_rows:
        from_core = lnk_row.getValue("From_Core")
        to_core = lnk_row.getValue("To_Core")

        # Get and avg CAVs for the core pair
        x_cav = arcpy.SearchCursor(
            core_lyr,
            where_clause="{} = {}".format(lm_env.COREFN, from_core),
            fields="norm_cav").next().getValue("norm_cav")
        y_cav = arcpy.SearchCursor(
            core_lyr,
            where_clause="{} = {}".format(lm_env.COREFN, to_core),
            fields="norm_cav").next().getValue("norm_cav")

        avg_cav = (x_cav + y_cav) / 2

        # Calc weighted sum
        csp = (
            (lm_env.CLOSEWEIGHT * lnk_row.getValue("Rel_Close")) +
            (lm_env.PERMWEIGHT * lnk_row.getValue("Rel_Perm")) +
            (lm_env.CAVWEIGHT * avg_cav))

        # Get ECIV for the core pair
        if lm_env.COREPAIRSTABLE_IN and lm_env.ECIVFIELD:
            neciv = arcpy.SearchCursor(
                lm_env.COREPAIRSTABLE_IN,
                where_clause="({0}={2} AND {1}={3}) "
                "OR ({1}={2} AND {0}={3})".format(
                    lm_env.FROMCOREFIELD, lm_env.TOCOREFIELD,
                    from_core, to_core),
                fields="neciv").next().getValue("neciv")
            csp += lm_env.ECIVWEIGHT * neciv

        # Increment weighted sum with Climate Gradient
        if lm_env.CCERAST_IN:
            csp += (lnk_row.getValue("Clim_Lnk_Priority") *
                    lm_env.CEDWEIGHT)
        lnk_row.setValue("CSP", csp)
        lnk_rows.updateRow(lnk_row)
    del lnk_rows

    # Normalize CSP values
    cps_norm_fld = "CSP_Norm"
    normalize_field(lcp_lines, "CSP", cps_norm_fld)

    # If user requests, flag low quality corridors not to use
    if lm_env.CPSNORM_CUTOFF:
        csp_trim_fld = "CSP_Norm_Trim"
        expression = "set_bool(!{}!, {})".format(cps_norm_fld,
                                                 lm_env.CPSNORM_CUTOFF)
        codeblock = (
            """def set_bool(cps_norm, cps_cutoff):
                   if cps_norm <= cps_cutoff:
                       return 0
                   else:
                       return 1""")

        check_add_field(lcp_lines, csp_trim_fld, "SHORT")
        arcpy.CalculateField_management(lcp_lines, csp_trim_fld, expression,
                                        "PYTHON_9.3", codeblock)


def chk_cav_wts():
    """Check weights used in CAV calculation."""
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


def calc_cav(core_lyr):
    """Calculate Core Area Value (CAV) and its components for each core."""
    lm_util.gprint("Calculating Core Area Value (CAV) and its components for "
                   "each core")
    chk_cav_wts()

    # Check/add fields
    for field in ("mean_res", "norm_res", "area", "norm_size", "perimeter",
                  "ap_ratio", "norm_ratio", "cav", "norm_cav", "cclim_env",
                  "fclim_env", "ocav", "nocav"):
        check_add_field(lm_env.COREFC, field, "DOUBLE")

    if not check_add_field(lm_env.COREFC, "ecav", "DOUBLE"):
        if lm_env.ECAVWEIGHT > 0:
            lm_util.gprint("Warning: ECAVWEIGHT > 0 but no ecav field in  "
                           "Cores feature class")
        arcpy.CalculateField_management(lm_env.COREFC, "ecav", "0",
                                       "PYTHON_9.3")
    check_add_field(lm_env.COREFC, "necav", "DOUBLE")

    # Current flow centrality (CFC, CF_Central) is copied from
    # Centrality Mapper
    if not check_add_field(lm_env.COREFC, "CF_Central", "DOUBLE"):
        # Default to 0s
        arcpy.CalculateField_management(lm_env.COREFC, "CF_Central", "0",
                                        "PYTHON_9.3")
    if lm_env.CFCWEIGHT > 0:
        # Copy values from Centrality Mapper output
        # (core_centrality.gdb.project_Cores) if available
        centrality_cores = os.path.join(lm_env.CORECENTRALITYGDB,
                                        lm_env.PREFIX + "_Cores")
        if arcpy.Exists(centrality_cores):
            arcpy.AddJoin_management(core_lyr, lm_env.COREFN,
                                     centrality_cores, lm_env.COREFN)
            arcpy.CalculateField_management(
                core_lyr, lm_env.CORENAME + ".CF_Central",
                "!" + lm_env.PREFIX + "_Cores.CF_Central!", "PYTHON_9.3")
            arcpy.RemoveJoin_management(core_lyr)
        # Ensure cores have at least one non-0 value for CFC (could have been
        # copied above or set earlier)
        max_val = value_range(lm_env.COREFC, "CF_Central")[1]
        if max_val is None or max_val == 0:
            raise AppError(
                "ERROR: A Current Flow Centrality Weight (CFCWEIGHT) was "
                "provided but no Current Flow Centrality (CF_Central) "
                "values are available. Please run Centrality Mapper on "
                "this project, then run Linkage Priority.")

    check_add_field(lm_env.COREFC, "ncfc", "DOUBLE")

    # Calc mean resistance
    core_mean(lm_env.RESRAST_IN, core_lyr, "mean_res")

    # Calc area, perimeter and ratio
    arcpy.CalculateField_management(core_lyr, "area", "!SHAPE.AREA!",
                                    "PYTHON_9.3")
    arcpy.CalculateField_management(core_lyr, "perimeter", "!SHAPE.LENGTH!",
                                    "PYTHON_9.3")
    arcpy.CalculateField_management(core_lyr, "ap_ratio",
                                    "!area! / !perimeter!", "PYTHON_9.3")

    # Normalize CAV inputs
    normalize_field(core_lyr, "mean_res", "norm_res", lm_env.RESNORMETH,
                    True)
    normalize_field(core_lyr, "area", "norm_size", lm_env.SIZENORMETH)
    normalize_field(core_lyr, "ap_ratio", "norm_ratio", lm_env.APNORMETH)
    normalize_field(core_lyr, "ecav", "necav", lm_env.ECAVNORMETH)
    normalize_field(core_lyr, "CF_Central", "ncfc", lm_env.CFCNORMETH)

    # Calc OCAV
    if lm_env.OCAVRAST_IN:
        # Get max and min
        lm_util.build_stats(lm_env.OCAVRAST_IN)
        result = arcpy.GetRasterProperties_management(lm_env.OCAVRAST_IN,
                                                      "MAXIMUM")
        max_ocav = float(result.getOutput(0))
        result = arcpy.GetRasterProperties_management(lm_env.OCAVRAST_IN,
                                                      "MINIMUM")
        min_ocav = float(result.getOutput(0))
        # Calc score range normalization on input
        ocav_raster = ((arcpy.sa.Raster(lm_env.OCAVRAST_IN) - min_ocav)
                       / (max_ocav - min_ocav))
        # Calc aerial mean ocav for each core
        core_mean(ocav_raster, core_lyr, "ocav")
        normalize_field(core_lyr, "ocav", "nocav", NM_SCORE)

        # Calc CAV
        arcpy.CalculateField_management(
            core_lyr, "cav",
            "(!norm_res! * " + str(lm_env.RESWEIGHT) + ") + (!norm_size! * " +
            str(lm_env.SIZEWEIGHT) + ") + (!norm_ratio! * " +
            str(lm_env.APWEIGHT) + ") + (!necav! * " + str(lm_env.ECAVWEIGHT)
            + ") + (!ncfc! * " + str(lm_env.CFCWEIGHT) + ") + (!nocav! * " +
            str(lm_env.OCAVWEIGHT) + ")", "PYTHON_9.3")

    else:
        # Calc CAV
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
    lm_util.gprint("Calculating permeability for each LCP line")
    normalize_field(lcp_lines, "cwd_to_Path_Length_Ratio", "Rel_Perm",
                    lm_env.RELPERMNORMETH, True)


def add_output_path(in_str):
    """Append LinkMap GDB path to inputted value."""
    return os.path.join(lm_env.OUTPUTGDB,
                        '_'.join([lm_env.PREFIX, in_str]))


def make_core_lyr():
    """Create feature layer from cores feature class.

    Raise error if core feature class is not found.
    """
    if not arcpy.Exists(lm_env.COREFC):
        raise AppError("Cores feature class not found. Please rerun "
                       "Linkage Pathways tool")
    return arcpy.MakeFeatureLayer_management(lm_env.COREFC, "core_lyr")


def get_lcp_fc():
    """Get LCP feature class. Raise error if not found."""
    lcp_lines = os.path.join(lm_env.LINKMAPGDB, lm_env.PREFIX + "_LCPs")
    if not arcpy.Exists(lcp_lines):
        raise AppError("LCP feature class not found. Please rerun "
                       "Linkage Pathways tool")
    if arcpy.GetCount_management(lcp_lines) == 0:
        raise AppError("LCP feature class contains no records")
    return lcp_lines


def create_run_gdbs():
    """Create scratch and if necessary intermediate GDB."""
    if not os.path.isdir(lm_env.SCRATCHDIR):
        os.makedirs(lm_env.SCRATCHDIR)

    if not arcpy.Exists(lm_env.SCRATCHGDB):
        arcpy.CreateFileGDB_management(
            lm_env.SCRATCHDIR, os.path.basename(lm_env.SCRATCHGDB))


def chk_lnk_tbls():
    """Check that LM finished with steps 3 and 5."""
    if (not os.path.isfile(os.path.join(lm_env.DATAPASSDIR,
                                        "linkTable_s3.csv"))
            or not os.path.isfile(os.path.join(lm_env.DATAPASSDIR,
                                               "linkTable_s5.csv"))):
        raise AppError("ERROR: Project directory must contain a successful "
                       "Linkage Mapper run with Steps 3 and 5.")


def run_analysis():
    """Run main Linkage Priority analysis."""
    lm_util.gprint("Retreiving outputs from Linkage Pathways model run""")
    lcp_lines = get_lcp_fc()
    core_lyr = make_core_lyr()

    chk_lnk_tbls()
    create_run_gdbs()

    calc_permeability(lcp_lines)
    calc_closeness(lcp_lines)
    calc_cav(core_lyr)

    # Calculate Corridor Specific Value and Blended Priority raster
    if lm_env.CALCCSPBP in ([lm_env.CALC_CSP, lm_env.CALC_CSPBP]):
        calc_csp(lcp_lines, core_lyr)
        if lm_env.CALCCSPBP == lm_env.CALC_CSPBP:
            calc_blended_priority(lcp_lines)

    # Save a copy of Cores as the "Output for ModelBuilder Precondition"
    if lm_env.OUTPUTFORMODELBUILDER:
        arcpy.CopyFeatures_management(lm_env.COREFC,
                                      lm_env.OUTPUTFORMODELBUILDER)


def read_lm_params(proj_dir):
    """Read Linkage Pathways input parameters from log file."""
    # Get log file for last LM run
    spath = os.path.join(proj_dir, "run_history", "log",
                         "*_Linkage Mapper.txt")
    entries = sorted(glob.glob(spath), key=os.path.getctime, reverse=True)
    last_lm_log = next(iter(entries or []), None)
    if not last_lm_log:
        raise AppError("ERROR: Log file for last Linkage Mapper run not "
                       "found. Please ensure Linkage Mapper is run "
                       "for this project before running Linkage Priority.")

    # Read parameters section from file and turn into tuple for passing
    parms = ""
    with open(last_lm_log) as log_file:
        for line in log_file:
            if line[0:13] == "Parameters:\t[":
                parms = line[13:len(line) - 3].replace("\\\\", "\\")
                break
    if parms == "":
        raise AppError("ERROR: Log file for last Linkage Mapper run does not "
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
    stime = lm_util.start_time()

    if argv is None:
        argv = sys.argv  # Get parameters from ArcGIS tool dialog
    try:
        get_lm_params(argv)
        lm_env.configure(lm_env.TOOL_LP, argv)
        lm_util.gprint("\nLinkage Priority Version " + lm_env.releaseNum)
        lm_util.check_project_dir()
        log_setup()
        run_analysis()
    except AppError as err:
        lm_util.gprint(err.message)
        lm_util.write_log(err.message)
    except arcpy.ExecuteError:
        msg = arcpy.GetMessages(2)
        arcpy.AddError(msg)
        lm_util.write_log(msg)
        exc_traceback = sys.exc_info()[2]
        lm_util.gprint(
            "Traceback (most recent call last):\n" +
            "".join(traceback.format_tb(exc_traceback)[:-1]))
    except Exception:
        exc_value, exc_traceback = sys.exc_info()[1:]
        arcpy.AddError(exc_value)
        lm_util.gprint(
            "Traceback (most recent call last):\n" +
            "".join(traceback.format_tb(exc_traceback)))
    finally:
        if not lm_env.KEEPINTERMEDIATE:
            lm_util.delete_dir(lm_env.SCRATCHDIR)
        arcpy.CheckInExtension("Spatial")
        lm_util.run_time(stime)
        lm_util.close_log_file()


if __name__ == "__main__":
    main()
