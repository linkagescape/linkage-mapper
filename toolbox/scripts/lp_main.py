"""Linkage Priority main module."""
# Authors: John Gallo and Randal Greene 2017


import os
import sys
import traceback
from datetime import datetime

import arcpy
from arcpy.sa import *

from lp_config import lp_env
from lm_config import tool_env as lm_env
import lm_util


_SCRIPT_NAME = "lp_main.py"
TFORMAT = "%m/%d/%y %H:%M:%S"


def delete_datasets_in_workspace():
    """Delete all datasets in workspace."""
    datasets = arcpy.ListDatasets()
    for dataset in datasets:
        delete_dataset(dataset)


def delete_dataset(dataset):
    """Delete one dataset."""
    try:
        arcpy.Delete_management(dataset)
    except arcpy.ExecuteError:
        arcpy.AddWarning("Error deleting scratch/intermediate/temporary dataset %s. Program will continue." % dataset)


def delete_arcpy_temp_datasets():
    """Delete datasets left behind by arcpy."""
    if arcpy.Exists(os.path.join(lm_env.SCRATCHDIR, "scratch.gdb")):
        # datasets in scratch gdb (including core_resistance_stats and old hangovers)
        arcpy.env.workspace = os.path.join(lm_env.SCRATCHDIR, "scratch.gdb")
        delete_datasets_in_workspace()
    # datasets in intermediate gdb (including old hangovers)
    if not lp_env.KEEPINTERMEDIATE:
        if arcpy.Exists(os.path.join(lm_env.SCRATCHDIR, "intermediate.gdb")):
            arcpy.env.workspace = os.path.join(lm_env.SCRATCHDIR, "intermediate.gdb")
            delete_datasets_in_workspace()
    # datasets in current OS directory
    arcpy.env.workspace = os.getcwd()
    delete_datasets_in_workspace()


def print_runtime(stime):
    """Print process time when running from script."""
    etime = datetime.now()
    rtime = etime - stime
    hours, minutes = ((rtime.days * 24 + rtime.seconds // 3600),
                      (rtime.seconds // 60) % 60)
    print "End time: %s" % etime.strftime(TFORMAT)
    print "Elapsed time: %s hrs %s mins" % (hours, minutes)


def main(argv=None):
    """Main function for Linkage Priority tool."""
    start_time = datetime.now()
    print "Start time: %s" % start_time.strftime(TFORMAT)

    if argv is None:
        # get parameters passed from tool dialog, if any
        argv = sys.argv
    try:
        # preparation steps
        lp_env.configure(argv)
        check_lp_project_dir()
        check_out_sa_license()
        arc_wksp_setup()
        config_lm()
        log_setup()
        # primary analysis
        run_analysis()
    except arcpy.ExecuteError:
        # arcpy error
        msg = arcpy.GetMessages(2)
        arcpy.AddError(msg)
        lm_util.write_log(msg)
        exc_traceback = sys.exc_info()[2]
        lm_util.gprint("Traceback (most recent call last):\n" + "".join(traceback.format_tb(exc_traceback)[:-1]))
    except Exception:
        # other exceptions
        exc_value, exc_traceback = sys.exc_info()[1:]
        arcpy.AddError(exc_value)
        lm_util.gprint("Traceback (most recent call last):\n" + "".join(traceback.format_tb(exc_traceback)))
    finally:
        delete_arcpy_temp_datasets()
        arcpy.CheckInExtension("Spatial")
        print_runtime(start_time)
        lm_util.close_log_file()


def check_lp_project_dir():
    """Checks to make sure path name is not too long. Long path names can cause problems with ESRI grids."""
    if len(lp_env.PROJDIR) > 140:
        msg = ('ERROR: Project directory "' + lp_env.PROJDIR + '" is too deep.  Please choose a shallow directory' +
               r'(something like "C:\PUMA").')
        raise Exception(msg)

    if ("-" in lp_env.PROJDIR or " " in lp_env.PROJDIR or
            "." in lp_env.PROJDIR):
        msg = ("ERROR: Project directory cannot contain spaces, dashes, or special characters.")
        raise Exception(msg)


def check_out_sa_license():
    """Check out the ArcGIS Spatial Analyst extension license."""
    if arcpy.CheckExtension("Spatial") == "Available":
        arcpy.CheckOutExtension("Spatial")
    else:
        msg = ("ERROR: Spatial Analyst extension not available.")
        raise Exception(msg)


def arc_wksp_setup():
    """Setup ArcPy workspace."""
    arcpy.env.overwriteOutput = True


def config_lm():
    """Configure Linkage Mapper."""
    # get log file for last LM run
    folder = os.path.join(lp_env.PROJDIR, "run_history", "log") # lm_env.LOGDIR
    if not os.path.exists(folder):
        raise Exception("ERROR: Log file for last Linkage Mapper run not found. Please ensure Linkage Mapper is run " +
                        "for this project before running Linkage Priority.")
    entries = (os.path.join(folder, filename) for filename in os.listdir(folder))
    entries = ((os.stat(filepath), filepath) for filepath in entries)
    entries = ((stats.st_ctime, filepath) for stats, filepath in entries)
    last_lm_log = ""
    for cdate, filepath in sorted(entries, reverse=True):
        if "_Linkage Mapper" in filepath:
            last_lm_log = filepath
            break
    if last_lm_log == "":
        raise Exception("ERROR: Log file for last Linkage Mapper run not found")

    # read parameters section from file and turn into tuple for passing
    parms = ""
    with open(last_lm_log) as log_file:
        for line in log_file:
            if line[0:13] == "Parameters:\t[":
                parms = line[13:len(line) - 3]
                break
    lm_arg = tuple(parms.replace("'", "").split(", "))
    if parms == "":
        raise Exception("ERROR: Log file for last Linkage Mapper run does not contain a Parameters line")

    # use LM config in addition to LP config by passing in parms from last LM run
    try:
        lm_env.configure(lm_env.TOOL_LP, lm_arg)
    except Exception:
        raise Exception("ERROR: Parameters from log file for last Linkage Mapper run could not be used. " +
                        "Please Ensure the last run was completed with the current version of Linkage Mapper.")

    lm_util.gprint("\nLinkage Priority Version " + lm_env.releaseNum)


def log_setup():
    """Set up Linkage Mapper logging."""
    lm_env.logFilePath = lm_util.create_log_file(lm_env.MESSAGEDIR, lm_env.TOOL, lp_env.PARAMS)
    if lp_env.LPCUSTSETTINGS_IN:
        lm_util.write_log("Linkage Priority (LP) settings from " + lp_env.LPCUSTSETTINGS_IN + ":")
    else:
        lm_util.write_log("Linkage Priority (LP) settings from lp_settings.py:")
    lm_util.write_log("RELPERMNORMETH: %s" % (lp_env.RELPERMNORMETH))
    lm_util.write_log("RELCLOSENORMETH: %s" % (lp_env.RELCLOSENORMETH))
    lm_util.write_log("CALCLP: %s" % (lp_env.CALCLP))
    lm_util.write_log("NORMCORRNORMETH: %s" % (lp_env.NORMCORRNORMETH))
    lm_util.write_log("RESNORMETH: %s" % (lp_env.RESNORMETH))
    lm_util.write_log("SIZENORMETH: %s" % (lp_env.SIZENORMETH))
    lm_util.write_log("APNORMETH: %s" % (lp_env.APNORMETH))
    lm_util.write_log("ECAVNORMETH: %s" % (lp_env.ECAVNORMETH))
    lm_util.write_log("MINCPV: %s" % (lp_env.MINCPV))
    lm_util.write_log("NORMALIZERCI: %s" % (lp_env.NORMALIZERCI))
    lm_util.write_log("TRUNCNORMETH: %s" % (lp_env.TRUNCNORMETH))
    lm_util.write_log("CALCBP: %s" % (lp_env.CALCBP))
    lm_util.write_log("NORMALIZELP: %s" % (lp_env.NORMALIZELP))
    lm_util.write_log("NORMALIZEBP: %s" % (lp_env.NORMALIZEBP))
    lm_util.write_log("KEEPINTERMEDIATE: %s" % (lp_env.KEEPINTERMEDIATE))
    lm_util.write_log("MAXCSPWEIGHT: %s" % (lp_env.MAXCSPWEIGHT))
    lm_util.write_log("MEANCSPWEIGHT: %s" % (lp_env.MEANCSPWEIGHT))
    lm_util.write_log("")


def check_add_field(feature_class, field_name, data_type):
    """Check if field exists, and if not then add."""
    exists = False
    field_names = [field.name for field in arcpy.ListFields(feature_class)]
    if field_name in field_names:
        exists = True
    if not exists:
        arcpy.AddField_management(feature_class, field_name, data_type)
    return exists


def normalize_field(in_table, in_field, out_field, normalization_method, invert=False):
    """Normalize values in in_field into out_field using score range or max score method, with optional inversion."""
    check_add_field(in_table, out_field, "DOUBLE")
    min_val = arcpy.SearchCursor(in_table, "", "", "", in_field + " A").next().getValue(in_field)
    max_val = arcpy.SearchCursor(in_table, "", "", "", in_field + " D").next().getValue(in_field)
    if max_val > 0:
        try:
            if normalization_method == 0:
                # 0 to 1 score range normalization
                if invert:
                    arcpy.CalculateField_management(in_table, out_field,
                                                    "(" + str(max_val) + " - !" + in_field + "!) / " + str(max_val - min_val),
                                                    "PYTHON_9.3")
                else:
                    arcpy.CalculateField_management(in_table, out_field,
                                                    "(!" + in_field  + "! - " + str(min_val) + ") / " + str(max_val - min_val),
                                                    "PYTHON_9.3")
            else:
                # max_val score normalization
                if invert:
                    arcpy.CalculateField_management(in_table, out_field,
                                                    "(" + str(max_val + min_val) + " - !" + in_field + "!) / " +
                                                    str(max_val), "PYTHON_9.3")
                else:
                    arcpy.CalculateField_management(in_table, out_field, "!" + in_field + "! / " + str(max_val), "PYTHON_9.3")
        except Exception:
            # other exception - assume it was caused by situation described in exception message
            exc_value, exc_traceback = sys.exc_info()[1:]
            arcpy.AddError(exc_value)
            lm_util.gprint("Traceback (most recent call last):\n" + "".join(traceback.format_tb(exc_traceback)))
            raise Exception("ERROR! MOST LIKELY CAUSE: One or more core areas are smaller than a pixel in the Resistance " +
                            "layer and/or Raster Analysis Cell Size environment setting. Try enlarging small core areas, " +
                            "resampling the Resistance layer or adjusting the Cell Size environment setting.")
    else:
        # set all 0s
        arcpy.CalculateField_management(in_table, out_field, "0", "PYTHON_9.3")


def normalize_raster(in_raster, normalization_method, invert=False):
    """Normalize values in in_raster using score range or max score method, with optional inversion."""
    lm_util.build_stats(in_raster)
    result = arcpy.GetRasterProperties_management(in_raster, "MINIMUM")
    min_val = float(result.getOutput(0))
    result = arcpy.GetRasterProperties_management(in_raster, "MAXIMUM")
    max_val = float(result.getOutput(0))
    if max_val > 0:
        if normalization_method == 0:
            # 0-1 score range normalization
            if invert:
                return (max_val - in_raster) / (max_val - min_val)
            else:
                return (in_raster - min_val) / (max_val - min_val)
        else:
            # max score normalization
            if invert:
                return (max_val + min_val - in_raster) / max_val
            else:
                return in_raster / max_val
    else:
        # set all 0s
        return in_raster * 0


def calc_permeability(lcp_lines):
    """Calculate raw and relative permeability for each Least Cost Path."""
    # raw
    lm_util.gprint("Calculating raw permeability for each LCP line")
    check_add_field(lcp_lines, "Raw_Perm", "DOUBLE")
    arcpy.CalculateField_management(lcp_lines, "Raw_Perm", "!LCP_Length! / !CW_Dist!", "PYTHON_9.3")

    # relative
    lm_util.gprint("Calculating relative permeability for each LCP line")
    normalize_field(lcp_lines, "Raw_Perm", "Rel_Perm", lp_env.RELPERMNORMETH)


def calc_closeness(lcp_lines):
    """Calculate relative closeness for each Least Cost Path."""
    lm_util.gprint("Calculating relative closeness for each LCP line")
    normalize_field(lcp_lines, "LCP_Length", "Rel_Close", lp_env.RELCLOSENORMETH, True)


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
        if not os.path.exists(os.path.join(lm_env.DATAPASSDIR, "nlcc", "nlc" + nlc_str)):
            break
        arcpy.env.workspace = os.path.join(lm_env.DATAPASSDIR, "nlcc", "nlc" + nlc_str)
        # process each corridor raster in folder
        for input_raster in arcpy.ListRasters():
            # max score normalization with inversion
            inv_norm_raster = normalize_raster(Raster(input_raster), lp_env.NORMCORRNORMETH, True)
            if not os.path.exists(os.path.join(lm_env.DATAPASSDIR, "nlcc", "nlc" + nlc_str, "inv_norm")):
                arcpy.CreateFolder_management(os.path.join(lm_env.DATAPASSDIR, "nlcc", "nlc" + nlc_str), "inv_norm")
            inv_norm_raster.save(os.path.join(lm_env.DATAPASSDIR, "nlcc", "nlc" + nlc_str, "inv_norm", input_raster))
        nlc_idx += 1
    arcpy.env.workspace = prev_ws


def cav():
    """Calculate Core Area Value (CAV) and its components for each core."""
    lm_util.gprint("Calculating Core Area Value (CAV) and its components for each core")
    arcpy.MakeFeatureLayer_management(lp_env.COREFC, "core_lyr")

    # check weights and warn if issues
    if lp_env.OCAVRAST_IN:
        if lp_env.RESWEIGHT + lp_env.SIZEWEIGHT + lp_env.APWEIGHT + lp_env.ECAVWEIGHT + lp_env.CFCWEIGHT +\
                lp_env.OCAVWEIGHT <> 1.0:
            lm_util.gprint("Warning: RESWEIGHT + SIZEWEIGHT + APWEIGHT + ECAVWEIGHT + CFCWEIGHT + OCAVWEIGHT <> 1.0")
    else:
        if lp_env.RESWEIGHT + lp_env.SIZEWEIGHT + lp_env.APWEIGHT + lp_env.ECAVWEIGHT + lp_env.CFCWEIGHT <> 1.0:
            lm_util.gprint("Warning: RESWEIGHT + SIZEWEIGHT + APWEIGHT + ECAVWEIGHT + CFCWEIGHT <> 1.0")
    if lp_env.OCAVWEIGHT > 0 and not lp_env.OCAVRAST_IN:
        lm_util.gprint("Warning: OCAVWEIGHT > 0 but no OCAV raster input provided")
    if lp_env.OCAVWEIGHT == 0 and lp_env.OCAVRAST_IN:
        lm_util.gprint("Warning: OCAV raster input provided, but OCAVWEIGHT = 0")

    # check/add fields
    check_add_field(lp_env.COREFC, "mean_res", "DOUBLE")
    check_add_field(lp_env.COREFC, "norm_res", "DOUBLE")
    check_add_field(lp_env.COREFC, "area", "DOUBLE")
    check_add_field(lp_env.COREFC, "norm_size", "DOUBLE")
    check_add_field(lp_env.COREFC, "perimeter", "DOUBLE")
    check_add_field(lp_env.COREFC, "ap_ratio", "DOUBLE")
    check_add_field(lp_env.COREFC, "norm_ratio", "DOUBLE")
    check_add_field(lp_env.COREFC, "cav", "DOUBLE")
    check_add_field(lp_env.COREFC, "norm_cav", "DOUBLE")
    check_add_field(lp_env.COREFC, "clim_env", "DOUBLE")
    check_add_field(lp_env.COREFC, "nclim_env", "DOUBLE")
    check_add_field(lp_env.COREFC, "fut_clim", "DOUBLE")
    check_add_field(lp_env.COREFC, "nfut_clim", "DOUBLE")
    check_add_field(lp_env.COREFC, "ocav", "DOUBLE")
    check_add_field(lp_env.COREFC, "nocav", "DOUBLE")
    if not check_add_field(lp_env.COREFC, "ecav", "DOUBLE"):
        if lp_env.ECAVWEIGHT > 0:
            lm_util.gprint("Warning: ECAVWEIGHT > 0 but no ecav field in Cores feature class")
        arcpy.CalculateField_management(lp_env.COREFC, "ecav", "0")
    check_add_field(lp_env.COREFC, "necav", "DOUBLE")

    # current flow centrality (CFC, CF_Central) is copied from Centrality Mapper
    if not check_add_field(lp_env.COREFC, "CF_Central", "DOUBLE"):
        # default to 0s
        arcpy.CalculateField_management(lp_env.COREFC, "CF_Central", "0")
    if lp_env.CFCWEIGHT > 0:
        # copy values from Centrality Mapper output (core_centrality.gdb.project_Cores) if available
        centrality_cores = os.path.join(lm_env.CORECENTRALITYGDB, lm_env.PREFIX + "_Cores")
        if arcpy.Exists(centrality_cores):
            arcpy.AddJoin_management("core_lyr", lp_env.COREFN, centrality_cores, lp_env.COREFN)
            arcpy.CalculateField_management("core_lyr", lp_env.CORENAME + ".CF_Central",
                                            "[" + lm_env.PREFIX + "_Cores.CF_Central]")
            arcpy.RemoveJoin_management("core_lyr")
        # ensure cores have at least one non-0 value for CFC (could have been copied above or set earlier)
        max_val = arcpy.SearchCursor(lm_env.COREFC, "", "", "", "CF_Central D").next().getValue("CF_Central")
        if max_val is None or max_val == 0:
            msg = ("ERROR: A Current Flow Centrality Weight (CFCWEIGHT) was provided but no Current Flow Centrality " +
                   "(CF_Central) values are available. Please run Centrality Mapper on this project, then run " +
                   "Linkage Priority.")
            raise Exception(msg)
    check_add_field(lp_env.COREFC, "ncfc", "DOUBLE")

    # calc mean resistance
    stats_table = ZonalStatisticsAsTable(lp_env.COREFC, lp_env.COREFN, lp_env.RESRAST_IN,
                                        os.path.join(lm_env.SCRATCHDIR, "scratch.gdb", "core_resistance_stats"))
    arcpy.AddJoin_management("core_lyr", lp_env.COREFN, stats_table, lp_env.COREFN)
    arcpy.CalculateField_management("core_lyr", lp_env.CORENAME + ".mean_res", "[core_resistance_stats.MEAN]")
    arcpy.RemoveJoin_management("core_lyr")

    # calc area, perimeter and ratio
    arcpy.CalculateField_management("core_lyr", "area", "!SHAPE.AREA!", "PYTHON_9.3")
    arcpy.CalculateField_management("core_lyr", "perimeter", "!SHAPE.LENGTH!", "PYTHON_9.3")
    arcpy.CalculateField_management("core_lyr", "ap_ratio", "!area! / !perimeter!", "PYTHON_9.3")

    # normalize CAV inputs
    # resistance - invert
    normalize_field("core_lyr", "mean_res", "norm_res", lp_env.RESNORMETH, True)

    # size
    normalize_field("core_lyr", "area", "norm_size", lp_env.SIZENORMETH)

    # area/perimeter ratio
    normalize_field("core_lyr", "ap_ratio", "norm_ratio", lp_env.APNORMETH)

    # ecav
    normalize_field("core_lyr", "ecav", "necav", lp_env.ECAVNORMETH)

    # cfc
    normalize_field("core_lyr", "CF_Central", "ncfc", lp_env.CFCNORMETH)

    # calc OCAV
    if lp_env.OCAVRAST_IN:
        # get max and min
        lm_util.build_stats(lp_env.OCAVRAST_IN)
        result = arcpy.GetRasterProperties_management(lp_env.OCAVRAST_IN, "MAXIMUM")
        max_ocav = float(result.getOutput(0))
        result = arcpy.GetRasterProperties_management(lp_env.OCAVRAST_IN, "MINIMUM")
        min_ocav = float(result.getOutput(0))
        # calc score range normalization on input
        ocav_raster = (Raster(lp_env.OCAVRAST_IN) - min_ocav) / (max_ocav - min_ocav)
        # calc aerial mean ocav for each core
        ocav_table = ZonalStatisticsAsTable(lp_env.COREFC, lp_env.COREFN, ocav_raster,
                                            os.path.join(lm_env.SCRATCHDIR, "scratch.gdb", "core_ocav_stats"))
        arcpy.AddJoin_management("core_lyr", lp_env.COREFN, ocav_table, lp_env.COREFN)
        arcpy.CalculateField_management("core_lyr", lp_env.CORENAME + ".ocav", "[core_ocav_stats.MEAN]")
        arcpy.RemoveJoin_management("core_lyr")
        # calc score range normalization on output
        normalize_field("core_lyr", "ocav", "nocav", 0)
        # calc CAV
        arcpy.CalculateField_management("core_lyr", "cav",
                                        "(!norm_res! * " + str(lp_env.RESWEIGHT) + ") + (!norm_size! * " +
                                        str(lp_env.SIZEWEIGHT) + ") + (!norm_ratio! * " + str(lp_env.APWEIGHT) +
                                        ") + (!necav! * " + str(lp_env.ECAVWEIGHT) + ") + (!ncfc! * " +
                                        str(lp_env.CFCWEIGHT) + ") + (!nocav! * " + str(lp_env.OCAVWEIGHT) + ")",
                                        "PYTHON_9.3")

    else:
        # calc CAV
        arcpy.CalculateField_management("core_lyr", "cav",
                                        "(!norm_res! * " + str(lp_env.RESWEIGHT) + ") + (!norm_size! * " +
                                        str(lp_env.SIZEWEIGHT) + ") + (!norm_ratio! * " + str(lp_env.APWEIGHT) +
                                        ") + (!necav! * " + str(lp_env.ECAVWEIGHT) + ") + (!ncfc! * " +
                                        str(lp_env.CFCWEIGHT) + ")", "PYTHON_9.3")

    # normalize CAV with score range normalization
    normalize_field("core_lyr", "cav", "norm_cav", 0)


def clim_env():
    """Calculate Climate Envelope(s) for each core."""
    lm_util.gprint("Calculating Current Climate Envelope (clim_env) for each core")

    # calc score range normalization on current climate envelope
    cce_raster = normalize_raster(Raster(lp_env.CCERAST_IN), 0)

    # calc aerial mean climate envelope for each core
    cce_table = ZonalStatisticsAsTable(lp_env.COREFC, lp_env.COREFN, cce_raster,
                                      os.path.join(lm_env.SCRATCHDIR, "scratch.gdb", "core_clim_env_stats"))
    arcpy.AddJoin_management("core_lyr", lp_env.COREFN, cce_table, lp_env.COREFN)
    arcpy.CalculateField_management("core_lyr", lp_env.CORENAME + ".clim_env", "[core_clim_env_stats.MEAN]")
    arcpy.RemoveJoin_management("core_lyr")

    # score range normalize resulting values
    normalize_field("core_lyr", "clim_env", "nclim_env", 0)

    # future climate
    if lp_env.FCERAST_IN:
        lm_util.gprint("Calculating Future Climate Envelope (fut_clim) for each core")

        # calc score range normalization on future climate envelope
        fce_raster = normalize_raster(Raster(lp_env.FCERAST_IN), 0)

        # calc aerial mean climate envelope for each core
        fce_table = ZonalStatisticsAsTable(lp_env.COREFC, lp_env.COREFN, fce_raster,
                                          os.path.join(lm_env.SCRATCHDIR, "scratch.gdb", "core_fclim_env_stats"))
        arcpy.AddJoin_management("core_lyr", lp_env.COREFN, fce_table, lp_env.COREFN)
        arcpy.CalculateField_management("core_lyr", lp_env.CORENAME + ".fut_clim", "[core_fclim_env_stats.MEAN]")
        arcpy.RemoveJoin_management("core_lyr")

        # score range normalize resulting values
        normalize_field("core_lyr", "fut_clim", "nfut_clim", 0)


def eciv():
    """Normalize Expert Corridor Importance Value (ECIV) for each corridor."""
    if lp_env.COREPAIRSTABLE_IN and lp_env.ECIVFIELD:
        lm_util.gprint("Normalizing Expert Corridor Importance Value (ECIV) for each corridor")

        # calc score range normalization
        normalize_field(lp_env.COREPAIRSTABLE_IN, lp_env.ECIVFIELD, "neciv", 0)


def csp(sum_rasters, count_non_null_cells_rasters, max_rasters, lcp_lines):
    """Calculate Corridor Specific Priority (CSP) for each corridor."""
    lm_util.gprint("Calculating Corridor Specific Priority (CSP) for each corridor")
    prev_ws = arcpy.env.workspace

    # check weights
    if lp_env.CCERAST_IN:
        if lp_env.CLOSEWEIGHT + lp_env.PERMWEIGHT + lp_env.CAVWEIGHT + lp_env.ECIVWEIGHT + lp_env.CEDWEIGHT <> 1.0:
            lm_util.gprint("Warning: CLOSEWEIGHT + PERMWEIGHT + CAVWEIGHT + ECIVWEIGHT + CEDWEIGHT <> 1.0")
    else:
        if lp_env.CLOSEWEIGHT + lp_env.PERMWEIGHT + lp_env.CAVWEIGHT + lp_env.ECIVWEIGHT <> 1.0:
            lm_util.gprint("Warning: CLOSEWEIGHT + PERMWEIGHT + CAVWEIGHT + ECIVWEIGHT <> 1.0")
    if lp_env.CEDWEIGHT > 0 and not lp_env.CCERAST_IN:
        lm_util.gprint("Warning: CEDWEIGHT > 0, but no Current Climate Envelope raster input provided")
    if lp_env.CEDWEIGHT == 0 and lp_env.CCERAST_IN:
        lm_util.gprint("Warning: Current Climate Envelope raster input provided, but CEDWEIGHT = 0")
    if lp_env.ECIVWEIGHT > 0 and ((not lp_env.COREPAIRSTABLE_IN) or (not lp_env.ECIVFIELD)):
        lm_util.gprint("Warning: ECIVWEIGHT > 0, but no Expert Corridor Importance Value field provided")
    if lp_env.ECIVWEIGHT == 0 and lp_env.COREPAIRSTABLE_IN and lp_env.ECIVFIELD:
        lm_util.gprint("Warning: Expert Corridor Importance Value field provided, but ECIVWEIGHT = 0")

    # could be multiple folders
    nlc_idx = 0
    while True:
        nlc_str = ""
        if nlc_idx > 0:
            nlc_str = str(nlc_idx)
        if not os.path.exists(os.path.join(lm_env.DATAPASSDIR, "nlcc", "nlc" + nlc_str, "inv_norm")):
            break
        arcpy.env.workspace = (os.path.join(lm_env.DATAPASSDIR, "nlcc", "nlc" + nlc_str, "inv_norm"))
        csp_rasters = []
        count_rasters = []
        # process each corridor raster in folder
        for inRaster in arcpy.ListRasters():
            # check for max 1
            result = arcpy.GetRasterProperties_management(inRaster, "MAXIMUM")
            max_in = float(result.getOutput(0))
            if max_in <> 1.0:
                lm_util.gprint("Warning: maximum " + inRaster + " value <> 1.0")

            # get cores from raster name
            name_parts = inRaster.partition("_")
            from_core = name_parts[0]
            to_core = name_parts[2]

            # check for corresponding link
            links = arcpy.SearchCursor(lcp_lines,
                                       "(From_Core = " + from_core + " AND To_Core = " + to_core +
                                       ") OR (From_Core = " + to_core + " AND To_Core = " + from_core + ")")
            if links:
                link = links.next()
                # get and avg CAVs for the core pair
                x_cav = arcpy.SearchCursor("core_lyr", lp_env.COREFN + " = " + from_core, "", "",
                                           "").next().getValue("norm_cav")
                y_cav = arcpy.SearchCursor("core_lyr", lp_env.COREFN + " = " + to_core, "", "",
                                           "").next().getValue("norm_cav")
                avg_cav = (x_cav + y_cav) / 2

                # get ECIV for the core pair
                neciv = 0
                if lp_env.COREPAIRSTABLE_IN and lp_env.ECIVFIELD:
                    neciv = arcpy.SearchCursor(lp_env.COREPAIRSTABLE_IN,
                                               "(" + lp_env.FROMCOREFIELD + " = " + from_core + " AND " +
                                               lp_env.TOCOREFIELD + " = " + to_core + ") OR" + "(" +
                                               lp_env.TOCOREFIELD + " = " + from_core + " AND " +
                                               lp_env.FROMCOREFIELD + " = " + to_core +
                                               ")").next().getValue("neciv")

                # get and difference climate envelopes for the core pair to create CED
                if lp_env.CCERAST_IN:
                    x_clim_env = arcpy.SearchCursor("core_lyr", lp_env.COREFN + " = " + from_core, "", "",
                                                    "").next().getValue("nclim_env")
                    y_clim_env = arcpy.SearchCursor("core_lyr", lp_env.COREFN + " = " + to_core, "", "",
                                                    "").next().getValue("nclim_env")
                    if lp_env.FCERAST_IN:
                        # climate envelope difference is calculated relative to future climate envelope
                        # cooler core uses future
                        if (x_clim_env > y_clim_env) and lp_env.HIGHERCE_COOLER:
                            y_clim_env = arcpy.SearchCursor("core_lyr", lp_env.COREFN + " = " + to_core, "", "",
                                                            "").next().getValue("nfut_clim")
                        else:
                            x_clim_env = arcpy.SearchCursor("core_lyr", lp_env.COREFN + " = " + from_core, "", "",
                                                            "").next().getValue("nfut_clim")
                    diff_clim_env = abs(x_clim_env - y_clim_env)

                    # calc weighted sum
                    output_raster = ((lp_env.CLOSEWEIGHT * link.getValue("Rel_Close")) +
                                    (lp_env.PERMWEIGHT * link.getValue("Rel_Perm")) +
                                    (0.0001 * Raster(inRaster)) + (lp_env.CAVWEIGHT * avg_cav) +
                                    (lp_env.ECIVWEIGHT * neciv) + (lp_env.CEDWEIGHT * diff_clim_env))
                else:
                    # calc weighted sum
                    output_raster = ((lp_env.CLOSEWEIGHT * link.getValue("Rel_Close")) +
                                    (lp_env.PERMWEIGHT * link.getValue("Rel_Perm")) +
                                    (0.0001 * Raster(inRaster)) + (lp_env.CAVWEIGHT * avg_cav) +
                                    (lp_env.ECIVWEIGHT * neciv))
                if lp_env.KEEPINTERMEDIATE:
                    # also save a copy for debugging purposes
                    arcpy.CopyRaster_management(output_raster,
                                                os.path.join(lm_env.SCRATCHDIR, "intermediate.gdb",
                                                             lm_env.PREFIX + "_CSP_" + from_core + "_" + to_core),
                                                None, None, None, None, None, "32_BIT_FLOAT")

                # get max and min
                lm_util.build_stats(output_raster)
                result = arcpy.GetRasterProperties_management(output_raster, "MAXIMUM")
                max_csp = float(result.getOutput(0))
                result = arcpy.GetRasterProperties_management(output_raster, "MINIMUM")
                min_csp = float(result.getOutput(0))

                # determine threshold value
                diff_csp = max_csp - min_csp
                thres_csp = max_csp - (diff_csp * lp_env.PROPCSPKEEP)

                # apply threshold
                con_raster = Con(output_raster, output_raster, "#", "VALUE > " + str(thres_csp))
                is_null_raster = IsNull(con_raster)
                count_raster = EqualTo(is_null_raster, 0)
                count_rasters.append(count_raster)
                csp_rasters.append(con_raster)
                if lp_env.KEEPINTERMEDIATE:
                    # also save a copy for debugging purposes
                    arcpy.CopyRaster_management(con_raster,
                                                os.path.join(lm_env.SCRATCHDIR, "intermediate.gdb",
                                                             lm_env.PREFIX + "_CSP_TOP_" + from_core + "_" + to_core),
                                                None, None, None, None, None, "32_BIT_FLOAT")
                del link, links

        # perform intermediate calculations on CSPs leading toward CPV
        arcpy.env.scratchWorkspace = lm_env.SCRATCHGDB
        sum_raster = CellStatistics(csp_rasters, "SUM", "DATA")
        sum_rasters.append(sum_raster)
        count_non_null_cells = CellStatistics(count_rasters, "SUM", "DATA")
        count_non_null_cells_rasters.append(count_non_null_cells)
        max_raster = CellStatistics(csp_rasters, "MAXIMUM", "DATA")
        max_rasters.append(max_raster)
        nlc_idx += 1

    arcpy.env.workspace = prev_ws


def cpv(sum_rasters, count_non_null_cells_rasters, max_rasters, cpv_raster):
    """Combine CSPs using Max and Mean to create overall Corridor Priority Value (CPV)."""
    lm_util.gprint("Combining CSPs using Max and Mean to create overall Corridor Priority Value (CPV)")
    sum_all_raster = CellStatistics(sum_rasters, "SUM", "DATA")
    count_all_non_null_cells_rasters = CellStatistics(count_non_null_cells_rasters, "SUM", "DATA")
    max_all_raster = CellStatistics(max_rasters, "MAXIMUM", "DATA")
    cpv_raster_tmp = (max_all_raster * lp_env.MAXCSPWEIGHT) + ((sum_all_raster / count_all_non_null_cells_rasters)
                                                               * lp_env.MEANCSPWEIGHT)
    cpv_raster_tmp.save(cpv_raster)
    if not lp_env.KEEPINTERMEDIATE:
        # clean-up CSP input rasters
        # could be multiple nlc folders
        nlc_idx = 0
        prev_ws = arcpy.env.workspace
        while True:
            nlc_str = ""
            if nlc_idx > 0:
                nlc_str = str(nlc_idx)
            if not os.path.exists(os.path.join(lm_env.DATAPASSDIR, "nlcc", "nlc" + nlc_str)):
                break
            arcpy.env.workspace = os.path.join(lm_env.DATAPASSDIR, "nlcc", "nlc" + nlc_str, "inv_norm")
            for inRaster in arcpy.ListRasters():
                lm_util.delete_data(inRaster)
            nlc_idx += 1
        arcpy.env.workspace = prev_ws


def rci(cpv_raster, rci_raster):
    """Clip CPV to the MINCPV and renormalize to create relative corridor importance (RCI)."""
    lm_util.gprint("Calculating overall Relative Corridor Importance (RCI)")
    tmp_raster = ExtractByAttributes(cpv_raster, " VALUE >= " + str(lp_env.MINCPV))
    if lp_env.NORMALIZERCI:
        # calc score range normalization
        rci_raster_tmp = normalize_raster(tmp_raster, 0)
    else:
        rci_raster_tmp = tmp_raster
    arcpy.CopyRaster_management(rci_raster_tmp, rci_raster, None, None, None, None, None, "32_BIT_FLOAT")


def linkage_priority(rci_raster, trunc_raster, lp_raster):
    """Clip RCI to extent of truncated raster (LP)."""
    lm_util.gprint("Calculating overall Linkage Priority")
    # mask RCI based on extent of truncated raster (which is based on CWDTHRESH)
    tmp_raster2 = ExtractByMask(rci_raster, trunc_raster)
    if lp_env.NORMALIZELP:
        # calc score range normalization
        lp_raster_tmp = normalize_raster(tmp_raster2, 0)
    else:
        lp_raster_tmp = tmp_raster2
    arcpy.CopyRaster_management(lp_raster_tmp, lp_raster, None, None, None, None, None, "32_BIT_FLOAT")


def norm_trunc(trunc_raster, norm_trunc_raster):
    """Invert and normalize truncated raster (NORMTRUNC)."""
    lm_util.gprint("Inverting and normalizing truncated raster (NORMTRUNC)")
    norm_trunc_raster_tmp = normalize_raster(Raster(trunc_raster), lp_env.TRUNCNORMETH, True)
    arcpy.CopyRaster_management(norm_trunc_raster_tmp, norm_trunc_raster, None, None, None, None, None, "32_BIT_FLOAT")


def blended_priority(norm_trunc_raster, lp_raster, bp_raster):
    """Calculate overall Blended Priority."""
    lm_util.gprint("Calculating overall Blended Priority")
    tmp_raster3 = (lp_env.TRUNCWEIGHT * Raster(norm_trunc_raster)) + (lp_env.LPWEIGHT * Raster(lp_raster))
    if lp_env.NORMALIZEBP:
        # calc score range normalization
        bp_raster_tmp = normalize_raster(tmp_raster3, 0)
    else:
        bp_raster_tmp = tmp_raster3
    arcpy.CopyRaster_management(bp_raster_tmp, bp_raster, None, None, None, None, None, "32_BIT_FLOAT")


def run_analysis():
    """Run main Linkage Priority analysis."""
    lm_util.gprint("Checking inputs")

    # check that LM finished with steps 3 and 5
    if not os.path.exists(os.path.join(lm_env.DATAPASSDIR, "linkTable_s3.csv")) or\
            not os.path.exists(os.path.join(lm_env.DATAPASSDIR, "linkTable_s5.csv")):
        msg = ("ERROR: Project directory must contain a successful Linkage Mapper run with Steps 3 and 5.")
        raise Exception(msg)

    # check/create gdb for scratch
    if not os.path.isdir(lm_env.SCRATCHDIR):
        os.makedirs(lm_env.SCRATCHDIR)
    if not arcpy.Exists(os.path.join(lm_env.SCRATCHDIR, "scratch.gdb")):
        arcpy.CreateFileGDB_management(lm_env.SCRATCHDIR, "scratch.gdb")
    arcpy.env.scratchWorkspace = os.path.join(lm_env.SCRATCHDIR, "scratch.gdb")

    # check/create gdb for intermediate
    if lp_env.KEEPINTERMEDIATE:
        if not arcpy.Exists(os.path.join(lm_env.SCRATCHDIR, "intermediate.gdb")):
            arcpy.CreateFileGDB_management(lm_env.SCRATCHDIR, "intermediate.gdb")

    # set key dataset locations
    lcp_lines = os.path.join(lm_env.OUTPUTDIR, "link_maps.gdb", lm_env.PREFIX + "_LCPs")
    cpv_raster = os.path.join(lm_env.OUTPUTGDB, lm_env.PREFIX + "_CPV")
    rci_raster = os.path.join(lm_env.OUTPUTGDB, lm_env.PREFIX + "_RCI")
    cutoff_text = str(lm_env.CWDTHRESH)
    if cutoff_text[-6:] == "000000":
        cutoff_text = cutoff_text[0:-6]+"m"
    elif cutoff_text[-3:] == "000":
        cutoff_text = cutoff_text[0:-3]+"k"
    trunc_raster = (lm_env.OUTPUTGDB + "\\" + lm_env.PREFIX + "_corridors_truncated_at_" + cutoff_text)
    norm_trunc_raster = os.path.join(lm_env.OUTPUTGDB, lm_env.PREFIX + "_NORMTRUNC")
    lp_raster = os.path.join(lm_env.OUTPUTGDB, lm_env.PREFIX + "_linkage_priority")
    bp_raster = os.path.join(lm_env.OUTPUTGDB, lm_env.PREFIX + "_blended_priority")

    # calc permeability
    calc_permeability(lcp_lines)

    # calc relative closeness
    calc_closeness(lcp_lines)

    # invert and normalize each corridor
    inv_norm()

    # calc Core Area Value (CAV) and its components for each core
    cav()

    # normalize Expert Corridor Importance Value (ECIV)
    eciv()

    # calc climate envelope
    if lp_env.CCERAST_IN:
        clim_env()

    # calc Corridor Specific Priority (CSP)
    prev_ws = arcpy.env.workspace
    sum_rasters = []
    count_non_null_cells_rasters = []
    max_rasters = []
    csp(sum_rasters, count_non_null_cells_rasters, max_rasters, lcp_lines)

    # calc Corridor Priority Value (CPV)
    cpv(sum_rasters, count_non_null_cells_rasters, max_rasters, cpv_raster)
    arcpy.env.workspace = prev_ws

    # calc Relative Corridor Importance (RCI)
    rci(cpv_raster, rci_raster)

    # calc Linkage Priority (LP)
    if lp_env.CALCLP:
        linkage_priority(rci_raster, trunc_raster, lp_raster)

    # calc Blended Priority (BP)
    if lp_env.CALCBP:
        if not lm_env.WRITETRUNCRASTER:
            msg = "When CALCBP = True, set WRITETRUNCRASTER = True in Linkage Mapper"
            lm_util.raise_error(msg)
        if not lp_env.CALCLP:
            msg = "When CALCBP = True, set CALCLP = True"
            lm_util.raise_error(msg)
        norm_trunc(trunc_raster, norm_trunc_raster)
        blended_priority(norm_trunc_raster, lp_raster, bp_raster)

    # save a copy of Cores as the "Output for ModelBuilder Precondition"
    if lp_env.OUTPUTFORMODELBUILDER:
        arcpy.CopyFeatures_management(lp_env.COREFC, lp_env.OUTPUTFORMODELBUILDER)


if __name__ == "__main__":
    # execution starts here
    # use parameters passed from ArcGIS
    main()
    # # hard code parameters for debugging (0th is the working folder, which is not explicit in the ArcGIS tool dialog)
    # args = [r"C:\GIS\LandAdvisor\LinkageMapper\LinkageMapper1_0_9_3\toolbox\scripts",
    #         r"C:\GIS\LandAdvisor\LinkageMapper\LinkageMapper1_0_9_3\MojaveTest\rg2d",
    #         r"C:\GIS\LandAdvisor\LinkageMapper\LinkageMapper1_0_9_3\MojaveTest\core_areas\singlecores_dissolved.shp",
    #         "CID",
    #         r"C:\GIS\LandAdvisor\LinkageMapper\LinkageMapper1_0_9_3\MojaveTest\anthro\human_mod.tif",
    #         "#",
    #         "0.33",
    #         "0.33",
    #         "0.34",
    #         "0.0",
    #         "0.0",
    #         "0.0",
    #         "#",
    #         "#",
    #         "#",
    #         "#",
    #         "#",
    #         "#",
    #         "false",
    #         "0.33",
    #         "0.33",
    #         "0.34",
    #         "0.0",
    #         "0.0",
    #         "0.008",
    #         "0.5",
    #         "0.5",
    #         "#",
    #         "#"]
    # main(args)
