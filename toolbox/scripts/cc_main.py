#!/usr/bin/env python2.6
# Author: Darren Kavanagh

"""Climate Linkage Mapper

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
import traceback
from datetime import datetime, timedelta
# import subprocess
# import pickle

import arcpy
import arcpy.sa as sa

from cc_config import cc_env
import cc_util
import lm_master
from lm_config import tool_env as lm_env
import lm_util

_SCRIPT_NAME = "cc_main.py"
# __version__ = "0.2.4"

FR_COL = "From_Core"
TO_COL = "To_Core"


def main(argv=None):
    """Main function for Climate Linkage Mapper tool"""
    tformat = "%m/%d/%y %H:%M:%S"
    stime = datetime.now()

    print "Start time: %s" % stime.strftime(tformat)  # Redundant in Arc
  
    zonal_tbl = "zstats.dbf"

    if argv is None:
        argv = sys.argv
    try:
        cc_env.configure(argv)        
        gisdbase = os.path.join(cc_env.proj_dir, "gwksp")
        cc_util.check_cc_project_dir()
        cc_util.remove_grass_wkspc(gisdbase)
        cc_util.add_grass_path(cc_env.gisbase)
        
        # Make sure no dll conflict with grass and ArcGIS
        gdal_check('Startup')
        
        import cc_grass_cwd  # Cannot import until configured

        # Check out the ArcGIS Spatial Analyst extension license
        os.environ['ESRI_SOFTWARE_CLASS'] = 'Professional'
        if arcpy.CheckExtension("Spatial") == "Available":
            arcpy.CheckOutExtension("Spatial")
        else:
            raise

        # Setup workspace
        arcpy.env.overwriteOutput = True
        arcpy.env.cellSize = "MAXOF"  # Setting to default. For batch runs.
        arcpy.env.workspace = cc_util.mk_proj_dir(cc_env.out_dir)
        arcpy.env.scratchWorkspace = cc_util.mk_proj_dir(cc_env.tmp_dir)
        # lm_outdir = cc_util.mk_proj_dir("lm_out")
        lm_outdir = cc_env.proj_dir

        # Configure Linkage Mapper
        lm_arg = (_SCRIPT_NAME, lm_outdir, cc_env.prj_core_fc, cc_env.core_fld,
                  cc_env.prj_resist_rast, "false", "false", "#", "#", "true",
                  "false", cc_env.prune_network, cc_env.max_nn, cc_env.nn_unit,
                  cc_env.keep_constelations, "true", "#", "#", "#")
        lm_env.configure(lm_env.TOOL_CC, lm_arg)
        lm_util.gprint('\nClimate Linkage Mapper Version ' + lm_env.releaseNum)
        
        lm_util.create_dir(lm_env.DATAPASSDIR)

        # Set up logging
        lm_util.create_dir(lm_env.LOGDIR)
        lm_util.create_dir(lm_env.MESSAGEDIR)
        lm_env.logFilePath=lm_util.create_log_file(lm_env.MESSAGEDIR, 
                                            lm_env.TOOL, lm_env.PARAMS) 
                                           
        # Clip inputs and create project area raster
        cc_copy_inputs()

        # Get zonal statistics for cores and climate
        lm_util.gprint("\nCALCULATING ZONAL STATISTICS FROM CLIMATE RASTER")
        climate_stats = arcpy.sa.ZonalStatisticsAsTable(cc_env.prj_core_fc,
            cc_env.core_fld, cc_env.prj_climate_rast, zonal_tbl, "DATA", "ALL")

        # Create core pairings table and limit based upon climate threshold
        core_pairings = create_pair_tbl(climate_stats)

        # Generate link table, calculate CWD and run Linkage Mapper
        if int(arcpy.GetCount_management(core_pairings).getOutput(0)) == 0:
            arcpy.AddWarning("\nNo core pairs within climate threshold. "
                             "Program will end")
        else:
            # Process pairings and generate link table
            grass_cores = process_pairings(core_pairings)
            if not grass_cores:
                arcpy.AddWarning("\nNo core pairs within Euclidean distances. "
                                 "Progam program will end")
            else:
                # Create CWD using Grass
                cc_grass_cwd.grass_cwd(grass_cores)

                # Alternative ways to call Grass module
                # pfile_nm = os.path.join(cc_env.out_dir, "tmp.pck")
                # with open(pfile_nm, "wb") as pfile:
                    # pickle.dump((argv, lm_arg, grass_cores), pfile)
                # mod_path = cc_env.code_dir
                # cargs = ("import sys; sys.path.append('" + mod_path + "'); "
                         # "import cc_grass_cwd; "
                         # "cc_grass_cwd.grass_cwd('" + pfile_nm + "')")
                # cargs = cargs.replace("\\", "\\\\")
                # mybat = os.path.join(mod_path, "cc_rgrass")
                # cmd = '""%s" "%s""' % (mybat, cargs)
                # os.system(cmd)
                # output = subprocess.Popen(['python', '-c', cargs],
                                          # stdout=subprocess.PIPE,
                                          # shell=True).stdout.read()

                # Run Linkage Mapper
                lm_util.gprint("\nRUNNING LINKAGE MAPPER "
                                 "TO CREATE CLIMATE CORRIDORS")
                lm_master.lm_master()

    except arcpy.ExecuteError:
        arcpy.AddError(arcpy.GetMessages(2))
        exc_traceback = sys.exc_info()[2]
        lm_util.gprint("Traceback (most recent call last):\n" +
                         "".join(traceback.format_tb(exc_traceback)[:-1]))

    except Exception:
        exc_value, exc_traceback = sys.exc_info()[1:]
        arcpy.AddError(exc_value)
        lm_util.gprint("Traceback (most recent call last):\n" +
                         "".join(traceback.format_tb(exc_traceback)))
    finally:
        cc_util.delete_feature(cc_env.prj_climate_rast)
        cc_util.delete_feature(cc_env.prj_resist_rast)
        if cc_env.prj_resist_rast <> cc_env.prj_area_rast:
            cc_util.delete_feature(cc_env.prj_area_rast)
        cc_util.delete_feature(cc_env.prj_core_fc)  # Keeping for reruns
        if cc_env.simplify_cores:            
            cc_util.delete_feature(cc_env.core_simp)
        cc_util.delete_feature(cc_env.tmp_dir)
        
        # Delete files left behind by arcpy
        cpath = os.getcwd()
        cc_util.delete_feature(os.path.join(cpath, ".prj"))
        cc_util.delete_feature(os.path.join(cpath, "info"))
        
        arcpy.CheckInExtension("Spatial")

        # Print process time when running from script
        etime = datetime.now()
        rtime = etime - stime
        hours, minutes = ((rtime.days * 24 + rtime.seconds // 3600),
                       (rtime.seconds // 60) % 60)
        print "End time: %s" % etime.strftime(tformat)
        print "Elapsed time: %s hrs %s mins" % (hours, minutes)


def cc_copy_inputs():
    """Clip Climate Linkage Mapper inputs to smallest extent"""
    ext_poly = os.path.join(cc_env.out_dir,"ext_poly.shp")  # Extent polygon
    try:
        lm_util.gprint("\nCOPYING LAYERS AND, IF NECESSARY, REDUCING EXTENT")

        climate_extent = arcpy.Raster(cc_env.climate_rast).extent
        resist_extent = arcpy.Raster(cc_env.resist_rast).extent
        xmin = max(climate_extent.XMin, resist_extent.XMin)
        ymin = max(climate_extent.YMin, resist_extent.YMin)
        xmax = min(climate_extent.XMax, resist_extent.XMax)
        ymax = min(climate_extent.YMax, resist_extent.YMax)
        if cc_env.resist_rast is not None:       
            # Set to minimum extent if resistance raster was given
            arcpy.env.extent = arcpy.Extent(xmin, ymin, xmax, ymax)
            arcpy.CopyRaster_management(cc_env.resist_rast, cc_env.prj_resist_rast)

        arcpy.CopyRaster_management(cc_env.climate_rast,
                                    cc_env.prj_climate_rast)

        arcpy.env.extent = None

        ## Create project area raster 
        ## Replaced with extent coords code below
        # proj_area_rast = sa.Con(sa.IsNull(cc_env.prj_climate_rast),
                                # sa.Int(cc_env.prj_climate_rast), 1)
        # proj_area_rast.save(cc_env.prj_area_rast)

        
        # # Clip core feature class
        
        # This invokes gdal- replace with polygon based on extent coords below
        # arcpy.RasterToPolygon_conversion(proj_area_rast, ext_poly,
                                         # "NO_SIMPLIFY", "VALUE")
        
        # Create array of boundary points                    
        array = arcpy.Array()
        pnt = arcpy.Point(xmin,ymin)
        array.add(pnt)
        pnt = arcpy.Point(xmax,ymin)
        array.add(pnt)
        pnt = arcpy.Point(xmax,ymax)
        array.add(pnt)
        pnt = arcpy.Point(xmin,ymax)
        array.add(pnt)
        # Add in the first point of the array again to close polygon boundary
        array.add(array.getObject(0))
        # Create a polygon geometry object using the array object
        ext_feat = arcpy.Polygon(array)
        arcpy.CopyFeatures_management(ext_feat, ext_poly)
        # Clip core feature class
        arcpy.Clip_analysis(cc_env.core_fc, ext_poly, cc_env.prj_core_fc)                                          
                            
    except Exception:
        raise
    finally:
        cc_util.delete_feature(ext_poly)



def create_pair_tbl(climate_stats):
    """Create core pair table and limit to climate threshold """
    cpair_tbl = pair_cores("corepairs.dbf")
    if int(arcpy.GetCount_management(cpair_tbl).getOutput(0)) > 0:
        limit_cores(cpair_tbl, climate_stats)
    return cpair_tbl


def pair_cores(cpair_tbl):
    """Create table with all possible core to core combinations"""
    srow, srows, outputrow, irows = None, None, None, None

    try:
        lm_util.gprint("\nCREATING CORE PAIRINGS TABLE")
        arcpy.CreateTable_management(cc_env.out_dir, cpair_tbl, "", "")
        arcpy.AddField_management(cpair_tbl, FR_COL, "Short", "", "",
                                  "", "", "NON_NULLABLE")
        arcpy.AddField_management(cpair_tbl, TO_COL, "Short", "", "",
                                  "", "", "NON_NULLABLE")
        arcpy.DeleteField_management(cpair_tbl, "Field1")

        srows = arcpy.SearchCursor(cc_env.prj_core_fc, "", "",
                                   cc_env.core_fld, cc_env.core_fld + " A")

        cores_list = [srow.getValue(cc_env.core_fld) for srow in srows]
        cores_product = list(itertools.combinations(cores_list, 2))

        lm_util.gprint("There are " + str(len(cores_list)) + " unique "
                         "cores and " + str(len(cores_product)) + " pairings")

        irows = arcpy.InsertCursor(cpair_tbl)
        for nrow in cores_product:
            outputrow = irows.newRow()
            outputrow.setValue(FR_COL, int(nrow[0]))
            outputrow.setValue(TO_COL, int(nrow[1]))
            irows.insertRow(outputrow)

        return cpair_tbl

    except Exception:
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
    """Limit core pairs based upon climate threshold"""
    pair_vw = "dist_tbvw"
    stats_vw = "stats_tbvw"
    core_id = cc_env.core_fld.upper()

    try:
        lm_util.gprint("\nLIMITING CORE PAIRS BASED UPON CLIMATE "
                         "THRESHOLD")

        arcpy.MakeTableView_management(pair_tbl, pair_vw)
        arcpy.MakeTableView_management(stats_tbl, stats_vw)

        # Add basic stats to distance table
        lm_util.gprint("Joining zonal statistics to pairings table")
        add_stats(stats_vw, core_id, "fr", pair_vw, TO_COL)
        add_stats(stats_vw, core_id, "to", pair_vw, FR_COL)

        # Calculate difference of 2 std
        lm_util.gprint("Calculating difference of 2 std")
        diffu_2std = "diffu_2std"
        arcpy.AddField_management(pair_vw, diffu_2std, "Float", "", "",
                                  "", "", "NULLABLE")
        arcpy.CalculateField_management(pair_vw, diffu_2std,
                                        "abs(!frumin2std! - !toumin2std!)",
                                        "PYTHON")

        # Filter distance table based on inputed threshold and delete rows
        lm_util.gprint("Filtering table based on threshold")
        diffu2std_fld = arcpy.AddFieldDelimiters(pair_vw, diffu_2std)
        expression = diffu2std_fld + " <= " + str(cc_env.climate_threshold)
        arcpy.SelectLayerByAttribute_management(pair_vw, "NEW_SELECTION",
                                                expression)
        rows_del = int(arcpy.GetCount_management(pair_vw).getOutput(0))
        if rows_del > 0:
            arcpy.DeleteRows_management(pair_vw)
        lm_util.gprint(str(rows_del) + " rows deleted")

    except Exception:
        raise
    finally:
        cc_util.delete_feature(stats_vw)
        cc_util.delete_feature(pair_vw)


def add_stats(stats_vw, core_id, fld_pre, table_vw, join_col):
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


def process_pairings(pairings):
    """Limit core pairings based on distance inputs and create linkage table

    Requires ArcInfo license.

    """
    lm_util.gprint("\nLIMITING CORE PAIRS BASED ON INPUTED DISTANCES AND "
                     "GENERATING LINK TABLE")
    # Simplify cores based on booolean in config
    if cc_env.simplify_cores:
        corefc = simplify_corefc()
    else:
        corefc = cc_env.prj_core_fc
    core_pairs, frm_cores = pairs_from_list(pairings)
    # Create link table
    core_list = create_lnk_tbl(corefc, core_pairs, frm_cores)
    return sorted(core_list)


def pairs_from_list(pairings):
    """Get list of core pairings and 'from cores'"""
    frm_cores = set()
    core_pairs = []
    srows = arcpy.SearchCursor(pairings, "", "", FR_COL + "; " + TO_COL)
    for srow in srows:
        from_core = srow.getValue(FR_COL)
        to_core = str(srow.getValue(TO_COL))
        frm_cores.add(from_core)
        core_pairs.append([str(from_core), to_core])
    frm_cores = map(str, sorted(frm_cores))
    return core_pairs, frm_cores


def create_lnk_tbl(corefc, core_pairs, frm_cores):
    """Create link table file and limit based on near table results"""
    fcore_vw = "fcore_vw"
    tcore_vw = "tcore_vw"
    ntblid_fn = "FID"
    fmid_fn = "IN_FID"
    toid_fn = "NEAR_FID"
    ndist_fn = "NEAR_DIST"
    jtocore_fn = cc_env.core_fld[:8] + "_1"  # dbf field length
    near_tbl = os.path.join(cc_env.out_dir, "neartbl.dbf")
    link_file = os.path.join(lm_env.DATAPASSDIR, "linkTable_s2.csv")

    link_tbl, srow, srows = None, None, None

    try:
        link_tbl = open(link_file, 'wb')
        writer = csv.writer(link_tbl, delimiter=',')
        headings = ["# link", "coreId1", "coreId2", "cluster1", "cluster2",
                    "linkType", "eucDist", "lcDist", "eucAdj", "cwdAdj"]
        writer.writerow(headings)

        core_list = set()
        no_cores = str(len(frm_cores))
        i = 1

        coreid_fld = arcpy.AddFieldDelimiters(corefc, cc_env.core_fld)

        for core_no, frm_core in enumerate(frm_cores):
            # From cores
            expression = coreid_fld + " = " + frm_core
            arcpy.MakeFeatureLayer_management(corefc, fcore_vw, expression)

            # To cores
            to_cores_lst = [x[1] for x in core_pairs if frm_core == x[0]]
            to_cores = ', '.join(to_cores_lst)
            expression = coreid_fld + " in (" + to_cores + ")"
            arcpy.MakeFeatureLayer_management(corefc, tcore_vw, expression)
            lm_util.gprint("Calculating Euclidean distance/s from Core "
                             + frm_core + " to " + str(len(to_cores_lst))
                             + " other cores" + " (" + str(core_no + 1) + "/"
                             + no_cores + ")")

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

    except Exception:
        raise
    finally:
        cc_util.delete_feature(near_tbl)
        cc_util.delete_feature(os.path.splitext(corefc)[0] + "_Pnt.shp")
        if link_tbl:
            link_tbl.close()
        if srow:
            del srow
        if srows:
            del srows

    return core_list


def simplify_corefc():
    """Simplify core feature class"""
    lm_util.gprint("Simplifying polygons to speed up core pair "
                     "distance calculations")
    corefc = cc_env.core_simp
    climate_rast = arcpy.Raster(cc_env.prj_climate_rast)
    tolerance = climate_rast.meanCellHeight / 3
    arcpy.cartography.SimplifyPolygon(cc_env.prj_core_fc, corefc,
        "POINT_REMOVE", tolerance, "#", "NO_CHECK")
    return corefc

def gdal_check(msg):
    # Code to check GDAL dlls and system path
    import subprocess
    gdal = subprocess.Popen("where gdal*", stdout=subprocess.PIPE,
                            shell=True).stdout.read()
    # lm_util.gprint("\nGDAL DLL/s at " + msg + ': ' + gdal)
    gdalList = gdal.split('\n')                        
    if 'arcgis' in gdalList[1].lower():
        lm_util.gprint("\nGDAL DLL/s at " + msg + ': ' + gdal)
        arcpy.AddWarning("It looks like there is a conflict between ArcGIS")
        arcpy.AddWarning("and GRASS. This might be caused by conflicts with ")
        arcpy.AddWarning("pre-loaded ArcGIS extensions like Geostatistical") 
        arcpy.AddWarning("Analyst.")
        arcpy.AddWarning("\nPlease DISABLE ANY EXTENSIONS YOU ARE NOT USING") 
        arcpy.AddWarning("(Click on Customize >> Extensions) and try again.")
        arcpy.AddWarning("\nIf that doesn't work you can try closing Arc and ")
        arcpy.AddWarning("instead run the tool using the 'CC Run Script.py' ")
        arcpy.AddWarning("python script.  This script can be found in the ")
        arcpy.AddWarning("'scripts' directory, located where the Linkage") 
        arcpy.AddWarning("Mapper toolbox is installed.\n")
        raise Exception("ArcGIS-GRASS GDAL DLL conflict")
    
if __name__ == "__main__":
    main()
