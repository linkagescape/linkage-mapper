"""Restore for maximum return of investment (ROI).

Script to iteratively run linkage mapper and barrier mapper tools restoring
for max ROI.
"""

import sys
import os
import shutil
import time

import numpy as npy
import arcpy

import lm_master
import barrier_master
from lm_config import tool_env as cfg
import lm_util as lu
from lm_util import gprint


_SCRIPT_NAME = os.path.basename(__file__)

arcpy.CheckOutExtension("spatial")


def main():
    """Iterate over LM, BM, and restoration tasks."""
    # USER SETTINGS ######################################################

    # Restoration Settings
    # ALL input data must be in the same projection
    start_time = time.clock()

    # Set to True to restore highest ROI. Set to False to restore strongest
    # barrier
    restore_max_roi = False

    # Resistance value of restored habitat.  Must be 1 or greater.
    restored_resistance_val = 1

    # No spaces or special chars in paths or gdb names
    restoration_data_gdb = ("C:\\barrierClassAnalysis\\"
                            "RestorationINPUTS_July2013.gdb")

    # No spaces in path, avoid using dropbox or network drive
    # Project directories will be created in this (iter1, iter2...) as will an
    # output geodatabase
    output_dir = "C:\\barrierClassAnalysis\\output"

    # Resistance raster. Should be in input GDB
    resistance_ras = "URWA_resis"
    # Core area feature class. Should be in input GDB 'URWA_HCAs_Doug_Grant'
    core_fc = 'URWA_HCAs_Doug_Grant'

    core_fn = 'HCA_ID'  # Core area field name

    radius = 450  # Restoration radius in meters
    iterations = 13  # Number of restorations to perform

    # If less than this proportion of ag in circle, don't consider restoring
    # circle
    min_ag_threshold = 0.75

    # Don't consider barriers below this improvement score (average improvement
    # per meter diameter restored)
    min_improvement_val = 0

    # Average per-m2 parcel cost per pixel. Snapped to resistance raster.
    parcel_cost_ras = 'DougGrantParcelCost_m2_projected_90m'

    # Right now this is just a raster with all pixels set to 0.113174
    restoration_cost_ras = 'restCostPer_m2'

    ag_ras = "ARESmaskp_projected"  # 1=Ag, 0=Not Ag

    # Some restorations benefit multiple corridors.
    # 'Maximum' takes the greatest improvement across core area pairs
    # 'Sum' adds improvement scores acreoss all pairs.
    barrier_combine_method = 'Maximum'

    # Use cwd_thresh = None for no threshold. Use cwd_thresh = X to not
    # consider restorations more than X map units away from each core area.
    cwd_thresh = None

    # END USER SETTINGS ######################################################

    try:
        # Setup path and create directories
        gprint('Hey! Make sure everything is in the same projection!\n')
        gprint('Setting up paths and creating directories')
        sys.path.append('..\\toolbox\\scripts')
        res_ras = os.path.join(restoration_data_gdb, resistance_ras)
        core_fc_path = os.path.join(restoration_data_gdb, core_fc)

        # Set up a NEW output gdb (leave previous ones on drive)
        i = None
        for i in range(1, 200):
            output_gdb = 'restorationOutput' + str(i) + '.gdb'
            if not arcpy.Exists(os.path.join(output_dir, output_gdb)):
                break
            gprint('Previous output GDB ' + output_gdb + ' exists.  '
                   'Delete to save disk space.')
        arcpy.CreateFileGDB_management(output_dir, output_gdb)
        output_gdb = os.path.join(output_dir, output_gdb)
        log_file = os.path.join(output_gdb,
                                'Iterate Barriers' + str(i) + '.py')

        # Write a copy of this file to output dir as a record of settings
        shutil.copyfile(__file__, log_file)

        arcpy.env.cellSize = res_ras
        arcpy.env.extent = res_ras
        arcpy.env.snapRaster = res_ras
        arcpy.env.overwriteOutput = True
        arcpy.env.scratchWorkspace = output_gdb
        arcpy.env.workspace = output_gdb

        spatialref = arcpy.Describe(res_ras).spatialReference
        mapunits = spatialref.linearUnitName
        gprint('Cell size = ' + str(arcpy.env.cellSize) + ' ' + mapunits + 's')

        # Calculate fraction of ag within radius of each pixel
        gprint('Calculating purchase cost, fraction of ag, etc within radius '
               'of each pixel.')
        ag_ras = os.path.join(restoration_data_gdb, ag_ras)
        in_neighborhood = arcpy.sa.NbrCircle(radius, "MAP")
        arcpy.env.extent = ag_ras
        out_focal_stats = arcpy.sa.FocalStatistics(
            ag_ras, in_neighborhood, "MEAN", "NODATA")
        proportion_ag_ras = os.path.join(output_gdb, 'proportionAgRas')
        out_focal_stats.save(proportion_ag_ras)
        arcpy.env.extent = res_ras

        # Calculate purchase cost of circles
        parcel_cost_ras = os.path.join(restoration_data_gdb,
                                       parcel_cost_ras)
        arcpy.env.extent = parcel_cost_ras
        out_focal_stats = arcpy.sa.FocalStatistics(
            parcel_cost_ras, in_neighborhood, "MEAN", "DATA")
        cost_focal_stats_ras = os.path.join(output_gdb, 'cost_focal_stats_ras')
        out_focal_stats.save(cost_focal_stats_ras)
        arcpy.env.extent = res_ras

        circle_area = float(npy.pi * radius * radius)
        outras = arcpy.sa.Raster(cost_focal_stats_ras) * circle_area
        purch_cost_ras = os.path.join(output_gdb, 'purchaseCostRaster')
        outras.save(purch_cost_ras)
        lu.delete_data(cost_focal_stats_ras)

        restoration_cost_ras = os.path.join(restoration_data_gdb,
                                            restoration_cost_ras)
        outras = (arcpy.sa.Raster(purch_cost_ras)
                  + (arcpy.sa.Raster(restoration_cost_ras) *
                     radius * radius * npy.pi))
        total_cost_ras = os.path.join(output_gdb, 'totalCostRaster')
        outras.save(total_cost_ras)

        # Create mask to remove areas without cost data
        arcpy.env.extent = total_cost_ras
        cost_mask_ras = os.path.join(output_gdb, 'costMaskRaster')
        cost_thresh = 0
        out_con = arcpy.sa.Con((arcpy.sa.Raster(total_cost_ras) >
                                float(cost_thresh)), 1)
        out_con.save(cost_mask_ras)
        arcpy.env.extent = res_ras

        # Create mask to remove areas below ag threshold
        out_con = arcpy.sa.Con((arcpy.sa.Raster(proportion_ag_ras) >
                                float(min_ag_threshold)), 1)
        ag_mask_ras = os.path.join(output_gdb, 'agMaskRaster')
        out_con.save(ag_mask_ras)

        do_step_1 = 'true'
        do_step_2 = 'true'
        do_step_5 = 'false'
        all_restored_areas_ras = ''

        for cur_iter in range(1, iterations + 1):
            start_time1 = time.clock()

            # Some env settings get changed by linkage mapper and must be
            # reset here
            arcpy.env.cellSize = res_ras
            arcpy.env.extent = res_ras
            arcpy.env.snapRaster = res_ras
            arcpy.env.scratchWorkspace = output_gdb
            arcpy.env.workspace = output_gdb

            lu.dashline(1)
            gprint('Running iteration number ' + str(cur_iter))
            proj_dir = os.path.join(output_dir, 'iter' + str(cur_iter)
                                    + 'Proj')
            lu.create_dir(output_dir)
            lu.delete_dir(proj_dir)
            lu.create_dir(proj_dir)
            if cur_iter > 1:  # Copy previous s2 linktable to new project dir
                datapass_dir = os.path.join(proj_dir, 'datapass')
                lu.create_dir(datapass_dir)
                proj_dir1 = os.path.join(output_dir, 'iter1Proj')
                datapass_dir_iter1 = os.path.join(proj_dir1, 'datapass')
                s2_link_tbl_iter1 = os.path.join(datapass_dir_iter1,
                                                 'linkTable_s2.csv')
                s2_link_tbl = os.path.join(datapass_dir, 'linkTable_s2.csv')
                shutil.copyfile(s2_link_tbl_iter1, s2_link_tbl)

            # Run Linkage Mapper

            # Copy distances text file from earlier LM run to the output
            # directory- speeds things up!
            dist_file = os.path.join(output_dir, core_fc + '_dists.txt')

            if not os.path.exists(dist_file):
                if cur_iter == 1:
                    gprint('Will calculate distance file.')
                    dist_file = '#'
                else:
                    proj_dir1 = os.path.join(output_dir, 'iter1Proj')
                    dist_file1 = os.path.join(proj_dir1, core_fc
                                              + '_dists.txt')
                    # Put a copy here for future runs
                    shutil.copyfile(dist_file1, dist_file)

            arcpy.env.scratchWorkspace = output_gdb
            arcpy.env.workspace = output_gdb

            argv = ('lm_master.py', proj_dir, core_fc_path, core_fn, res_ras,
                    do_step_1, do_step_2, 'Cost-Weighted & Euclidean',
                    dist_file, 'true', 'true', 'false', '4', 'Cost-Weighted',
                    'true', do_step_5, 'true', '200000', '10000', '#', '#',
                    '#', '#')
            gprint('Running ' + str(argv))
            cfg.lm_configured = False
            lm_master.lm_master(argv)
            do_step_1 = 'false'  # Can skip for future iterations
            do_step_2 = 'false'  # Can skip for future iterations
            do_step_5 = 'false'  # Skipping for future iterations

            start_radius = str(radius)
            end_radius = str(radius)
            radius_step = '0'
            save_radius_ras = 'false'
            write_pct_ras = 'false'

            argv = ('barrier_master.py', proj_dir, res_ras, start_radius,
                    end_radius, radius_step, barrier_combine_method,
                    save_radius_ras, write_pct_ras, cwd_thresh)
            gprint('Running ' + str(argv))
            barrier_master.bar_master(argv)

            # Some env settings get changed by linkage mapper and must be
            # reset here
            arcpy.env.cellSize = res_ras
            arcpy.env.extent = res_ras
            arcpy.env.snapRaster = res_ras
            arcpy.env.scratchWorkspace = output_gdb
            arcpy.env.workspace = output_gdb

            gprint('Finding restoration circles with max barrier score / ROI')
            # Find points with max ROI
            prefix = os.path.basename(proj_dir)
            if barrier_combine_method == 'Sum':
                sum_suffix = 'Sum'
            else:
                sum_suffix = ''
            barrier_fn = (prefix + "_BarrierCenters" + sum_suffix + "_Rad" +
                          str(radius))
            barrier_ras = os.path.join(proj_dir, 'output', 'barriers.gdb',
                                       barrier_fn)
            if not arcpy.Exists(barrier_ras):
                msg = ('Error: cannot find barrier output: ' + barrier_ras)
                lu.raise_error(msg)

            if cur_iter > 1:
                gprint('Creating mask for previously restored areas')
                in_neighborhood = arcpy.sa.NbrCircle(radius, "MAP")
                arcpy.env.extent = all_restored_areas_ras
                out_focal_stats = arcpy.sa.FocalStatistics(
                    all_restored_areas_ras,
                    in_neighborhood, "MEAN", "DATA")
                all_restored_focal_ras = os.path.join(
                    output_gdb, 'allRestFocRas_iter' + str(cur_iter))

                # Anything > 0 would include a restored area
                out_focal_stats.save(all_restored_focal_ras)
                arcpy.env.extent = res_ras
                rest_mask_ras = os.path.join(
                    output_gdb, 'restMaskRaster_iter' + str(cur_iter))
                minval = 0
                out_con = arcpy.sa.Con(
                    (arcpy.sa.Raster(all_restored_focal_ras) ==
                     float(minval)),
                    1)
                out_con.save(rest_mask_ras)

            # Candidate areas have not been restored, have cost data, meet
            # minimum improvement score criteria, and have enough ag in them
            candidate_barrier_ras = os.path.join(
                output_gdb, 'candidateBarrierRaster' + '_iter' + str(cur_iter))
            if cur_iter > 1:
                gprint('Creating candidate restoration raster using barrier '
                       'results, previous restorations, and selection '
                       'criteria')

                # ROI scores will be in terms of total improvement
                # (= score * diameter)
                out_calc = (arcpy.sa.Raster(cost_mask_ras) *
                            arcpy.sa.Raster(ag_mask_ras) *
                            arcpy.sa.Raster(barrier_ras) *
                            arcpy.sa.Raster(rest_mask_ras) *
                            (radius * 2))
            else:
                out_calc = (arcpy.sa.Raster(cost_mask_ras) *
                            arcpy.sa.Raster(ag_mask_ras) *
                            arcpy.sa.Raster(barrier_ras) *
                            radius * 2)

            min_barrier_score = min_improvement_val * radius * 2
            if restored_resistance_val != 1:
                out_calc_2 = (out_calc -
                              (2 * radius * (restored_resistance_val - 1)))
                out_con = arcpy.sa.Con(
                    (out_calc_2 >= float(min_barrier_score)), out_calc_2)
            else:
                out_con = arcpy.sa.Con((out_calc >= float(min_barrier_score)),
                                       out_calc)
            out_con.save(candidate_barrier_ras)
            lu.build_stats(candidate_barrier_ras)

            purchase_roi_ras = os.path.join(output_gdb, 'purchaseRoiRaster'
                                            + '_iter' + str(cur_iter))
            out_calc = (arcpy.sa.Raster(candidate_barrier_ras) /
                        arcpy.sa.Raster(purch_cost_ras))
            out_calc.save(purchase_roi_ras)
            lu.build_stats(purchase_roi_ras)

            total_roi_ras = os.path.join(output_gdb, 'purchaseRestRoiRaster'
                                         + '_iter' + str(cur_iter))
            out_calc = (arcpy.sa.Raster(candidate_barrier_ras) /
                        arcpy.sa.Raster(total_cost_ras))
            out_calc.save(total_roi_ras)
            lu.build_stats(total_roi_ras)

            max_barrier = arcpy.GetRasterProperties_management(
                candidate_barrier_ras, "MAXIMUM")
            gprint('Maximum barrier improvement score: '
                   + str(max_barrier.getOutput(0)))
            if max_barrier < 0:
                arcpy.AddWarning("\nNo barriers found that meet CWD or Ag "
                                 "threshold criteria.")

            max_purch_roi = arcpy.GetRasterProperties_management(
                purchase_roi_ras, "MAXIMUM")
            gprint('Maximum purchase ROI score: '
                   + str(max_purch_roi.getOutput(0)))

            max_roi = arcpy.GetRasterProperties_management(total_roi_ras,
                                                           "MAXIMUM")
            gprint('Maximum total ROI score: ' + str(max_roi.getOutput(0)))

            if restore_max_roi:
                out_point = os.path.join(output_gdb, 'maxRoiPoint'
                                         + '_iter' + str(cur_iter))
                gprint('Choosing circle with maximum ROI to restore')
                out_con = arcpy.sa.Con(
                    (arcpy.sa.Raster(total_roi_ras)
                     >= float(max_roi.getOutput(0))),
                    total_roi_ras)
                max_roi_ras = os.path.join(output_gdb, 'max_roi_ras')
                out_con.save(max_roi_ras)
                # Save max ROI to point
                try:
                    arcpy.RasterToPoint_conversion(max_roi_ras, out_point)
                except Exception:
                    msg = ('Error: it looks like there are no viable '
                           'restoration candidates.')
                    lu.raise_error(msg)

            else:  # Restoring strongest barrier instead
                out_point = os.path.join(output_gdb, 'maxBarrierPoint'
                                         + '_iter' + str(cur_iter))
                gprint('Choosing circle with maximum BARRIER IMPROVEMENT SCORE'
                       ' to restore')
                out_con = arcpy.sa.Con(
                    (arcpy.sa.Raster(candidate_barrier_ras) >=
                     float(max_barrier.getOutput(0))),
                    candidate_barrier_ras)
                max_barrier_ras = os.path.join(output_gdb, 'maxBarrierRaster')
                out_con.save(max_barrier_ras)
                # Save max barrier to point
                try:
                    arcpy.RasterToPoint_conversion(max_barrier_ras, out_point)
                except Exception:
                    msg = ('Error: it looks like there are no viable '
                           'restoration candidates.')
                    lu.raise_error(msg)

            gprint('Done evaluating candidate restorations')
            result = int(arcpy.GetCount_management(out_point).getOutput(0))
            if result > 1:
                # Would be better to retain point with max barrier score when
                # we have multiple points with same ROI
                arcpy.AddWarning('Deleting points with identical '
                                 'ROI/improvement score values')

                arcpy.DeleteIdentical_management(out_point, "grid_code",
                                                 0.1, 0.1)

            arcpy.sa.ExtractMultiValuesToPoints(
                out_point,
                [
                    [candidate_barrier_ras, "barrierScore"],
                    [purch_cost_ras, "purchCost"],
                    [total_cost_ras, "totalCost"],
                    [purchase_roi_ras, "purchaseROI"],
                    [total_roi_ras, "totalROI"]
                ],
                "NONE")

            arcpy.AddField_management(out_point, "restorationNumber", "SHORT")
            arcpy.CalculateField_management(out_point, "restorationNumber",
                                            cur_iter)
            arcpy.AddField_management(out_point, "radius", "DOUBLE")
            arcpy.CalculateField_management(out_point, "radius", radius)
            arcpy.AddField_management(out_point, "barrierScore_per_m",
                                      "DOUBLE")
            arcpy.CalculateField_management(
                out_point, "barrierScore_per_m",
                "(float(!barrierScore!) / (!radius! * 2))",
                "PYTHON")

            gprint('\nCreating restoration circles')
            if restore_max_roi:
                circle_fc = os.path.join(output_gdb, 'maxRoiCircle'
                                         + '_iter' + str(cur_iter))
            else:
                circle_fc = os.path.join(output_gdb, 'maxBarrierCircle'
                                         + '_iter' + str(cur_iter))
            arcpy.Buffer_analysis(out_point, circle_fc, radius)
            gprint('Rasterizing restoration circles')
            if restore_max_roi:
                circle_ras = os.path.join(output_gdb, 'maxRoicircle_ras'
                                          + '_iter' + str(cur_iter))
            else:
                circle_ras = os.path.join(output_gdb, 'maxBarrierCircleRas'
                                          + '_iter' + str(cur_iter))
            arcpy.FeatureToRaster_conversion(circle_fc, 'totalROI', circle_ras,
                                             arcpy.env.cellSize)

            # restore raster
            gprint('Digitally restoring resistance raster')
            res_ras_restored = os.path.join(output_gdb, 'resRastRestored'
                                            + '_iter' + str(cur_iter))
            out_con = arcpy.sa.Con(arcpy.sa.IsNull(circle_ras), res_ras,
                                   restored_resistance_val)
            out_con.save(res_ras_restored)

            all_restored_areas_ras = os.path.join(
                output_gdb, 'allRestoredAreas_iter' + str(cur_iter))
            prev_restored_areas_ras = os.path.join(
                output_gdb, 'allRestoredAreas_iter' + str(cur_iter-1))
            if cur_iter == 1:
                out_con = arcpy.sa.Con(arcpy.sa.IsNull(circle_ras), 0, 1)
            else:
                # Add this restoration to areas restored
                out_con = arcpy.sa.Con(arcpy.sa.IsNull(circle_ras),
                                       prev_restored_areas_ras, 1)
            out_con.save(all_restored_areas_ras)

            lu.delete_data(circle_ras)

            # Use for next iteration resistance raster
            res_ras = res_ras_restored

            # Add circle into feature class with all circles
            if restore_max_roi:
                all_circles_fc = os.path.join(output_gdb, "allCirclesMaxROI")
            else:
                all_circles_fc = os.path.join(output_gdb,
                                              "allCirclesMaxBarriers")
            if cur_iter == 1:
                arcpy.CopyFeatures_management(circle_fc, all_circles_fc)
            else:
                arcpy.Append_management(circle_fc, all_circles_fc, "TEST")
            gprint('Finished iteration #' + str(cur_iter))
            start_time1 = lu.elapsed_time(start_time1)

        gprint('\nDone with iterations.')
        start_time = lu.elapsed_time(start_time)
        gprint('Outputs saved in: ' + output_gdb)
        gprint('Back up your project directories if you want to save '
               'corridor/barrier results.')

    # Return GEOPROCESSING specific errors
    except arcpy.ExecuteError:
        lu.dashline(1)
        gprint('****Iteration script failed. Details follow.****')
        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except Exception:
        lu.dashline(1)
        gprint('****Iteration script failed. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)


if __name__ == "__main__":
    main()
