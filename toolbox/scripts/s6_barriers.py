# Author: Brad McRae

"""Detect influential barriers given CWD calculations from Step 3.

Reguired Software:
ArcGIS Desktop 10.3+ or ArcGIS Pro with Spatial Analyst extension
Numpy
"""

from os import path
import time

import numpy as npy
import arcpy

from lm_config import tool_env as cfg
import lm_util as lu
from lm_util import gprint
from lm_retry_decorator import Retry


_SCRIPT_NAME = "s6_Barriers.py"

# Writing to tifs allows long filenames, needed for large radius values.
# Virtually no speed penalty in this case based on tests with large dataset.
TIF = '.tif'

# If search area overlaps a core, set output to NoData, which ends up as 0 in
# mosaic xxx
SET_CORES_TO_NULL = False


def step6_calc_barriers():
    """Detect influential barriers given CWD calculations from Step 3."""
    try:
        arcpy.CheckOutExtension("spatial")
        lu.dashline(0)
        gprint('Running script ' + _SCRIPT_NAME)

        if cfg.BARRIER_CWD_THRESH is not None:
            lu.dashline(1)
            gprint('Invoking CWD Threshold of ' + str(cfg.BARRIER_CWD_THRESH)
                   + ' map units.')

        if cfg.SUM_BARRIERS:
            sum_suffix = '_Sum'
            cfg.BARRIERBASEDIR = cfg.BARRIERBASEDIR + sum_suffix
            base_name, extension = path.splitext(cfg.BARRIERGDB)
            cfg.BARRIERGDB = base_name + sum_suffix + extension

            gprint('\nBarrier scores will be SUMMED across core pairs.')
        else:
            sum_suffix = ''

        if not arcpy.Exists(cfg.BARRIERGDB):
            # Create output geodatabase
            arcpy.CreateFileGDB_management(cfg.OUTPUTDIR,
                                           path.basename(cfg.BARRIERGDB))

        start_radius = int(cfg.STARTRADIUS)
        end_radius = int(cfg.ENDRADIUS)
        radius_step = int(cfg.RADIUSSTEP)
        if radius_step == 0:
            end_radius = start_radius  # Calculate at just one radius value
            radius_step = 1

        link_table_file = lu.get_prev_step_link_table(step=6)
        arcpy.env.workspace = cfg.SCRATCHDIR
        arcpy.env.scratchWorkspace = cfg.ARCSCRATCHDIR
        prefix = path.basename(cfg.PROJECTDIR)
        # For speed:
        arcpy.env.pyramid = "NONE"
        arcpy.env.rasterStatistics = "NONE"

        # set the analysis extent and cell size to that of the resistance
        # surface
        arcpy.env.extent = cfg.RESRAST
        arcpy.env.cellSize = cfg.RESRAST
        arcpy.env.snapRaster = cfg.RESRAST
        spatialref = arcpy.Describe(cfg.RESRAST).spatialReference
        map_units = (str(spatialref.linearUnitName)).lower()
        if len(map_units) > 1 and map_units[-1] != 's':
            map_units = map_units + 's'

        if (float(arcpy.env.cellSize) > start_radius
                or start_radius > end_radius):
            msg = ('Error: minimum detection radius must be greater than '
                   'cell size (' + str(arcpy.env.cellSize) +
                   ') \nand less than or equal to maximum detection radius.')
            lu.raise_error(msg)

        link_table = lu.load_link_table(link_table_file)
        num_links = link_table.shape[0]
        num_corridor_links = lu.report_links(link_table)
        if num_corridor_links == 0:
            lu.dashline(1)
            msg = '\nThere are no linkages. Bailing.'
            lu.raise_error(msg)

        # set up directories for barrier and barrier mosaic grids
        gprint("Creating intermediate output folder: " + cfg.BARRIERBASEDIR)
        lu.delete_dir(cfg.BARRIERBASEDIR)
        lu.create_dir(cfg.BARRIERBASEDIR)
        arcpy.CreateFolder_management(cfg.BARRIERBASEDIR, cfg.BARRIERDIR_NM)
        cbarrierdir = path.join(cfg.BARRIERBASEDIR, cfg.BARRIERDIR_NM)

        cores_to_process = npy.unique(
            link_table[:, cfg.LTB_CORE1:cfg.LTB_CORE2 + 1])
        max_core_num = max(cores_to_process)

        # Set up focal directories.
        # To keep there from being > 100 grids in any one directory,
        # outputs are written to:
        # barrier\focalX_ for cores 1-99 at radius X
        # barrier\focalX_1 for cores 100-199
        # etc.
        lu.dashline(0)

        for radius in range(start_radius, end_radius + 1, radius_step):
            core1path = lu.get_focal_path(1, radius)
            path1 = path.split(core1path)[0]
            path2, dir2 = path.split(path1)
            arcpy.CreateFolder_management(path.dirname(path2),
                                          path.basename(path2))
            arcpy.CreateFolder_management(path.dirname(path1),
                                          path.basename(path1))

            if max_core_num > 99:
                gprint('Creating subdirectories for ' + str(radius) + ' ' +
                       str(map_units) + ' radius analysis scale.')
                focal_dir_base_name = dir2

                cp100 = cores_to_process.astype('int32') / 100
                ind = npy.where(cp100 > 0)
                dir_nums = npy.unique(cp100[ind])
                for dir_num in dir_nums:
                    focal_dir = focal_dir_base_name + str(dir_num)
                    gprint('...' + focal_dir)
                    arcpy.CreateFolder_management(path2, focal_dir)

        # Create resistance raster with filled-in Nodata values for later use
        arcpy.env.extent = cfg.RESRAST
        resist_fill_ras = path.join(cfg.SCRATCHDIR, "resist_fill")
        output = arcpy.sa.Con(arcpy.sa.IsNull(cfg.RESRAST), 1000000000,
                              arcpy.sa.Raster(cfg.RESRAST) - 1)
        output.save(resist_fill_ras)

        core_list = link_table[:, cfg.LTB_CORE1:cfg.LTB_CORE2 + 1]
        core_list = npy.sort(core_list)

        # Loop through each search radius to calculate barriers in each link
        rad_id = 0  # Keep track of no of radii processed - used for temp dir
        for radius in range(start_radius, end_radius + 1, radius_step):
            rad_id = rad_id + 1
            link_table_tmp = link_table.copy()

            @Retry(10)
            # Can't pass vars in and modify them.
            def do_radius_loop():
                """Do radius loop."""
                link_table = link_table_tmp.copy()
                start_time = time.clock()
                link_loop = 0
                pct_done = 0
                gprint('\nMapping barriers at a radius of ' + str(radius) +
                       ' ' + str(map_units))
                if cfg.SUM_BARRIERS:
                    gprint('using SUM method')
                else:
                    gprint('using MAXIMUM method')
                if num_corridor_links > 1:
                    gprint('0 percent done')
                last_mosaic_ras = None
                last_mosaic_ras_pct = None
                for x in range(0, num_links):
                    pct_done = lu.report_pct_done(
                        link_loop, num_corridor_links, pct_done)
                    if ((link_table[x, cfg.LTB_LINKTYPE] > 0) and
                            (link_table[x, cfg.LTB_LINKTYPE] < 1000)):
                        link_loop = link_loop + 1
                        # source and target cores
                        corex = int(core_list[x, 0])
                        corey = int(core_list[x, 1])

                        # Get cwd rasters for source and target cores
                        cwd_ras1 = lu.get_cwd_path(corex)
                        cwd_ras2 = lu.get_cwd_path(corey)

                        # Mask out areas above CWD threshold
                        cwd_tmp1 = None
                        cwd_tmp2 = None
                        if cfg.BARRIER_CWD_THRESH is not None:
                            if x == 1:
                                lu.dashline(1)
                                gprint('  Using CWD threshold of '
                                       + str(cfg.BARRIER_CWD_THRESH)
                                       + ' map units.')
                            arcpy.env.extent = cfg.RESRAST
                            arcpy.env.cellSize = cfg.RESRAST
                            arcpy.env.snapRaster = cfg.RESRAST
                            cwd_tmp1 = path.join(cfg.SCRATCHDIR,
                                                 "tmp" + str(corex))
                            out_con = arcpy.sa.Con(
                                cwd_ras1 < float(cfg.BARRIER_CWD_THRESH),
                                cwd_ras1)
                            out_con.save(cwd_tmp1)
                            cwd_ras1 = cwd_tmp1
                            cwd_tmp2 = path.join(cfg.SCRATCHDIR,
                                                 "tmp" + str(corey))
                            out_con = arcpy.sa.Con(
                                cwd_ras2 < float(cfg.BARRIER_CWD_THRESH),
                                cwd_ras2)
                            out_con.save(cwd_tmp2)
                            cwd_ras2 = cwd_tmp2

                        focal_ras1 = lu.get_focal_path(corex, radius)
                        focal_ras2 = lu.get_focal_path(corey, radius)

                        link = lu.get_links_from_core_pairs(link_table,
                                                            corex, corey)
                        lc_dist = float(link_table[link, cfg.LTB_CWDIST])

                        # Detect barriers at radius using neighborhood stats
                        # Create the Neighborhood Object
                        inner_radius = radius - 1
                        outer_radius = radius

                        dia = 2 * radius
                        in_neighborhood = ("ANNULUS " + str(inner_radius)
                                           + " " + str(outer_radius) + " MAP")

                        @Retry(10)
                        def exec_focal():
                            """Execute focal statistics."""
                            if not path.exists(focal_ras1):
                                arcpy.env.extent = cwd_ras1
                                out_focal_stats = arcpy.sa.FocalStatistics(
                                    cwd_ras1, in_neighborhood,
                                    "MINIMUM", "DATA")
                                if SET_CORES_TO_NULL:
                                    # Set areas overlapping cores to NoData xxx
                                    out_focal_stats2 = arcpy.sa.Con(
                                        out_focal_stats > 0, out_focal_stats)
                                    out_focal_stats2.save(focal_ras1)
                                else:
                                    out_focal_stats.save(focal_ras1)
                                arcpy.env.extent = cfg.RESRAST

                            if not path.exists(focal_ras2):
                                arcpy.env.extent = cwd_ras2
                                out_focal_stats = arcpy.sa.FocalStatistics(
                                    cwd_ras2, in_neighborhood,
                                    "MINIMUM", "DATA")
                                if SET_CORES_TO_NULL:
                                    # Set areas overlapping cores to NoData xxx
                                    out_focal_stats2 = arcpy.sa.Con(
                                        out_focal_stats > 0, out_focal_stats)
                                    out_focal_stats2.save(focal_ras2)
                                else:
                                    out_focal_stats.save(focal_ras2)
                                arcpy.env.extent = cfg.RESRAST
                        exec_focal()

                        lu.delete_data(cwd_tmp1)
                        lu.delete_data(cwd_tmp2)

                        barrier_ras = path.join(
                            cbarrierdir, "b" + str(radius) + "_" + str(corex)
                            + "_" + str(corey)+'.tif')

                        # Need to set nulls to 0,
                        # also create trim rasters as we go
                        if cfg.SUM_BARRIERS:
                            out_ras = ((lc_dist - arcpy.sa.Raster(focal_ras1) -
                                        arcpy.sa.Raster(focal_ras2) - dia)
                                       / dia)
                            out_con = arcpy.sa.Con(arcpy.sa.IsNull(out_ras),
                                                   0, out_ras)
                            out_con2 = arcpy.sa.Con(out_con < 0, 0, out_con)
                            out_con2.save(barrier_ras)

                            # Execute FocalStatistics to fill out search radii
                            in_neighborhood = ("CIRCLE " + str(outer_radius)
                                               + " MAP")
                            fill_ras = path.join(
                                cbarrierdir, "b" + str(radius) + "_"
                                + str(corex) + "_" + str(corey) + "_fill.tif")
                            out_focal_stats = arcpy.sa.FocalStatistics(
                                barrier_ras, in_neighborhood,
                                "MAXIMUM", "DATA")
                            out_focal_stats.save(fill_ras)

                            if cfg.WRITE_TRIM_RASTERS:
                                trm_ras = path.join(
                                    cbarrierdir, "b" + str(radius) + "_"
                                    + str(corex) + "_" + str(corey)
                                    + "_trim.tif")
                                ras_list = [fill_ras, resist_fill_ras]
                                out_cell_statistics = arcpy.sa.CellStatistics(
                                    ras_list, "MINIMUM")
                                out_cell_statistics.save(trm_ras)

                        else:

                            @Retry(10)
                            def clac_ben():
                                """Calculate potential benefit.

                                Calculate potential benefit per map unit
                                restored.
                                """
                                out_ras = (
                                    (lc_dist - arcpy.sa.Raster(focal_ras1)
                                     - arcpy.sa.Raster(focal_ras2) - dia)
                                    / dia)
                                out_ras.save(barrier_ras)
                            clac_ben()

                        if cfg.WRITE_PCT_RASTERS:
                            # Calculate % potential benefit per unit restored
                            barrier_ras_pct = path.join(
                                cbarrierdir, "b" + str(radius) + "_"
                                + str(corex) + "_" + str(corey)
                                + '_pct.tif')

                            @Retry(10)
                            def calc_ben_pct():
                                """Calc benefit percentage."""
                                outras = (100 * (arcpy.sa.Raster(barrier_ras)
                                                 / lc_dist))
                                outras.save(barrier_ras_pct)
                            calc_ben_pct()

                        # Mosaic barrier results across core area pairs
                        mosaic_dir = path.join(cfg.SCRATCHDIR, 'mos'
                                               + str(rad_id) + '_'
                                               + str(x + 1))
                        lu.create_dir(mosaic_dir)

                        mos_fn = 'mos_temp'
                        tmp_mosaic_ras = path.join(mosaic_dir, mos_fn)
                        tmp_mosaic_ras_trim = path.join(mosaic_dir,
                                                        'mos_temp_trm')
                        arcpy.env.workspace = mosaic_dir
                        if link_loop == 1:
                            last_mosaic_ras_trim = None
                            # For first grid copy rather than mosaic
                            arcpy.CopyRaster_management(barrier_ras,
                                                        tmp_mosaic_ras)
                            if cfg.SUM_BARRIERS and cfg.WRITE_TRIM_RASTERS:
                                arcpy.CopyRaster_management(
                                    trm_ras, tmp_mosaic_ras_trim)
                        else:
                            if cfg.SUM_BARRIERS:
                                out_con = arcpy.sa.Con(
                                    arcpy.sa.Raster(barrier_ras) < 0,
                                    last_mosaic_ras,
                                    arcpy.sa.Raster(barrier_ras)
                                    + arcpy.sa.Raster(last_mosaic_ras))
                                out_con.save(tmp_mosaic_ras)
                                if cfg.WRITE_TRIM_RASTERS:
                                    out_con = arcpy.sa.Con(
                                        arcpy.sa.Raster(trm_ras) < 0,
                                        last_mosaic_ras_trim,
                                        arcpy.sa.Raster(trm_ras)
                                        + arcpy.sa.Raster(last_mosaic_ras_trim)
                                        )
                                    out_con.save(tmp_mosaic_ras_trim)

                            else:
                                in_rasters = (";".join([barrier_ras,
                                                        last_mosaic_ras]))

                                @Retry(10)
                                def mosaic_to_new():
                                    """Mosaic to new raster."""
                                    arcpy.MosaicToNewRaster_management(
                                        input_rasters=in_rasters,
                                        output_location=mosaic_dir,
                                        raster_dataset_name_with_extension\
                                        =mos_fn,
                                        pixel_type="32_BIT_FLOAT",
                                        cellsize=arcpy.env.cellSize,
                                        number_of_bands="1",
                                        mosaic_method="MAXIMUM",
                                        mosaic_colormap_mode="MATCH")
                                mosaic_to_new()

                        if link_loop > 1:  # Clean up from previous loop
                            lu.delete_data(last_mosaic_ras)
                            last_mosaic_dir = path.dirname(last_mosaic_ras)
                            lu.clean_out_workspace(last_mosaic_dir)
                            lu.delete_dir(last_mosaic_dir)

                        last_mosaic_ras = tmp_mosaic_ras
                        if cfg.WRITE_TRIM_RASTERS:
                            last_mosaic_ras_trim = tmp_mosaic_ras_trim
                        if cfg.WRITE_PCT_RASTERS:
                            mos_pct_fn = 'mos_temp_pct'
                            mosaic_dir_pct = path.join(cfg.SCRATCHDIR, 'mosP'
                                                       + str(rad_id) + '_'
                                                       + str(x+1))
                            lu.create_dir(mosaic_dir_pct)
                            tmp_mosaic_ras_pct = path.join(mosaic_dir_pct,
                                                           mos_pct_fn)
                            if link_loop == 1:
                                # If this is the first grid then copy
                                # rather than mosaic
                                if cfg.SUM_BARRIERS:
                                    out_con = arcpy.sa.Con(
                                        arcpy.sa.Raster(barrier_ras_pct)
                                        < 0, 0,
                                        arcpy.sa.Con(arcpy.sa.IsNull
                                                     (barrier_ras_pct),
                                                     0, barrier_ras_pct))
                                    out_con.save(tmp_mosaic_ras_pct)
                                else:
                                    arcpy.CopyRaster_management(
                                        barrier_ras_pct, tmp_mosaic_ras_pct)

                            else:
                                if cfg.SUM_BARRIERS:

                                    @Retry(10)
                                    def sum_barriers():
                                        """Sum barriers."""
                                        out_con = arcpy.sa.Con(
                                            arcpy.sa.Raster(barrier_ras_pct)
                                            < 0,
                                            last_mosaic_ras_pct,
                                            arcpy.sa.Raster(barrier_ras_pct)
                                            + arcpy.sa.Raster(
                                                last_mosaic_ras_pct))
                                        out_con.save(tmp_mosaic_ras_pct)
                                    sum_barriers()
                                else:
                                    in_rasters = (";".join([barrier_ras_pct,
                                                  last_mosaic_ras_pct]))

                                    @Retry(10)
                                    def max_barriers():
                                        """Get max barriers."""
                                        arcpy.MosaicToNewRaster_management(
                                            input_rasters=in_rasters,
                                            output_location=mosaic_dir_pct,
                                            raster_dataset_name_with_extension
                                            =mos_pct_fn,
                                            pixel_type="32_BIT_FLOAT",
                                            cellsize=arcpy.env.cellSize,
                                            number_of_bands="1",
                                            mosaic_method="MAXIMUM",
                                            mosaic_colormap_mode="MATCH")
                                    max_barriers()

                            if link_loop > 1:  # Clean up from previous loop
                                lu.delete_data(last_mosaic_ras_pct)
                                last_mosaic_dir_pct = path.dirname(
                                    last_mosaic_ras_pct)
                                lu.clean_out_workspace(last_mosaic_dir_pct)
                                lu.delete_dir(last_mosaic_dir_pct)

                            last_mosaic_ras_pct = tmp_mosaic_ras_pct

                        if not cfg.SAVEBARRIERRASTERS:
                            lu.delete_data(barrier_ras)
                            if cfg.WRITE_PCT_RASTERS:
                                lu.delete_data(barrier_ras_pct)
                            if cfg.WRITE_TRIM_RASTERS:
                                lu.delete_data(trm_ras)

                        # Temporarily disable links in linktable -
                        # don't want to mosaic them twice
                        for y in range(x + 1, num_links):
                            corex1 = int(core_list[y, 0])
                            corey1 = int(core_list[y, 1])
                            if corex1 == corex and corey1 == corey:
                                link_table[y, cfg.LTB_LINKTYPE] = (
                                    link_table[y, cfg.LTB_LINKTYPE] + 1000)
                            elif corex1 == corey and corey1 == corex:
                                link_table[y, cfg.LTB_LINKTYPE] = (
                                    link_table[y, cfg.LTB_LINKTYPE] + 1000)

                if num_corridor_links > 1 and pct_done < 100:
                    gprint('100 percent done')
                gprint('Summarizing barrier data for search radius.')
                # Rows that were temporarily disabled
                rows = npy.where(link_table[:, cfg.LTB_LINKTYPE] > 1000)
                link_table[rows, cfg.LTB_LINKTYPE] = (
                    link_table[rows, cfg.LTB_LINKTYPE] - 1000)
                # -----------------------------------------------------------------
                # Set negative values to null or zero and write geodatabase.
                mosaic_fn = (prefix + "_BarrierCenters" + sum_suffix + "_Rad" +
                             str(radius))
                mosaic_ras = path.join(cfg.BARRIERGDB, mosaic_fn)
                arcpy.env.extent = cfg.RESRAST

                out_set_null = arcpy.sa.SetNull(tmp_mosaic_ras,
                                                tmp_mosaic_ras,
                                                "VALUE < 0")  # xxx orig
                out_set_null.save(mosaic_ras)

                lu.delete_data(tmp_mosaic_ras)

                if cfg.SUM_BARRIERS and cfg.WRITE_TRIM_RASTERS:
                    mosaic_fn = (prefix + "_BarrierCircles_RBMin" + sum_suffix
                                 + "_Rad" + str(radius))
                    mosaic_ras_trim = path.join(cfg.BARRIERGDB, mosaic_fn)
                    arcpy.CopyRaster_management(tmp_mosaic_ras_trim,
                                                mosaic_ras_trim)
                    lu.delete_data(tmp_mosaic_ras)

                if cfg.WRITE_PCT_RASTERS:
                    # Do same for percent raster
                    mosaic_pct_fn = (prefix + "_BarrierCenters_Pct"
                                     + sum_suffix + "_Rad" + str(radius))
                    arcpy.env.extent = cfg.RESRAST
                    out_set_null = arcpy.sa.SetNull(tmp_mosaic_ras_pct,
                                                    tmp_mosaic_ras_pct,
                                                    "VALUE < 0")
                    mosaic_ras_pct = path.join(cfg.BARRIERGDB, mosaic_pct_fn)
                    out_set_null.save(mosaic_ras_pct)
                    lu.delete_data(tmp_mosaic_ras_pct)

                # 'Grow out' maximum restoration gain to
                # neighborhood size for display
                in_neighborhood = "CIRCLE " + str(outer_radius) + " MAP"
                # Execute FocalStatistics
                fill_ras_fn = "barriers_fill" + str(outer_radius) + TIF
                fill_ras = path.join(cfg.BARRIERBASEDIR, fill_ras_fn)
                out_focal_stats = arcpy.sa.FocalStatistics(
                    mosaic_ras, in_neighborhood, "MAXIMUM", "DATA")
                out_focal_stats.save(fill_ras)

                if cfg.WRITE_PCT_RASTERS:
                    # Do same for percent raster
                    fill_ras_pct_fn = (
                        "barriers_fill_pct" + str(outer_radius) + TIF)
                    fill_ras_pct = path.join(cfg.BARRIERBASEDIR,
                                             fill_ras_pct_fn)
                    out_focal_stats = arcpy.sa.FocalStatistics(
                        mosaic_ras_pct, in_neighborhood, "MAXIMUM", "DATA")
                    out_focal_stats.save(fill_ras_pct)

                # Place copies of filled rasters in output geodatabase
                arcpy.env.workspace = cfg.BARRIERGDB
                fill_ras_fn = (prefix + "_BarrrierCircles" + sum_suffix
                               + "_Rad" + str(outer_radius))
                arcpy.CopyRaster_management(fill_ras, fill_ras_fn)
                if cfg.WRITE_PCT_RASTERS:
                    fill_ras_pct_fn = (prefix + "_BarrrierCircles_Pct"
                                       + sum_suffix + "_Rad"
                                       + str(outer_radius))
                    arcpy.CopyRaster_management(fill_ras_pct,
                                                fill_ras_pct_fn)

                if not cfg.SUM_BARRIERS and cfg.WRITE_TRIM_RASTERS:
                    # Create pared-down version of filled raster- remove pixels
                    # that don't need restoring by allowing a pixel to only
                    # contribute its resistance value to restoration gain
                    out_ras_fn = "barriers_trm" + str(outer_radius) + TIF
                    out_ras = path.join(cfg.BARRIERBASEDIR, out_ras_fn)
                    ras_list = [fill_ras, resist_fill_ras]
                    out_cell_statistics = arcpy.sa.CellStatistics(ras_list,
                                                                  "MINIMUM")
                    out_cell_statistics.save(out_ras)

                    # SECOND ROUND TO CLIP BY DATA VALUES IN BARRIER RASTER
                    out_ras_2fn = ("barriers_trm" + sum_suffix
                                   + str(outer_radius) + "_2" + TIF)
                    out_ras2 = path.join(cfg.BARRIERBASEDIR, out_ras_2fn)
                    output = arcpy.sa.Con(arcpy.sa.IsNull(fill_ras),
                                          fill_ras, out_ras)
                    output.save(out_ras2)
                    out_ras_fn = (prefix + "_BarrierCircles_RBMin"
                                  + sum_suffix + "_Rad"
                                  + str(outer_radius))
                    arcpy.CopyRaster_management(out_ras2, out_ras_fn)
                start_time = lu.elapsed_time(start_time)

            # Call the above function
            do_radius_loop()

        # Combine rasters across radii
        gprint('\nCreating summary rasters...')
        if start_radius != end_radius:
            radii_suffix = ('_Rad' + str(int(start_radius)) + 'To'
                            + str(int(end_radius)) + 'Step'
                            + str(int(radius_step)))
            mosaic_fn = "bar_radii"
            mosaic_pct_fn = "bar_radii_pct"
            arcpy.env.workspace = cfg.BARRIERBASEDIR
            for radius in range(start_radius, end_radius + 1, radius_step):
                # Fixme: run speed test with gdb mosaicking above and here
                radius_fn = (prefix + "_BarrierCenters" + sum_suffix + "_Rad"
                             + str(radius))
                radius_ras = path.join(cfg.BARRIERGDB, radius_fn)

                if radius == start_radius:
                    # If this is the first grid then copy rather than mosaic
                    arcpy.CopyRaster_management(radius_ras, mosaic_fn)
                else:
                    mosaic_ras = path.join(cfg.BARRIERBASEDIR, mosaic_fn)
                    arcpy.Mosaic_management(radius_ras, mosaic_ras,
                                            "MAXIMUM", "MATCH")

                if cfg.WRITE_PCT_RASTERS:
                    radius_pct_fn = (prefix + "_BarrierCenters_Pct"
                                     + sum_suffix + "_Rad" + str(radius))
                    radius_ras_pct = path.join(cfg.BARRIERGDB, radius_pct_fn)

                    if radius == start_radius:
                        # If this is the first grid then copy rather than
                        # mosaic
                        arcpy.CopyRaster_management(radius_ras_pct,
                                                    mosaic_pct_fn)
                    else:
                        mosaic_ras_pct = path.join(cfg.BARRIERBASEDIR,
                                                   mosaic_pct_fn)
                        arcpy.Mosaic_management(radius_ras_pct,
                                                mosaic_ras_pct,
                                                "MAXIMUM", "MATCH")

            # Copy results to output geodatabase
            arcpy.env.workspace = cfg.BARRIERGDB
            mosaic_fn = prefix + "_BarrierCenters" + sum_suffix + radii_suffix
            arcpy.CopyRaster_management(mosaic_ras, mosaic_fn)

            if cfg.WRITE_PCT_RASTERS:
                mosaic_pct_fn = (prefix + "_BarrierCenters_Pct" + sum_suffix +
                                 radii_suffix)
                arcpy.CopyRaster_management(mosaic_ras_pct, mosaic_pct_fn)

            # GROWN OUT rasters
            fill_mosaic_fn = "barriers_radii_fill" + TIF
            fill_mosaic_pct_fn = "barriers_radii_fill_pct" + TIF
            fill_mosaic_ras = path.join(cfg.BARRIERBASEDIR, fill_mosaic_fn)
            trim_mosaic_ras_pct = path.join(cfg.BARRIERBASEDIR,
                                            fill_mosaic_pct_fn)

            arcpy.env.workspace = cfg.BARRIERBASEDIR
            for radius in range(start_radius, end_radius + 1, radius_step):
                radius_fn = "barriers_fill" + str(radius) + TIF
                # fixme- do this when only a single radius too
                radius_ras = path.join(cfg.BARRIERBASEDIR, radius_fn)
                if radius == start_radius:
                    # If this is the first grid then copy rather than mosaic
                    arcpy.CopyRaster_management(radius_ras, fill_mosaic_fn)
                else:
                    arcpy.Mosaic_management(radius_ras, fill_mosaic_ras,
                                            "MAXIMUM", "MATCH")

                if cfg.WRITE_PCT_RASTERS:
                    radius_pct_fn = "barriers_fill_pct" + str(radius) + TIF
                    # fixme- do this when only a single radius too
                    radius_ras_pct = path.join(cfg.BARRIERBASEDIR,
                                               radius_pct_fn)
                    if radius == start_radius:
                        # For first grid copy rather than mosaic
                        arcpy.CopyRaster_management(radius_ras_pct,
                                                    fill_mosaic_pct_fn)
                    else:
                        arcpy.Mosaic_management(radius_ras_pct,
                                                trim_mosaic_ras_pct,
                                                "MAXIMUM", "MATCH")

            # Copy result to output geodatabase
            arcpy.env.workspace = cfg.BARRIERGDB
            fill_mosaic_fn = (prefix + "_BarrierCircles" + sum_suffix
                              + radii_suffix)
            arcpy.CopyRaster_management(fill_mosaic_ras, fill_mosaic_fn)
            if cfg.WRITE_PCT_RASTERS:
                fill_mosaic_pct_fn = (prefix + "_BarrierCircles_Pct"
                                      + sum_suffix + radii_suffix)
                arcpy.CopyRaster_management(trim_mosaic_ras_pct,
                                            fill_mosaic_pct_fn)

            # GROWN OUT AND TRIMMED rasters (Can't do percent)
            if cfg.WRITE_TRIM_RASTERS:
                trim_mosaic_fn = "bar_radii_trm"
                arcpy.env.workspace = cfg.BARRIERBASEDIR
                trim_mosaic_ras = path.join(cfg.BARRIERBASEDIR, trim_mosaic_fn)
                for radius in range(start_radius, end_radius + 1, radius_step):
                    radius_fn = (prefix + "_BarrierCircles_RBMin" + sum_suffix
                                 + "_Rad" + str(radius))
                    # fixme- do this when only a single radius too
                    radius_ras = path.join(cfg.BARRIERGDB, radius_fn)

                    if radius == start_radius:
                        # For first grid copy rather than mosaic
                        arcpy.CopyRaster_management(radius_ras, trim_mosaic_fn)
                    else:
                        arcpy.Mosaic_management(radius_ras, trim_mosaic_ras,
                                                "MAXIMUM", "MATCH")
                # Copy result to output geodatabase
                arcpy.env.workspace = cfg.BARRIERGDB
                trim_mosaic_fn = (prefix + "_BarrierCircles_RBMin" + sum_suffix
                                  + radii_suffix)
                arcpy.CopyRaster_management(trim_mosaic_ras, trim_mosaic_fn)

        if not cfg.SAVE_RADIUS_RASTERS:
            arcpy.env.workspace = cfg.BARRIERGDB
            rasters = arcpy.ListRasters()
            for raster in rasters:
                if 'rad' in raster.lower() and 'step' not in raster.lower():
                    lu.delete_data(raster)

        arcpy.env.workspace = cfg.BARRIERGDB
        rasters = arcpy.ListRasters()
        for raster in rasters:
            gprint('\nBuilding output statistics and pyramids\n'
                   'for raster ' + raster)
            lu.build_stats(raster)

        # Clean up temporary files and directories
        if not cfg.SAVEBARRIERRASTERS:
            lu.delete_dir(cbarrierdir)
            lu.delete_dir(cfg.BARRIERBASEDIR)

        if not cfg.SAVEFOCALRASTERS:
            for radius in range(start_radius, end_radius + 1, radius_step):
                core1path = lu.get_focal_path(1, radius)
                path1 = path.split(core1path)[0]
                path2 = path.split(path1)[0]
                lu.delete_dir(path2)

    # Return GEOPROCESSING specific errors
    except arcpy.ExecuteError:
        lu.dashline(1)
        gprint('****Failed in step 6. Details follow.****')
        lu.exit_with_geoproc_error(_SCRIPT_NAME)

    # Return any PYTHON or system specific errors
    except Exception:
        lu.dashline(1)
        gprint('****Failed in step 6. Details follow.****')
        lu.exit_with_python_error(_SCRIPT_NAME)

    return
