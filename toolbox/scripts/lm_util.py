# Authors: Brad McRae and Darren Kavanagh

"""Contains functions called by linkage mapper and barrier mapper scripts."""

import os
import sys
import subprocess
from datetime import datetime as dt
import time
import traceback
import ConfigParser
import shutil
import gc
import ctypes
import locale
from lm_retry_decorator import Retry


import numpy as npy
import arcpy

from lm_config import tool_env as cfg
try:
    test = cfg.releaseNum
except Exception:
    cfg.releaseNum = 'unknown'

_SCRIPT_NAME = "lm_util.py"


def cwd_cutoff_str(cutoff):
    """Convert CDW cutoff to text and abbreviate if possible."""
    cutoff_text = str(cutoff)
    if cutoff_text[-6:] == "000000":
        cutoff_text = cutoff_text[0:-6] + "m"
    elif cutoff_text[-3:] == "000":
        cutoff_text = cutoff_text[0:-3] + "k"
    return cutoff_text


def get_linktable_row(linkid, linktable):
    """Returns the linkTable row index for a given link ID"""
    try:
        # Most likely.  linkTables tend to be in order with no skipped links.
        if linktable[linkid - 1, cfg.LTB_LINKID] == linkid:
            linktablerow = linkid - 1
            return linktablerow
        else:
            numLinks = linktable.shape[0]
            for linktablerow in range(0, numLinks):
                if linktable[linktablerow, cfg.LTB_LINKID] == linkid:
                    return linktablerow
        return -1  # Not found
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)


def get_link_type_desc(linktypecode):
    """For a linkType code returns description to attribute link and link maps.

    NOTE: These must map to LT codes in lm_config (eg LT_CPLK)

    """
    if linktypecode < 0:  # These are dropped links
        activelink = '0'
        if linktypecode == -1:
            linktypedesc = '"Not_nearest_neighbors"'
        elif linktypecode == -11:
            linktypedesc = '"Too_long_Euclidean_distance"'
        elif linktypecode == -12:
            linktypedesc = '"Too_long_CWD"'
        elif linktypecode == -13:
            linktypedesc = '"Too_short_Euclidean_distance"'
        elif linktypecode == -14:
            linktypedesc = '"Too_short_CWD"'
        elif linktypecode == -15:
            linktypedesc = '"Intermediate_core"'
        elif linktypecode == -100:
            linktypedesc = '"User_removed"'
        else:
            linktypedesc = '"Unknown_inactive"'
    else:  # These are active links
        activelink = '1'
        if linktypecode == 1:
            linktypedesc = '"Connects_cores"'
        elif linktypecode == 10:
            linktypedesc = '"Connects_ nearest_neighbors"'
        elif linktypecode == 20:
            linktypedesc = '"Connects_constellations"'
        elif linktypecode == 100:
            linktypedesc = '"User_retained"'
        else:
            linktypedesc = '"Unknown_active"'

    return activelink, linktypedesc


def get_links_from_core_pairs(linktable, firstCore, secondCore):
    """Given two cores, finds their matching row in the link table"""
    try:
        rows = npy.zeros((0), dtype="int32")
        numLinks = linktable.shape[0]
        for link in range(0, numLinks):
            corex = int(linktable[link, cfg.LTB_CORE1])
            corey = int(linktable[link, cfg.LTB_CORE2])
            if int(corex) == int(firstCore) and int(corey) == int(secondCore):
                rows = npy.append(rows, link)
            elif (int(corex) == int(secondCore) and int(corey) ==
                  int(firstCore)):
                rows = npy.append(rows, link)
        return rows
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)


def drop_links(linktable, maxeud, mineud, maxcwd, mincwd,
               DISABLE_LEAST_COST_NO_VAL):
    """Inactivates links that fail to meet min or max length criteria"""
    try:
        numLinks = linktable.shape[0]
        numDroppedLinks = 0
        coreList = linktable[:, cfg.LTB_CORE1:cfg.LTB_CORE2 + 1]
        coreList = npy.sort(coreList)
        if DISABLE_LEAST_COST_NO_VAL:
            for x in range(0, numLinks):
                linkid = str(int(linktable[x, cfg.LTB_LINKID]))
                if linktable[x, cfg.LTB_CWDIST] == -1:
                    #Check only enabled corridor links
                    if (linktable[x, cfg.LTB_LINKTYPE] > 0):
                        corex = str(int(linktable[x, cfg.LTB_CORE1]))
                        corey = str(int(linktable[x, cfg.LTB_CORE2]))
                        gprint(
                            "The least-cost corridor between " + str(corex) +
                            " and " + str(corey) + " (link #" + linkid + ") "
                            "has an unknown length in cost distance units. "
                            "This means it is longer than the max "
                            "cost-weighted distance specified in a previous "
                            "step OR it passes through NODATA cells "
                            "and will be dropped.\n")
                        # Disable link
                        linktable[x, cfg.LTB_LINKTYPE] = cfg.LT_TLLC
                        numDroppedLinks = numDroppedLinks + 1

        # Check for corridors that are too long in Euclidean or cost-weighted
        # distance
        if (maxeud is not None) or (maxcwd is not None):
            for x in range(0, numLinks):
                linkid = str(int(linktable[x, cfg.LTB_LINKID]))
                if maxeud is not None:
                    if linktable[x, cfg.LTB_EUCDIST] > maxeud:
                        # Check only enabled corridor links
                        if (linktable[x, cfg.LTB_LINKTYPE] > 0) and (
                            linktable[x, cfg.LTB_LINKTYPE] != cfg.LT_KEEP):
                            corex = str(int(coreList[x, 0]))
                            corey = str(int(coreList[x, 1]))
                            gprint("Link #" + linkid +
                                          " connecting cores " + str(corex) +
                                          " and " + str(corey) + " is  " +
                                          str(linktable[x, cfg.LTB_EUCDIST]) +
                                          " units long- too long in "
                                          "Euclidean distance.")
                            # Disable link
                            linktable[x, cfg.LTB_LINKTYPE] = cfg.LT_TLEC
                            numDroppedLinks = numDroppedLinks + 1

                if maxcwd is not None:
                    # Check only enabled corridor links
                    if (linktable[x, cfg.LTB_LINKTYPE] > 0) and (
                        linktable[x, cfg.LTB_LINKTYPE] != cfg.LT_KEEP):
                        if (linktable[x, cfg.LTB_CWDIST] > maxcwd):
                            corex = str(int(linktable[x, cfg.LTB_CORE1]))
                            corey = str(int(linktable[x, cfg.LTB_CORE2]))
                            gprint(
                                "Link #" + linkid + " connecting cores " +
                                str(corex) + " and " + str(corey) +
                                " is " +
                                str(linktable[x, cfg.LTB_CWDIST]) +
                                " units long- too long in cost-distance "
                                "units.")
                            #  Disable link
                            linktable[x, cfg.LTB_LINKTYPE] = cfg.LT_TLLC
                            numDroppedLinks = numDroppedLinks + 1

        if mineud is not None or mincwd is not None:
            for x in range(0, numLinks):
                linkid = str(int(linktable[x, cfg.LTB_LINKID]))
                # Check only enabled corridor links
                if (linktable[x, cfg.LTB_LINKTYPE] > 0):
                    if mineud is not None:
                        if (linktable[x, cfg.LTB_EUCDIST] < mineud) and (
                            linktable[x, cfg.LTB_LINKTYPE] != cfg.LT_KEEP):
                            corex = str(int(coreList[x, 0]))
                            corey = str(int(coreList[x, 1]))
                            gprint(
                                "Link #" + linkid + " connecting cores " +
                                str(corex) + " and " + str(corey) + " is "
                                "only " + str(linktable[x, cfg.LTB_EUCDIST]) +
                                " units long- too short in Euclidean "
                                "distance.")
                            # Disable link
                            linktable[x, cfg.LTB_LINKTYPE] = cfg.LT_TSEC
                            numDroppedLinks = numDroppedLinks + 1

                    if mincwd is not None:
                        if ((linktable[x, cfg.LTB_CWDIST] < mincwd) and
                            (linktable[x, cfg.LTB_CWDIST]) != -1) and (
                            linktable[x, cfg.LTB_LINKTYPE] != cfg.LT_KEEP):
                            if (linktable[x, cfg.LTB_LINKTYPE] > 0):
                                corex = str(int(linktable[x, cfg.LTB_CORE1]))
                                corey = str(int(linktable[x, cfg.LTB_CORE2]))
                                gprint(
                                    "Link #" + linkid + " connecting cores " +
                                    str(corex) + " and " + str(corey) +
                                    " is only " +
                                    str(linktable[x, cfg.LTB_CWDIST]) +
                                    " units long- too short in cost distance "
                                    "units.")
                                # Disable link
                                linktable[x, cfg.LTB_LINKTYPE] = cfg.LT_TSLC
                                numDroppedLinks = numDroppedLinks + 1
        return linktable, numDroppedLinks
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)



def get_core_list(coreFC, coreFN):
    """Returns a list of core area IDs from polygon file"""
    try:

        # This returns a list with length equal to number of shapes, not unique
        # core id's.

        # Get the number of core shapes
        arcpy.env.extent = arcpy.Describe(coreFC).Extent

        # Get core data into numpy array
        coreList = npy.zeros((1, 2))
        cur = arcpy.SearchCursor(coreFC)
        row = cur.next()
        i = 0
        while row:
            if i > 0:
                coreList = npy.append(coreList,  npy.zeros((1, 2)), axis=0)
            coreList[i, 0] = row.getValue(coreFN)
            coreList[i, 1] = row.getValue(coreFN)
            row = cur.next()
            i = i + 1

        del cur, row
        coreCount=coreList.shape[0]
        if coreCount < 2:
            dashline(1)
            msg =('\nERROR: Less than two core areas detected. This can '
                   '\nhappen if you have selected a core area in ArcMap. It '
                   '\ncan also happen when resistance and core area maps are '
                   '\nmissing spatial reference data or are in different '
                   '\nprojections. '
                   '\nBailing because there is nothing to connect.')
            raise_error(msg)

        arcpy.env.extent = "MAXOF"  # For downstream operations
        return coreList

    except arcpy.ExecuteError:
        exit_with_geoproc_error(_SCRIPT_NAME)
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)



def get_core_targets(core, linktable):
    """Returns a list of other core areas the core area is connected to."""
    try:
        targetList = npy.zeros((len(linktable), 2), dtype="int32")
        # possible targets. 1st column = cfg.LTB_CORE1
        targetList[:, 0] = linktable[:, cfg.LTB_CORE1]
        # possible targets. 2nd column = cfg.LTB_CORE2
        targetList[:, 1] = linktable[:, cfg.LTB_CORE2]
        # Copy of cfg.LTB_LINKTYPE column
        validPair = linktable[:, cfg.LTB_LINKTYPE]
        validPair = npy.where((validPair == cfg.LT_KEEP), cfg.LT_CORR,
                               validPair)
        validPair = npy.where((validPair == cfg.LT_CORR), 1, 0)  # map corridor.
        targetList[:, 0] = npy.multiply(targetList[:, 0], validPair)
        targetList[:, 1] = npy.multiply(targetList[:, 1], validPair)

        rows, cols = npy.where(targetList == int(core))
        targetList = targetList[rows, 1 - cols]
        targetList = npy.unique(npy.asarray(targetList))
    except arcpy.ExecuteError:
        exit_with_geoproc_error(_SCRIPT_NAME)
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)
    return targetList


def start_time():
    """Print and return start time."""
    start_time = dt.now()
    print "Start time: {}".format(start_time.strftime("%m/%d/%y %H:%M:%S"))
    return start_time


def s2hhmmss(s):
    """Convert seconds to hours, minutes and seconds."""
    m, s = divmod(int(s), 60)
    h, m = divmod(m, 60)
    return h, m, s


def run_time(stime):
    """Print program execution time when running from script."""
    etime = dt.now()
    hours, minutes, seconds = s2hhmmss((etime - stime).total_seconds())
    print "End time: {}".format(etime.strftime("%m/%d/%y %H:%M:%S"))
    print "Execution time: {} hrs {} mins {} seconds".format(
        hours, minutes, seconds)


def elapsed_time(start_time):
    """Print elapsed time given a start time and return a new start time."""
    now = time.clock()
    hours, minutes, seconds = s2hhmmss(now - start_time)
    msg = "Task took:"
    if minutes == 0:
        gprint("{} {} seconds.\n".format(msg, seconds))
    elif hours == 0:
        gprint("{} {} minutes and {} seconds.\n".format(
            msg, minutes, seconds))
    else:
        gprint("{} {} hours and {} minutes and {} seconds.\n".format(
            msg, hours, minutes, seconds))
    return now


def report_pct_done(current, goal, last):
    """Reports percent done"""
    try:
        goal = float(goal)
        pctDone = ((float(current) / goal) * 100)
        pctDone = 10 * (npy.floor(pctDone/10))
        if pctDone - last >= 10:
            gprint(str(int(pctDone)) + " percent done")
            return 10*int((npy.floor(pctDone/10)))
        else:
            return last
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)


def report_links(linktable):
    """Prints number of links in a link table"""
    try:
        numLinks = linktable.shape[0]
        gprint('There are ' + str(numLinks) + ' links in the '
                          'table.')
        linkTypes = linktable[:, cfg.LTB_LINKTYPE]
        numCorridorLinks = sum(linkTypes == cfg.LT_CORR) + sum(linkTypes ==
                           cfg.LT_NNC) + sum(linkTypes == cfg.LT_KEEP)
        numComponentLinks = sum(linkTypes == cfg.LT_CLU)
        if numComponentLinks > 0:
            gprint('This includes ' + str(numCorridorLinks) +
                              ' potential corridor links and ' +
                          str(numComponentLinks) + ' component links.')
        elif numCorridorLinks > 0:
            gprint('This includes ' + str(numCorridorLinks) +
                              ' potential corridor links.')
        else:
            numCorridorLinks = 0
            gprint('\n***NOTE: There are NO corridors to map!')
            dashline(2)

    except arcpy.ExecuteError:
        exit_with_geoproc_error(_SCRIPT_NAME)
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)
    return numCorridorLinks + numComponentLinks


def build_stats(raster):
    """Builds statistics and pyramids for output rasters"""
    try:
        arcpy.CalculateStatistics_management(raster, "1", "1", "#")
    except Exception:
        gprint('Statistics failed. They can still be calculated manually.')
    try:
        arcpy.BuildPyramids_management(raster)
    except Exception:
        gprint('Pyramids failed. They can still be built manually.')
    return



############################################################################
## Adjacency and allocation functions ##########################
############################################################################

def get_adj_using_shift_method(alloc):
    """Returns table listing adjacent core areas using a shift method.

    The method involves shifting the allocation grid one pixel and then looking
    for pixels with different allocations across shifted grids.

    """
    cellSize = arcpy.Describe(alloc).MeanCellHeight
    arcpy.env.cellSize = cellSize

    posShift = arcpy.env.cellSize
    negShift = -1 * float(arcpy.env.cellSize)

    arcpy.env.workspace = cfg.SCRATCHDIR

    gprint('Calculating adjacencies crossing allocation boundaries...')
    start_time = time.clock()
    arcpy.Shift_management(alloc, "alloc_r", posShift, "0")

    alloc_r = "alloc_r"
    adjTable_r = get_allocs_from_shift(arcpy.env.workspace, alloc, alloc_r)
    arcpy.Shift_management(alloc, "alloc_ul", negShift, posShift)

    alloc_ul = "alloc_ul"
    adjTable_ul = get_allocs_from_shift(arcpy.env.workspace, alloc, alloc_ul)
    arcpy.Shift_management(alloc, "alloc_ur", posShift, posShift)

    alloc_ur = "alloc_ur"
    adjTable_ur = get_allocs_from_shift(arcpy.env.workspace, alloc, alloc_ur)
    arcpy.Shift_management(alloc, "alloc_u", "0", posShift)

    alloc_u = "alloc_u"
    adjTable_u = get_allocs_from_shift(arcpy.env.workspace, alloc, alloc_u)
    start_time = elapsed_time(start_time)

    adjTable = combine_adjacency_tables(adjTable_r, adjTable_u, adjTable_ur,
                                        adjTable_ul)

    return adjTable


def combine_adjacency_tables(adjTable_r, adjTable_u, adjTable_ur, adjTable_ul):
    """Combines tables describing whether core areas are adjacent based on
    allocation zones that touch on horizontal, vertical, and diagonal axes
    """
    try:
        adjTable = npy.append(adjTable_r, adjTable_u, axis=0)
        adjTable = npy.append(adjTable, adjTable_ur, axis=0)
        adjTable = npy.append(adjTable, adjTable_ul, axis=0)

        pairs = npy.sort(adjTable[:, 0:2])
        adjTable[:, 0:2] = pairs

        # sort by 1st core Id then by 2nd core Id
        ind = npy.lexsort((adjTable[:, 1], adjTable[:, 0]))
        adjTable = adjTable[ind]

        numDists = len(adjTable)
        x = 1
        while x < numDists:
            if (adjTable[x, 0] == adjTable[x - 1, 0] and
                adjTable[x, 1] == adjTable[x - 1, 1]):
                adjTable[x - 1, 0] = 0  # mark for deletion
            x = x + 1

        if numDists > 0:
            delRows = npy.asarray(npy.where(adjTable[:, 0] == 0))
            delRowsVector = npy.zeros((delRows.shape[1]), dtype="int32")
            delRowsVector[:] = delRows[0, :]
            adjTable = delete_row(adjTable, delRowsVector)
            del delRows
            del delRowsVector

        return adjTable
    except arcpy.ExecuteError:
        exit_with_geoproc_error(_SCRIPT_NAME)
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)


def get_allocs_from_shift(workspace, alloc, alloc_sh):
    """Returns a table of adjacent allocation zones using grid shift method"""
    try:
        arcpy.env.scratchWorkspace = cfg.ARCSCRATCHDIR
        arcpy.env.workspace = workspace
        combine_ras = os.path.join(arcpy.env.workspace, "combine")
        count = 0
        statement = ('comb_ras = arcpy.sa.Combine([alloc, alloc_sh]); '
                     'comb_ras.save(combine_ras)')
        while True:
            try:
                exec(statement)
            except Exception:
                count, tryAgain = retry_arc_error(count, statement)
                if not tryAgain:
                    exec(statement)
            else:
                break
        allocLookupTable = get_alloc_lookup_table(arcpy.env.workspace,
                                                  combine_ras)
        # Overwrite setting does not work for Combine, so delete raster
        delete_data(combine_ras)
        return allocLookupTable[:, 1:3]

    except arcpy.ExecuteError:
        exit_with_geoproc_error(_SCRIPT_NAME)
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)


# Try this in script 4 too... and anything else with minras...
def get_alloc_lookup_table(workspace, combine_ras):
    """Returns a table of adjacent allocation zones.

    Requires a raster with allocation zone attributes.

    """
    try:
        fldlist = arcpy.ListFields(combine_ras)

        valFld = fldlist[1].name
        allocFld = fldlist[3].name
        allocFld_sh = fldlist[4].name

        allocLookupTable = npy.zeros((0, 3), dtype="int32")
        appendRow = npy.zeros((1, 3), dtype="int32")

        rows = arcpy.SearchCursor(combine_ras)
        row = rows.next()
        while row:
            alloc = row.getValue(allocFld)
            alloc_sh = row.getValue(allocFld_sh)
            if alloc != alloc_sh:
                appendRow[0, 0] = row.getValue(valFld)
                appendRow[0, 1] = alloc
                appendRow[0, 2] = alloc_sh
                allocLookupTable = npy.append(allocLookupTable, appendRow,
                                              axis=0)
            row = rows.next()
        del row
        del rows

        return allocLookupTable
    except arcpy.ExecuteError:
        exit_with_geoproc_error(_SCRIPT_NAME)
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)


############################################################################
## Bounding Circle Functions ##########################
############################################################################
def get_centroids(shapefile, field):
    """Returns centroids of features"""
    try:
        pointArray = npy.zeros((0, 3), dtype="float32")
        xyCumArray = npy.zeros((0, 3), dtype="float32")
        xyArray = npy.zeros((1, 3), dtype="float32")
        rows = arcpy.SearchCursor(shapefile)
        row = rows.next()
        while row:
            feat = row.shape
            center = feat.centroid
            center = str(center)
            xy = center.split(" ")
            if "," in xy[0] or "," in xy[1]:
                msg = ('ERROR: It appears that your region settings are not in '
                        'USA format (decimal commas are used instead of decimal '
                        'points). '
                    'Please change your region settings in Windows to English (USA) or '
                    'another convention that uses decimal points.  You may need to'
                    'modify the coordinate system of your input files as well.')
                raise_error(msg)
            xyArray[0, 0] = float(xy[0])
            xyArray[0, 1] = float(xy[1])
            value = row.getValue(field)
            xyArray[0, 2] = int(value)
            xyCumArray = npy.append(xyCumArray, xyArray, axis=0)
            row = rows.next()
        del row, rows
        pointArray = npy.append(pointArray, xyCumArray, axis=0)

        return pointArray

    except arcpy.ExecuteError:
        exit_with_geoproc_error(_SCRIPT_NAME)
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)


def get_bounding_circle_data(extentBoxList, corex, corey, distbuff):
    """Returns centroid and radius of circles bounding the extent boxes."""
    try:
        circlePointData = npy.zeros((1, 5), dtype='float32')
        numBoxes = extentBoxList.shape[0]
        # doing loop because can't get where to work correctly on this list
        for line in range(numBoxes):
            if extentBoxList[line, 0] == corex:
                corexData = extentBoxList[line, :]
            if extentBoxList[line, 0] == corey:
                coreyData = extentBoxList[line, :]
                break

        x_ulx = corexData[1]
        x_lrx = corexData[2]
        x_uly = corexData[3]
        x_lry = corexData[4]
        y_ulx = coreyData[1]
        y_lrx = coreyData[2]
        y_uly = coreyData[3]
        y_lry = coreyData[4]

        xmin = min(x_ulx, y_ulx)
        ymin = min(x_lry, y_lry)
        xmax = max(x_lrx, y_lrx)
        ymax = max(x_uly, y_uly)

        centX = xmin + (xmax - xmin) / 2
        centY = ymin + (ymax - ymin) / 2
        radius = (npy.sqrt(pow((xmax - xmin) / 2, 2) +
                  pow((ymax - ymin) / 2, 2)))
        if distbuff != 0:
            radius = radius + int(distbuff)
        circlePointData[0, :] = [centX, centY, corex, corey, radius]
    except arcpy.ExecuteError:
        exit_with_geoproc_error(_SCRIPT_NAME)
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)

    return circlePointData


def get_box_data(field_val, extent):
    """Create Numpy array with bounding box data."""
    box_data = npy.zeros((1, 5), dtype='float32')
    box_data[0, :] = [field_val, extent.XMin, extent.XMax, extent.YMax,
                      extent.YMin]
    return box_data


def get_sel_ext_box_coords(feature, field_name, field_val):
    """Get coordinates of bounding box for a unique feature."""
    shp_field = arcpy.Describe(feature).shapeFieldName
    search_row = arcpy.SearchCursor(
        feature,
        where_clause="{} = {}".format(field_name, field_val),
        fields=shp_field).next()
    extent = search_row.getValue(shp_field).extent
    del search_row
    return get_box_data(field_val, extent)


def get_ext_box_coords(feature):
    """Get coordinates of bounding box that contains all features."""
    return get_box_data(1, arcpy.Describe(feature).extent)


def make_points(workspace, pointArray, outFC):
    """Creates a shapefile with points specified by coordinates in pointArray

       outFC is just the filename, not path.
       pointArray is x,y,corex,corey,radius

    """
    try:
        wkspbefore = arcpy.env.workspace
        arcpy.env.workspace = workspace
        delete_data(outFC)
        arcpy.CreateFeatureclass_management(workspace, outFC, "POINT")
        if pointArray.shape[1] > 3:
            arcpy.AddField_management(outFC, "corex", "LONG")
            arcpy.AddField_management(outFC, "corey", "LONG")
            arcpy.AddField_management(outFC, "radius", "DOUBLE")
            arcpy.AddField_management(outFC, "cores_x_y", "TEXT")
        else:
            arcpy.AddField_management(outFC, "XCoord", "DOUBLE")
            arcpy.AddField_management(outFC, "YCoord", "DOUBLE")
            arcpy.AddField_management(outFC, cfg.COREFN, "LONG")
        rows = arcpy.InsertCursor(outFC)

        numPoints = pointArray.shape[0]
        for i in range(numPoints):
            point = arcpy.Point()
            point.ID = i
            point.X = float(pointArray[i, 0])
            point.Y = float(pointArray[i, 1])
            row = rows.newRow()
            row.shape = point
            row.setValue("ID", i)
            if pointArray.shape[1] > 3:
                row.setValue("corex", int(pointArray[i, 2]))
                row.setValue("corey", int(pointArray[i, 3]))
                row.setValue("cores_x_y", str(int(pointArray[i, 2])) + '_' +
                             str(int(pointArray[i, 3])))
                row.setValue("radius", float(pointArray[i, 4]))
            else:
                row.setValue("XCoord", float(pointArray[i, 0]))
                row.setValue("YCoord", float(pointArray[i, 1]))
                row.setValue(cfg.COREFN, float(pointArray[i, 2]))

            rows.insertRow(row)
            del row
            del point
        del rows
        arcpy.env.workspace = wkspbefore

    except arcpy.ExecuteError:
        exit_with_geoproc_error(_SCRIPT_NAME)
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)
    return


############################################################################
## LCP Shapefile Functions #################################################
############################################################################

def create_lcp_shapefile(ws,linktable, sourceCore, targetCore, lcpLoop):
    """Creates lcp shapefile.

    Shows locations of least-cost path lines attributed with corridor
    info/status.

    """
    try:
        arcpy.env.workspace = ws
        rows = get_links_from_core_pairs(linktable, sourceCore, targetCore)
        link = rows[0]

        lcpline = os.path.join(ws, "lcpline.shp")
        lcpRas = os.path.join(ws, "lcp")

        arcpy.RasterToPolyline_conversion(lcpRas, lcpline, "NODATA", "",
                                           "NO_SIMPLIFY")

        lcplineDslv = os.path.join(ws, "lcplineDslv.shp")
        arcpy.Dissolve_management(lcpline, lcplineDslv)

        arcpy.AddField_management(lcplineDslv, "Link_ID", "LONG", "5")
        arcpy.CalculateField_management(lcplineDslv, "Link_ID",
                                     int(linktable[link, cfg.LTB_LINKID]),
                                     "PYTHON_9.3")

        linktypecode = linktable[link, cfg.LTB_LINKTYPE]
        activelink, linktypedesc = get_link_type_desc(linktypecode)

        arcpy.AddField_management(lcplineDslv, "Active", "SHORT")
        arcpy.CalculateField_management(lcplineDslv, "Active", activelink,
                                     "PYTHON_9.3")

        arcpy.AddField_management(lcplineDslv, "Link_Info", "TEXT")
        arcpy.CalculateField_management(lcplineDslv, "Link_Info",
                                         linktypedesc, "PYTHON_9.3")

        arcpy.AddField_management(lcplineDslv, "From_Core", "LONG", "5")
        arcpy.CalculateField_management(lcplineDslv, "From_Core",
                                         int(sourceCore), "PYTHON_9.3")
        arcpy.AddField_management(lcplineDslv, "To_Core", "LONG", "5")
        arcpy.CalculateField_management(lcplineDslv, "To_Core",
                                         int(targetCore), "PYTHON_9.3")

        arcpy.AddField_management(lcplineDslv, "Euc_Dist", "DOUBLE", "10",
                                   "2")
        arcpy.CalculateField_management(lcplineDslv, "Euc_Dist",
                                     linktable[link, cfg.LTB_EUCDIST],
                                     "PYTHON_9.3")

        arcpy.AddField_management(lcplineDslv, "CW_Dist", "DOUBLE", "10", "2")
        arcpy.CalculateField_management(lcplineDslv, "CW_Dist",
                                     linktable[link, cfg.LTB_CWDIST],
                                     "PYTHON_9.3")
        arcpy.AddField_management(lcplineDslv, "LCP_Length", "DOUBLE", "10",
                                   "2")
        rows = arcpy.UpdateCursor(lcplineDslv)
        row = rows.next()
        while row:
            feat = row.shape
            lcpLength = int(feat.length)
            row.setValue("LCP_Length", lcpLength)
            rows.updateRow(row)
            row = rows.next()
        del row, rows

        try:
            distRatio1 = (float(linktable[link, cfg.LTB_CWDIST])
                        / float(linktable[link, cfg.LTB_EUCDIST]))
        except ZeroDivisionError:
            distRatio1 = -1

        arcpy.AddField_management(lcplineDslv, "cwd2Euc_R", "DOUBLE", "10",
                                   "2")
        arcpy.CalculateField_management(lcplineDslv, "cwd2Euc_R", distRatio1,
                                     "PYTHON_9.3")

        try:
            distRatio2 = (float(linktable[link, cfg.LTB_CWDIST])
                        / float(lcpLength))
        except ZeroDivisionError:
            distRatio2 = -1

        arcpy.AddField_management(lcplineDslv, "cwd2Path_R", "DOUBLE", "10",
                                   "2")
        arcpy.CalculateField_management(lcplineDslv, "cwd2Path_R", distRatio2,
                                     "PYTHON_9.3")

        lcpLoop = lcpLoop + 1
        lcpShapefile = os.path.join(cfg.DATAPASSDIR, "lcpLines_s3.shp")
        if lcpLoop == 1:
            if arcpy.Exists(lcpShapefile):
                try:
                    arcpy.Delete(lcpShapefile)

                except Exception:
                    dashline(1)
                    msg = ('ERROR: Could not remove LCP shapefile ' +
                           lcpShapefile + '. Was it open in ArcMap?\n You may '
                           'need to re-start ArcMap to release the file lock.')
                    raise_error(msg)

            arcpy.Copy_management(lcplineDslv, lcpShapefile)
        else:
            arcpy.Append_management(lcplineDslv, lcpShapefile, "TEST")

        return lcpLoop

    except arcpy.ExecuteError:
        exit_with_geoproc_error(_SCRIPT_NAME)
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)


def get_lcp_shapefile(lastStep, thisStep):
    """Returns path of lcp shapefile generated by previous step.

    """
    try:
        if thisStep == 5:
            oldLcpShapefile = os.path.join(
                cfg.DATAPASSDIR, "lcpLines_s" + str(lastStep) + ".shp")
            # If last step wasn't step 4 then must be step 3
            if not arcpy.Exists(oldLcpShapefile):
                # step 3
                oldLcpShapefile = os.path.join(
                    cfg.DATAPASSDIR, "lcpLines_s" + str(lastStep - 1) +
                    ".shp")

        elif thisStep > 5:
                # Steps not necessarily in order.  Look for highest step
                # number.
                lastStep = 8
                while True:
                    oldLcpShapefile = os.path.join(
                        cfg.DATAPASSDIR, "lcpLines_s" + str(lastStep) + ".shp")
                    if arcpy.Exists(oldLcpShapefile):
                        break
                    else:
                        lastStep = lastStep - 1
                        if lastStep == 2: # No previous shapefile found
                            dashline(1)
                            msg = ('ERROR: Could not find LCP shapefile from a '
                                    '\nstep previous to step '+ str(thisStep)
                                     +'.')
                            raise_error(msg)

        else:
            oldLcpShapefile = os.path.join(
                cfg.DATAPASSDIR, "lcpLines_s" + str(lastStep) + ".shp")
        return oldLcpShapefile

    except arcpy.ExecuteError:
        exit_with_geoproc_error(_SCRIPT_NAME)
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)


def update_lcp_shapefile(linktable, lastStep, thisStep):
    """Updates lcp shapefiles with new link information/status"""

    try:
        lcpShapefile = os.path.join(cfg.DATAPASSDIR, "lcpLines_s" +
                                    str(thisStep) + ".shp")

        numLinks = linktable.shape[0]
        # g1' g2' THEN c1 c2

        if thisStep > 5:
            linkTableTemp = linktable
        else:
            extraCols = npy.zeros((numLinks, 3), dtype="float64")
            linkTableTemp = npy.append(linktable, extraCols, axis=1)
            del extraCols
            linkTableTemp[:, cfg.LTB_LCPLEN] = -1
            linkTableTemp[:, cfg.LTB_CWDEUCR] = -1
            linkTableTemp[:, cfg.LTB_CWDPATHR] = -1


        if lastStep != thisStep:
            oldLcpShapefile = get_lcp_shapefile(lastStep, thisStep)
            delete_data(lcpShapefile)
            arcpy.Copy_management(oldLcpShapefile, lcpShapefile)
            #if thisStep > 5:
            arcpy.AddField_management(lcpShapefile, "Eff_Resist", "FLOAT") ###
            arcpy.AddField_management(lcpShapefile, "cwd2EffR_r", "FLOAT")
            arcpy.AddField_management(lcpShapefile, "CF_Central", "FLOAT") ###
        rows = arcpy.UpdateCursor(lcpShapefile)
        row = rows.next()
        line = 0
        while row:
            linkid = row.getValue("Link_ID")
            linktypecode = linktable[linkid - 1, cfg.LTB_LINKTYPE]
            activelink, linktypedesc = get_link_type_desc(linktypecode)
            row.setValue("Link_Info", linktypedesc)
            row.setValue("Active", activelink)
            if thisStep > 5:
                current = linkTableTemp[linkid - 1, cfg.LTB_CURRENT]
                effResist = linkTableTemp[linkid - 1, cfg.LTB_EFFRESIST]
                CWDTORRatio = linkTableTemp[linkid - 1, cfg.LTB_CWDTORR]
                #fixme: linkid - 1 assumes linktable ordered
                row.setValue("Eff_Resist", effResist)
                row.setValue("cwd2EffR_r",CWDTORRatio)
                row.setValue("CF_Central", current)
            else:
                row.setValue("Eff_Resist", -1)
                row.setValue("cwd2EffR_r",-1)
                row.setValue("CF_Central", -1)
            rows.updateRow(row)

            linktablerow = get_linktable_row(linkid, linkTableTemp)
            linkTableTemp[linktablerow, cfg.LTB_LCPLEN] = row.getValue(
                "LCP_Length")
            linkTableTemp[linktablerow, cfg.LTB_CWDEUCR] = row.getValue(
                "cwd2Euc_R")
            linkTableTemp[linktablerow, cfg.LTB_CWDPATHR] = row.getValue(
                "cwd2Path_R")
            row = rows.next()
            line = line + 1
        # delete cursor and row points to remove locks on the data
        del row, rows

        outputLcpShapefile = os.path.join(cfg.OUTPUTDIR, cfg.PREFIX +
                                        "_lcpLines_s" + str(thisStep) + ".shp")
        if arcpy.Exists(outputLcpShapefile):
            try:
                arcpy.Delete_management(outputLcpShapefile)
            except Exception:
                dashline(1)
                msg = ('ERROR: Could not remove lcp shapefile from output '
                       'directory: ' + outputLcpShapefile +
                       '. Is it open in ArcMap?\n You may '
                       'need to re-start ArcMap to release the file lock.')
                raise_error(msg)

        arcpy.Copy_management(lcpShapefile, outputLcpShapefile)

        return linkTableTemp

    except arcpy.ExecuteError:
        exit_with_geoproc_error(_SCRIPT_NAME)
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)


############################################################################
## Graph Functions ########################################################
############################################################################
def delete_row(A, delrow):
    """Deletes rows from a matrix

    From gapdt.py by Viral Shah

    """
    try:
        m = A.shape[0]
        n = A.shape[1]
        keeprows = npy.delete(npy.arange(0, m), delrow)
        keepcols = npy.arange(0, n)
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)
    return A[keeprows][:, keepcols]


def delete_col(A, delcol):
    """Deletes columns from a matrix

    From gapdt.py by Viral Shah

    """
    try:
        m = A.shape[0]
        n = A.shape[1]
        keeprows = npy.arange(0, m)
        keepcols = npy.delete(npy.arange(0, n), delcol)
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)
    return A[keeprows][:, keepcols]


def components_no_sparse(G):
    """Returns components of a graph while avoiding use of sparse matrices

    From gapdt.py by Viral Shah

    """
    try:
        U, V = npy.where(G)
        n = G.shape[0]
        D = npy.arange(0, n, dtype='int32')
        star = npy.zeros(n, 'int32')
        all_stars = False
        while True:
            star = check_stars(D, star)
            D = conditional_hooking(D, star, U, V)
            star = check_stars(D, star)
            D = unconditional_hooking(D, star, U, V)

            # pointer jumping
            D = D[D]

            if all_stars:
                return relabel(D, 1)

            if sum(star) == n:
                all_stars = True
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)


def relabel(oldlabel, offset=0): # same as gapdt
    """Utility for components code

    From gapdt.py by Viral Shah

    """
    newlabel = npy.zeros(npy.size(oldlabel), dtype='int32')
    s = npy.sort(oldlabel)
    perm = npy.argsort(oldlabel)
    f = npy.where(npy.diff(npy.concatenate(([s[0] - 1], s))))
    newlabel[f] = 1
    newlabel = npy.cumsum(newlabel)
    newlabel[perm] = npy.copy(newlabel)
    return newlabel - 1 + offset


def conditional_hooking (D, star, u, v):
    """Utility for components code (updated Jan 2012)

    From gapdt.py by Viral Shah

    """
    Du = D[u]
    Dv = D[v]

    hook = npy.where ((star[u] == 1) & (Du > Dv))
    D[Du[hook]] = Dv[hook]

    return D

def unconditional_hooking (D, star, u, v):
    """Utility for components code (updated Jan 2012)

    From gapdt.py by Viral Shah

    """

    Du = D[u]
    Dv = D[v]

    hook = npy.where((star[u] == 1) & (Du != Dv))
    D[Du[hook]] = Dv[hook]

    return D

def check_stars(D, star): # same as gapdt
    """Utility for components code

    From gapdt.py by Viral Shah

    """
    star[:] = 1
    notstars = npy.where(D != D[D])
    star[notstars] = 0
    star[D[D[notstars]]] = 0
    star = star[D]
    return star


############################################################################
## Input Functions ########################################################
############################################################################
def load_link_table(linkTableFile):
    """Reads link table created by previous step """
    try:
        linkTable1 = npy.loadtxt(linkTableFile, dtype='Float64',
                             comments='#', delimiter=',')
        if len(linkTable1) == linkTable1.size:  # Just one connection
            linktable = npy.zeros((1, len(linkTable1)), dtype='Float64')
            linktable[:, 0:len(linkTable1)] = linkTable1[0:len(linkTable1)]
        else:
            linktable = linkTable1
        return linktable
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)


############################################################################
## Output Functions ########################################################
############################################################################
def gprint(string):
    if string[0:7] == "Warning":
        arcpy.AddWarning(string)
    else:
        arcpy.AddMessage(string)
    try:
        if cfg.LOGMESSAGES:
            write_log(string)
    except Exception:
        pass

def create_log_file(messageDir, toolName, inParameters):
    ft = tuple(time.localtime())
    timeNow = time.ctime()
    fileName = ('%s_%s_%s_%s%s_%s.txt' % (ft[0], ft[1], ft[2], ft[3], ft[4], toolName))
    filePath = os.path.join(messageDir,fileName)
    try:
        logFile=open(filePath,'a')
    except Exception:
        logFile=open(filePath,'w')
    if inParameters is not None:
        logFile.write('*'*70 + '\n')
        logFile.write('Linkage Mapper log file: %s \n\n' % (toolName))
        logFile.write('Start time:\t%s \n' % (timeNow))
        logFile.write('Parameters:\t%s \n\n' % (inParameters))
    logFile.close()
    dashline()
    gprint('A record of run settings and messages can be found in your '
           'log directory:')
    gprint(cfg.MESSAGEDIR)
    dashline(2)

    return filePath

def write_log(string):
    try:
        logFile=open(cfg.logFilePath,'a')
    except Exception:
        logFile=open(cfg.logFilePath,'w')
    try:
        #Sometimes int objects returned for arc failures so need str below
        logFile.write(str(string) + '\n')
    except IOError:
        pass
    finally:
        logFile.close()


def write_custom_to_log(settings_file):
    """Write custom settings to log file."""
    write_log("Custom settings from {}:".format(settings_file))
    with open(settings_file) as custfile:
        for line in custfile.readlines():
            if not line.startswith('#') and "=" in line:
                write_log(line[:line.find("#")].replace(' =', ':'))
    write_log("")


def close_log_file():
    timeNow = time.ctime()
    try:
        write_log('\nStop time:\t\t%s \n\n' % (timeNow))
    except Exception:
        pass


def write_link_table(linktable, outlinkTableFile, *inLinkTableFile):
    """Writes link tables to pass link data between steps """
    try:

        numLinks = linktable.shape[0]
        outFile = open(outlinkTableFile, "w")

        if linktable.shape[1] == 10:
            outFile.write("# link,coreId1,coreId2,cluster1,cluster2,linkType,"
                           "eucDist,lcDist,eucAdj,cwdAdj\n")

            for x in range(0, numLinks):
                for y in range(0, 9):
                    outFile.write(str(linktable[x, y]) + ",")
                outFile.write(str(linktable[x, 9]))
                outFile.write("\n")

        elif linktable.shape[1] == 13:
            outFile.write("#link,coreId1,coreId2,cluster1,cluster2,linkType,"
                           "eucDist,lcDist,eucAdj,cwdAdj,lcpLength,"
                           "cwdToEucRatio,cwdToPathRatio\n")

            for x in range(0, numLinks):
                for y in range(0, 12):
                    outFile.write(str(linktable[x, y]) + ",")
                outFile.write(str(linktable[x, 12]))
                outFile.write("\n")
        else: #Called from pinch point or centrality tool, has 16 columns
            outFile.write("#link,coreId1,coreId2,cluster1,cluster2,linkType,"
                           "eucDist,lcDist,eucAdj,cwdAdj,lcpLength,"
                           "cwdToEucRatio,cwdToPathRatio,Eff_Resist,"
                           "CWDTORRatio,CF_Centrality\n")

            for x in range(0, numLinks):
                for y in range(0, 15):
                    outFile.write(str(linktable[x, y]) + ",")
                outFile.write(str(linktable[x, 15]))
                outFile.write("\n")


        outFile.write("# Linkage Mapper Version " + cfg.releaseNum)
        outFile.write("\n# ---Run Settings---")
        outFile.write("\n# Project Directory: " + cfg.PROJECTDIR)
        if cfg.TOOL == 'linkage_mapper':
            outFile.write("\n# Core Area Feature Class: " + cfg.COREFC)

            outFile.write("\n# Core Area Field Name: " + cfg.COREFN)
            outFile.write("\n# Resistance Raster: " + cfg.RESRAST_IN)
            outFile.write("\n# Step 1 - Identify Adjacent Core Areas: " +
                           str(cfg.STEP1))
            outFile.write("\n# Step 1 Adjacency Method Includes Cost-Weighted "
                           "Distance: " + str(cfg.S1ADJMETH_CW))
            outFile.write("\n# Step 1 Adjacency Method Includes Euclidean "
                           "Distance: " + str(cfg.S1ADJMETH_EU))
            outFile.write("\n# Step 2 - Construct a Network of Core Areas: " +
                           str(cfg.STEP2))
            outFile.write("\n# Conefor Distances Text File: " +
                          str(cfg.S2EUCDISTFILE))
            outFile.write("\n# Network Adjacency Method Includes Cost-Weighted "
                           "Distance: " + str(cfg.S2ADJMETH_CW))
            outFile.write("\n# Network Adjacency Method Includes Euclidean "
                           "Distance: " + str(cfg.S2ADJMETH_EU))
            outFile.write("\n# Step 3 - Calculate Cost-Weighted Distances and "
                           "Least-Cost Paths: " + str(cfg.STEP3))
            outFile.write("\n# Drop Corridors that Intersect Core Areas: "
                           + str(cfg.S3DROPLCCS))
            outFile.write("\n# Step 4 - Refine Network: " + str(cfg.STEP4))
            if cfg.IGNORES4MAXNN:
                outFile.write("\n# Option A - Number of Connected Nearest Neighbors: Unlimimted")
            else:
                outFile.write("\n# Option A - Number of Connected Nearest Neighbors: "
                               + str(cfg.S4MAXNN))
            outFile.write("\n# Option B - Nearest Neighbor Measurement Unit is "
                           "Cost-Weighted Distance: " + str(cfg.S4DISTTYPE_CW))
            outFile.write("\n# Option C - Connect Neighboring Constellations : "
                           + str(cfg.S4CONNECT))
            outFile.write("\n# Step 5 - Calculate Normalize and Mosaic "
                           "Corridors: " + str(cfg.STEP5))
            outFile.write("\n# Bounding Circles Buffer Distance: "
                           + str(cfg.BUFFERDIST))
            outFile.write("\n# Maximum Cost-Weighted Corridor Distance: "
                           + str(cfg.MAXCOSTDIST))
            outFile.write("\n# Maximum Euclidean Corridor Distance: "
                           + str(cfg.MAXEUCDIST))
            outFile.write("\n# Minimum Cost-Weighted Corridor Distance: "
                           + str(cfg.MINCOSTDIST))
            outFile.write("\n# Minimum Euclidean Corridor Distance: "
                           + str(cfg.MINEUCDIST))


        elif cfg.TOOL == 'pinchpoint_mapper':
            outFile.write("\n# Pinchpoints Analyzed: "
                           + str(cfg.DOPINCH))
            outFile.write("\n# Resistance Raster: " + cfg.RESRAST)
            outFile.write("\n# CWD Cutoff Distance: "
                           + str(cfg.CWDCUTOFF))
            outFile.write("\n# Resistance Raster ValuesS quared for "
                            "Circuitscape Analyses: "+ str(cfg.SQUARERESISTANCES))
            outFile.write("\n# Network Centrality Analyzed: "
                          + str(cfg.DOCENTRALITY))
            outFile.write("\n# Core Area Feature Class: " + cfg.COREFC)

            outFile.write("\n# Core Area Field Name: " + cfg.COREFN)
            FN = str(inLinkTableFile)
            FN = FN.replace("('", "")
            FN = FN.replace("',)", "")
            outFile.write("\n# Link Data Derived from: " + FN)


        outFile.close()
    except arcpy.ExecuteError:
        exit_with_geoproc_error(_SCRIPT_NAME)
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)
    return


def write_adj_file(outcsvfile, adjTable):
    """Outputs adjacent core areas to pass adjacency info between steps"""
    outfile = open(outcsvfile, "w")
    outfile.write("#Edge" + "," + str(cfg.COREFN) + "," + str(cfg.COREFN) +
                  "_1" + "\n")
    for x in range(0, len(adjTable)):
        outfile.write(str(x) + "," + str(adjTable[x, 0]) + "," +
                      str(adjTable[x, 1]) + "\n")
    outfile.close()


def write_link_maps(linkTableFile, step):
    """Writes stick maps (aka link maps)

    These are vector line files showing links connecting core area pairs. Links
    contain attributes about corridor characteristics.

    """
    try:
        arcpy.env.workspace = cfg.OUTPUTDIR
        linktable = load_link_table(linkTableFile)

        numLinks = linktable.shape[0]
        arcpy.toolbox = "management"

        coresForLinework = "cores_for_linework.shp"

        # Preferred method to get geometric center
        pointArray = get_centroids(cfg.COREFC, cfg.COREFN)

        make_points(arcpy.env.workspace, pointArray, coresForLinework)
        numLinks = linktable.shape[0]

        coreLinks = linktable

        # create coreCoords array, with geographic centers of cores
        coreCoords = npy.zeros(pointArray.shape, dtype='float64')
        coreCoords[:, 0] = pointArray[:, 2]
        coreCoords[:, 1] = pointArray[:, 0]
        coreCoords[:, 2] = pointArray[:, 1]

        # Create linkCoords array
        linkCoords = npy.zeros((len(coreLinks), 13))
        linkCoords[:, 0] = coreLinks[:, cfg.LTB_LINKID]
        linkCoords[:, 1:3] = npy.sort(
            coreLinks[:, cfg.LTB_CORE1:cfg.LTB_CORE2 + 1])
        linkCoords[:, 3:5] = coreLinks[:, cfg.LTB_EUCDIST:cfg.LTB_CWDIST + 1]
        linkCoords[:, 9] = coreLinks[:, cfg.LTB_LINKTYPE]
        if step > 5:
            linkCoords[:, 10] = coreLinks[:, cfg.LTB_EFFRESIST]
            linkCoords[:, 11] = coreLinks[:, cfg.LTB_CWDTORR]
            linkCoords[:, 12] = coreLinks[:, cfg.LTB_CURRENT]
        else:
            linkCoords[:, 10] = -1
            linkCoords[:, 11] = -1
            linkCoords[:, 12] = -1

        if len(coreLinks) > 0:
            ind = npy.lexsort((linkCoords[:, 2], linkCoords[:, 1]))
            linkCoords = linkCoords[ind]

        # Get core coordinates into linkCoords
        for i in range(0, len(linkCoords)):
            grp1 = linkCoords[i, 1]
            grp2 = linkCoords[i, 2]

            for core in range(0, len(coreCoords)):
                if coreCoords[core, 0] == grp1:
                    linkCoords[i, 5] = coreCoords[core, 1]
                    linkCoords[i, 6] = coreCoords[core, 2]
                elif coreCoords[core, 0] == grp2:
                    linkCoords[i, 7] = coreCoords[core, 1]
                    linkCoords[i, 8] = coreCoords[core, 2]

        if len(coreLinks) > 0:
            ind = npy.argsort((linkCoords[:, 0]))  # Sort by LTB_LINKID
            linkCoords = linkCoords[ind]

        coreLinksShapefile = cfg.PREFIX + '_sticks_s' + str(step) + '.shp'

        file = os.path.join(arcpy.env.workspace,coreLinksShapefile)
        if arcpy.Exists(file):
            try:
                arcpy.Delete_management(file)
            except Exception:
                dashline(1)
                msg = ('ERROR: Could not remove shapefile ' +
                       coreLinksShapefile + '. Was it open in ArcMap?\n You may '
                       'need to re-start ArcMap to release the file lock.')
                raise_error(msg)

        # make coreLinks.shp using linkCoords table
        # will contain linework between each pair of connected cores
        arcpy.CreateFeatureclass_management(
            arcpy.env.workspace, coreLinksShapefile, "POLYLINE")

        # ADD ATTRIBUTES
        arcpy.AddField_management(coreLinksShapefile, "Link_ID", "LONG")
        arcpy.AddField_management(coreLinksShapefile, "Active", "SHORT")
        arcpy.AddField_management(coreLinksShapefile, "Link_Info", "TEXT")
        arcpy.AddField_management(coreLinksShapefile, "From_Core", "LONG")
        arcpy.AddField_management(coreLinksShapefile, "To_Core", "LONG")
        arcpy.AddField_management(coreLinksShapefile, "Euc_Dist", "FLOAT")
        arcpy.AddField_management(coreLinksShapefile, "CW_Dist", "FLOAT")
        arcpy.AddField_management(coreLinksShapefile, "cwd2Euc_R", "FLOAT")
        arcpy.AddField_management(coreLinksShapefile, "Eff_Resist", "FLOAT")
        arcpy.AddField_management(coreLinksShapefile, "cwd2EffR_r", "FLOAT")
        arcpy.AddField_management(coreLinksShapefile, "CF_Central", "FLOAT")
        #Create an Array and Point object.
        lineArray = arcpy.Array()
        pnt = arcpy.Point()

        # linkCoords indices:
        numLinks = len(linkCoords)

        #Open a cursor to insert rows into the shapefile.
        cur = arcpy.InsertCursor(coreLinksShapefile)

        ##Loop through each record in linkCoords table
        for i in range(0, numLinks):

            #Set the X and Y coordinates for origin vertex.
            pnt.x = linkCoords[i, 5]
            pnt.y = linkCoords[i, 6]
            #Insert it into the line array
            lineArray.add(pnt)

            #Set the X and Y coordinates for destination vertex
            pnt.x = linkCoords[i, 7]
            pnt.y = linkCoords[i, 8]
            #Insert it into the line array
            lineArray.add(pnt)

            #Insert the new poly into the feature class.
            feature = cur.newRow()
            feature.shape = lineArray
            cur.insertRow(feature)

            lineArray.removeAll()
        del cur

        #Add attribute data to link shapefile
        rows = arcpy.UpdateCursor(coreLinksShapefile)
        row = rows.next()
        line = 0
        while row:
            # linkCoords indices
            row.setValue("Link_ID", linkCoords[line, 0])
            if linkCoords[line, 9] == 2:
                row.setValue("Link_Info", "Group_Pair")
            linktypecode = linkCoords[line, 9]
            activelink, linktypedesc = get_link_type_desc(linktypecode)
            row.setValue("Active", activelink)
            row.setValue("Link_Info", linktypedesc)

            row.setValue("From_Core", linkCoords[line, 1])
            row.setValue("To_Core", linkCoords[line, 2])
            row.setValue("Euc_Dist", linkCoords[line, 3])
            row.setValue("CW_Dist", linkCoords[line, 4])
            if linkCoords[line, 4] <= 0 or linkCoords[line, 3] <= 0:
                row.setValue("cwd2Euc_R", -1)
            else:
                row.setValue("cwd2Euc_R", linkCoords[line, 4] /
                             linkCoords[line, 3])
            row.setValue("Eff_Resist", linkCoords[line, 10])
            row.setValue("cwd2EffR_r", linkCoords[line, 11])
            row.setValue("CF_Central", linkCoords[line, 12])

            rows.updateRow(row)
            row = rows.next()
            line = line + 1

        del row, rows

        #clean up temp files
        delete_data(coresForLinework)

        return

    except arcpy.ExecuteError:
        exit_with_geoproc_error(_SCRIPT_NAME)
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)


def set_dataframe_sr():
    """Sets data frame spatial reference to input core area projection.

       Differing spatial reference can cause problems in step 1.
    """
    try:
        sr = arcpy.Describe(cfg.COREFC).spatialReference
        gprint('Setting data frame spatial reference to that of '
                'core area feature class.')
    except Exception:
        try:
            sr = arcpy.Describe(cfg.RESRAST).spatialReference
        except Exception:
            return
    try:
        mxd = arcpy.mapping.MapDocument("current")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        df.spatialReference = sr
    except Exception:
        pass


def create_dir(lmfolder):
    """Creates folder if it doesn't exist."""
    if not os.path.exists(lmfolder):
        arcpy.CreateFolder_management(os.path.dirname(lmfolder),
                                       os.path.basename(lmfolder))


def move_old_results():
    """Updates project directory structure to new version.

    """
    try:
        oldFolder = cfg.CWDBASEDIR_OLD
        newFolder = cfg.CWDBASEDIR
        move_results_folder(oldFolder, newFolder)

        oldFolder = cfg.ADJACENCYDIR_OLD
        newFolder = cfg.ADJACENCYDIR
        move_results_folder(oldFolder, newFolder)

        oldFolder = cfg.LCCBASEDIR_OLD
        newFolder = cfg.LCCBASEDIR
        move_results_folder(oldFolder, newFolder)

        oldFolder = cfg.LOGDIR_OLD
        newFolder = cfg.LOGDIR
        move_results_folder(oldFolder, newFolder)

        oldFolder = cfg.MESSAGEDIR_OLD
        newFolder = cfg.MESSAGEDIR
        move_results_folder(oldFolder, newFolder)

    except Exception:
        exit_with_python_error(_SCRIPT_NAME)


def move_results_folder(oldFolder, newFolder):
    try:
        if (os.path.exists(oldFolder) and not os.path.exists(newFolder)):
            os.rename(oldFolder, newFolder)
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)


def delete_file(filename):
    """Delete file from disk."""
    try:
        os.remove(filename)
    except OSError:
        pass


def delete_dir(dir_path):
    """Delete directory from disk ignoring errors."""
    if os.path.isdir(dir_path):
        shutil.rmtree(dir_path, ignore_errors=True)


def delete_data(item):
    """Delete data from disk using ArcPy."""
    try:
        arcpy.Delete_management(item)
    except arcpy.ExecuteError:
        pass


def clean_out_workspace(ws):
    """Delete all datasets in an ArcPy workspace."""
    if arcpy.Exists(ws):
        cur_ws = arcpy.env.workspace
        arcpy.env.workspace = ws
        datasets = arcpy.ListDatasets()
        for data in datasets:
                delete_data(data)
        arcpy.env.workspace = cur_ws


def make_raster_paths(no_rast, base_dir, sub_dir):
    r"""Set up raster directories to insure < 100 grids in any one directory.

    Outputs are written to: base\sub for rasters 1-99, base\sub1 for rasters
    100-199, etc.
    """
    delete_dir(base_dir)
    try:
        os.makedirs(os.path.join(base_dir, sub_dir))
        for dir_no in range(1, (no_rast / 100) + 1):
            os.mkdir(os.path.join(base_dir,
                                  ''.join([sub_dir, str(dir_no)])))
    except OSError:
        exit_with_python_error(_SCRIPT_NAME)


def rast_path(count, base_dir, sub_dir):
    """Return the path for the raster corresponding to its count."""
    dir_count = count / 100
    if dir_count > 0:
        rast_path = os.path.join(base_dir,
                                 ''.join([sub_dir, str(dir_count)]))
    else:
        rast_path = os.path.join(base_dir, sub_dir)
    return rast_path


def get_cwd_path(core):
    """Return the path for the cwd raster corresponding to a core area."""
    dir_path = rast_path(core, cfg.CWDBASEDIR, cfg.CWDSUBDIR_NM)
    fname = "cwd_{}".format(core)
    return os.path.join(dir_path, fname)


def get_focal_path(core,radius):
    """Returns the path for the focal raster corresponding to a core area """
    dirCount = int(core / 100)
    focalDir1 = cfg.FOCALSUBDIR1_NM + str(radius)
    if dirCount > 0:
        return os.path.join(cfg.BARRIERBASEDIR, focalDir1,
                         cfg.FOCALSUBDIR2_NM + str(dirCount),
                         cfg.FOCALGRID_NM + str(core))
    else:
        return os.path.join(cfg.BARRIERBASEDIR, focalDir1,
                         cfg.FOCALSUBDIR2_NM, cfg.FOCALGRID_NM + str(core))


def check_project_dir():
    """Checks to make sure path name is not too long.

    Long path names can cause problems with ESRI grids.
    """
    if len(cfg.PROJECTDIR) > 140:
        msg = ('ERROR: Project directory "' + cfg.PROJECTDIR +
               '" is too deep.  Please choose a shallow directory'
               '(something like "C:\PUMA").')
        raise_error(msg)

    if "-" in cfg.PROJECTDIR or " " in cfg.PROJECTDIR or "." in cfg.PROJECTDIR:
        msg = ('ERROR: Project directory cannot contain spaces, dashes, or '
                'special characters.')
        raise_error(msg)
    head=cfg.PROJECTDIR
    for i in range(1,100):
        if len(head) < 4: # We've gotten to the base of the tree
            break
        head,tail=os.path.split(head)
        if tail[0].isdigit():
            msg = ('ERROR: No directory names in project directory path can start with a number or '
                    'else Arc will crash. Please change name of "' + tail + '" or choose a new directory.')
            raise_error(msg)
    return

def get_prev_step_link_table(step):
    """Returns the name of the link table created by the previous step"""
    try:
        prevStep = step - 1
        if step > 5:
            prevStepLinkTable = os.path.join(cfg.DATAPASSDIR,
                                         'linkTable_s5_plus.csv')
            gprint('\nLooking for ' + prevStepLinkTable)

            if os.path.exists(prevStepLinkTable):
                return prevStepLinkTable
            prevStepLinkTable = os.path.join(cfg.DATAPASSDIR,
                                             'linkTable_s5.csv')
            gprint('\nLooking for ' + prevStepLinkTable)

            if os.path.exists(prevStepLinkTable):
                return prevStepLinkTable

        if step > 4:
            prevStepLinkTable = os.path.join(cfg.DATAPASSDIR,
                                             'linkTable_s4.csv')
            gprint('\nLooking for ' + prevStepLinkTable)

            if os.path.exists(prevStepLinkTable):
                return prevStepLinkTable
            else:
                prevStep = 3  # Can skip step 4

        prevStepLinkTable = os.path.join(cfg.DATAPASSDIR, 'linkTable_s' +
                                         str(prevStep) + '.csv')
        gprint('\nLooking for ' + cfg.DATAPASSDIR +
                          '\linkTable_s' + str(prevStep) + '.csv')
        if os.path.exists(prevStepLinkTable):
            return prevStepLinkTable
        else:
            msg = ('\nERROR: Could not find a linktable from step previous to '
                   'step #' + str(step) + ' in datapass directory.  See above '
                   'for valid linktable files.')
            raise_error(msg)

    except arcpy.ExecuteError:
        exit_with_geoproc_error(_SCRIPT_NAME)
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)


def get_this_step_link_table(step):
    """Returns name of link table to write for current step"""
    try:
        if step > 5:
            filename = os.path.join(cfg.DATAPASSDIR, 'linkTable_s5_plus.csv')

        else:
            filename = os.path.join(cfg.DATAPASSDIR, 'linkTable_s' + str(step)
                             + '.csv')
        return filename

    except arcpy.ExecuteError:
        exit_with_geoproc_error(_SCRIPT_NAME)
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)


def clean_up_link_tables(step):
    """Remove link tables from previous runs."""
    try:
        filename = os.path.join(cfg.DATAPASSDIR, 'linkTable_s7_s8.csv')
        delete_file(filename)
        filename = os.path.join(cfg.DATAPASSDIR, 'linkTable_s5_plus.csv')
        delete_file(filename)
        for stepNum in range(step, 9):
            filename = os.path.join(cfg.DATAPASSDIR, 'linkTable_s' +
                                    str(stepNum) + '.csv')
            delete_file(filename)
            lcpFC = os.path.join(cfg.DATAPASSDIR,'lcpLines_s' +
                                    str(stepNum) + '.shp')
            delete_data(lcpFC)
        filename = os.path.join(cfg.OUTPUTDIR, 'linkTable_final.csv')
        delete_file(filename)
        filename = os.path.join(cfg.OUTPUTDIR, cfg.PREFIX + '_linkTable_final.csv')
        delete_file(filename)
        filename = os.path.join(cfg.OUTPUTDIR, 'linkTable_s5.csv')
        delete_file(filename)
        filename = os.path.join(cfg.OUTPUTDIR, cfg.PREFIX + 'linkTable_s5.csv')
        delete_file(filename)
    except arcpy.ExecuteError:
        exit_with_geoproc_error(_SCRIPT_NAME)
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)


def copy_final_link_maps(step):
    """Copies final link maps from datapass to the output directory"""
    try:
        PREFIX = cfg.PREFIX
        coreLinksShapefile = os.path.join(cfg.OUTPUTDIR, PREFIX + '_sticks_s'
                                          + str(step) + '.shp')
        lcpShapefile = os.path.join(cfg.DATAPASSDIR, 'lcpLines_s' +
                                    str(step) + '.shp')

        if not arcpy.Exists(cfg.LINKMAPGDB):
            arcpy.CreateFileGDB_management(os.path.dirname(cfg.LINKMAPGDB),
                             os.path.basename(cfg.LINKMAPGDB))

        if not arcpy.Exists(cfg.LOGLINKMAPGDB):
            arcpy.CreateFileGDB_management(os.path.dirname(cfg.LOGLINKMAPGDB),
                             os.path.basename(cfg.LOGLINKMAPGDB))

        if arcpy.Exists(coreLinksShapefile):
            arcpy.MakeFeatureLayer_management(coreLinksShapefile, "flinks")


            field = "Active"
            expression = field + " = " + str(1)
            arcpy.SelectLayerByAttribute_management("flinks", "NEW_SELECTION",
                                          expression)

            activeLinksShapefile = os.path.join(cfg.LINKMAPGDB,
                                               PREFIX + '_Sticks')
            arcpy.CopyFeatures_management("flinks", activeLinksShapefile)

            rename_fields(activeLinksShapefile)

            expression = field + " = " + str(0)
            arcpy.SelectLayerByAttribute_management("flinks", "NEW_SELECTION",
                                          expression)
            inActiveLinksShapefile = os.path.join(cfg.LINKMAPGDB,
                                                  PREFIX + '_Inactive_Sticks')
            arcpy.CopyFeatures_management("flinks", inActiveLinksShapefile)


        if arcpy.Exists(lcpShapefile):
            arcpy.MakeFeatureLayer_management(lcpShapefile, "flcp")
            field = "Active"
            expression = field + " = " + str(1)
            arcpy.SelectLayerByAttribute_management("flcp", "NEW_SELECTION", expression)

            activeLcpShapefile = os.path.join(cfg.LINKMAPGDB,
                                              PREFIX + '_LCPs')

            arcpy.CopyFeatures_management("flcp", activeLcpShapefile)
            rename_fields(activeLcpShapefile)

            expression = field + " = " + str(0)
            arcpy.SelectLayerByAttribute_management("flcp", "NEW_SELECTION", expression)
            inActiveLcpShapefile = os.path.join(cfg.LINKMAPGDB,
                                               PREFIX + '_Inactive_LCPs')
            arcpy.CopyFeatures_management("flcp", inActiveLcpShapefile)

        # Move stick and lcp maps for each step to log directory to reduce
        # clutter in output
        for i in range(2, 9):
            # delete log file from pre 0.7.7 versions
            oldLogLinkFile = os.path.join(cfg.LOGDIR, PREFIX + '_sticks_s'
                                        + str(i) + '.shp')
            delete_data(oldLogLinkFile)

            oldLinkFile = os.path.join(cfg.OUTPUTDIR, PREFIX + '_sticks_s'
                                        + str(i) + '.shp')
            logLinkFile = os.path.join(cfg.LOGLINKMAPGDB, PREFIX + '_sticks_s'
                                        + str(i))
            if arcpy.Exists(oldLinkFile):
                try:
                    move_map(oldLinkFile, logLinkFile)
                except Exception:
                    pass

            # delete log file from pre 0.7.7 versions
            oldLogLcpShapeFile = os.path.join(cfg.LOGDIR, PREFIX + '_lcpLines_s' +
                                           str(i) + '.shp')
            delete_data(oldLogLcpShapeFile)
            oldLcpShapeFile = os.path.join(cfg.OUTPUTDIR, PREFIX + '_lcpLines_s'
                                           + str(i) + '.shp')
            logLcpShapeFile = os.path.join(cfg.LOGLINKMAPGDB, PREFIX + '_lcpLines_s' +
                                           str(i))
            if arcpy.Exists(oldLcpShapeFile):
                try:
                    move_map(oldLcpShapeFile, logLcpShapeFile)
                except Exception:
                    pass
        return
    except arcpy.ExecuteError:
        exit_with_geoproc_error(_SCRIPT_NAME)
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)


def move_map(oldMap, newMap):
    """Moves a map to a new location """
    if arcpy.Exists(oldMap):
        delete_data(newMap)
        try:
            arcpy.CopyFeatures_management(oldMap, newMap)
            delete_data(oldMap)
        except Exception:
            pass
    return


@Retry(10)
def call_circuitscape(cspath, outConfigFile):
    """Call Circuitscape."""
    mem_flag = False
    fail_flag = False
    gprint('     Calling Circuitscape:')
    proc = subprocess.Popen([cspath, outConfigFile],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            shell=True)
    while proc.poll() is None:
        output = proc.stdout.readline()

        if 'Traceback' in output:
            gprint("\nCircuitscape failed.")
            fail_flag = True
            if 'memory' in output:
                mem_flag = True
        if ('Processing' not in output and 'laplacian' not in output
                and 'node_map' not in output
                and (('--' in output) or ('sec' in output)
                     or (fail_flag is True))):
            gprint("      " + output.replace("\r\n", ""))

    # Catch any output lost if process closes too quickly
    output = proc.communicate()[0]
    for line in output.split('\r\n'):
        if 'Traceback' in line:
            gprint("\nCircuitscape failed.")
            if 'valid sources' in output.lower():
                gprint('Corridors may be too narrow. Try upping your CWD '
                       'cutoff distance.')
            if 'memory' in line:
                mem_flag = True
        if ('Processing' not in line and 'laplacian' not in line
                and 'node_map' not in line
                and (('--' in line) or ('sec' in line)
                     or (fail_flag is True))):
            gprint("      " + str(line))
    return mem_flag


def rename_fields(FC):
    try:
        arcpy.AddField_management(FC, "cw_to_Euc_Dist_Ratio", "FLOAT")
        arcpy.CalculateField_management(FC, "cw_to_Euc_Dist_Ratio","!cwd2Euc_R!", "PYTHON")
        arcpy.DeleteField_management(FC, "cwd2Euc_R")

        fieldList = arcpy.ListFields(FC)
        for field in fieldList:
            if str(field.name) == "cwd2Path_R":
                arcpy.AddField_management(FC, "cwd_to_Path_Length_Ratio", "FLOAT")
                arcpy.CalculateField_management(FC, "cwd_to_Path_Length_Ratio","!cwd2Path_R!", "PYTHON")
                arcpy.DeleteField_management(FC, "cwd2Path_R")

        arcpy.AddField_management(FC, "Current_Flow_Centrality", "FLOAT")
        arcpy.CalculateField_management(FC, "Current_Flow_Centrality","!CF_Central!", "PYTHON")
        arcpy.DeleteField_management(FC, "CF_Central")

        arcpy.AddField_management(FC, "Effective_Resistance", "FLOAT")
        arcpy.CalculateField_management(FC, "Effective_Resistance","!Eff_Resist!", "PYTHON")
        arcpy.DeleteField_management(FC, "Eff_Resist")

        arcpy.AddField_management(FC, "cwd_to_Eff_Resist_Ratio", "FLOAT")
        arcpy.CalculateField_management(FC, "cwd_to_Eff_Resist_Ratio","!cwd2EffR_r!", "PYTHON")
        arcpy.DeleteField_management(FC, "cwd2EffR_r")

    except arcpy.ExecuteError:
        exit_with_geoproc_error(_SCRIPT_NAME)
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)

############################################################################
##Error Checking and Handling Functions ####################################
############################################################################
def print_arcgis_failures(statement, failures):
    """ Reports ArcGIS call that's failing and decides whether to restart
        iteration.

    """
    dashline(1)
    gprint('***Problem encountered executing statement:')
    gprint('"' + statement + '"')

    print_warnings()

    failures = failures + 1
    return failures


def print_drive_warning():
    gprint('\n********************************************************')
    drive, depth, realpath = get_dir_depth(cfg.PROJECTDIR)
    if drive.lower() != 'c' or depth > 3 or 'dropbox' in realpath.lower():
        gprint('NOTE: ArcGIS errors are more likely when writing to remote '
            'drives or deep file structures. We recommend shallow '
            'project directories on local drives, like C:\puma. '
            'Errors may also result from conflicts with anti-virus '
            'software (known problems with AVG). We have also seen '
            'conflicts when writing to synchronized folders like DROPBOX.\n\n'
            'Note also that Linkage Mapper tools often work best when run '
            'from ArcCatalog instead of ArcMap. \n ')
    else:
        gprint('NOTE: Linkage Mapper tools often work best when run  '
            'from ArcCatalog instead of ArcMap. \n')
    localdict = locale.localeconv()
    if localdict['decimal_point'] != '.':
        msg = ('ERROR: It looks like decimals are indicated by commas instead of decimal points.\n'
            'Try changing your Windows Regional and Language settings to English (United States)\n'
            'or another convention that uses decimal points.\n')
        raise_error(msg)


def get_dir_depth(dir):
    import string
    realpath = os.path.normpath(dir)
    drive = realpath[0]
    depth = 0
    for i in range(0, len(realpath)):
        if realpath[i] == os.path.sep:
            depth = depth + 1
    return drive, depth, realpath


def check_cores(FC,FN):
    """Checks for positive integer core IDs with appropriate naming."""
    try:
        if not cfg.STEP1:
            try:
                adjList = npy.loadtxt(cfg.EUCADJFILE, dtype='string',
                                         comments = "x", delimiter=',')
                prevCoreFN =adjList[0,1]
            except Exception: # If file not found
                prevCoreFN = cfg.COREFN
            if prevCoreFN != cfg.COREFN:
                msg = ('\nError: Core field name must be the same as used in '
                        'Linkage Mapper step 1 ("' + prevCoreFN + '"). '
                        '\nPlease make sure you are using the same core area '
                        'file and resistance raster as well.')
                raise_error(msg)

        invalidFNs = ['fid','id','oid','shape']
        if FN.lower() in invalidFNs:
            dashline(1)
            msg = ('ERROR: Core area field name "ID", "FID", "OID", and "Shape" '
                   'are reserved for ArcGIS. Please choose another field- '
                    'must be a positive integer.')
            raise_error(msg)

        cid_field = arcpy.ListFields(FC, FN)[0]
        if cid_field.type not in ('Integer', 'SmallInteger'):
            dashline(1)
            raise_error('ERROR: Core area field must be in Integer format.')

        coreList = get_core_list(FC,FN)
        if npy.amin(coreList) < 1:
            dashline(1)
            msg = ('ERROR: Core area field must contain only positive integers. ')
            raise_error(msg)

        rows = arcpy.SearchCursor(FC)
        row = rows.next()
        feat = row.shape
        center = str(feat.centroid)
        xy = center.split(" ")
        if "," in xy[0]:
            msg = ('ERROR: It appears that your region settings are not in '
                  'USA format (decimal commas are used instead of decimal '
                  'points). '
                  'Please change your region settings in Windows to English (USA) or '
                  'another convention that uses decimal points.  You may need to'
                  'to modify the coordinate system of your input files as well.')
            raise_error(msg)

    except arcpy.ExecuteError:
        exit_with_geoproc_error(_SCRIPT_NAME)
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)


def retry_arc_error(count, statement):
    """Re-tries ArcGIS calls in case of server problems or other 'hiccups'."""
    try:
        if count < 5:
            count = count + 1
            sleepTime = 20*count

            arcpy.AddWarning('-------------------------------------------------')
            arcpy.AddWarning('Failed to execute ' + statement + ' on try '
                              '#' + str(count) + '.\n')

            print_warnings()

            arcpy.AddWarning("Will try again. ")
            arcpy.AddWarning('---------TRYING AGAIN IN ' +
                                   str(int(sleepTime)) + ' SECONDS---------\n')
            snooze(sleepTime)
            return count, True

        elif count < 7:
            sleepTime = 300
            count = count + 1
            arcpy.AddWarning('Failed to execute ' + statement + ' on try #' +
                        str(count) + '.\n Could be an ArcGIS hiccup.  Trying'
                        ' again in 5 minutes.\n')
            snooze(sleepTime)

            return count, True

        else:
            sleepTime = 300
            count = count + 1
            arcpy.AddWarning('Failed to execute ' + statement + ' on try #' +
                        str(count) + '.\n Could be an ArcGIS hiccup.  Trying'
                        ' one last time in 5 minutes.\n')
            snooze(sleepTime)

            return count, False
    except Exception:
        exit_with_python_error(_SCRIPT_NAME)


def print_warnings():
    tb = sys.exc_info()[2]  # get the traceback object
    # tbinfo contains the error's line number and the code
    tbinfo = traceback.format_tb(tb)[0]
    line = tbinfo.split(", ")[1]
    filename = tbinfo.split(", ")[0]
    filename = filename.rsplit("File ")[1]

    if arcpy.GetMaxSeverity > 1:
        msg = ("The following ArcGIS error is being reported "
                    "on line " + line + " of " + filename + ":")
        arcpy.AddWarning(msg)
        write_log(msg)
        arcpy.AddWarning(arcpy.GetMessages(2))
        write_log(arcpy.GetMessages(2))
        print_drive_warning()

    else:
        msg = ("The following error is being reported at "
                        + line + " of " + filename + ":")
        err = traceback.format_exc().splitlines()[-1]
        arcpy.AddWarning(msg)
        arcpy.AddWarning(err + '\n')
        write_log(msg)
        write_log(err)


def snooze(sleepTime):
    for i in range(1,int(sleepTime)+1):
        time.sleep(1)
        # Dummy operations to give user ability to cancel:
        installD = arcpy.GetInstallInfo("desktop")



def exit_with_geoproc_error(filename):
    """Handle geoprocessor errors and provide details to user"""
    dashline()
    tb = sys.exc_info()[2]  # get the traceback object
    # tbinfo contains the error's line number and the code
    tbinfo = traceback.format_tb(tb)[0]
    line = tbinfo.split(", ")[1]
    msg = ("Geoprocessing error on **" + line + "** of " + filename + " "
                "in Linkage Mapper Version " + str(cfg.releaseNum) + ":")
    arcpy.AddError(msg)
    write_log(msg) #xxx
    dashline(1)
    msg=arcpy.GetMessages(2)
    arcpy.AddError(arcpy.GetMessages(2))
    write_log(msg)
    dashline()
    print_drive_warning()
    close_log_file()
    exit(1)


def exit_with_python_error(filename):
    """Handle python errors and provide details to user"""
    dashline()
    tb = sys.exc_info()[2]  # get the traceback object
    # tbinfo contains the error's line number and the code
    tbinfo = traceback.format_tb(tb)[0]
    line = tbinfo.split(", ")[1]

    err = traceback.format_exc().splitlines()[-1]
    msg = ("Python error on **" + line + "** of " + filename + " "
                "in Linkage Mapper Version " + str(cfg.releaseNum) + ":")
    arcpy.AddError(msg)
    arcpy.AddError(err)
    write_log(msg)
    write_log(err)
    close_log_file()
    exit(1)


def raise_error(msg):
    arcpy.AddError(msg)
    write_log(msg)
    close_log_file()
    exit(1)


def dashline(lspace=0):
    """Output dashed line in tool output dialog.

       0 = No empty line
       1 = Empty line before
       2 = Empty line after

    """
    if lspace == 1:
        gprint(' ')
    gprint('---------------------------------')
    if lspace == 2:
        gprint(' ')


############################################################################
## Circuitscape Functions ##################################################
############################################################################
@Retry(5)
def setCircuitscapeOptions():
    """Sets default options for calling Circuitscape.

    """
    options = {}
    options['data_type']='raster'
    options['version']='unknown'
    options['low_memory_mode']=False
    options['scenario']='pairwise'
    options['habitat_file']='(Browse for a habitat map file)'
    options['habitat_map_is_resistances']=True
    options['point_file']=('(Browse for file with '
                          'locations of focal points or areas)')
    options['point_file_contains_polygons']=True
    options['connect_four_neighbors_only']=False
    options['connect_using_avg_resistances']=True
    options['use_polygons']=False
    options['polygon_file']='(Browse for a short-circuit region file)'
    options['source_file']='(Browse for a current source file)'
    options['ground_file']='(Browse for a ground point file)'
    options['ground_file_is_resistances']=True
    options['use_unit_currents']=False
    options['use_direct_grounds']=False
    options['remove_src_or_gnd']='not entered'
    options['output_file']='(Choose a base name for output files)'
    options['write_cur_maps']=True
    options['write_cum_cur_map_only']=True
    options['log_transform_maps']=False
    options['write_volt_maps']=False
    options['solver']='cg+amg'
    options['compress_grids']=False
    options['print_timings']=False
    options['use_mask']=False
    options['mask_file']='None'
    options['use_included_pairs']=False
    options['included_pairs_file']='None'
    options['use_variable_source_strengths']=False
    options['variable_source_file']='None'
    options['write_max_cur_maps']=False
    options['set_focal_node_currents_to_zero']=True

    return options

def writeCircuitscapeConfigFile(configFile, options):
    """Creates a configuration file for calling Circuitscape.

    """
    config = ConfigParser.ConfigParser()

    sections={}
    section='Version'
    sections['version']=section

    section='Connection scheme for raster habitat data'
    sections['connect_four_neighbors_only']=section
    sections['connect_using_avg_resistances']=section

    section='Short circuit regions (aka polygons)'
    sections['use_polygons']=section
    sections['polygon_file']=section

    section='Options for advanced mode'
    sections['source_file']=section
    sections['ground_file']=section
    sections['ground_file_is_resistances']=section
    sections['use_unit_currents']=section
    sections['use_direct_grounds']=section
    sections['remove_src_or_gnd']=section

    section='Calculation options'
    sections['solver']=section
    sections['print_timings']=section
    sections['low_memory_mode']=section

    section='Output options'
    sections['output_file']=section
    sections['write_cur_maps']=section
    sections['write_cum_cur_map_only']=section
    sections['log_transform_maps']=section
    sections['write_volt_maps']=section
    sections['compress_grids']=section
    sections['write_max_cur_maps']=section
    sections['set_focal_node_currents_to_zero']=section

    section='Mask file'
    sections['use_mask']=section
    sections['mask_file']=section

    section='Options for pairwise and one-to-all and all-to-one modes'
    sections['use_included_pairs']=section
    sections['included_pairs_file']=section
    sections['point_file']=section
    sections['point_file_contains_polygons']=section

    section='Options for one-to-all and all-to-one modes'
    sections['use_variable_source_strengths']=section
    sections['variable_source_file']=section

    section='Habitat raster or graph'
    sections['habitat_file']=section
    sections['habitat_map_is_resistances']=section

    section="Circuitscape mode"
    sections['scenario']=section
    sections['data_type']=section

    if options['ground_file_is_resistances']=='not entered':
        options['ground_file_is_resistances'] = False
    if options['point_file_contains_polygons']=='not entered':
        options['point_file_contains_polygons'] = False

    for option in sections:
        try:
            config.add_section(sections[option])
        except Exception:
            pass
    for option in sections:
        config.set(sections[option], option, options[option])

    f = open(configFile, 'w')
    config.write(f)
    f.close()

def warn(string):
    arcpy.AddWarning(string)
    try:
        if cfg.LOGMESSAGES:
            write_log(string)
    except Exception:
        pass


class MEMORYSTATUSEX(ctypes.Structure):
    _fields_ = [("dwLength", ctypes.c_uint),
                ("dwMemoryLoad", ctypes.c_uint),
                ("ullTotalPhys", ctypes.c_ulonglong),
                ("ullAvailPhys", ctypes.c_ulonglong),
                ("ullTotalPageFile", ctypes.c_ulonglong),
                ("ullAvailPageFile", ctypes.c_ulonglong),
                ("ullTotalVirtual", ctypes.c_ulonglong),
                ("ullAvailVirtual", ctypes.c_ulonglong),
                ("sullAvailExtendedVirtual", ctypes.c_ulonglong),]

    def __init__(self):
        # have to initialize this to the size of MEMORYSTATUSEX
        self.dwLength = 2*4 + 7*8     # size = 2 ints, 7 longs
        return super(MEMORYSTATUSEX, self).__init__()

def get_mem():
    stat = MEMORYSTATUSEX()
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(stat))
    totMem = float(int(10 * float(stat.ullTotalPhys)/1073741824))/10
    availMem = float(int(10 * float(stat.ullAvailPhys)/1073741824))/10
    return totMem, availMem
