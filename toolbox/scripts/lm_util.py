#!/usr/bin/env python2.5

"""Contains functions called by linkage mapper and barrier mapper scripts."""  

__filename__ = "lm_util.py"
__version__ = "0.6.5"

import os
import sys
import time
import traceback
import ConfigParser
import shutil

import numpy as npy
import arcgisscripting

from lm_config import Config as Cfg

gp = Cfg.gp
gprint = gp.addmessage

def get_linktable_row(linkid, linktable):
    """Returns the linkTable row index for a given link ID"""
    try:
        # Most likely.  linkTables tend to be in order with no skipped links.
        if linktable[linkid - 1, Cfg.LTB_LINKID] == linkid:
            linktablerow = linkid - 1
            return linktablerow
        else:
            numLinks = linktable.shape[0]
            for linktablerow in range(0, numLinks):
                if linktable[linktablerow, Cfg.LTB_LINKID] == linkid:
                    return linktablerow
        return -1  # Not found
    except:
        raise_python_error(__filename__)


def get_link_type_desc(linktypecode):
    """For a linkType code returns description to attribute link and link maps.

    NOTE: These must map to LT codes in lm_config (eg LT_CPLK)

    """
    gp.addmessage('linktype', str(linktypecode))
    if linktypecode < 0:  # These are dropped links
        activelink = '0'
        if linktypecode == -1:
            linktypedesc = '"Not_nearest_neighbors"'
#        elif linktypecode == -2:
#            linktypedesc = '"Not_2nd_nearest_neighbors"'
#        elif linktypecode == -3:
#            linkTypeCodes='"Not_3rd_nearest_neighbors"'
#        elif linktypecode == -4:
#            linkTypeCodes='"Not_4th_nearest_neighbors"'
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
       # elif linktypecode == 11:
           # linktypedesc = '"1st_nearest_neighbors"'
       # elif linktypecode == 12:
           # linktypedesc = '"2nd_nearest_neighbors"'
       # elif linktypecode == 13:
           # linktypedesc = '"3rd_nearest_neighbors"'
       # elif linktypecode == 14:
           # linktypedesc = '"4th_nearest_neighbors"'
        elif linktypecode == 20:
            linktypedesc = '"Connects_constellations"'
        else:
            linktypedesc = '"Unknown_active"'
        #Add user retained?
    return activelink, linktypedesc


def get_links_from_core_pairs(linktable, firstCore, secondCore):
    """Given two cores, finds their matching row in the link table"""
    try:
        rows = npy.zeros((0), dtype="int32")
        numLinks = linktable.shape[0]
        for link in range(0, numLinks):
            corex = int(linktable[link, Cfg.LTB_CORE1])
            corey = int(linktable[link, Cfg.LTB_CORE2])
            if int(corex) == int(firstCore) and int(corey) == int(secondCore):
                rows = npy.append(rows, link)
            elif (int(corex) == int(secondCore) and int(corey) ==
                  int(firstCore)):
                rows = npy.append(rows, link)
        return rows
    except:
        raise_python_error(__filename__)


def drop_links(linktable, maxeud, mineud, maxcwd, mincwd,
               disableLeastCostNoVal):
    """Inactivates links that fail to meet min or max length criteria"""
    try:
        # dashline(1)
        numLinks = linktable.shape[0]
        numDroppedLinks = 0
        coreList = linktable[:, Cfg.LTB_CORE1:Cfg.LTB_CORE2 + 1]
        coreList = npy.sort(coreList)
        if disableLeastCostNoVal:
            for x in range(0, numLinks):
                linkid = str(int(linktable[x, Cfg.LTB_LINKID]))
                if linktable[x, Cfg.LTB_CWDIST] == -1:
                    #Check only enabled corridor links
                    if (linktable[x, Cfg.LTB_LINKTYPE] > 0):
                        corex = str(int(linktable[x, Cfg.LTB_CORE1]))
                        corey = str(int(linktable[x, Cfg.LTB_CORE2]))
                        gp.addmessage(
                            "The least-cost corridor between " + str(corex) +
                            " and " + str(corey) + " (link #" + linkid + ") "
                            "has an unknown length in cost distance units. "
                            "This means it is longer than the max "
                            "cost-weighted distance specified in the 'Calc "
                            "CWDs' script OR it passes through NODATA cells "
                            "and will be dropped.\n")
                        # Disable link
                        linktable[x, Cfg.LTB_LINKTYPE] = Cfg.LT_TLLC
                        numDroppedLinks = numDroppedLinks + 1

        # Check for corridors that are too long in Euclidean or cost-weighted
        # distance
        if maxeud is not None or maxcwd is not None:
            for x in range(0, numLinks):
                linkid = str(int(linktable[x, Cfg.LTB_LINKID]))
                if maxeud is not None:
                    if linktable[x, Cfg.LTB_EUCDIST] > maxeud:
                        # Check only enabled corridor links
                        if (linktable[x, Cfg.LTB_LINKTYPE] > 0):
                            corex = str(int(coreList[x, 0]))
                            corey = str(int(coreList[x, 1]))
                            gp.addmessage("Link #" + linkid +
                                          " connecting cores " + str(corex) +
                                          " and " + str(corey) + " is  " +
                                          str(linktable[x, Cfg.LTB_EUCDIST]) +
                                          " units long- too long in "
                                          "Euclidean distance.")
                            # Disable link
                            linktable[x, Cfg.LTB_LINKTYPE] = Cfg.LT_TLEC
                            numDroppedLinks = numDroppedLinks + 1

                    if maxcwd is not None:
                        # Check only enabled corridor links
                        if (linktable[x, Cfg.LTB_LINKTYPE] > 0):
                            # Check for -1 Cfg.LTB_CWDIST
                            if (linktable[x, Cfg.LTB_CWDIST] > maxcwd):
                                corex = str(int(linktable[x, Cfg.LTB_CORE1]))
                                corey = str(int(linktable[x, Cfg.LTB_CORE2]))
                                gp.addmessage(
                                    "Link #" + linkid + " connecting cores " +
                                    str(corex) + " and " + str(corey) +
                                    " is " +
                                    str(linktable[x, Cfg.LTB_CWDIST]) +
                                    " units long- too long in cost-distance "
                                    "units.")
                                #  Disable link
                                linktable[x, Cfg.LTB_LINKTYPE] = Cfg.LT_TLLC
                                numDroppedLinks = numDroppedLinks + 1

        if mineud is not None or mincwd is not None:
            for x in range(0, numLinks):
                linkid = str(int(linktable[x, Cfg.LTB_LINKID]))
                # Check only enabled corridor links
                if (linktable[x, Cfg.LTB_LINKTYPE] > 0):
                    if mineud is not None:
                        if linktable[x, Cfg.LTB_EUCDIST] < mineud:
                            corex = str(int(coreList[x, 0]))
                            corey = str(int(coreList[x, 1]))
                            gp.addmessage(
                                "Link #" + linkid + " connecting cores " +
                                str(corex) + " and " + str(corey) + " is "
                                "only " + str(linktable[x, Cfg.LTB_EUCDIST]) +
                                " units long- too short in Euclidean "
                                "distance.")
                            # Disable link
                            linktable[x, Cfg.LTB_LINKTYPE] = Cfg.LT_TSEC
                            numDroppedLinks = numDroppedLinks + 1

                    if mincwd is not None:
                        if ((linktable[x, Cfg.LTB_CWDIST] < mincwd) and
                            (linktable[x, Cfg.LTB_CWDIST]) != -1):
                            if (linktable[x, Cfg.LTB_LINKTYPE] > 0):
                                corex = str(int(linktable[x, Cfg.LTB_CORE1]))
                                corey = str(int(linktable[x, Cfg.LTB_CORE2]))
                                gp.addmessage(
                                    "Link #" + linkid + " connecting cores " +
                                    str(corex) + " and " + str(corey) +
                                    " is only " +
                                    str(linktable[x, Cfg.LTB_CWDIST]) +
                                    " units long- too short in cost distance "
                                    "units.")
                                # Disable link
                                linktable[x, Cfg.LTB_LINKTYPE] = Cfg.LT_TSLC
                                numDroppedLinks = numDroppedLinks + 1
        return linktable, numDroppedLinks
    except:
        raise_python_error(__filename__)


def get_zonal_minimum(dbfFile):
    """Finds the minimum value in a table generated by zonal statistics"""
    try:
        rows = gp.searchcursor(dbfFile)
        row = rows.Next()
        if row is not None:
            coreMin = row.Min
            while row:
                if coreMin > row.Min:
                    coreMin = row.Min
                row = rows.next()
        else:
            coreMin = None
        del row
        del rows
        return coreMin
    except arcgisscripting.ExecuteError:
        raise_geoproc_error(__filename__)
    except:
        raise_python_error(__filename__)


def get_core_list():
    """Returns a list of core area IDs from polygon file"""
    try:
        # Get the number of cores
        # FIXME: I think this returns number of shapes, not number of unique
        # cores.
        coreCount = int(gp.GetCount_management(Cfg.COREFC).GetOutput(0))
        # Get core data into numpy array
        coreList = npy.zeros((coreCount, 2))
        cur = gp.SearchCursor(Cfg.COREFC)
        row = cur.Next()
        i = 0
        while row:
            coreList[i, 0] = row.GetValue(Cfg.COREFN)
            coreList[i, 1] = row.GetValue(Cfg.COREFN)
            row = cur.Next()
            i = i + 1
        del cur, row
    except arcgisscripting.ExecuteError:
        raise_geoproc_error(__filename__)
    except:
        raise_python_error(__filename__)
    return coreList


def get_core_targets(core, linktable):
    """Returns a list of other core areas the core area is connected to."""
    try:
        targetList = npy.zeros((len(linktable), 2), dtype="int32")
        # possible targets. 1st column = Cfg.LTB_CORE1
        targetList[:, 0] = linktable[:, Cfg.LTB_CORE1]
        # possible targets. 2nd column = Cfg.LTB_CORE2
        targetList[:, 1] = linktable[:, Cfg.LTB_CORE2]
        # Copy of Cfg.LTB_LINKTYPE column
        validPair = linktable[:, Cfg.LTB_LINKTYPE]
        validPair = npy.where(validPair == Cfg.LT_CORR, 1, 0)  # map corridor.
        targetList[:, 0] = npy.multiply(targetList[:, 0], validPair)
        targetList[:, 1] = npy.multiply(targetList[:, 1], validPair)

        rows, cols = npy.where(targetList == int(core))
        targetList = targetList[rows, 1 - cols]
        targetList = npy.unique(npy.asarray(targetList))
    except arcgisscripting.ExecuteError:
        raise_geoproc_error(__filename__)
    except:
        raise_python_error(__filename__)
    return targetList


def elapsed_time(start_time):
    """Returns elapsed time given a start time"""
    try:
        now = time.clock()
        elapsed = now - start_time
        secs = int(elapsed)
        mins = int(elapsed / 60)
        hours = int(mins / 60)
        mins = mins - hours * 60
        secs = secs - mins * 60 - hours * 3600
        if mins == 0:
            gp.addmessage('That took ' + str(secs) + ' seconds.\n')
        elif hours == 0:
            gp.addmessage('That took ' + str(mins) + ' minutes and ' +
                              str(secs) + ' seconds.\n')
        else:
            gp.addmessage('That took ' + str(hours) + ' hours ' +
                              str(mins) + ' minutes and ' + str(secs) +
                              ' seconds.\n')
        return now
    except:
        raise_python_error(__filename__)


def report_pct_done(current, goal, last):
    """Reports percent done"""
    try:
        goal = float(goal) 
        pctDone = ((float(current) / goal) * 100)
        pctDone = 10 * (npy.floor(pctDone/10))
        if pctDone - last >= 10:
            gp.addmessage(str(int(pctDone)) + " percent done")        
            return 10*int((npy.floor(pctDone/10)))
        else:
            return last
    except:
        raise_python_error(__filename__)


def report_links(linktable):
    """Prints number of links in a link table"""
    try:
        numLinks = linktable.shape[0]
        gp.addmessage('There are ' + str(numLinks) + ' links in the '
                          'table.')
        linkTypes = linktable[:, Cfg.LTB_LINKTYPE]
        numCorridorLinks = sum(linkTypes == Cfg.LT_CORR) + sum(linkTypes ==
                                                               Cfg.LT_NNC)
        numComponentLinks = sum(linkTypes == Cfg.LT_CLU)
        if numComponentLinks > 0:
            gp.addmessage('This includes ' + str(numCorridorLinks) +
                              ' potential corridor links and ' +
                          str(numComponentLinks) + ' component links.')
        elif numCorridorLinks > 0:
            gp.addmessage('This includes ' + str(numCorridorLinks) +
                              ' potential corridor links.')
        else:
            numCorridorLinks = 0
            gp.addmessage('\n***NOTE: There are NO corridors to map!')
            dashline(2)

    except arcgisscripting.ExecuteError:
        raise_geoproc_error(__filename__)
    except:
        raise_python_error(__filename__)
    return numCorridorLinks + numComponentLinks


############################################################################
## Adjacency and allocation functions ##########################
############################################################################

def get_adj_using_shift_method(alloc):
    """Returns table listing adjacent core areas using a shift method.

    The method involves shifting the allocation grid one pixel and then looking
    for pixels with different allocations across shifted grids.

    """
    cellSize = gp.Describe(alloc).MeanCellHeight
    gp.CellSize = cellSize

    posShift = gp.CellSize
    negShift = -1 * float(gp.CellSize)

    gp.workspace = Cfg.SCRATCHDIR

    gp.addmessage('Calculating adjacencies crossing horizontal allocation '
                  'boundaries...')
    start_time = time.clock()
    gp.Shift_management(alloc, "alloc_r", posShift, "0")

    alloc_r = "alloc_r"
    adjTable_r = get_allocs_from_shift(gp.workspace, alloc, alloc_r)
    start_time = elapsed_time(start_time)

    gp.addmessage('Calculating adjacencies crossing upper-left diagonal '
                      'allocation boundaries...')
    gp.Shift_management(alloc, "alloc_ul", negShift, posShift)

    alloc_ul = "alloc_ul"
    adjTable_ul = get_allocs_from_shift(gp.workspace, alloc, alloc_ul)
    start_time = elapsed_time(start_time)

    gp.addmessage('Calculating adjacencies crossing upper-right diagonal '
                      'allocation boundaries...')
    gp.Shift_management(alloc, "alloc_ur", posShift, posShift)

    alloc_ur = "alloc_ur"
    adjTable_ur = get_allocs_from_shift(gp.workspace, alloc, alloc_ur)
    start_time = elapsed_time(start_time)

    gp.addmessage('Calculating adjacencies crossing vertical allocation '
                      'boundaries...')
    gp.Shift_management(alloc, "alloc_u", "0", posShift)

    alloc_u = "alloc_u"
    adjTable_u = get_allocs_from_shift(gp.workspace, alloc, alloc_u)
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
    except arcgisscripting.ExecuteError:
        raise_geoproc_error(__filename__)
    except:
        raise_python_error(__filename__)


def get_allocs_from_shift(workspace, alloc, alloc_sh):
    """Returns a table of adjacent allocation zones using grid shift method"""
    try:
        combine_ras = os.path.join(gp.workspace, "combine")
        count = 0
        statement = ('gp.SingleOutputMapAlgebra_sa("combine(" + alloc + '
                     '", " + alloc_sh + ")", combine_ras, alloc, alloc_sh)')
        while True:
            try:
                exec statement
            except:
                count, tryAgain = hiccup_test(count, statement)
                if not tryAgain:
                    exec statement
            else:
                break
        allocLookupTable = get_alloc_lookup_table(gp.workspace,
                                                  combine_ras)
        return allocLookupTable[:, 1:3]

    except arcgisscripting.ExecuteError:
        raise_geoproc_error(__filename__)
    except:
        raise_python_error(__filename__)


# Try this in script 4 too... and anything else with minras...
def get_alloc_lookup_table(workspace, combine_ras):
    """Returns a table of adjacent allocation zones.

    Requires a raster with allocation zone attributes.

    """
    try:
        fldlist = gp.listfields(combine_ras)

        valFld = fldlist[1].name
        allocFld = fldlist[3].name
        allocFld_sh = fldlist[4].name

        allocLookupTable = npy.zeros((0, 3), dtype="int32")
        appendRow = npy.zeros((1, 3), dtype="int32")

        rows = gp.searchcursor(combine_ras)
        row = rows.next()
        while row:
            alloc = row.getvalue(allocFld)
            alloc_sh = row.getvalue(allocFld_sh)
            if alloc != alloc_sh:
                appendRow[0, 0] = row.getvalue(valFld)
                appendRow[0, 1] = alloc
                appendRow[0, 2] = alloc_sh
                allocLookupTable = npy.append(allocLookupTable, appendRow,
                                              axis=0)
            row = rows.next()
        del row
        del rows

        return allocLookupTable
    except arcgisscripting.ExecuteError:
        raise_geoproc_error(__filename__)
    except:
        raise_python_error(__filename__)


############################################################################
## Bounding Circle Functions ##########################
############################################################################
def new_extent(fc, field, value):
    """Returns the maximum area extent of features where field == value"""
    try:
        shapeFieldName = gp.describe(fc).shapefieldname

        # searchRows = gp.searchcursor(fc, " "'"core_ID"'" = " +
        #                                  str(currentID))
        searchRows = gp.searchcursor(fc, field + ' = ' + str(value))
        searchRow = searchRows.next()
        # get the 1st features extent
        extentObj = searchRow.getvalue(shapeFieldName).extent
        xMin = extentObj.xmin
        yMin = extentObj.ymin
        xMax = extentObj.xmax
        yMax = extentObj.yMax
        searchRow = searchRows.next()  # now move on to the other features
        while searchRow:
            extentObj = searchRow.getvalue(shapeFieldName).extent
            if extentObj.xmin < xMin:
                xMin = extentObj.xmin
            if extentObj.ymin < yMin:
                yMin = extentObj.ymin
            if extentObj.xmax > xMax:
                xMax = extentObj.xmax
            if extentObj.ymax > yMax:
                yMax = extentObj.ymax
            searchRow = searchRows.next()
        del searchRow
        del searchRows
    except arcgisscripting.ExecuteError:
        raise_geoproc_error(__filename__)
    except:
        raise_python_error(__filename__)

    # return [xMin,yMin,xMax,yMax]
    # "%s %s %s %s" % tuple(lstExt)
    # strExtent = (str(xMin) + ' ' + str(yMin) + ' ' + str(xMax) + ' ' +
    #              str(yMax))
    # strExtent
    return  str(xMin), str(yMin), str(xMax), str(yMax)


def get_centroids(shapefile, field):
    """Returns centroids of features"""
    try:
        pointArray = npy.zeros((0, 3), dtype="float32")
        xyCumArray = npy.zeros((0, 3), dtype="float32")
        xyArray = npy.zeros((1, 3), dtype="float32")
        rows = gp.SearchCursor(shapefile)
        row = rows.Next()
        while row:
            feat = row.shape
            center = feat.Centroid
            center = str(center)
            xy = center.split(" ")
            xyArray[0, 0] = float(xy[0])
            xyArray[0, 1] = float(xy[1])
            value = row.GetValue(field)
            xyArray[0, 2] = int(value)
            xyCumArray = npy.append(xyCumArray, xyArray, axis=0)
            row = rows.Next()
        del row, rows
        pointArray = npy.append(pointArray, xyCumArray, axis=0)

        return pointArray

    except arcgisscripting.ExecuteError:
        raise_geoproc_error(__filename__)
    except:
        raise_python_error(__filename__)


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
    except arcgisscripting.ExecuteError:
        raise_geoproc_error(__filename__)
    except:
        raise_python_error(__filename__)

    return circlePointData


def get_extent_box_coords(fieldValue=None):
    """Get coordinates of bounding box that contains selected features"""
    try:
        # get all features, not just npy.where Cfg.COREFN = fieldValue
        if fieldValue is None:
            fieldValue = 1
            desc = gp.Describe
            extent = desc(Cfg.FCORES).extent
            lr = extent.lowerright
            ul = extent.upperleft
            ulx = ul.x
            uly = ul.y
            lrx = lr.x
            lry = lr.y
        else:
            ulx, lry, lrx, uly = new_extent(Cfg.FCORES, Cfg.COREFN, fieldValue)

        ulx = float(ulx)
        lrx = float(lrx)
        uly = float(uly)
        lry = float(lry)
        boxData = npy.zeros((1, 5), dtype='float32')
        boxData[0, :] = [fieldValue, ulx, lrx, uly, lry]
    except arcgisscripting.ExecuteError:
        raise_geoproc_error(__filename__)
    except:
        raise_python_error(__filename__)

    return boxData


def make_points(workspace, pointArray, outFC):
    """Creates a shapefile with points specified by coordinates in pointArray

       outFC is just the filename, not path.
       pointArray is x,y,corex,corey,radius

    """
    try:
        wkspbefore = gp.workspace
        gp.workspace = workspace
        if gp.exists(outFC):
            gp.delete_management(outFC)
        gp.CreateFeatureclass_management(workspace, outFC, "POINT")
        #for field in fieldArray:
        if pointArray.shape[1] > 3:
            gp.addfield(outFC, "corex", "SHORT")
            gp.addfield(outFC, "corey", "SHORT")
            gp.addfield(outFC, "radius", "DOUBLE")
            gp.addfield(outFC, "cores_x_y", "TEXT")
        else:
            gp.addfield(outFC, "XCoord", "DOUBLE")
            gp.addfield(outFC, "YCoord", "DOUBLE")
            gp.addfield(outFC, Cfg.COREFN, "SHORT")
        rows = gp.InsertCursor(outFC)

        numPoints = pointArray.shape[0]
        for i in range(numPoints):
            point = gp.CreateObject("Point")
            point.ID = i
            point.X = float(pointArray[i, 0])
            point.Y = float(pointArray[i, 1])
            row = rows.NewRow()
            row.shape = point
            row.SetValue("ID", i)
            if pointArray.shape[1] > 3:
                row.SetValue("corex", int(pointArray[i, 2]))
                row.SetValue("corey", int(pointArray[i, 3]))
                row.SetValue("cores_x_y", str(int(pointArray[i, 2])) + '_' +
                             str(int(pointArray[i, 3])))
                row.SetValue("radius", float(pointArray[i, 4]))
            else:
                row.SetValue("XCoord", float(pointArray[i, 0]))
                row.SetValue("YCoord", float(pointArray[i, 1]))
                row.SetValue(Cfg.COREFN, float(pointArray[i, 2]))

            rows.InsertRow(row)
            del row
            del point
        del rows
        gp.workspace = wkspbefore

    except arcgisscripting.ExecuteError:
        raise_geoproc_error(__filename__)
    except:
        raise_python_error(__filename__)
    return


############################################################################
## LCP Shapefile Functions #################################################
############################################################################

def create_lcp_shapefile(linktable, sourceCore, targetCore, lcpLoop):
    """Creates lcp shapefile.

    Shows locations of least-cost path lines attributed with corridor
    info/status.

    """
    try:
        rows = get_links_from_core_pairs(linktable, sourceCore, targetCore)
        link = rows[0]

        lcpline = os.path.join(Cfg.SCRATCHDIR, "lcpline.shp")
        gp.RasterToPolyline_conversion("lcp", lcpline, "NODATA", "",
                                           "NO_SIMPLIFY")

        lcplineDslv = os.path.join(Cfg.SCRATCHDIR, "lcplineDslv.shp")
        gp.Dissolve_management(lcpline, lcplineDslv)

        gp.AddField_management(lcplineDslv, "Link_ID", "SHORT", "5")
        gp.CalculateField_management(lcplineDslv, "Link_ID",
                                         int(linktable[link, Cfg.LTB_LINKID]))

        linktypecode = linktable[link, Cfg.LTB_LINKTYPE]
        activelink, linktypedesc = get_link_type_desc(linktypecode)

        gp.AddField_management(lcplineDslv, "Active", "SHORT")
        gp.CalculateField_management(lcplineDslv, "Active", activelink)

        gp.AddField_management(lcplineDslv, "Link_Info", "TEXT")
        gp.CalculateField_management(lcplineDslv, "Link_Info",
                                         linktypedesc)

        gp.AddField_management(lcplineDslv, "From_Core", "SHORT", "5")
        gp.CalculateField_management(lcplineDslv, "From_Core",
                                         int(sourceCore))
        gp.AddField_management(lcplineDslv, "To_Core", "SHORT", "5")
        gp.CalculateField_management(lcplineDslv, "To_Core",
                                         int(targetCore))

        gp.AddField_management(lcplineDslv, "Euc_Dist", "DOUBLE", "10",
                                   "2")
        gp.CalculateField_management(lcplineDslv, "Euc_Dist",
                                         linktable[link, Cfg.LTB_EUCDIST])

        gp.AddField_management(lcplineDslv, "CW_Dist", "DOUBLE", "10", "2")
        gp.CalculateField_management(lcplineDslv, "CW_Dist",
                                         linktable[link, Cfg.LTB_CWDIST])
        gp.AddField_management(lcplineDslv, "LCP_Length", "DOUBLE", "10",
                                   "2")
        rows = gp.UpdateCursor(lcplineDslv)
        row = rows.Next()
        while row:
            feat = row.shape
            lcpLength = int(feat.length)
            row.SetValue("LCP_Length", lcpLength)
            rows.UpdateRow(row)
            row = rows.Next()
        del row, rows

        try:
            distRatio1 = (float(linktable[link, Cfg.LTB_CWDIST])
                        / float(linktable[link, Cfg.LTB_EUCDIST]))
        except ZeroDivisionError:
            distRatio1 = -1
            
        gp.AddField_management(lcplineDslv, "cwd2Euc_R", "DOUBLE", "10",
                                   "2")
        gp.CalculateField_management(lcplineDslv, "cwd2Euc_R", distRatio1)

        try:
            distRatio2 = (float(linktable[link, Cfg.LTB_CWDIST]) 
                        / float(lcpLength))
        except ZeroDivisionError:
            distRatio2 = -1
            
        gp.AddField_management(lcplineDslv, "cwd2Path_R", "DOUBLE", "10",
                                   "2")
        gp.CalculateField_management(lcplineDslv, "cwd2Path_R", distRatio2)

        lcpLoop = lcpLoop + 1
        lcpShapefile = os.path.join(Cfg.DATAPASSDIR, "lcpLines_s3.shp")
        gp.RefreshCatalog(Cfg.DATAPASSDIR)
        if lcpLoop == 1:
            if gp.Exists(lcpShapefile):
                try:
                    gp.Delete(lcpShapefile)

                except:
                    dashline(1)
                    msg = ('ERROR: Could not remove LCP shapefile ' +
                           lcpShapefile + '. Was it open in ArcMap?\n You may '
                           'need to re-start ArcMap to release the file lock.')
                    gp.AddError(msg)
                    exit(1)

            gp.copy_management(lcplineDslv, lcpShapefile)
        else:
            gp.Append_management(lcplineDslv, lcpShapefile, "TEST")
 
        return lcpLoop

    except arcgisscripting.ExecuteError:
        raise_geoproc_error(__filename__)
    except:
        raise_python_error(__filename__)


def update_lcp_shapefile(linktable, lastStep, thisStep):
    """Updates lcp shapefiles with new link information/status"""
    try:
        lcpShapefile = os.path.join(Cfg.DATAPASSDIR, "lcpLines_s" +
                                    str(thisStep) + ".shp")

        numLinks = linktable.shape[0]
        # g1' g2' THEN c1 c2

        if thisStep > 5:
            linkTableTemp = linktable
            #linkTableTemp[:, Cfg.LTB_CURRENT] = linkTableTemp[:, Cfg.LTB_LCPLEN]      
        else:    
            extraCols = npy.zeros((numLinks, 3), dtype="float64")
            linkTableTemp = npy.append(linktable, extraCols, axis=1)
            del extraCols
            linkTableTemp[:, Cfg.LTB_LCPLEN] = -1
            linkTableTemp[:, Cfg.LTB_CWDEUCR] = -1
            linkTableTemp[:, Cfg.LTB_CWDPATHR] = -1
        

        if lastStep != thisStep:
            if thisStep == 5:
                oldLcpShapefile = os.path.join(
                    Cfg.DATAPASSDIR, "lcpLines_s" + str(lastStep) + ".shp")
                # If last step wasn't step 4 then must be step 3
                if not gp.exists(oldLcpShapefile):
                    # step 3
                    oldLcpShapefile = os.path.join(
                        Cfg.DATAPASSDIR, "lcpLines_s" + str(lastStep - 1) +
                        ".shp")
            
            elif thisStep > 5:
                oldLcpShapefile = os.path.join(
                    Cfg.DATAPASSDIR, "lcpLines_s" + str(lastStep) + ".shp")
                if not gp.exists(oldLcpShapefile):
                    # step 4
                    oldLcpShapefile = os.path.join(
                        Cfg.DATAPASSDIR, "lcpLines_s" + str(lastStep - 1) +
                        ".shp")
                    if not gp.exists(oldLcpShapefile):
                        # step 3
                        oldLcpShapefile = os.path.join(
                            Cfg.DATAPASSDIR, "lcpLines_s" + str(lastStep - 1) +
                            ".shp")
            else:
                oldLcpShapefile = os.path.join(
                    Cfg.DATAPASSDIR, "lcpLines_s" + str(lastStep) + ".shp")

            if gp.exists(lcpShapefile):
                try:
                    gp.delete_management(lcpShapefile)
                except:
                    dashline(1)
                    msg = ('ERROR: Could not remove LCP shapefile ' +
                           lcpShapefile + '. Is it open in ArcMap?\n You may '
                           'need to re-start ArcMap to release the file lock.')
                    gp.AddError(msg)
                    exit(1)
            gp.copy_management(oldLcpShapefile, lcpShapefile)
            if thisStep > 5:
                gp.AddField_management(lcpShapefile, "Eff_Resist", "FLOAT") ### 
                gp.AddField_management(lcpShapefile, "cwd2EffR_r", "FLOAT")
                gp.AddField_management(lcpShapefile, "CF_Central", "FLOAT") ### 
        rows = gp.UpdateCursor(lcpShapefile)
        row = rows.Next()
        line = 0
        while row:
            linkid = row.getvalue("Link_ID")
            linktypecode = linktable[linkid - 1, Cfg.LTB_LINKTYPE]
            activelink, linktypedesc = get_link_type_desc(linktypecode)
            row.SetValue("Link_Info", linktypedesc)
            row.SetValue("Active", activelink)
            if thisStep > 5:
                current = linkTableTemp[linkid - 1, Cfg.LTB_CURRENT] 
                effResist = linkTableTemp[linkid - 1, Cfg.LTB_EFFRESIST] 
                CWDTORRatio = linkTableTemp[linkid - 1, Cfg.LTB_CWDTORR]
                #fixme: linkid - 1 assumes linktable ordered
                row.SetValue("Eff_Resist", effResist)
                row.SetValue("cwd2EffR_r",CWDTORRatio)
                row.SetValue("CF_Central", current)
            rows.UpdateRow(row)

            linktablerow = get_linktable_row(linkid, linkTableTemp)
            linkTableTemp[linktablerow, Cfg.LTB_LCPLEN] = row.getvalue(
                "LCP_Length")
            linkTableTemp[linktablerow, Cfg.LTB_CWDEUCR] = row.getvalue(
                "cwd2Euc_R")
            linkTableTemp[linktablerow, Cfg.LTB_CWDPATHR] = row.getvalue(
                "cwd2Path_R")
            row = rows.Next()
            line = line + 1
        # delete cursor and row points to remove locks on the data
        del row, rows


        PREFIX = Cfg.PREFIX
        outputLcpShapefile = os.path.join(Cfg.OUTPUTDIR, PREFIX + "_lcpLines_s" +
                                          str(thisStep) + ".shp")
        if gp.exists(outputLcpShapefile):
            try:
                gp.delete_management(outputLcpShapefile)
            except:
                dashline(1)
                msg = ('ERROR: Could not remove lcp shapefile from output '
                       'directory: ' + outputLcpShapefile +
                       '. Is it open in ArcMap?\n You may '
                       'need to re-start ArcMap to release the file lock.')
                gp.AddError(msg)
                exit(1)
        gp.copy_management(lcpShapefile, outputLcpShapefile)

        return linkTableTemp

    except arcgisscripting.ExecuteError:
        raise_geoproc_error(__filename__)
    except:
        raise_python_error(__filename__)


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
    except:
        raise_python_error(__filename__)
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
    except:
        raise_python_error(__filename__)
    return A[keeprows][:, keepcols]


def delete_row_col(A, delrow, delcol):
    """Deletes rows and columns from a matrix

    From gapdt.py by Viral Shah

    """
    try:
        m = A.shape[0]
        n = A.shape[1]

        keeprows = npy.delete(npy.arange(0, m), delrow)
        keepcols = npy.delete(npy.arange(0, n), delcol)
    except:
        raise_python_error(__filename__)
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
    except:
        raise_python_error(__filename__)


def relabel(oldlabel, offset=0):
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


def conditional_hooking(D, star, u, v):
    """Utility for components code

    From gapdt.py by Viral Shah

    """
    Du = D[u]
    Dv = D[v]

    hook = npy.where((Du == D[Du]) & (Dv < Du))
    Du = Du[hook]
    Dv = Dv[hook]

    D[Du] = Dv
    return D


def unconditional_hooking(D, star, u, v):
    """Utility for components code

    From gapdt.py by Viral Shah

    """
    Du = D[u]
    Dv = D[v]

    hook = npy.where((star[u] == 1) & (Du != Dv))
    D[Du[hook]] = Dv[hook]

    return D


def check_stars(D, star):
    """Utility for components code

    From gapdt.py by Viral Shah

    """
    star[:] = 1
    notstars = npy.where(D != D[D])
    star[notstars] = 0
    star[D[D[notstars]]] = 0
    star = star[D]
    return star


def pointer_jumping(D):
    """Utility for components code

    From gapdt.py by Viral Shah

    """
    n = D.size
    Dold = npy.zeros(n, dtype='int32')

    while any(Dold != D):
        Dold = D
        D = D[D]
    return D


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
    except:
        raise_python_error(__filename__)


############################################################################
## Output Functions ########################################################
############################################################################
def write_link_table(linktable, outlinkTableFile, *inLinkTableFile):
    """Writes link tables to pass link data between steps """
    try:

        numLinks = linktable.shape[0]
        outFile = open(outlinkTableFile, "w")
        # gprint('linktable')
        # gprint(str(linktable.shape[1]))
        # gprint(str(linktable.astype('int32')))
        
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
                           "cwdToEucRatio,cwdToPathRatio,Eff_Resist,CWDTORRatio,CF_Centrality\n")

            for x in range(0, numLinks):
                for y in range(0, 15):
                    outFile.write(str(linktable[x, y]) + ",")
                outFile.write(str(linktable[x, 15]))
                outFile.write("\n")


        outFile.write("# Linkage Mapper Version " + __version__)
        outFile.write("\n# ---Run Settings---")
        outFile.write("\n# Project Directory: " + Cfg.PROJECTDIR)
        if Cfg.TOOL == 'linkage_mapper': 
            outFile.write("\n# Core Area Feature Class: " + Cfg.COREFC)

            outFile.write("\n# Core Area Field Name: " + Cfg.COREFN)
            outFile.write("\n# Resistance Raster: " + Cfg.RESRAST_IN)
            outFile.write("\n# Step 1 - Identify Adjacent Core Areas: " +
                           str(Cfg.STEP1))
            outFile.write("\n# Step 1 Adjacency Method Includes Cost-Weighted "
                           "Distance: " + str(Cfg.S1ADJMETH_CW))
            outFile.write("\n# Step 1 Adjacency Method Includes Euclidean "
                           "Distance: " + str(Cfg.S1ADJMETH_EU))
            outFile.write("\n# Step 2 - Construct a Network of Core Areas: " +
                           str(Cfg.STEP2))
            outFile.write("\n# Conefor Distances Text File: " + 
                          str(Cfg.S2EUCDISTFILE))
            outFile.write("\n# Network Adjacency Method Includes Cost-Weighted "
                           "Distance: " + str(Cfg.S2ADJMETH_CW))
            outFile.write("\n# Network Adjacency Method Includes Euclidean "
                           "Distance: " + str(Cfg.S2ADJMETH_EU))
            outFile.write("\n# Step 3 - Calculate Cost-Weighted Distances and "
                           "Least-Cost Paths: " + str(Cfg.STEP3))
            outFile.write("\n# Drop Corridors that Intersect Core Areas: "
                           + Cfg.S3DROPLCCS)
            outFile.write("\n# Step 4 - Refine Network" + str(Cfg.STEP4))
            outFile.write("\n# Option A - Number of Connected Nearest Neighbors: "
                           + str(Cfg.S4MAXNN))
            outFile.write("\n# Option B - Nearest Neighbor Measurement Unit is "
                           "Cost-Weighted Distance: " + str(Cfg.S4DISTTYPE_CW))
            outFile.write("\n# Option C - Connect Neighboring Constellations : "
                           + str(Cfg.S4CONNECT))
            outFile.write("\n# Step 5 - Calculate Normalize and Mosaic "
                           "Corridors: " + str(Cfg.STEP5))
            outFile.write("\n# Bounding Circles Buffer Distance: "
                           + str(Cfg.BUFFERDIST))
            outFile.write("\n# Maximum Cost-Weighted Corridor Distance: "
                           + str(Cfg.MAXCOSTDIST))
            outFile.write("\n# Maximum Euclidean Corridor Distance: "
                           + str(Cfg.MAXEUCDIST))
            outFile.write("\n# Minimum Cost-Weighted Corridor Distance: "
                           + str(Cfg.MINCOSTDIST))
            outFile.write("\n# Minimum Euclidean Corridor Distance: "
                           + str(Cfg.MINEUCDIST))


        elif Cfg.TOOL == 'pinchpoint_mapper':
            outFile.write("\n# Pinchpoints Analyzed: "
                           + str(Cfg.DOPINCH))
            outFile.write("\n# Resistance Raster: " + Cfg.RESRAST)
            outFile.write("\n# CWD Cutoff Distance: "
                           + str(Cfg.CWDCUTOFF))
            outFile.write("\n# Resistance Raster ValuesS quared for "
                            "Circuitscape Analyses: "+ str(Cfg.SQUARERESISTANCES))
            outFile.write("\n# Network Centrality Analyzed: "
                          + str(Cfg.DOCENTRALITY))
            outFile.write("\n# Core Area Feature Class: " + Cfg.COREFC)

            outFile.write("\n# Core Area Field Name: " + Cfg.COREFN)
            FN = str(inLinkTableFile)
            FN = FN.replace("('", "")
            FN = FN.replace("',)", "")
            outFile.write("\n# Link Data Derived from: " + FN)
                           
                           
        outFile.close()
    except arcgisscripting.ExecuteError:
        raise_geoproc_error(__filename__)
    except:
        raise_python_error(__filename__)
    return


def write_adj_file(outcsvfile, adjTable):
    """Outputs adjacent core areas to pass adjacency info between steps"""
    outfile = open(outcsvfile, "w")
    outfile.write("#Edge" + "," + str(Cfg.COREFN) + "," + str(Cfg.COREFN) +
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
        gp.OverwriteOutput = True
        gp.workspace = Cfg.OUTPUTDIR
        gp.RefreshCatalog(Cfg.OUTPUTDIR)
        linktable = load_link_table(linkTableFile)
        
        #linktable = npy.loadtxt(linkTableFile, dtype ='Float64', comments='#',
        #                    delimiter=',')
        numLinks = linktable.shape[0]
        gp.toolbox = "management"

        coresForLinework = "cores_for_linework.shp"

        # Preferred method to get geometric center
        pointArray = get_centroids(Cfg.COREFC, Cfg.COREFN)

        make_points(gp.workspace, pointArray, coresForLinework)
        numLinks = linktable.shape[0]
        # rows,cols = npy.where(
        #   linktable[:, Cfg.LTB_LINKTYPE:Cfg.LTB_LINKTYPE + 1] == Cfg.LT_CORR)

        coreLinks = linktable
               
        # create coreCoords array, with geographic centers of cores
        coreCoords = npy.zeros(pointArray.shape, dtype='float64')
        coreCoords[:, 0] = pointArray[:, 2]
        coreCoords[:, 1] = pointArray[:, 0]
        coreCoords[:, 2] = pointArray[:, 1]

        # Create linkCoords array
        linkCoords = npy.zeros((len(coreLinks), 13))
        linkCoords[:, 0] = coreLinks[:, Cfg.LTB_LINKID]
        linkCoords[:, 1:3] = npy.sort(
            coreLinks[:, Cfg.LTB_CORE1:Cfg.LTB_CORE2 + 1])
        linkCoords[:, 3:5] = coreLinks[:, Cfg.LTB_EUCDIST:Cfg.LTB_CWDIST + 1]
        linkCoords[:, 9] = coreLinks[:, Cfg.LTB_LINKTYPE]
        if step > 5:
            linkCoords[:, 10] = coreLinks[:, Cfg.LTB_EFFRESIST]
            linkCoords[:, 11] = coreLinks[:, Cfg.LTB_CWDTORR]
            linkCoords[:, 12] = coreLinks[:, Cfg.LTB_CURRENT]
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

        coreLinksShapefile = Cfg.PREFIX + '_sticks_s' + str(step) + '.shp'

        file = os.path.join(gp.workspace,coreLinksShapefile)
        if gp.Exists(file):
            try:
                gp.Delete(file)

            except:
                dashline(1)
                msg = ('ERROR: Could not remove shapefile ' +
                       coreLinksShapefile + '. Was it open in ArcMap?\n You may '
                       'need to re-start ArcMap to release the file lock.')
                gp.AddError(msg)
                exit(1)

        # make coreLinks.shp using linkCoords table
        # will contain linework between each pair of connected cores
        gp.CreateFeatureclass(gp.workspace, coreLinksShapefile,
                                  "POLYLINE")
        
        # ADD ATTRIBUTES
        gp.AddField_management(coreLinksShapefile, "Link_ID", "SHORT")
        gp.AddField_management(coreLinksShapefile, "Active", "SHORT")
        gp.AddField_management(coreLinksShapefile, "Link_Info", "TEXT")
        gp.AddField_management(coreLinksShapefile, "From_Core", "SHORT")
        gp.AddField_management(coreLinksShapefile, "To_Core", "SHORT")
        gp.AddField_management(coreLinksShapefile, "Euc_Dist", "FLOAT")
        gp.AddField_management(coreLinksShapefile, "CW_Dist", "FLOAT")
        gp.AddField_management(coreLinksShapefile, "cwd2Euc_R", "FLOAT")
        gp.AddField_management(coreLinksShapefile, "Eff_Resist", "FLOAT")
        gp.AddField_management(coreLinksShapefile, "cwd2EffR_r", "FLOAT")
        gp.AddField_management(coreLinksShapefile, "CF_Central", "FLOAT")
        #Create an Array and Point object.
        lineArray = gp.CreateObject("Array")
        pnt = gp.CreateObject("Point")

        # linkCoords indices:
        numLinks = len(linkCoords)

        #Open a cursor to insert rows into the shapefile.
        cur = gp.InsertCursor(coreLinksShapefile)

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
            feature = cur.NewRow()
            feature.shape = lineArray
            cur.InsertRow(feature)

            lineArray.RemoveAll()
        del cur

        #Add attribute data to link shapefile
        rows = gp.UpdateCursor(coreLinksShapefile)
        row = rows.Next()
        line = 0
        while row:
            # linkCoords indices
            row.SetValue("Link_ID", linkCoords[line, 0])
            if linkCoords[line, 9] == 2:
                row.SetValue("Link_Info", "Group_Pair")
            linktypecode = linkCoords[line, 9]
            activelink, linktypedesc = get_link_type_desc(linktypecode)
            row.SetValue("Active", activelink)
            row.SetValue("Link_Info", linktypedesc)

            row.SetValue("From_Core", linkCoords[line, 1])
            row.SetValue("To_Core", linkCoords[line, 2])
            row.SetValue("Euc_Dist", linkCoords[line, 3])
            row.SetValue("CW_Dist", linkCoords[line, 4])
            if linkCoords[line, 4] <= 0 or linkCoords[line, 3] <= 0:
                row.SetValue("cwd2Euc_R", -1)
            else:
                row.SetValue("cwd2Euc_R", linkCoords[line, 4] /
                             linkCoords[line, 3])
            row.SetValue("Eff_Resist", linkCoords[line, 10])
            row.SetValue("cwd2EffR_r", linkCoords[line, 11])
            row.SetValue("CF_Central", linkCoords[line, 12])
            
            rows.UpdateRow(row)
            row = rows.Next()
            line = line + 1

        del row, rows

        #clean up temp files
        gp.delete_management(coresForLinework)

        return

    except arcgisscripting.ExecuteError:
        raise_geoproc_error(__filename__)
    except:
        raise_python_error(__filename__)


############################################################################
## File and Path Management Functions ######################################
############################################################################


def move_old_results():
    try:
        oldFolder = Cfg.CWDBASEDIR_OLD
        newFolder = Cfg.CWDBASEDIR
        move_results_folder(oldFolder, newFolder)
            
        oldFolder = Cfg.ADJACENCYDIR_OLD
        newFolder = Cfg.ADJACENCYDIR
        move_results_folder(oldFolder, newFolder)            

        oldFolder = Cfg.LCCBASEDIR_OLD
        newFolder = Cfg.LCCBASEDIR
        move_results_folder(oldFolder, newFolder) 

        oldFolder = Cfg.BARRIERBASEDIR_OLD
        newFolder = Cfg.BARRIERBASEDIR
        move_results_folder(oldFolder, newFolder) 

        oldFolder = Cfg.CIRCUITBASEDIR_OLD
        newFolder = Cfg.CIRCUITBASEDIR
        move_results_folder(oldFolder, newFolder) 
        
        oldFolder = Cfg.CENTRALITYBASEDIR_OLD 
        newFolder = Cfg.CENTRALITYBASEDIR 
        move_results_folder(oldFolder, newFolder) 

    except:
        raise_python_error(__filename__)
        
        
def move_results_folder(oldFolder, newFolder):        
    try:
        if (os.path.exists(oldFolder) and not os.path.exists(newFolder)):
            os.rename(oldFolder, newFolder)
            #shutil.copytree(oldFolder,newFolder)
            #shutil.rmtree(oldFolder)
    except:
        raise_python_error(__filename__)
        

def delete_file(file):
    try:
        os.remove(file)
    except:
        pass
    return

    
def delete_dir(dir):
    try:
        gp.RefreshCatalog(dir)
        shutil.rmtree(dir)       
    except:
        # In case rmtree was unsuccessful due to lock on data
        try:
            if gp.Exists(dir):
                gp.RefreshCatalog(dir)
                gp.delete_management(dir)            
        except:
            pass
    return

    
def delete_data(dataset):
    try:
        gp.delete_management(dataset)
    except:
        pass


def get_cwd_path(core):
    """Returns the path for the cwd raster corresponding to a core area """
    dirCount = int(core / 100)
    if dirCount > 0:
        return os.path.join(Cfg.CWDBASEDIR, Cfg.CWDSUBDIR_NM + str(dirCount),
                         "cwd_" + str(core))
    else:
        return os.path.join(Cfg.CWDBASEDIR, Cfg.CWDSUBDIR_NM, "cwd_"
                         + str(core))

                         
def get_focal_path(core,radius):
    """Returns the path for the focal raster corresponding to a core area """
    dirCount = int(core / 100)
    focalDir1 = Cfg.FOCALSUBDIR1_NM + str(radius)
    if dirCount > 0:
        return os.path.join(Cfg.BARRIERBASEDIR, focalDir1,
                         Cfg.FOCALSUBDIR2_NM + str(dirCount), 
                         Cfg.FOCALGRID_NM + str(core))
    else:
        return os.path.join(Cfg.BARRIERBASEDIR, focalDir1,
                         Cfg.FOCALSUBDIR2_NM, Cfg.FOCALGRID_NM + str(core))
                         
                                                  
def check_project_dir():
    """Checks to make sure path name is not too long.

    Long path names can cause problems with ESRI grids.

    """
    if len(Cfg.PROJECTDIR) > 100:
        msg = ('ERROR: Project directory "' + Cfg.PROJECTDIR +
               '" is too deep.  Please choose a shallow directory'
               '(something like "C:\ANBO").')
        gp.AddError(msg)
        gp.AddMessage(gp.GetMessages(2))
        exit(1)
    return


def get_prev_step_link_table(step):
    """Returns the name of the link table created by the previous step"""
    try:
        prevStep = step - 1

        
        if (step == 7) or (step == 8):
            if step == 7:
                prevStepLinkTable = os.path.join(Cfg.DATAPASSDIR,
                                             'linkTable_s8.csv')
                gp.addmessage('\nLooking for ' + prevStepLinkTable)

                if os.path.exists(prevStepLinkTable):
                    return prevStepLinkTable
            else:
                prevStepLinkTable = os.path.join(Cfg.DATAPASSDIR,
                                             'linkTable_s7.csv')
                gp.addmessage('\nLooking for ' + prevStepLinkTable)

                if os.path.exists(prevStepLinkTable):
                    return prevStepLinkTable

            prevStepLinkTable = os.path.join(Cfg.DATAPASSDIR,
                                             'linkTable_s5.csv')
            gp.addmessage('\nLooking for ' + prevStepLinkTable)

            if os.path.exists(prevStepLinkTable):
                return prevStepLinkTable

            prevStepLinkTable = os.path.join(Cfg.DATAPASSDIR,
                                             'linkTable_s4.csv')
            gp.addmessage('\nLooking for ' + prevStepLinkTable)

            if os.path.exists(prevStepLinkTable):
                return prevStepLinkTable
            else:
                prevStep = 3  # Can skip steps 4 & 5
        
        
        if step == 6:
            prevStepLinkTable = os.path.join(Cfg.DATAPASSDIR,
                                             'linkTable_s5.csv')
            gp.addmessage('\nLooking for ' + prevStepLinkTable)

            if os.path.exists(prevStepLinkTable):
                return prevStepLinkTable

            prevStepLinkTable = os.path.join(Cfg.DATAPASSDIR,
                                             'linkTable_s4.csv')
            gp.addmessage('\nLooking for ' + prevStepLinkTable)

            if os.path.exists(prevStepLinkTable):
                return prevStepLinkTable
            else:
                prevStep = 3  # Can skip steps 4 & 5
                
        if step == 5:
            prevStepLinkTable = os.path.join(Cfg.DATAPASSDIR,
                                             'linkTable_s4.csv')
            gp.addmessage('\nLooking for ' + prevStepLinkTable)

            if os.path.exists(prevStepLinkTable):
                return prevStepLinkTable
            else:
                prevStep = 3  # Can skip step 4

        prevStepLinkTable = os.path.join(Cfg.DATAPASSDIR, 'linkTable_s' +
                                         str(prevStep) + '.csv')
        gp.addmessage('\nLooking for ' + Cfg.DATAPASSDIR +
                          '\linkTable_s' + str(prevStep) + '.csv')
        if os.path.exists(prevStepLinkTable):
            return prevStepLinkTable
        else:
            msg = ('\nERROR: Could not find a linktable from step previous to '
                   'step #' + str(step) + ' in datapass directory.  See above '
                   'for valid linktable files.')
            gp.AddError(msg)
            exit(1)

    except arcgisscripting.ExecuteError:
        raise_geoproc_error(__filename__)
    except:
        raise_python_error(__filename__)


def get_this_step_link_table(step):
    """Returns name of link table to write for current step"""
    try:
        filename = os.path.join(Cfg.DATAPASSDIR, 'linkTable_s' + str(step)
                             + '.csv')
        return filename

    except arcgisscripting.ExecuteError:
        raise_geoproc_error(__filename__)
    except:
        raise_python_error(__filename__)


def clean_up_link_tables(step):
    """Remove link tables from previous runs."""
    try:
        filename = os.path.join(Cfg.DATAPASSDIR, 'linkTable_s7_s8.csv')
        if os.path.isfile(filename):
            os.remove(filename)
        for stepNum in range(step, 9):
            filename = os.path.join(Cfg.DATAPASSDIR, 'linkTable_s' +
                                    str(stepNum) + '.csv')
            if os.path.isfile(filename):
                os.remove(filename)

        filename = os.path.join(Cfg.OUTPUTDIR, 'linkTable_final.csv')
        if os.path.isfile(filename):
            os.remove(filename)
        filename = os.path.join(Cfg.OUTPUTDIR, Cfg.PREFIX + '_linkTable_final.csv')
        if os.path.isfile(filename):
            os.remove(filename)

        filename = os.path.join(Cfg.OUTPUTDIR, 'linkTable_s5.csv')
        if os.path.isfile(filename):
            os.remove(filename)
        filename = os.path.join(Cfg.OUTPUTDIR, Cfg.PREFIX + 'linkTable_s5.csv')
        if os.path.isfile(filename):
            os.remove(filename)
            
    except arcgisscripting.ExecuteError:
        raise_geoproc_error(__filename__)
    except:
        raise_python_error(__filename__)


def copy_final_link_maps(step):
    """Copies final link maps from datapass to the output directory"""
    try:
        PREFIX = Cfg.PREFIX
        coreLinksShapefile = os.path.join(Cfg.OUTPUTDIR, PREFIX + '_sticks_s'
                                          + str(step) + '.shp')
        lcpShapefile = os.path.join(Cfg.DATAPASSDIR, 'lcpLines_s' +
                                    str(step) + '.shp')

        if not gp.exists(Cfg.LINKMAPGDB):
            gp.createfilegdb(Cfg.OUTPUTDIR, os.path.basename(Cfg.LINKMAPGDB))
        if gp.exists(coreLinksShapefile):
            gp.MakeFeatureLayer(coreLinksShapefile, "flinks")
            field = "Active"
            expression = field + " = " + str(1)
            gp.selectlayerbyattribute("flinks", "NEW_SELECTION",
                                          expression)

            activeLinksShapefile = os.path.join(Cfg.LINKMAPGDB,
                                               PREFIX + '_Sticks')
            gp.CopyFeatures_management("flinks", activeLinksShapefile)

            expression = field + " = " + str(0)
            gp.selectlayerbyattribute("flinks", "NEW_SELECTION",
                                          expression)
            inActiveLinksShapefile = os.path.join(Cfg.LINKMAPGDB,
                                                  PREFIX + '_Inactive_Sticks')
            gp.CopyFeatures_management("flinks", inActiveLinksShapefile)


        if gp.exists(lcpShapefile):
            gp.MakeFeatureLayer(lcpShapefile, "flcp")
            field = "Active"
            expression = field + " = " + str(1)
            gp.selectlayerbyattribute("flcp", "NEW_SELECTION", expression)

            activeLcpShapefile = os.path.join(Cfg.LINKMAPGDB, 
                                              PREFIX + '_LCPs')
            gp.CopyFeatures_management("flcp", activeLcpShapefile)

            expression = field + " = " + str(0)
            gp.selectlayerbyattribute("flcp", "NEW_SELECTION", expression)
            inActiveLcpShapefile = os.path.join(Cfg.LINKMAPGDB,
                                               PREFIX + '_Inactive_LCPs')
            gp.CopyFeatures_management("flcp", inActiveLcpShapefile)

        # Move stick and lcp maps for each step to log directory to reduce
        # clutter in output
        for i in range(2, 9):
            oldLinkFile = os.path.join(Cfg.OUTPUTDIR, PREFIX + '_sticks_s' 
                                        + str(i) + '.shp')
            logLinkFile = os.path.join(Cfg.LOGDIR, PREFIX + '_sticks_s' 
                                        + str(i) + '.shp')
            if gp.exists(oldLinkFile):
                try:
                    move_map(oldLinkFile, logLinkFile)
                except:
                    pass
            oldLcpShapeFile = os.path.join(Cfg.OUTPUTDIR, PREFIX + '_lcpLines_s'
                                           + str(i) + '.shp')
            logLcpShapeFile = os.path.join(Cfg.LOGDIR, PREFIX + '_lcpLines_s' +
                                           str(i) + '.shp')
            if gp.exists(oldLcpShapeFile):
                try:
                    move_map(oldLcpShapeFile, logLcpShapeFile)
                except:
                    pass
        return
    except arcgisscripting.ExecuteError:
        raise_geoproc_error(__filename__)
    except:
        raise_python_error(__filename__)


def move_map(oldMap, newMap):
    """Moves a map to a new location """
    if gp.exists(oldMap):
        if gp.exists(newMap):
            try:
                gp.delete_management(newMap)
            except:
                pass
        try:
            gp.CopyFeatures_management(oldMap, newMap)
            gp.delete_management(oldMap)
        except:
            pass
    return


############################################################################
##Error Checking and Handling Functions ####################################
############################################################################
def print_conefor_warning():
    """Warns that some links have no euclidean distances in conefor file."""
    gp.addmessage('\nWARNING: At least one potential link was dropped '
                      'because')
    gp.addmessage('there was no Euclidean distance value in the input '
                      'Euclidean')
    gp.addmessage('distance file from Conefor extension.\n')
    gp.addmessage('This may just mean that there were core areas that were'
                      ' adjacent')
    gp.addmessage('but were farther apart than the optional maximum '
                      'distance used ')
    gp.addmessage('when running Conefor.  But it can also mean that '
                      'distances  were')
    gp.addmessage('calculated using a different core area shapefile or the'
                      ' wrong field')
    gp.addmessage('in the same core area shapefile.\n')


def check_steps():
    """Check to make sure there are no skipped steps in a sequence of chosen
    steps (except step 4 which is optional)

    """
    skipStep = False
    if Cfg.STEP1 and not Cfg.STEP2 and Cfg.STEP3:
        skipStep = True
    if Cfg.STEP2 and not Cfg.STEP3 and (Cfg.STEP4 or Cfg.STEP5):
        skipStep = True
    if skipStep:
        try:
            dashline(1)
            msg = ("Error: You can start or stop at different steps, but you "
                   "can't SKIP any except for step 4.\n")
            gp.AddError(msg)
            exit(0)
        except:
            raise_python_error(__filename__)
    return


def check_cores():
    """Checks for positive integer core IDs with appropriate naming."""
    try:
        invalidFNs = ['fid','id','oid','shape']
        if Cfg.COREFN.lower() in invalidFNs:
            dashline(1)
            msg = ('ERROR: Core area field name "ID", "FID", "OID", and "Shape" are reserved '
                    'for ArcGIS. Please choose another field- must be a '
                    'positive integer.')
            gp.AddError(msg)
            exit(1)      

        fieldList = gp.ListFields(Cfg.COREFC)
        for field in fieldList:
            if str(field.Name) == Cfg.COREFN:
                FT = str(field.Type)
                if (FT != 'SmallInteger' and FT != 'SHORT' and FT != 'Integer' 
                    and FT != 'LONG'):
                    dashline(1)
                    msg = ('ERROR: Core area field must be in Short Integer '
                            'format.')
                    gp.AddError(msg)
                    exit(1)                   

        coreList = get_core_list()
        # test = coreList - coreList.astype(int)
        # if npy.any(test) == True:
            # dashline(1)
            # msg = ('ERROR: Core area field must be in integer format. ')
            # gp.AddError(msg)
            # exit(1)

        if npy.amin(coreList) < 1:
            dashline(1)
            msg = ('ERROR: Core area field must contain only positive integers. ')
            gp.AddError(msg)
            exit(1)

    except arcgisscripting.ExecuteError:
        raise_geoproc_error(__filename__)
    except:
        raise_python_error(__filename__)
    
    
    
    
def hiccup_test(count, statement):
    """Re-tries ArcGIS calls in case of server problems or 'other hiccups'."""
    try:
        if count < 10:
            sleepTime = 10 * count
            count = count + 1
            dashline(1)
            if count == 1:
                gp.addmessage('Failed to execute ' + statement + ' on try '
                                  '#' + str(count) + '.\n Could be an ArcGIS '
                                  'hiccup.')
                # dashline(2)
                gp.addmessage("Here's the error being reported: ")

                for msg in range(0, gp.MessageCount):
                    if gp.GetSeverity(msg) == 2:
                        gp.AddReturnMessage(msg)
                    print gp.AddReturnMessage(msg)
                    # dashline(2)
            else:
                gp.addmessage('Failed again executing ' + statement +
                                  ' on try #' + str(count) +
                                  '.\n Could be an ArcGIS hiccup- scroll up '
                                  'for error description.')

                gp.addmessage('---------Trying again in ' +
                                  str(sleepTime) + ' seconds---------\n')
            time.sleep(sleepTime)
            return count, True
        else:
            sleepTime = 60
            count = count + 1
            dashline(1)
            gp.addmessage('Failed to execute ' + statement + ' on try #' +
                          str(count) + '.\n Could be an ArcGIS hiccup.  Trying'
                          'again in 1 minute.\n')
            time.sleep(sleepTime)
            return count, True
    except:
        raise_python_error(__filename__)


def raise_geoproc_error(filename):
    """Handle geoprocessor errors and provide details to user"""
    dashline(1)
    tb = sys.exc_info()[2]  # get the traceback object
    # tbinfo contains the error's line number and the code
    tbinfo = traceback.format_tb(tb)[0]
    line = tbinfo.split(", ")[1]

    gp.AddError("Geoprocessing error on **" + line + "** of " + filename + " "
                "in Linkage Mapper Version " + str(__version__) + ":")

    dashline(1)
    for msg in range(0, gp.MessageCount):
        if gp.GetSeverity(msg) == 2:
            gp.AddReturnMessage(msg)
        # dashline(2)
        print gp.AddReturnMessage(msg)
        # dashline(2)
    exit(0)


def raise_python_error(filename):
    """Handle python errors and provide details to user"""
    dashline(1)
    tb = sys.exc_info()[2]  # get the traceback object
    # tbinfo contains the error's line number and the code
    tbinfo = traceback.format_tb(tb)[0]
    line = tbinfo.split(", ")[1]

    err = traceback.format_exc().splitlines()[-1]

    gp.AddError("Python error on **" + line + "** of " + filename + " "
                "in Linkage Mapper Version " + str(__version__) + ":")
    gp.AddError(err)

    # dashline(2)
    exit(0)


def dashline(lspace=0):
    """Output dashed line in tool output dialog.

       0 = No empty line
       1 = Empty line before
       2 = Empty line after

    """
    if lspace == 1:
        gp.addmessage('\n')
    gp.addmessage('---------------------------------')
    if lspace == 2:
        gp.addmessage('\n')

        
############################################################################
## Config File #############################################################
############################################################################
def set_lm_options():
    """Set default options

    """

    options = {}
    options['lcc_cutoff'] = None
    return options
        
        
############################################################################
## Circuitscape Functions ##################################################
############################################################################
def setCircuitscapeOptions():
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
    options['print_timings']=True
    options['use_mask']=False
    options['mask_file']='None' 
    options['use_included_pairs']=False
    options['included_pairs_file']='None' 
    options['use_variable_source_strengths']=False
    options['variable_source_file']='None' 

    return options

def writeCircuitscapeConfigFile(configFile, options):
   
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
        except:
            pass
    for option in sections:
        config.set(sections[option], option, options[option])

    f = open(configFile, 'w')
    config.write(f)
    f.close()
        
