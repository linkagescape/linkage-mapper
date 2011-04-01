#!/usr/bin/env python2.5

"""Linkage Mapper configuration module

Assigns input parameter from ToolBox to variables and reads defaults from
linkage_mapper.ini

"""

import os.path as path
import sys

import arcgisscripting

def str2bool(pstr):
        return pstr == 'true'

def setadjmeth(inparam):
    if inparam == "Cost-Weighted":
        meth_cw = True
        meth_eu = False
    elif inparam == "Euclidean":
        meth_cw = False
        meth_eu = True
    else:
        meth_cw = True
        meth_eu = True
    return meth_cw, meth_eu

def nullfloat(innum):
    if innum == '#':
        nfloat = None
    else:
        nfloat=float(innum)
    return nfloat

class Config():
    PROJECTDIR = sys.argv[1]            # Project directory
    COREFC = sys.argv[2]                # Core area feature class
    COREFN = sys.argv[3]                # Core area field name
    RESRAST = sys.argv[4]               # Resistance raster

    STEP1 = str2bool(sys.argv[5])
    S1ADJMETH_CW, S1ADJMETH_EU = setadjmeth(sys.argv[6])

    STEP2 = str2bool(sys.argv[7])
    S2EUCDISTFILE = sys.argv[8]
    S2ADJMETH_CW, S2ADJMETH_EU = setadjmeth(sys.argv[9])

    STEP3 = str2bool(sys.argv[10])
    S3DROPLCCS = sys.argv[11]             # Drop LCC's with intermediate cores

    STEP4 = str2bool(sys.argv[12])
    S4MAXNN = int(sys.argv[13])           # No of connected nearest neighbors
    # Nearest neighbor measurement unit
    S4DISTTYPE_CW, S4DISTTYPE_EU  = setadjmeth(sys.argv[14])
    S4CONNECT = str2bool(sys.argv[15])

    STEP5 = str2bool(sys.argv[16])

    BUFFERDIST = nullfloat(sys.argv[17])
    MAXCOSTDIST = nullfloat(sys.argv[18])
    MAXEUCDIST = nullfloat(sys.argv[19])

    MINCOSTDIST = None
    MINEUCDIST = None

    OUTPUTDIR = path.join(PROJECTDIR, "output")
    SCRATCHDIR = path.join(PROJECTDIR, "scratch")
    LOGDIR = path.join(PROJECTDIR, "log")
    DATAPASSDIR = path.join(PROJECTDIR, "datapass")
    DATAPASSARCHDIR = path.join(PROJECTDIR, "datapass_archive")
    ADJACENCYDIR = path.join(PROJECTDIR, "adj")
    CWDBASEDIR = path.join(PROJECTDIR, "cwd")
    CWDSUBDIR_NM = "cw"
    LCCBASEDIR = path.join(PROJECTDIR, "nlcc")
    LCCNLCDIR_NM = "nlc"
    LCCMOSAICDIR = path.join(LCCBASEDIR, "mosaic")

    OUTPUTGDB = path.join(OUTPUTDIR, "linkages.gdb")

    # This is how "wide" corridors will be (measured in cost-weighted
    # distances) in a truncated raster
    CWDTHRESH = 100000
    if MAXCOSTDIST is None:
        TMAXCWDIST = None
    else:
        TMAXCWDIST = MAXCOSTDIST + CWDTHRESH  # This will limit cw calcs

    SAVENORMLCCS = False  # Set to True to save individual normalized LCC grids
    FCORES = "fcores"
    BNDCIRCEN = "boundingCircleCenter.shp"
    BNDCIR = "boundingCircle.shp"

    # Link table column numbers
    LTB_LINKID = 0  # Link ID
    LTB_CORE1 = 1  # Core ID of 1st core area link connects
    LTB_CORE2 = 2  # Core ID of 2nd core area link connects
    LTB_CLUST1 = 3  # Component ID of 1st core area link connects
    LTB_CLUST2 = 4  # Component ID of 2nd core area link connects
    LTB_LINKTYPE = 5
    LTB_EUCDIST = 6
    LTB_CWDIST = 7
    LTB_EUCADJ = 8
    LTB_CWDADJ = 9
    LTB_LCPLEN = 10
    LTB_CWDEUCR = 11
    LTB_CWDPATHR = 12

    # Linkage type values (NOTE- these map to codes in get_linktype_desc)
    LT_CPLK = -1  # Not_nearest N neighbors
    LT_TLEC = -11  # Too long Euclidean distance
    LT_TLLC = -12  # Too long Cost-Weighted distance
    LT_TSEC = -13  # Too short Euclidean distance
    LT_TSLC = -14  # Too short Cost-Weighted distance
    LT_INT = -15  # Intermediate corea area detected
    LT_CORR = 1  # Connects cores (corridor)
    LT_NNC = 10 # Corridor connecting N Nearest Neighbors
    LT_CLU = 20  # Connects_constellations (corridor)
    LT_NNCT = 30  # TEMP NN corridor links (s4), may be able to get rid of this
# 1 corridor
# 10 NN corridor
# 11 1st nn (future)
# 12 2nd nn etc (future)
# 20 constel
# 21 1st nn constel (future)
# 22 2nd nn constel etc (future)
# 30 temp saved nnconstel.  maybe able to get rid fo this

   
    gp = arcgisscripting.create(9.3)
    gp.CheckOutExtension("Spatial")
    gp.OverwriteOutput = True
    gp.SnapRaster = RESRAST