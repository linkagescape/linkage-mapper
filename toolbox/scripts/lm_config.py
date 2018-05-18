#!/usr/bin/env python2
# Authors: Brad McRae and Darren Kavanagh

"""Linkage Mapper configuration module.

Assigns input parameters from ToolBox to variables, and sets constants.

"""

import imp
import os.path as path

import arcgisscripting

import lm_version as ver
import lm_settings


GP_NULL = '#'


def str2bool(pstr):
    """Convert ESRI boolean string to Python boolean type."""
    return pstr == 'true'


def setadjmeth(inparam):
    """Return boolean variables for adjacency methods."""
    meth_cw = "cost" in inparam.lower()
    meth_eu = "euclid" in inparam.lower()
    return meth_cw, meth_eu


def nullfloat(innum):
    """Convert ESRI float or null to Python float."""
    if innum == GP_NULL:
        nfloat = None
    else:
        nfloat = float(innum)
        if nfloat == 0:
            nfloat = None
    return nfloat


def nullstring(arg_string):
    """Convert ESRI nullstring to Python null."""
    if arg_string == GP_NULL:
        arg_string = None
    return arg_string


def config_global(config, arg):
    """Configure global variables for all tools."""
    config.PARAMS = str(arg)  # Convert to string in case '\' exists
    config.releaseNum = ver.releaseNum
    config.LOGMESSAGES = True
    # File names, directory paths & folder names
    proj_dir = arg[1]
    config.PROJECTDIR = proj_dir  # Project directory
    config.SCRATCHDIR = path.join(proj_dir, "scratch")
    config.ARCSCRATCHDIR = path.join(config.SCRATCHDIR, "arcscratch")
    config.PREFIX = path.basename(proj_dir)
    config.DATAPASSDIR = path.join(proj_dir, "datapass")
    config.CWDADJFILE = path.join(config.DATAPASSDIR, "cwdAdj.csv")
    config.EUCADJFILE = path.join(config.DATAPASSDIR, "eucAdj.csv")
    config.OUTPUTDIR = path.join(proj_dir, "output")
    config.LOGDIR = path.join(proj_dir, "run_history")
    config.LOGDIR_OLD = path.join(proj_dir, "logFiles")
    config.logFile = None
    config.logFilePath = None
    config.logFileCopyPath = path.join(proj_dir, 'last_run_log.txt')
    config.MESSAGEDIR = path.join(config.LOGDIR, "log")
    config.MESSAGEDIR_OLD = path.join(config.LOGDIR, "Messages")
    config.ADJACENCYDIR = path.join(config.DATAPASSDIR, "adj")
    config.ADJACENCYDIR_OLD = path.join(proj_dir, "adj")
    config.CWDBASEDIR = path.join(config.DATAPASSDIR, "cwd")
    config.CWDBASEDIR_OLD = path.join(proj_dir, "cwd")
    config.CWDSUBDIR_NM = "cw"
    config.LCCBASEDIR = path.join(config.DATAPASSDIR, "nlcc")
    config.LCCBASEDIR_OLD = path.join(proj_dir, "nlcc")
    config.LCCNLCDIR_NM = "nlc"
    config.LCCMOSAICDIR = path.join(config.LCCBASEDIR, "mosaic")
    config.MOSAICGDB = path.join(config.LCCMOSAICDIR, "mosaic.gdb")
    config.FOCALSUBDIR1_NM = "focalr"
    config.FOCALSUBDIR2_NM = "f"
    config.FOCALGRID_NM = "focal"
    config.BARRIERBASEDIR = path.join(proj_dir, "barrier_tmp")
    config.BARRIERDIR_NM = "bar"
    config.BARRIERMOSAICDIR = "mosaic"
    config.CIRCUITBASEDIR = path.join(proj_dir, "pinchpt_tmp")
    config.CENTRALITYBASEDIR = path.join(proj_dir, "centrality_tmp")

    config.CIRCUITCONFIGDIR_NM = "config"
    config.CIRCUITOUTPUTDIR_NM = "output"

    # Save individual current maps from Circuitscape
    config.SAVECURRENTMAPS = False

    config.SAVECIRCUITDIR = False
    config.SAVE_TEMP_CIRCUIT_FILES = False

    # Write voltage maps from pinchpoint analysis
    config.WRITE_VOLT_MAPS = False

    config.FCORES = "fcores"

    config.OUTPUTGDB = path.join(config.OUTPUTDIR, "corridors.gdb")
    config.EXTRAGDB = path.join(config.OUTPUTDIR, "extra.gdb")
    config.OUTPUTGDB_OLD = path.join(config.OUTPUTDIR, "linkages.gdb")
    config.CWDGDB = path.join(config.OUTPUTDIR, "cwd.gdb")
    config.LINKMAPGDB = path.join(config.OUTPUTDIR, "link_maps.gdb")
    config.LOGLINKMAPGDB = path.join(config.LOGDIR, "link_maps.gdb")
    config.BARRIERGDB = path.join(config.OUTPUTDIR, "barriers.gdb")
    config.PINCHGDB = path.join(config.OUTPUTDIR, "pinchpoints.gdb")
    config.CORECENTRALITYGDB = path.join(config.OUTPUTDIR,
                                         "core_centrality.gdb")
    config.BNDCIRCEN = path.join(config.SCRATCHDIR,
                                 "boundingCircleCenter.shp")
    config.BNDCIRCENS = path.join(config.SCRATCHDIR,
                                  "boundingCircleCenters.shp")
    config.BNDCIR = path.join(config.SCRATCHDIR, "boundingCircle.shp")
    config.BNDCIRS = path.join(config.SCRATCHDIR, "boundingCircles.shp")
    config.BNDFC = "boundingFeature.shp"
    config.BOUNDRESIS = path.join(config.SCRATCHDIR, "boundResis")

    # Link table column numbers
    config.LTB_LINKID = 0  # Link ID
    config.LTB_CORE1 = 1  # Core ID of 1st core area link connects
    config.LTB_CORE2 = 2  # Core ID of 2nd core area link connects
    config.LTB_CLUST1 = 3  # Component ID of 1st core area link connects
    config.LTB_CLUST2 = 4  # Component ID of 2nd core area link connects
    config.LTB_LINKTYPE = 5
    config.LTB_EUCDIST = 6
    config.LTB_CWDIST = 7
    config.LTB_EUCADJ = 8
    config.LTB_CWDADJ = 9
    config.LTB_LCPLEN = 10
    config.LTB_CWDEUCR = 11
    config.LTB_CWDPATHR = 12
    config.LTB_EFFRESIST = 13
    config.LTB_CWDTORR = 14
    config.LTB_CURRENT = 15

    # Linkage type values (NOTE- these map to codes in get_linktype_desc)
    config.LT_CPLK = -1  # Not_nearest N neighbors
    config.LT_TLEC = -11  # Too long Euclidean distance
    config.LT_TLLC = -12  # Too long Cost-Weighted distance
    config.LT_TSEC = -13  # Too short Euclidean distance
    config.LT_TSLC = -14  # Too short Cost-Weighted distance
    config.LT_INT = -15  # Intermediate corea area detected
    config.LT_CORR = 1  # Connects cores (corridor)
    config.LT_NNC = 10  # Corridor connecting N Nearest Neighbors
    config.LT_CLU = 20  # Connects_constellations (corridor)
    config.LT_KEEP = 100  # user retained despite length
    # TEMP NN corridor links (s4), may be able to get rid of this
    config.LT_NNCT = 30

    # Temporary resistance raster copy to be created in master scripts
    config.RESRAST = path.join(config.SCRATCHDIR, 'resrast')


def config_lm(config, arg):
    """Configure global variables for Linkage Mapper."""
    config.CONNECTFRAGS = False
    config.COREFC = arg[2]  # Core area feature class
    splits = config.COREFC.split("\\")
    config.CORENAME = splits[len(splits) - 1].split(".")[0]
    config.COREFN = arg[3]  # Core area field name
    config.RESRAST_IN = arg[4]  # Resistance raster

    # Processing steps inputs
    config.STEP1 = str2bool(arg[5])

    config.S1ADJMETH_CW = True
    config.S1ADJMETH_EU = True

    config.STEP2 = str2bool(arg[6])
    config.S2ADJMETH_CW, config.S2ADJMETH_EU = setadjmeth(arg[7])
    config.S2EUCDISTFILE = nullstring(arg[8])

    config.STEP3 = str2bool(arg[9])
    # Drop LCC's passing through intermediate cores
    config.S3DROPLCCS = str2bool(arg[10])
    config.STEP4 = str2bool(arg[11])
    if arg[12] == "Unlimited":
        config.S4MAXNN = 99
        config.IGNORES4MAXNN = True
    else:
        config.S4MAXNN = int(arg[12])  # No of connected nearest neighbors
        config.IGNORES4MAXNN = False
    # NN Unit
    config.S4DISTTYPE_CW, config.S4DISTTYPE_EU = setadjmeth(arg[13])
    config.S4CONNECT = str2bool(arg[14])
    config.STEP5 = str2bool(arg[15])

    # Optional input parameters
    config.BUFFERDIST = nullfloat(arg[18])
    config.MAXCOSTDIST = nullfloat(arg[19])

    if config.S2EUCDISTFILE is not None:
        if config.S2EUCDISTFILE.lower() == 'cluster':
            # Custom code to consolidate nearby cores will be called
            config.S2EUCDISTFILE = None
            config.CONNECTFRAGS = True
            config.MAXCOSTDIST = None
            config.S1ADJMETH_CW = False
            config.S2ADJMETH_CW = False

    config.MAXEUCDIST = nullfloat(arg[20])

    # Optional parameters that only apply to 2.0.0+ toolbox,
    # for ArcGIS 10.x only
    if "10." in config.gp.GetInstallInfo('desktop')['Version']:
        config.WRITETRUNCRASTER = str2bool(arg[16])
        config.CWDTHRESH = int(arg[17])
        config.OUTPUTFORMODELBUILDER = nullstring(arg[21])
        config.LMCUSTSETTINGS = nullstring(arg[22])
    else:
        config.OUTPUTFORMODELBUILDER = None

        # These two settings are hardcoded based on the values
        # used in lm_settings, before they were daylighted in
        # the toolbox in LM 2.0.0+ for ArcGIS10.x
        config.WRITETRUNCRASTER = True
        config.CWDTHRESH = 200000

        config.LMCUSTSETTINGS = None

    if config.LMCUSTSETTINGS:
        cust_settings = (imp.load_source(
            config.LMCUSTSETTINGS.split(".")[0], config.LMCUSTSETTINGS))
        for setting in dir(cust_settings):
            if setting == setting.upper():
                setting_value = getattr(cust_settings, setting)
                setattr(config, setting, setting_value)
    else:
        for setting in dir(lm_settings):
            if setting == setting.upper():
                setting_value = getattr(lm_settings, setting)
                setattr(config, setting, setting_value)

    if config.MAXCOSTDIST is None:
        config.TMAXCWDIST = None
    elif config.CWDTHRESH is not None:
        config.TMAXCWDIST = None
    config.SCRATCHGDB = path.join(config.SCRATCHDIR, "scratch.gdb")
    config.CORERAS = path.join(config.SCRATCHDIR, "core_ras")
    return True


def config_barrier(config, arg):
    """Configure global variables for Barrier tool."""
    config.RESRAST_IN = arg[2]
    config.STARTRADIUS = arg[3]
    config.ENDRADIUS = arg[4]
    config.RADIUSSTEP = arg[5]
    config.BARRIER_METH = arg[6]
    config.SAVE_RADIUS_RASTERS = str2bool(arg[7])
    config.WRITE_PCT_RASTERS = str2bool(arg[8])
    if len(arg) > 9:
        config.BARRIER_CWD_THRESH = arg[9]
        if arg[9] == "#" or arg[9] == "" or arg[9] == 0 or arg[9] == "0":
            config.BARRIER_CWD_THRESH = None
    else:
        config.BARRIER_CWD_THRESH = None

    config.BARRIER_METH_MAX = 'max' in config.BARRIER_METH.lower()
    config.BARRIER_METH_SUM = 'sum' in config.BARRIER_METH.lower()

    if config.RADIUSSTEP == GP_NULL or config.ENDRADIUS == config.STARTRADIUS:
        config.RADIUSSTEP = 0
    if ((float(config.STARTRADIUS) + float(config.RADIUSSTEP))
            > float(config.ENDRADIUS)):
        config.RADIUSSTEP = 0
    if config.RADIUSSTEP == 0:
        config.SAVE_RADIUS_RASTERS = True

    # Save individual focal grids for barrier analysis
    config.SAVEFOCALRASTERS = False

    # Save temporary directories
    config.SAVEBARRIERDIR = False

    # Calculate minimum of resistance and improvement score
    config.WRITE_TRIM_RASTERS = False

    # Save individual barrier grids for each core area pair
    config.SAVEBARRIERRASTERS = False

    config.STEP1 = False


def config_climate(config, arg):
    """Configure global variables for Climate Corridor tool."""
    config.lm_configured = config_lm(config, arg)


def config_lp(config, arg):
    """Configure global variables for Linkage Priority tool."""
    config.lm_configured = config_lm(config, arg)


def config_circuitscape(config, arg):
    """Configure global variables for Circuitscape."""
    config.CSPATH = arg[-1]  # Path to Circuitscape
    config.COREFC = arg[2]
    config.COREFN = arg[3]

    if len(arg) == 5:
        config.DOCENTRALITY = True
        config.DOPINCH = False
        config.CWDCUTOFF = 0
        config.DO_ALLPAIRS = False
    else:
        config.DOPINCH = True
        config.DOCENTRALITY = False
        config.RESRAST_IN = arg[4]
        config.CWDCUTOFF = int(nullfloat(arg[5]))  # CDW cutoff distance
        config.SQUARERESISTANCES = str2bool(arg[6])  # Square resistance values

        # Do adjacent pair corridor pinchpoint calculations
        # using raster CWD maps
        config.DO_ADJACENTPAIRS = str2bool(arg[7])
        # Do all-pair current calculations using raster corridor map
        config.DO_ALLPAIRS = str2bool(arg[8])
        config.ALL_PAIR_CHOICE = arg[9]
        if config.DO_ALLPAIRS:
            if "pairwise" in config.ALL_PAIR_CHOICE.lower():
                config.ALL_PAIR_SCENARIO = 'pairwise'
            else:
                config.ALL_PAIR_SCENARIO = 'all-to-one'

    config.STEP1 = False

    config.SAVECENTRALITYDIR = False


class Configure(object):
    """Class container to hold global variables."""

    TOOL_LM = 'Linkage Mapper'
    TOOL_CC = 'Linkage Mapper Climate'
    TOOL_LP = 'Linkage Priority'
    TOOL_BM = 'Barrier mapper'
    TOOL_CS = 'Circuitscape'

    def __init__(self):
        """Initialize class and create single geoprocessor object."""
        self.gp = arcgisscripting.create(9.3)
        self.gp.CheckOutExtension("Spatial")
        self.gp.OverwriteOutput = True
        self.lm_configured = False

    def configure(self, tool, arg):
        """Assign variables for Configure class."""
        config_global(self, arg)

        if tool == Configure.TOOL_LM:
            self.lm_configured = config_lm(self, arg)
        elif tool == Configure.TOOL_CC:
            config_climate(self, arg)
        elif tool == Configure.TOOL_LP:
            config_lp(self, arg)
        elif tool == Configure.TOOL_BM:
            config_barrier(self, arg)
        elif tool == Configure.TOOL_CS:
            config_circuitscape(self, arg)
        else:
            raise RuntimeError('Undefined tool to configure')
        self.TOOL = tool


tool_env = Configure()  # Class instance that is use by tool modules
