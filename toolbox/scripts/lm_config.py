# Authors: Brad McRae and Darren Kavanagh

"""Linkage Mapper configuration module.

Assigns input parameters from ToolBox to variables, and sets constants.

"""

from os import path
import imp
import json

import arcpy

import lm_version as ver
import lm_util_config as util


def get_code_path():
    """Get full path to tool scripts folder."""
    return path.dirname(path.realpath(__file__))


def set_custom(settings, config_obj):
    """Add attributes for custom file to configuration object."""
    cust_settings = imp.load_source("set_mod", settings)
    for setting in dir(cust_settings):
        if setting == setting.upper():
            setting_value = getattr(cust_settings, setting)
            setattr(config_obj, setting, setting_value)


def setadjmeth(inparam):
    """Return boolean variables for adjacency methods."""
    meth_cw = "cost" in inparam.lower()
    meth_eu = "euclid" in inparam.lower()
    return meth_cw, meth_eu


def nullfloat(innum):
    """Convert ESRI float or null to Python float."""
    if innum == util.GP_NULL:
        nfloat = None
    else:
        nfloat = float(innum)
        if nfloat == 0:
            nfloat = None
    return nfloat


def config_global(config, arg):
    """Configure global variables for all tools."""
    config.releaseNum = ver.releaseNum
    config.LOGMESSAGES = True
    # File names, directory paths & folder names
    proj_dir = arg[1]
    config.PROJECTDIR = proj_dir  # Project directory
    config.SCRATCHDIR = path.join(proj_dir, "scratch")
    config.ARCSCRATCHDIR = path.join(config.SCRATCHDIR, "arcscratch")
    config.PREFIX = path.basename(proj_dir)
    config.DATAPASSDIR = path.join(proj_dir, "datapass")
    config.LM_PASSFILE = path.join(config.DATAPASSDIR, "lm_param.json")
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
    config.CORENAME = path.splitext(path.basename(config.COREFC))[0]
    config.COREFN = arg[3]  # Core area field name
    config.RESRAST_IN = arg[4]  # Resistance raster

    # Processing steps inputs
    config.STEP1 = util.str2bool(arg[5])

    config.S1ADJMETH_CW = True
    config.S1ADJMETH_EU = True

    config.STEP2 = util.str2bool(arg[6])
    config.S2ADJMETH_CW, config.S2ADJMETH_EU = setadjmeth(arg[7])
    config.S2EUCDISTFILE = util.nullstring(arg[8])

    config.STEP3 = util.str2bool(arg[9])
    # Drop LCC's passing through intermediate cores
    config.S3DROPLCCS = util.str2bool(arg[10])
    config.STEP4 = util.str2bool(arg[11])
    if arg[12] == "Unlimited":
        config.S4MAXNN = 99
        config.IGNORES4MAXNN = True
    else:
        config.S4MAXNN = int(arg[12])  # No of connected nearest neighbors
        config.IGNORES4MAXNN = False
    # NN Unit
    config.S4DISTTYPE_CW, config.S4DISTTYPE_EU = setadjmeth(arg[13])
    config.S4CONNECT = util.str2bool(arg[14])

    config.STEP5 = util.str2bool(arg[15])
    config.WRITETRUNCRASTER = util.str2bool(arg[16])
    config.CWDTHRESH = int(arg[17])

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

    config.OUTPUTFORMODELBUILDER = util.nullstring(arg[21])

    if arg[22] == util.GP_NULL:
        config.LMCUSTSETTINGS = path.join(get_code_path(),
                                          'lm_settings.py')
    else:
        config.LMCUSTSETTINGS = arg[22]

    set_custom(config.LMCUSTSETTINGS, config)

    if config.MAXCOSTDIST is None:
        config.TMAXCWDIST = None
    elif config.CWDTHRESH is not None:
        config.TMAXCWDIST = None

    config.CORERAS = path.join(config.SCRATCHDIR, "core_ras")


def config_barrier(config, arg):
    """Configure global variables for Barrier tool."""
    config.RESRAST_IN = arg[2]
    config.STARTRADIUS = arg[3]
    config.ENDRADIUS = arg[4]
    config.RADIUSSTEP = arg[5]
    config.BARRIER_METH = arg[6]
    config.SAVE_RADIUS_RASTERS = util.str2bool(arg[7])
    config.WRITE_PCT_RASTERS = util.str2bool(arg[8])
    if (len(arg) == 9 or
            (len(arg) > 9 and arg[9] in ["#", "", 0, "0"])):
        config.BARRIER_CWD_THRESH = None
    else:
        config.BARRIER_CWD_THRESH = arg[9]

    config.BARRIER_METH_MAX = 'max' in config.BARRIER_METH.lower()
    config.BARRIER_METH_SUM = 'sum' in config.BARRIER_METH.lower()

    if (config.RADIUSSTEP == util.GP_NULL
            or config.ENDRADIUS == config.STARTRADIUS):
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
    config_lm(config, arg)


def get_cwdthresh(lm_passfile):
    """Get CWDTHRESH from Linkage Pathways model run."""
    try:
        with open(lm_passfile, 'r') as params_file:
            settings = json.load(params_file)
    except FileNotFoundError:
        raise RuntimeError('Linkage Pathways parameters file not found. '
                           'Try re-running Linkage Pathways.')
    return settings['CWDTHRESH']


def config_lp(config, arg):
    """Configure global variables for Linkage Priority tool."""
    # Model Inputs
    # ------------
    config.COREFC = arg[2]
    config.COREFN = arg[3]
    config.RESRAST_IN = arg[4]

    # Core Area Value (CAV) Options
    # -----------------------------
    config.OCAVRAST_IN = util.nullstring(arg[5])
    config.RESWEIGHT = float(arg[6])
    config.SIZEWEIGHT = float(arg[7])
    config.APWEIGHT = float(arg[8])
    config.ECAVWEIGHT = float(arg[9])
    config.CFCWEIGHT = float(arg[10])
    config.OCAVWEIGHT = float(arg[11])

    # Corridor Specific Priority (CSP) Options
    # ----------------------------------------
    #  Expert Corridor Importance Value
    config.COREPAIRSTABLE_IN = util.nullstring(arg[12])
    config.FROMCOREFIELD = util.nullstring(arg[13])
    config.TOCOREFIELD = util.nullstring(arg[14])
    config.ECIVFIELD = util.nullstring(arg[15])

    # Climate Linkage Priority Value
    config.CCERAST_IN = util.nullstring(arg[16])
    config.FCERAST_IN = util.nullstring(arg[18])
    config.CANALOG_MIN = float(arg[19])
    config.CANALOG_MAX = float(arg[20])
    config.CANALOG_MINRMAX = float(arg[21])
    config.CANALOG_TARGET = float(arg[22])
    config.CANALOG_PIORITY = float(arg[23])
    config.CANALOG_WEIGHT = float(arg[24])
    config.CPREF_VALUE = float(arg[25])
    config.CPREF_MIN = float(arg[26])
    config.CPREF_MAX = float(arg[27])
    config.CPREF_WEIGHT = float(arg[28])

    # CSP Weights
    config.CLOSEWEIGHT = float(arg[29])
    config.PERMWEIGHT = float(arg[30])
    config.CAVWEIGHT = float(arg[31])
    config.ECIVWEIGHT = float(arg[32])
    config.CEDWEIGHT = float(arg[33])

    # CPS Trim Value
    config.CPSNORM_CUTOFF = nullfloat(arg[34])

    # Blended Priority Options
    # ------------------------
    config.TRUNCWEIGHT = float(arg[35])
    config.LPWEIGHT = float(arg[36])

    # Additional Options
    # ------------------
    config.OUTPUTFORMODELBUILDER = util.nullstring(arg[37])
    if arg[38] == util.GP_NULL:
        config.LPCUSTSETTINGS_IN = path.join(get_code_path(),
                                             'lp_settings.py')
    else:
        config.LPCUSTSETTINGS_IN = arg[38]

    # Settings from Linkage Pathways
    # ------------------------------
    config.CWDTHRESH = get_cwdthresh(config.LM_PASSFILE)

    #  Custom settings
    # ----------------
    set_custom(config.LPCUSTSETTINGS_IN, config)

    config.CALC_CSP = 1  # Calculate Corridor Specific Priority (CSP)
    config.CALC_CSPBP = 2  # Calculate CSP and Blended Priority

    # Misc
    # ----
    config.SCRATCHGDB = path.join(config.SCRATCHDIR, "scratch.gdb")
    config.CORENAME = path.splitext(path.basename(config.COREFC))[0]


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

        # CDW cutoff distance
        config.CWDCUTOFF = int(nullfloat(arg[5]))
        # Square resistance values
        config.SQUARERESISTANCES = util.str2bool(arg[6])

        # Do adjacent pair corridor pinchpoint calculations
        # using raster CWD maps
        config.DO_ADJACENTPAIRS = util.str2bool(arg[7])
        # Do all-pair current calculations using raster corridor map
        config.DO_ALLPAIRS = util.str2bool(arg[8])
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
    TOOL_CC = 'Climate Linkage Mapper'
    TOOL_LP = 'Linkage Priority'
    TOOL_BM = 'Barrier mapper'
    TOOL_CS = 'Circuitscape'

    def __init__(self):
        """Initialize class."""
        arcpy.CheckOutExtension("Spatial")
        arcpy.env.overwriteOutput = True
        self.TOOL = ''

    def configure(self, tool, arg):
        """Assign variables for Configure class."""
        config_global(self, arg)

        if tool == Configure.TOOL_LM:
            config_lm(self, arg)
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
