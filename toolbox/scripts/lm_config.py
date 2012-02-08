#!/usr/bin/env python2.
# Authors: Brad McRae and Darren Kavanagh

"""Linkage Mapper configuration module.

Assigns input parameters from ToolBox to variables, and sets constants.

"""

import os.path as path

import arcgisscripting

import lm_version as ver
import lm_settings 

GP_NULL = '#'
LINKAGE_MAPPER = 'linkage_mapper'
BARRIER_TOOL = 'barrier_mapper'
CIRCUITSCAPE = 'circuitscape'


def str2bool(pstr):
    """Convert ESRI boolean string to Python boolean type"""
    return pstr == 'true'


def setadjmeth(inparam):
    """Return boolean variables for distance methods"""
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
    """Convert ESRI float or null to Python float"""
    if innum == GP_NULL:
        nfloat = None
    else:
        nfloat = float(innum)
    return nfloat


def nullstring(arg_string):
    """Convert ESRI nullstring to Python null"""
    if arg_string == GP_NULL:
        arg_string = None
    return arg_string


def config_global(config, arg):    
    """Configure global variables for all tools"""    
    config.PARAMS = arg
    config.releaseNum = ver.releaseNum
    config.LOGMESSAGES = True
    # File names, directory paths & folder names
    proj_dir = arg[1]
    config.PROJECTDIR = proj_dir  # Project directory
    config.SCRATCHDIR = path.join(proj_dir, "scratch_cs")
    config.ARCSCRATCHDIR = path.join(config.SCRATCHDIR, "arcscratch")
    config.PREFIX = path.basename(proj_dir)
    config.DATAPASSDIR = path.join(proj_dir, "datapass")
    config.CWDADJFILE = path.join(config.DATAPASSDIR, "cwdAdj.csv")
    config.EUCADJFILE = path.join(config.DATAPASSDIR, "eucAdj.csv")
    config.OUTPUTDIR = path.join(proj_dir, "output")
    config.LOGDIR = path.join(proj_dir, "run_history")
    config.LOGDIR_OLD = path.join(proj_dir, "logFiles")
    config.logFile = None
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

    # Save individual focal grids for barrier analysis
    config.SAVEFOCALRASTERS = False

    config.SAVEBARRIERRASTERS = False  # Save individual barrier grids

    # Save individual current maps from Circuitscape
    config.SAVECURRENTMAPS = False

    config.SAVECIRCUITDIR = False
    config.SAVEBARRIERDIR = False
    config.SAVECENTRALITYDIR = False
    config.SAVECURRENTMAPS = False

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

    # 1 corridor
    # 10 NN corridor
    # 11 1st nn (future)
    # 12 2nd nn etc (future)
    # 20 constel
    # 21 1st nn constel (future)
    # 22 2nd nn constel etc (future)
    # 30 temp saved nnconstel.  maybe able to get rid of this

    #Temporary resistance raster copy to be created in master scripts
    config.RESRAST = path.join(config.SCRATCHDIR, 'resrast')


def config_lm(config, arg, scratch_dir):
    """ Configure global variables for Linkage Mapper"""
    config.TOOL = LINKAGE_MAPPER
    config.COREFC = arg[2]  # Core area feature class
    config.COREFN = arg[3]  # Core area field name
    config.RESRAST_IN = arg[4]  # Resistance raster

    # Processing steps inputs
    config.STEP1 = str2bool(arg[5])

    ### SETTING BOTH ADJ METHODS TO TRUE FOR S1 IN FUTURE RELEASES ###
    # config.S1ADJMETH_CW, config.S1ADJMETH_EU = setadjmeth(arg[6])

    config.S1ADJMETH_CW = True
    config.S1ADJMETH_EU = True

    config.STEP2 = str2bool(arg[6])
    config.S2ADJMETH_CW, config.S2ADJMETH_EU = setadjmeth(arg[7])
    config.S2EUCDISTFILE = nullstring(arg[8])
    config.STEP3 = str2bool(arg[9])
    # Drop LCC's passing through intermediate cores
    config.S3DROPLCCS = str2bool(arg[10])
    config.STEP4 = str2bool(arg[11])
    config.S4MAXNN = int(arg[12])  # No of connected nearest neighbors
    # NN Unit
    config.S4DISTTYPE_CW, config.S4DISTTYPE_EU = setadjmeth(arg[13])
    config.S4CONNECT = str2bool(arg[14])
    config.STEP5 = str2bool(arg[15])

    # Optional input parameters
    config.BUFFERDIST = nullfloat(arg[16])
    config.MAXCOSTDIST = nullfloat(arg[17])
    if config.MAXCOSTDIST == 0:
        config.MAXCOSTDIST = None
    config.MAXEUCDIST = nullfloat(arg[18])
    if config.MAXEUCDIST == 0:
        config.MAXEUCDIST = None
   
    for setting in dir(lm_settings):
        if setting == setting.upper():
            setting_value = getattr(lm_settings, setting)
            setattr(config, setting, setting_value)

    if config.MAXCOSTDIST is None:
        config.TMAXCWDIST = None
    elif config.CWDTHRESH is not None:
        config.TMAXCWDIST = None  # Max is disabled for now- see line below
        # Will limit cw calc
        # config.TMAXCWDIST = MAXCOSTDIST + CWDTHRESH

    config.SCRATCHGDB = path.join(scratch_dir, "scratch.gdb")
    # Permanent copy of core area FC made at step 1 and used in all steps
    # config.COREDIR = path.join(config.DATAPASSDIR, "corecopy")
    # config.COREFC = path.join(config.COREDIR, "core_copy.shp")
    # config.COREFN = "GRIDCODE"
    config.CORERAS = path.join(config.SCRATCHDIR, "core_ras")


def config_barrier(config, arg):
    """Configure global variables for Barrier tool"""
    config.TOOL = BARRIER_TOOL
    config.RESRAST_IN = arg[2]
    config.STARTRADIUS = arg[3]
    config.ENDRADIUS = arg[4]
    config.RADIUSSTEP = arg[5]
    if config.RADIUSSTEP == GP_NULL:
        config.RADIUSSTEP = 0
    config.STEP1 = False


def config_circuitscape(config, arg):
    """Configure global variables for Circuitscape"""
    config.TOOL = CIRCUITSCAPE
    config.COREFC = arg[2]
    config.COREFN = arg[3]
    config.DOCENTRALITY = str2bool(arg[4])
    config.DOPINCH = str2bool(arg[5])
    config.RESRAST_IN = arg[6]
    config.CWDCUTOFF = nullfloat(arg[7])  # CDW cutoff distance
    config.SQUARERESISTANCES = str2bool(arg[8])  # Square resistance values

    # Do adjacent pair corridor pinchpoint calculations using raster CWD maps
    config.DO_ADJACENTPAIRS = str2bool(arg[9])
    # Do all-pair current calculations using raster corridor map
    config.DO_ALLPAIRS = str2bool(arg[10])

    config.STEP1 = False


class Configure(object):
    """Class container to hold global variables"""
    def __init__(self):
        """Initialize class and create single geoprocessor object"""        
        self.gp = arcgisscripting.create(9.3)
        self.gp.CheckOutExtension("Spatial")
        self.gp.OverwriteOutput = True

    def configure(self, tool, arg):
        """Setup variables for Configure class"""
        config_global(self, arg)
        if tool == LINKAGE_MAPPER:
            config_lm(self, arg, self.SCRATCHDIR)
        elif tool == BARRIER_TOOL:
            config_barrier(self, arg)
        elif tool == CIRCUITSCAPE:
            config_circuitscape(self, arg)
        else:
            raise RuntimeError('Undefined tool to configure')


tool_env = Configure()  # Class instance that is use by tool modules
