# Authors: Darren Kavanagh and Brad McRae

"""Create instance of config object to store input parameters and settings.

Classes:
    ClimateConfig: Container to store tool settings.

"""

import os
import subprocess
import sys

import lm_util_config as util


class ClimateConfig(object):
    """Class container to hold Climate Tool global variables."""

    def __init__(self):
        """Init class (empty)."""
        pass

    def configure(self, arg):
        """Assign input parameters and set constants."""
        # Input parameter names
        self.PARAM_NAMES = [
            'proj_dir', 'core_lyr', 'core_fld', 'cli_lyr', 'resis_lyr',
            'grass_path',
            'min_eucdist', 'max_eucdist', 'cl_threshold', 'cl_cost',
            'pn_opt', 'pn_no_nn', 'pn_nn_unit', 'pn_connect']

        # Input files
        self.proj_dir = arg[1]  # Project directory
        self.core_fc = arg[2]  # Core area feature class
        self.core_fld = arg[3]  # Core area field name
        self.climate_rast = arg[4]  # Climate raster (+ path)
        self.resist_rast = util.nullstring(arg[5])  # Resist raster (+ path)

        # Setup GRASS environmental variables
        self.gisbase = os.environ['GISBASE'] = arg[6]  # GRASS path
        sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))
        self.gpath = (''.join([os.path.join(cc_env.gisbase, gpath) + os.pathsep
                               for gpath in [r'mysys\bin', 'bin', 'extrabin',
                                             'lib', r'etc\python', 'etc']]))

        # Overwrite default startup subprocess variables to insure the console
        # window is hidden. This is necessary as functions within the GRASS
        # scripting library open subprocesses with the console window open.
        subprocess.STARTUPINFO.dwFlags |= subprocess.STARTF_USESHOWWINDOW
        subprocess.STARTUPINFO.wShowWindow = subprocess.SW_HIDE

        # Tool settings
        self.min_euc_dist = float(arg[7])  # Min distance between core pairs
        self.max_euc_dist = float(arg[8])  # Max distance between core pairs
        self.climate_threshold = float(arg[9])  # Climate threshold
        self.climate_cost = float(arg[10])  # Cost incurred by climate variable

        # Prune network settings
        # For Linkage Mapper. No type conversions necessary
        self.prune_network = arg[11]
        self.max_nn = arg[12]  # No of connected nearest neighbors
        self.nn_unit = arg[13]  # NN Unit
        self.keep_constelations = arg[14]

        # Setup model global variables
        self.scratch_dir = os.path.join(self.proj_dir, "scratch")
        self.cc_gdb = os.path.join(self.scratch_dir, "cc.gdb")

        # Define project area files
        self.prj_core_fc = os.path.join(self.cc_gdb, "cc_cores")
        self.prj_resist_rast = os.path.join(self.cc_gdb, "cc_resist")
        self.prj_climate_rast = "climate"
        self.prj_core_rast = "cores"
        self.simplify_cores = True
        self.WRITETRUNCRASTER = "true"
        self.CWDTHRESH = "200000"
        self.OUTPUTFORMODELBUILDER = "#"
        self.LMCUSTSETTINGS = "#"


cc_env = ClimateConfig()
