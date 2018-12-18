# Authors: John Gallo and Randal Greene 2017

"""Linkage Priority configuration module."""

import os.path as path

import lm_util_config as util


class PriorityConfig(object):
    """Class container to hold Linkage Priority parameters and settings."""

    def __init__(self):
        """Init class (empty)."""
        pass

    def configure(self, arg):
        """Assign parameters and settings.

        Assign parameters and settings from passed arguments and advanced
        settings.
        """
        self.PARAMS = str(arg)  # Convert to string in case '\' exists

        # Model Inputs
        # ------------
        self.PROJDIR = arg[1]
        self.COREFC = arg[2]
        self.COREFN = arg[3]
        self.RESRAST_IN = arg[4]

        # Core Area Value (CAV) Options
        # -----------------------------
        self.OCAVRAST_IN = util.nullstring(arg[5])
        self.RESWEIGHT = float(arg[6])
        self.SIZEWEIGHT = float(arg[7])
        self.APWEIGHT = float(arg[8])
        self.ECAVWEIGHT = float(arg[9])
        self.CFCWEIGHT = float(arg[10])
        self.OCAVWEIGHT = float(arg[11])

        # Corridor Specific Priority (CSP) Options
        # ----------------------------------------
        #  Expert Corridor Importance Vale
        self.COREPAIRSTABLE_IN = util.nullstring(arg[12])
        self.FROMCOREFIELD = util.nullstring(arg[13])
        self.TOCOREFIELD = util.nullstring(arg[14])
        self.ECIVFIELD = util.nullstring(arg[15])

        # Climate Linkage Priority Value
        self.CCERAST_IN = util.nullstring(arg[16])
        self.FCERAST_IN = util.nullstring(arg[18])
        self.CANALOG_MIN = float(arg[19])
        self.CANALOG_MAX = float(arg[20])
        self.CANALOG_TARGET = float(arg[21])
        self.CANALOG_PIORITY = float(arg[22])
        self.CANALOG_WEIGHT = float(arg[23])
        self.CPREF_VALUE = float(arg[24])
        self.CPREF_MIN = float(arg[25])
        self.CPREF_MAX = float(arg[26])
        self.CPREF_WEIGHT = float(arg[27])

        # CSP Weights
        self.CLOSEWEIGHT = float(arg[28])
        self.PERMWEIGHT = float(arg[29])
        self.CAVWEIGHT = float(arg[30])
        self.ECIVWEIGHT = float(arg[31])
        self.CEDWEIGHT = float(arg[32])
        self.PROPCSPKEEP = float(arg[33])

        # Blended Priority Options
        # ------------------------
        self.TRUNCWEIGHT = float(arg[34])
        self.LPWEIGHT = float(arg[35])

        # Additional Options
        # ------------------
        self.OUTPUTFORMODELBUILDER = util.nullstring(arg[36])
        if arg[37] == util.GP_NULL:
            self.LPCUSTSETTINGS_IN = path.join(util.get_code_path(),
                                               'lp_settings.py')
        else:
            self.LPCUSTSETTINGS_IN = arg[37]

        # - - - - - - - - - - - - - - - - - -

        # core corename from feature class name
        splits = self.COREFC.split("\\")
        self.CORENAME = splits[len(splits) - 1].split(".")[0]

        util.set_custom(self.LPCUSTSETTINGS_IN, self)


lp_env = PriorityConfig()
