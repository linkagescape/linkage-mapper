"""Upgrade v1.x logfile to v2.0 logfile."""

import os
import sys
import glob
from ast import literal_eval

import arcpy


FILE_PREFIX = "_v2_Linkage"


class LogFileException(Exception):
    """Custom exception."""

    pass


def get_log_file(proj_dir):
    """Get last log file."""
    spath = os.path.join(proj_dir, "run_history", "log",
                         "*_Linkage Mapper.txt")
    entries = sorted(glob.glob(spath), key=os.path.getctime, reverse=True)
    last_lm_log = next(iter(entries or []), None)
    if not last_lm_log:
        raise LogFileException("Log file was not found.")
    if FILE_PREFIX in last_lm_log:
        raise LogFileException("Upgraded log file already exits.")
    return last_lm_log


def update_params(param_line):
    """Update parameter list from  v1.x to v2.0."""
    pline_heading = "Parameters:\t"

    if param_line[0:12] != pline_heading:
        raise LogFileException("Log file does not contain a Parameters line")

    lm_arg = literal_eval(param_line[12:len(param_line) - 1])

    lm_arg.append(lm_arg[17])  # MAXCOSTDIST (arg[18])
    lm_arg.append(lm_arg[18])  # MAXEUCDIST (arg[19])
    lm_arg.append('#')  # OUTPUTFORMODELBUILDER (arg[20])
    lm_arg.append('#')  # LMCUSTSETTINGS (arg[21])

    lm_arg[18] = lm_arg[16]  # BUFFERDIST
    lm_arg[17] = 200000  # CWDTHRESH
    lm_arg[16] = 'true'  # WRITETRUNCRASTER

    return "{}{}".format(pline_heading, lm_arg)


def main(argv=None):
    """Upgrade v1.x logfile to v2.0 logfile."""
    if argv is None:
        argv = sys.argv  # Get parameters from ArcGIS tool dialog

    proj_dir = argv[1]

    try:
        last_lm_log = get_log_file(proj_dir)
        new_log_file = last_lm_log.replace("_Linkage", FILE_PREFIX)

        lines = open(last_lm_log).read().splitlines()
        lines[4] = update_params(lines[4])
        open(new_log_file, 'w').write('\n'.join(lines))
    except LogFileException as err:
        arcpy.AddError(err)
    else:
        arcpy.AddMessage(
            "New Log file {} saved.".format(os.path.basename(new_log_file)))


if __name__ == "__main__":
    main()
