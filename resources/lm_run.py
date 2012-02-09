#!/usr/bin/env python2.5

"""Script to run Linkage Mapper"""
import sys


_SCRIPT_NAME = "lm_run.py"


def main():
    """Runs Linkage Mapper with coded inputs"""
    proj_dir = "V:\\demoProject"
    core_fc = "V:\\demoData\\Cores.shp"
    core_fl = "core_ID"
    resis_rast = "V:\\demoData\\resistances"
    # distance_file = "V:\\demoData\\distances_Cores.txt"
    
    arg = (_SCRIPT_NAME, proj_dir, core_fc, core_fl, resis_rast, "true", "true",
            "Cost-Weighted & Euclidean", "#",
            "true", "true", "false", "4", "Cost-Weighted", "false", "true",
            "100000", "100000", "100000")
    sys.argv = arg
    
    sys.path.append('../toolbox/scripts')
    import lm_master
    lm_master.lm_master()


if __name__ == "__main__":
    main()
