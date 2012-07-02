#!/usr/bin/env python2.6

"""Script to run Climate Corridor tool"""

import sys
import os.path as path


_SCRIPT_NAME = path.basename(__file__)


def main():
    """Runs Climate Corridor tool"""
    proj_dir = "W:\\ProjectAscEq"
    core_fc = "W:\\Base Test\\PROJECT\\in_data\\cores.shp"
    core_fl = "HCA_ID"
    climate_rast = "W:\\Base Test\\PROJECT\\in_data\\climate.img"
    #resis_rast = "#"
    resis_rast = "W:\\Base Test\\PROJECT\\in_data\\resist.img"
    #gisbase = "C:\\Program Files (x86)\\GRASS 6.4.2"
    #gisbase = "C:\\Program Files (x86)\\GRASS GIS 6.5.svn"
    gisbase = "C:\Program Files (x86)\GRASS GIS 7.0.svn"
    min_euc_dist = 2000
    max_euc_dist = 50000
    climate_threashold = 1
    climate_cost = 50000

    sys.path.append('..\\toolbox\\scripts')
    import cc_main

    argv = (_SCRIPT_NAME, proj_dir, core_fc, core_fl, climate_rast, resis_rast,
            gisbase, min_euc_dist, max_euc_dist, climate_threashold,
            climate_cost)

    cc_main.main(argv)


if __name__ == "__main__":
    main()
