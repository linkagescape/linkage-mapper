#!/usr/bin/env python2.5

"""Script to set system arguments for Linkage Mapper"""
import sys

sys.argv = ("lm_master.py", "C:/lm_test/demoProject", 
    "C:/lm_test/demoData/Cores.shp", "core_ID",
    "C:/lm_test/demoData/resistances", "true", "Cost-Weighted & Euclidean",
    "true", "Cost-Weighted & Euclidean", 
    "C:/lm_test/demoData/distances_Cores.txt", "true", "true", "false", "4", 
    "Cost-Weighted", "false", "true", "100000", "100000", "100000")