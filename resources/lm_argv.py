#!/usr/bin/env python2.5

"""Script to set system arguments for Linkage Mapper"""
import sys

sys.argv = ("lm_master", "project", "/data/Cores.shp", "core_ID",
    "/data/resistances", "true", "Cost-Weighted & Euclidean", "false",
    "C:/lm_test/project/distances_Cores.txt", "Cost-Weighted & Euclidean",
    "false", "false", "false", "4", "Cost-Weighted", "false", "false",
    "100000", "100000", "100000")