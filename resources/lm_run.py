#!/usr/bin/env python2.5

"""Script to run Linkage Mapper"""
import arcgisscripting


def main():
    """Runs Linkage Mapper with coded inputs"""
    gp = arcgisscripting.create(9.3)

    gp.AddToolbox("..\\toolbox\\Connectivity Mapping Tools.tbx")
    gp.Workspace = "C:\\lm_test"

    gp.LinkageMapper(
        "project", "/data/Cores.shp", "core_ID", "/data/resistances", "true",
        "Cost-Weighted & Euclidean", "true",
        "C:/lm_test/project/distances_Cores.txt", "Cost-Weighted & Euclidean",
        "true", "true", "false", "4", "Cost-Weighted", "false", "true",
        "100000", "100000", "100000")


if __name__ == "__main__":
    main()
