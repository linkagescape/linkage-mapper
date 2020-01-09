"""Script to run Linkage Priority tool.

Assumes Linkage Mapper scripts and demo data are in their default folders.

"""

import sys
import os


def in_params(demo_dir):
    """Define input paramaters."""
    gp_blank = "#"  # Geoprocessing's empty parameter value

    prj_dir = os.path.join(demo_dir,  # Folder containing demos
                           'demo_run',  # Container folder to hold model runs
                           'lm')  # Linkage Mapper model run folder

    return (
        [os.path.basename(__file__),  # Script Name ** Do not modify **
        prj_dir,  # Project Directory (1)
        os.path.join(demo_dir,
                     'demoData\\cores.shp'),  # Core Area Feature Class (2)
        "core_ID",  # Core Area Field Name (3)
        os.path.join(demo_dir,
                     'demoData\\resistances'),  # Resistance Raster (4)
        gp_blank,  # Other Core Area Value Raster (5)
        0.33,  # Resistance Weight in CAV Calculation (6)
        0.33,  # Size Weight in CAV Calculation (7)
        0.34,  # Area/Perimeter Weight in CAV Calculation (8)
        0,  # Expert Core Area Value Weight in CAV Calculation (9)
        0,  # Current Flow Centrality Weight in CAV Calculation (10)
        0,  # Other Core Area Value Weight in CAV Calculation (11)
        gp_blank,  # Core Pairs Table (12)
        gp_blank,  # From Core Field (13)
        gp_blank,  # To Core Field (14)
        gp_blank,  # Expert Corridor Importance Value Field (15)
        gp_blank,  # Current Climate Envelope Raster (16)
        "false",  # Modify the Advanced Climate Signature Parameters? (17)
        gp_blank,  # Future Climate Envelope Raster (18)
        0.5,  # Relative Priority of Minimum Climate Analog Ratio (19)
        0,  # Relative Priority of Maximum Climate Analog Ratio (20)
        1.15,  # Lowest Allowable Maximum Climate Analog Ratio (21)
        1,  # Targeted Climate Analog Ratio (22)
        1,  # Relative Priority of Achieving the Targeted Ratio (23)
        0.5,  # Climate Analog Linkage Priority Weight (24)
        351,  # Preferred Climate Signature Value for a Core (25)
        0.5,  # Relative Priority of Min Climate Pref Attainment Ratio (26)
        0,  # Relative Priority of Max Climate Pref Attainment Ratio (27)
        0,  # Climate Preference Linkage Priority Weight (28)
        0.33,  # Closeness Weight in CSP Calculation (29)
        0.33,  # Permeability Weight in CSP Calculation (30)
        0.34,  # Core Area Value Weight in CSP Calculation (31)
        0,  # Expert Corridor Importance Value Weight in CSP Calculation (32)
        0,  # Climate Envelope Difference Weight in CSP Calculation (33)
        0.6,  # CSP Threashold Value (34)
        0.5,  # Truncated Corridors Weight in Blended Priority Calc (35)
        0.5,  # Linkage Priority Weight in Blended Priority Calculation (36)
        gp_blank,  # Output for ModelBuilder Precondition (37)
        gp_blank]  # Custom Settings File (38)
        )


def main():
    """Set path and run model."""
    demo_path = (os.path.abspath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')))

    sys.path.append(os.path.join(demo_path, '..\\toolbox\\scripts'))
    import lp_main
    lp_main.main(in_params(demo_path))


if __name__ == "__main__":
    main()
