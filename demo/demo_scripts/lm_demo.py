"""Script to run Linkage Pathways Tool.

Assumes Linkage Mapper scripts and demo data are in their default folders.

"""

import sys
import os


def create_dir(in_dir):
    """Create directory if doesn't exist."""
    if not os.path.exists(in_dir):
        os.makedirs(in_dir)


def in_params(demo_dir):
    """Define input paramaters."""
    gp_blank = '#'  # Geoprocessing's empty parameter value

    prj_dir = os.path.join(demo_dir,  # Folder containing demos
                           'demo_run',  # Container folder to hold model runs
                           'lm')  # Linkage Mapper model run folder
    create_dir(prj_dir)

    return (
        [os.path.basename(__file__),  # Script Name ** Do not modify **
        prj_dir,  # Project Directory (1)
        os.path.join(demo_dir,
                     'demoData\\cores.shp'),  # Core Area Feature Class (2)
        "core_ID",  # Core Area Field Name (3)
        os.path.join(demo_dir,
                     'demoData\\resistances'),  # Resistance Raster (4)
        "true",  # Step 1 - Identify Adjacent Core Areas (5)
        "true",  # Step 2 - Construct a Network of Core Areas (6)
        "Cost-Weighted & Euclidean",  # Network Adjacency Method (7)
        gp_blank,  # Core Area Distances Text File (8)
        "true",  # Step 3 - Calculate CWD and LCP (9)
        "true",  # Drop Corridors that Intersect Core Areas (10)
        "true",  # Step 4 - Prune Network (11)
        4,  # Option A - Maxium Number of Connected NN (12)
        "Cost-Weighted",  # Option B - NN Measurement Unit (13)
        "true",  # Option C - Connect NC (14)
        "true",  # Step 5 - Calculate, Normalize and Mosaic Corridors (15)
        "true",  # Truncate Corridors (16)
        200000,  # CWD Threashold to Use in Truncating Corridors (17)
        10000,  # Bounding Circles Buffer Distance(18)
        100000,  # Maximum Cost-Weighted Corridor Distance (19)
        40000,  # Maximum Euclidean Corridor Distance (20)
        gp_blank,  # Output for ModelBuilder Precondition (21)
        gp_blank]  # Custom Settings File (22)
        )


def main():
    """Set path and run model."""
    demo_path = (os.path.abspath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')))

    sys.path.append(os.path.join(demo_path, '..\\toolbox\\scripts'))
    import lm_master
    lm_master.lm_master(in_params(demo_path))


if __name__ == "__main__":
    main()
