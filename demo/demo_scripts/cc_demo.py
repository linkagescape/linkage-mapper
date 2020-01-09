"""Script to run Climate Linkage Mapper tool.

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
    prj_dir = os.path.join(demo_dir,  # Folder containing demos
                           'demo_run',  # Container folder to hold model runs
                           'cc')  # Climate Linkage Mapper model run folder
    create_dir(prj_dir)

    return (
        [os.path.basename(__file__),  # Script Name ** Do not modify **
        prj_dir,  # Project Folder (1)
        os.path.join(demo_dir,
                     'demoData\\cc_cores.shp'),  # Core Area Feature Class (2)
        "HCA_ID",  # Core Area Field Name (3)
        os.path.join(demo_dir,
                     'demoData\\cc_climate.img'),  # Climate Raster (4)
        os.path.join(demo_dir,
                     'demoData\\resistances'),  # Resistance Raster (5)
        "C:\\Program Files\\GRASS GIS 7.4.0",  # GRASS GIS Folder (6)
        2000,  # Minium Distance Between Core Pairs (7)
        50000,  # Maxium Distance Between Core Pairs (8)
        1,  # Climate Threashold (9)
        50000,  # Climate Variable Cost (10)
        'true',  # Prune Network (11)
        4,  # Number of Connected Nearest Neighbors (12)
        "Cost-Weighted",  # Nearest Neighbor Measurement Unit (13)
        "true"]  # Connect Neighboring Constellations (14)
        )


def main():
    """Set path and run model."""
    demo_path = (os.path.abspath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')))

    sys.path.append(os.path.join(demo_path, '..\\toolbox\\scripts'))
    import cc_main
    cc_main.main(in_params(demo_path))


if __name__ == "__main__":
    main()
