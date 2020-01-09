"""Script to run Pinchpoint Mapper tool.

Assumes Linkage Mapper scripts are in their default folders.

"""

import sys
import os


def in_params(demo_dir):
    """Define input paramaters."""
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
        5000,  # CWD cutoff distance (5)
        "false",  # Square resistance values? (6)
        "true",  # Calculate adjacent par pinch points using Circuitscape (7)
        "true",  # Calculate raster centrality using Circuitscape (8)
        "All-to-one"]  # Circuitscape mode for raster centrality (9)
        )


def main():
    """Set path and run model."""
    demo_path = (os.path.abspath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')))

    sys.path.append(os.path.join(demo_path, '..\\toolbox\\scripts'))
    import circuitscape_master
    circuitscape_master.circuitscape_master(in_params(demo_path))


if __name__ == "__main__":
    main()
